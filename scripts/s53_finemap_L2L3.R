# -*- coding: utf-8 -*-
# s53_finemap_L2L3.R -- 任务A：L2 (chr10 YME1L1) / L3 (chr2 ANXA4) 的 susie_rss 精细定位
#
# 严格复用 s30_finemap_chr1.R 的方法与代码路径，仅将下列硬编码点参数化：
#   LOCUS / CHROM / LO37 / HI37 / GENES(基因体) / 输出文件名
#   以及 pair 表里 pip_mass_* 由「固定三基因」改为「按位点基因表动态生成」
#
# 方法（与 s30 一致）：
#   * 变异集与 harmonise 同 s28/s30（allele-aware key on pos38 + 排序等位对）
#   * LD = OneK1K 980 供者 PLINK bfile，--indep-pairwise 200 50 0.9，--r square（有符号 r）
#   * susieR::susie_rss(z, R, n, L=10, estimate_residual_variance=FALSE, check_prior=TRUE)
#   * outcome-only 精细定位用「bim × GWAS 交集」，丢弃回文变异（保守）
#
# Usage: Rscript s53_finemap_L2L3.R L2|L3 [prune_r2=0.9]

suppressPackageStartupMessages({
  library(data.table)
  library(susieR)
})

PROJ  <- "D:/endometriosis_project/11_sc_eqtl_mr_project"
PREP  <- file.path(PROJ, "00_data_raw/onek1k/coloc_prep")
TAB   <- file.path(PROJ, "tables")
WORK  <- file.path(PROJ, "00_data_raw/onek1k/finemap_work_L2L3")
PLINK <- "D:/endometriosis_project/_tools/plink.exe"
BFILE <- "D:/endometriosis_project/_onek1k_plink/plink_merged_980_donors"
dir.create(WORK, showWarnings = FALSE, recursive = TRUE)

args <- commandArgs(trailingOnly = TRUE)
LOCUS <- if (length(args) >= 1) args[1] else "L2"
PRUNE_R2 <- if (length(args) >= 2) as.numeric(args[2]) else 0.9
N_OVR <- if (length(args) >= 3) as.integer(args[3]) else NA_integer_
TAG <- if (PRUNE_R2 == 0.9) "" else sprintf("_r%s", gsub("\\.", "", as.character(PRUNE_R2)))
if (!is.na(N_OVR)) TAG <- paste0(TAG, sprintf("_n%d", N_OVR))

LOG <- sprintf("D:/endometriosis_project/_s53_finemap_%s.log", LOCUS)
con <- file(LOG, open = "wt", encoding = "UTF-8")
sink(con, split = TRUE); sink(con, type = "message")

# ---- 位点配置（基因体 = Ensembl GRCh37 REST /overlap/region，2026-09-28 实查）----
CFG <- list(
  L2 = list(chrom="10", lo37=26900000L, hi37=27900000L, target="YME1L1",
            genes=data.table(
              gene  = c("PDSS1","ABI1","ANKRD26","YME1L1","MASTL","ACBD5","PTCHD3","RAB18"),
              ensid = c("ENSG00000148459","ENSG00000136754","ENSG00000107890",
                        "ENSG00000136758","ENSG00000120539","ENSG00000107897",
                        "ENSG00000182077","ENSG00000099246"),
              start37 = c(26986588L,27035522L,27280843L,27399383L,27443753L,27484146L,27687116L,27793197L),
              end37   = c(27035727L,27150016L,27389421L,27444195L,27475853L,27531059L,27703297L,27831143L))),
  L3 = list(chrom="2", lo37=69300000L, hi37=70200000L, target="ANXA4",
            genes=data.table(
              gene  = c("ANTXR1","GFPT1","NFU1","AAK1","ANXA4","GMCL1","SNRNP27","MXD1","ASPRV1","PCBP1-AS1"),
              ensid = c("ENSG00000169604","ENSG00000198380","ENSG00000169599",
                        "ENSG00000115977","ENSG00000196975","ENSG00000087338",
                        "ENSG00000124380","ENSG00000059728","ENSG00000244617","ENSG00000179818"),
              start37 = c(69240310L,69546905L,69622882L,69688532L,69871557L,70056774L,70120692L,70124820L,70187226L,70189395L),
              end37   = c(69476459L,69614382L,69664760L,69901481L,70053596L,70108528L,70132707L,70170077L,70189397L,70315978L)))
)
if (!LOCUS %in% names(CFG)) stop("unknown locus")
C0 <- CFG[[LOCUS]]
CHROM <- C0$chrom; LO37 <- C0$lo37; HI37 <- C0$hi37; GENES <- C0$genes

cat("=== s53 finemap", LOCUS, "chr", CHROM, "@", format(Sys.time()), "===\n")
cat("prune r2 =", PRUNE_R2, " tag =", TAG, " window GRCh37:", LO37, "-", HI37, "\n\n")
print(GENES)

N_EQTL <- if (!is.na(N_OVR)) N_OVR else 980L

annot_gene <- function(p37) {
  out <- rep(NA_character_, length(p37))
  for (i in seq_len(nrow(GENES)))
    out[p37 >= GENES$start37[i] & p37 <= GENES$end37[i]] <- GENES$gene[i]
  need <- which(is.na(out))
  for (k in need) {
    x <- p37[k]
    d <- ifelse(x < GENES$start37, GENES$start37 - x,
                ifelse(x > GENES$end37, x - GENES$end37, 0L))
    j <- which.min(d)
    side <- if (x > GENES$end37[j]) "flank3" else "flank5"
    out[k] <- sprintf("%s_%s_%dkb", side, GENES$gene[j], round(d[j] / 1000))
  }
  out
}
gene_of <- function(lab) sub("_[0-9]+kb$", "", sub("^flank[35]_", "", lab))

ld_signed <- function(snps37, tag, pval = NULL) {
  f_in  <- file.path(WORK, sprintf("prune%s_%s.in", TAG, tag))
  f_out <- file.path(WORK, sprintf("prune%s_%s", TAG, tag))
  fwrite(data.table(V1 = snps37), f_in, col.names = FALSE, sep = "\t")
  cmd1 <- sprintf('"%s" --bfile "%s" --allow-no-sex --chr %s --extract "%s" --indep-pairwise 200 50 %s --out "%s"',
                  PLINK, BFILE, CHROM, f_in, PRUNE_R2, f_out)
  system(cmd1, ignore.stdout = TRUE, ignore.stderr = TRUE)
  keep <- file.path(WORK, sprintf("prune%s_%s.prune.in", TAG, tag))
  if (!file.exists(keep)) return(NULL)
  k <- fread(keep, header = FALSE)$V1
  if (length(k) > 1000 && !is.null(pval)) {
    dt <- data.table(snp = snps37, p = pval)
    k <- dt[snp %in% k & is.finite(p)][order(p)][1:1000]$snp
    k <- k[!is.na(k)]
  }
  if (length(k) < 2) return(NULL)
  f_k <- file.path(WORK, sprintf("keep%s_%s.txt", TAG, tag))
  fwrite(data.table(V1 = k), f_k, col.names = FALSE, sep = "\t")
  ld_o <- file.path(WORK, sprintf("ld%s_%s", TAG, tag))
  cmd2 <- sprintf('"%s" --bfile "%s" --allow-no-sex --chr %s --extract "%s" --write-snplist --r square --out "%s"',
                  PLINK, BFILE, CHROM, f_k, ld_o)
  system(cmd2, ignore.stdout = TRUE, ignore.stderr = TRUE)
  f_ld <- paste0(ld_o, ".ld"); f_sp <- paste0(ld_o, ".snplist")
  if (!file.exists(f_ld) || !file.exists(f_sp)) return(NULL)
  m <- as.matrix(fread(f_ld, header = FALSE))
  ids <- fread(f_sp, header = FALSE)$V1
  if (nrow(m) != length(ids) || nrow(m) < 2) return(NULL)
  dimnames(m) <- list(ids, ids)
  ok <- apply(m, 1, function(z) all(is.finite(z)))
  m <- m[ok, ok, drop = FALSE]
  if (nrow(m) < 2) return(NULL)
  m
}

harmonise <- function(e, b, gw, cell = NULL) {
  if (!is.null(cell)) e <- e[cell_type == cell]
  if (nrow(e) == 0) return(NULL)
  e <- merge(e, b, by.x = "snp_grch37", by.y = "snp", all.x = TRUE)
  e <- e[!is.na(a1) & !is.na(a2)]
  e[, pos38 := as.integer(sub("^[^:]+:", "", snp_grch38))]
  e <- e[!is.na(pos38)]
  e[, a1 := toupper(a1)]; e[, a2 := toupper(a2)]
  e[, allele_key := paste0(pos38, "_", pmin(a1, a2), "_", pmax(a1, a2))]
  g <- copy(gw)
  g[, ea := toupper(ea)]; g[, oa := toupper(oa)]
  g <- g[!is.na(pos38) & !is.na(beta) & !is.na(se) & se > 0]
  g[, allele_key := paste0(pos38, "_", pmin(ea, oa), "_", pmax(ea, oa))]
  g <- g[!duplicated(allele_key)]
  m <- merge(e, g, by = "allele_key", suffixes = c("_e", "_o"))
  if ("pos38_e" %in% names(m)) { setnames(m, "pos38_e", "pos38"); m[, pos38_o := NULL] }
  if (nrow(m) < 10) return(NULL)
  m <- m[(ea == a1 & oa == a2) | (ea == a2 & oa == a1)]
  if (nrow(m) < 10) return(NULL)
  m[, flip := (ea == a2) & (oa == a1)]
  m[, palindromic := paste0(pmin(a1, a2), pmax(a1, a2)) %in% c("AT", "CG")]
  m[, af_mism := abs(af - eaf)]; m[, af_flip := abs(af - (1 - eaf))]
  m <- m[!(palindromic == TRUE & pmin(af_mism, af_flip) > 0.05 & abs(af_mism - af_flip) < 0.05)]
  m[palindromic == TRUE & af_flip < af_mism, flip := !flip]
  m[, b2 := ifelse(flip, -beta, beta)]
  m[, maf_e := pmin(af, 1 - af)]; m[, maf_o := pmin(eaf, 1 - eaf)]
  m <- m[is.finite(maf_e) & is.finite(maf_o) & maf_e > 0 & maf_o > 0 &
           is.finite(slope) & is.finite(slope_se) & slope_se > 0 &
           is.finite(b2) & is.finite(se) & se > 0]
  m <- m[!is.na(rsid) & rsid != ""]
  if (nrow(m) < 10) return(NULL)
  m <- m[!duplicated(rsid)]
  m[, p37 := as.integer(sub("^[^:]+:", "", snp_grch37))]
  m[, N_out := n_cases + n_controls]
  setorder(m, pos38)
  m
}

run_susie <- function(beta, se, R, n, L = 10) {
  z <- beta / se
  z[!is.finite(z)] <- 0
  wmsg <- character(0)
  r <- tryCatch(
    withCallingHandlers(
      susie_rss(z = z, R = R, n = n, L = L,
                estimate_residual_variance = FALSE, check_prior = TRUE),
      warning = function(w) { wmsg <<- c(wmsg, conditionMessage(w)); invokeRestart("muffleWarning") }),
    error = function(e) list(.error = conditionMessage(e)))
  if (is.null(r$.error)) r$.warn <- paste(unique(wmsg), collapse = " | ") else r$.warn <- ""
  r
}
warn_flag <- function(r) {
  if (!is.null(r$.error)) return("error")
  w <- r$.warn
  if (is.null(w) || is.na(w) || w == "") return("ok")
  if (grepl("did not converge", w)) return("IBSS_NOT_CONVERGED")
  if (grepl("prior variance", w)) return("PRIOR_VAR_ZERO")
  paste0("warn:", substr(w, 1, 60))
}
mx <- function(x) { x <- x[is.finite(x)]; if (length(x) == 0) NA_real_ else max(x) }
argmax1 <- function(x) { x[!is.finite(x)] <- -Inf; if (length(x) == 0) NA_integer_ else which.max(x) }

# ---------------- load ----------------
disc_all <- fread(file.path(TAB, "15_discovery_significant_with_locus.csv"))
disc <- disc_all[locus_id == LOCUS]
mr31 <- fread(file.path(TAB, "31_final_MR_with_method.csv"))
setnames(mr31, old = c("p", "p_outcome"), new = c("p_mr", "p_outcome_mr"), skip_absent = TRUE)
disc <- merge(disc, mr31[method == "Wald ratio", .(locus_id, gene, cell_type,
                                                  instrument_grch38 = instrument,
                                                  p_mr, F_stat_mr = F_stat)],
              by = c("locus_id", "gene", "cell_type"), all.x = TRUE)
cat("\n", LOCUS, " pandas pairs =", nrow(disc), "\n"); print(disc[, .(gene, cell_type, instrument_grch38, F_stat_mr)])
instr37_map <- unique(fread(file.path(TAB, "01b_discovery_instruments_grch38.csv"))[
  , .(gene, cell_type, pos37, pos38, variant_id, variant_id_grch37)])
print(instr37_map[gene %in% GENES$gene][order(gene, cell_type)])

# bim_universe 只有 GRCh37 位置；GRCh38 映射借用目标基因的 eQTL 文件（同 s30 做法）
b  <- fread(file.path(PREP, sprintf("bim_universe_%s.csv", LOCUS)))
gw <- fread(file.path(PREP, sprintf("gwas_%s.csv", LOCUS)))
b  <- b[pos >= LO37 & pos <= HI37]
cat("bim universe in window =", nrow(b), " ; gwas rows =", nrow(gw), "\n")

OUT_V <- file.path(TAB, sprintf("52_finemap_variant%s_%s.csv", TAG, LOCUS))
OUT_P <- file.path(TAB, sprintf("52b_finemap_pair%s_%s.csv", TAG, LOCUS))
OUT_G <- file.path(TAB, sprintf("52c_finemap_pipmass%s_%s.csv", TAG, LOCUS))

# ============================== PART 1: outcome-only ==============================
cat("\n========== PART 1: outcome-only fine-mapping ==========\n")
g0 <- copy(gw)
g0[, ea := toupper(ea)]; g0[, oa := toupper(oa)]
g0 <- g0[!is.na(pos38) & !is.na(beta) & !is.na(se) & se > 0]
bb <- data.table(snp37 = b$snp, pos37 = as.integer(b$pos),
                 A1 = toupper(b$a1), A2 = toupper(b$a2))
eq_map <- unique(fread(file.path(PREP, sprintf("eqtl_%s_%s.csv", LOCUS, C0$target)))[
  , .(snp_grch37, snp_grch38)])
bb <- merge(bb, eq_map, by.x = "snp37", by.y = "snp_grch37", all.x = TRUE)
bb[, pos38_b := as.integer(sub("^[^:]+:", "", snp_grch38))]
bb <- bb[!is.na(pos38_b)]
bb[, allele_key := paste0(pos38_b, "_", pmin(A1, A2), "_", pmax(A1, A2))]
g0[, allele_key := paste0(pos38, "_", pmin(ea, oa), "_", pmax(ea, oa))]
g0 <- g0[!duplicated(allele_key)]
o1 <- merge(bb[, .(snp37, pos37, A1, A2, pos38 = pos38_b, allele_key)], g0,
            by = "allele_key", suffixes = c("_b", "_o"))
if ("pos38_b" %in% names(o1)) setnames(o1, "pos38_b", "pos38")
if ("pos38_o" %in% names(o1)) o1[, pos38_o := NULL]
if (!"pos38" %in% names(o1)) stop("pos38 column lost after merge")
o1 <- o1[(ea == A1 & oa == A2) | (ea == A2 & oa == A1)]
o1[, flip := (ea == A2) & (oa == A1)]
o1[, palindromic := paste0(pmin(A1, A2), pmax(A1, A2)) %in% c("AT", "CG")]
n_pal <- nrow(o1[palindromic == TRUE])
o1 <- o1[palindromic == FALSE]
o1[, bo := ifelse(flip, -beta, beta)]
o1[, p37 := pos37]
o1 <- o1[!duplicated(rsid)]
setorder(o1, pos38)
cat("outcome-only intersection =", nrow(o1), "(dropped palindromic =", n_pal, ")\n")

var_rows <- list(); pair_rows <- list(); gm_rows <- list()
LD1 <- ld_signed(o1$snp37, "outcome", o1$p)
if (is.null(LD1)) {
  cat("!! outcome LD unavailable\n")
} else {
  cat("outcome LD =", nrow(LD1), "x", ncol(LD1), "\n")
  ids <- rownames(LD1); oo <- o1[match(ids, snp37)]
  nn <- as.numeric(median(oo$n_cases + oo$n_controls))
  cat("outcome n (median) =", nn, "\n")
  dimnames(LD1) <- list(oo$rsid, oo$rsid)
  r <- run_susie(oo$bo, oo$se, LD1, nn, L = 10)
  cat("outcome converge flag:", warn_flag(r), "\n")
  if (!is.null(r$.error)) {
    cat("!! outcome susie_rss error:", r$.error, "\n")
  } else {
    pip <- susie_get_pip(r); cs <- r$sets$cs
    cat("outcome: niter =", r$niter, " nCS =", length(cs), " max PIP =", round(mx(pip), 4), "\n")
    oo[, PIP := pip]; oo[, gene_annot := annot_gene(p37)]; oo[, in_cs := ""]
    if (length(cs) > 0) for (j in seq_along(cs)) {
      idx <- cs[[j]]; nm <- names(cs)[j]
      oo[idx, in_cs := ifelse(in_cs == "", nm, paste(in_cs, nm, sep = ";"))]
      lead <- idx[argmax1(pip[idx])]
      cat(sprintf("  CS %s: size=%d lead=%s (GRCh37:%d, %s) PIP_lead=%.4f\n",
                  nm, length(idx), oo$rsid[lead], oo$p37[lead], oo$gene_annot[lead], pip[lead]))
    }
    var_rows[[1]] <- oo[, .(scope = "outcome_only", side = "outcome", rsid, snp37, pos37 = p37,
                            pos38, gene_annot, PIP, in_cs, beta = bo, se)]
    gm0 <- oo[, .(PIP_mass = sum(PIP, na.rm = TRUE), n_var = .N), by = .(gene_region = gene_of(gene_annot))]
    gm_rows[[1]] <- gm0[, .(scope = "outcome_only", side = "outcome", gene_region, PIP_mass, n_var)]
    pair_rows[[1]] <- data.table(
      scope = "outcome_only", nsnp = nrow(oo), niter = r$niter, n_cs = length(cs),
      converge = warn_flag(r), max_pip = mx(pip),
      argmax_rsid = oo$rsid[argmax1(pip)], argmax_pos37 = oo$p37[argmax1(pip)],
      argmax_gene = oo$gene_annot[argmax1(pip)], n_palindromic = n_pal,
      mr_instrument = NA_character_, mr_instrument_in_set = NA)
    cat("-- outcome PIP mass by gene region --\n"); print(gm_rows[[1]][order(-PIP_mass)])
  }
}

# ============================== PART 2: per pair ==============================
cat("\n========== PART 2: eQTL fine-mapping per pair ==========\n")
f_e_cache <- list()
for (i in seq_len(nrow(disc))) {
  g <- disc$gene[i]; ct <- disc$cell_type[i]
  tag <- sprintf("%s_%s", g, ct)
  cat(sprintf("\n[%02d/%02d] %s / %s\n", i, nrow(disc), g, ct))
  f_e <- file.path(PREP, sprintf("eqtl_%s_%s.csv", LOCUS, g))
  if (!file.exists(f_e)) { cat("  ! missing eqtl file\n"); next }
  if (is.null(f_e_cache[[g]])) f_e_cache[[g]] <- fread(f_e)
  m <- harmonise(f_e_cache[[g]], b, gw, cell = ct)
  if (is.null(m)) { cat("  ! harmonise failed\n"); next }
  cat("  harmonised nsnp =", nrow(m), "\n")
  LD <- ld_signed(m$snp_grch37, tag, m$pval_nominal)
  if (is.null(LD)) { cat("  ! LD unavailable\n"); next }
  ids <- rownames(LD); mm <- m[match(ids, snp_grch37)]
  cat("  LD =", nrow(LD), "x", ncol(LD), " ; N_out median =", round(median(mm$N_out)), "\n")
  dimnames(LD) <- list(mm$rsid, mm$rsid)
  N_out_med <- as.numeric(median(mm$N_out))
  re <- run_susie(mm$slope, mm$slope_se, LD, N_EQTL, L = 10)
  ro <- run_susie(mm$b2, mm$se, LD, N_out_med, L = 10)
  cat("  converge: eQTL =", warn_flag(re), " | outcome =", warn_flag(ro), "\n")
  rec <- function(r) {
    if (is.null(r$.error)) return(list(pip = susie_get_pip(r), cs = r$sets$cs, niter = r$niter, err = NA_character_))
    list(pip = rep(NA_real_, nrow(mm)), cs = NULL, niter = NA_integer_, err = r$.error)
  }
  E <- rec(re); O <- rec(ro)
  mm[, PIP_eqtl := E$pip]; mm[, PIP_outcome := O$pip]
  mm[, gene_annot := annot_gene(p37)]
  mm[, in_cs_eqtl := ""]; mm[, in_cs_outcome := ""]
  lead_e <- lead_o <- NA_integer_
  for (nm in names(E$cs)) { idx <- E$cs[[nm]]
    mm[idx, in_cs_eqtl := ifelse(in_cs_eqtl == "", nm, paste(in_cs_eqtl, nm, sep = ";"))]
    lead_e <- idx[argmax1(mm$PIP_eqtl[idx])] }
  for (nm in names(O$cs)) { idx <- O$cs[[nm]]
    mm[idx, in_cs_outcome := ifelse(in_cs_outcome == "", nm, paste(in_cs_outcome, nm, sep = ";"))]
    lead_o <- idx[argmax1(mm$PIP_outcome[idx])] }
  cat("  eQTL: nCS =", length(E$cs), " maxPIP =", round(mx(mm$PIP_eqtl), 4), "\n")
  cat("  out : nCS =", length(O$cs), " maxPIP =", round(mx(mm$PIP_outcome), 4), "\n")
  if (!is.na(lead_e)) cat(sprintf("    eQTL CS lead: %s GRCh37:%d %s PIP=%.4f\n",
      mm$rsid[lead_e], mm$p37[lead_e], mm$gene_annot[lead_e], mm$PIP_eqtl[lead_e]))
  if (!is.na(lead_o)) cat(sprintf("    out  CS lead: %s GRCh37:%d %s PIP=%.4f\n",
      mm$rsid[lead_o], mm$p37[lead_o], mm$gene_annot[lead_o], mm$PIP_outcome[lead_o]))

  instr38 <- disc$instrument_grch38[i]
  instr37 <- NA_integer_; instr37_id <- NA_character_; instr_in_set <- FALSE
  r2_inst_leadE <- NA_real_; r2_inst_leadO <- NA_real_
  im <- instr37_map[gene == g & cell_type == ct]
  if (nrow(im) >= 1) {
    instr37 <- as.integer(im$pos37[1]); instr37_id <- as.character(im$variant_id_grch37[1])
    instr_in_set <- instr37_id %in% mm$snp_grch37
    if (instr_in_set) {
      rid <- mm$rsid[match(instr37_id, mm$snp_grch37)]
      if (!is.na(lead_e)) r2_inst_leadE <- LD[rid, mm$rsid[lead_e]]^2
      if (!is.na(lead_o)) r2_inst_leadO <- LD[rid, mm$rsid[lead_o]]^2
    }
    cat(sprintf("    MR instrument %s -> %s (GRCh37:%d, %s) in_set=%s r2(inst,CS_eQTL)=%s r2(inst,CS_out)=%s\n",
                instr38, instr37_id, instr37, annot_gene(instr37), instr_in_set,
                ifelse(is.na(r2_inst_leadE), "NA", sprintf("%.3f", r2_inst_leadE)),
                ifelse(is.na(r2_inst_leadO), "NA", sprintf("%.3f", r2_inst_leadO))))
  }

  gm_e <- mm[, .(PIP_mass = sum(PIP_eqtl, na.rm = TRUE), n_var = .N), by = gene_of(gene_annot)]
  gm_o <- mm[, .(PIP_mass = sum(PIP_outcome, na.rm = TRUE), n_var = .N), by = gene_of(gene_annot)]
  setnames(gm_e, "gene_of", "gene_region"); setnames(gm_o, "gene_of", "gene_region")
  cat("  eQTL PIP mass:\n"); print(gm_e[order(-PIP_mass)])
  cat("  out  PIP mass:\n"); print(gm_o[order(-PIP_mass)])
  gm_rows[[length(gm_rows)+1]] <- gm_e[, .(scope = tag, side = "eqtl", gene_region, PIP_mass, n_var)]
  gm_rows[[length(gm_rows)+1]] <- gm_o[, .(scope = tag, side = "outcome", gene_region, PIP_mass, n_var)]

  # --- outcome-side palindromic sensitivity（与 s30_finemap_chr1.R 同口径）---
  nopal <- which(!mm$palindromic)
  cs_nopal <- NULL; conv_nopal <- NA_character_; mpip_nopal <- NA_real_
  lead_nopal_pos37 <- NA_integer_; lead_nopal_gene <- NA_character_
  if (length(nopal) >= 20) {
    rn <- run_susie(mm$b2[nopal], mm$se[nopal], LD[nopal, nopal, drop = FALSE], N_out_med, L = 10)
    conv_nopal <- warn_flag(rn)
    if (is.null(rn$.error)) {
      pn <- susie_get_pip(rn); cs_nopal <- rn$sets$cs; mpip_nopal <- mx(pn)
      if (length(cs_nopal) > 0) {
        i0 <- sort(unique(unlist(cs_nopal)))
        lead_nopal_pos37 <- mm$p37[nopal][i0][argmax1(pn[i0])]
        lead_nopal_gene  <- mm$gene_annot[nopal][i0][argmax1(pn[i0])]
      }
      cat(sprintf("    [sens] outcome 去 %d 个回文 SNP: nsnp=%d converge=%s nCS=%d maxPIP=%s CS_lead=%s\n",
                  sum(mm$palindromic), length(nopal), conv_nopal, length(cs_nopal),
                  ifelse(is.na(mpip_nopal), "NA", sprintf("%.4f", mpip_nopal)),
                  ifelse(is.na(lead_nopal_gene), "none",
                         sprintf("%d %s", lead_nopal_pos37, lead_nopal_gene))))
    } else cat("    [sens] outcome 去回文: error ", rn$.error, "\n")
  }

  var_rows[[length(var_rows)+1]] <- mm[, .(scope = tag, side = "eqtl", rsid, snp37 = snp_grch37,
      pos37 = p37, pos38, gene_annot, PIP = PIP_eqtl, in_cs = in_cs_eqtl, beta = slope, se = slope_se)]
  var_rows[[length(var_rows)+1]] <- mm[, .(scope = tag, side = "outcome", rsid, snp37 = snp_grch37,
      pos37 = p37, pos38, gene_annot, PIP = PIP_outcome, in_cs = in_cs_outcome, beta = b2, se = se)]

  pair_rows[[length(pair_rows)+1]] <- data.table(
    scope = tag, nsnp = nrow(mm), niter_eqtl = E$niter, niter_outcome = O$niter,
    converge_eqtl = warn_flag(re), converge_outcome = warn_flag(ro),
    n_cs_eqtl = length(E$cs), n_cs_outcome = length(O$cs),
    max_pip_eqtl = mx(mm$PIP_eqtl), max_pip_outcome = mx(mm$PIP_outcome),
    cs_lead_eqtl_rsid = if (is.na(lead_e)) NA_character_ else mm$rsid[lead_e],
    cs_lead_eqtl_pos37 = if (is.na(lead_e)) NA_integer_ else mm$p37[lead_e],
    cs_lead_eqtl_gene = if (is.na(lead_e)) NA_character_ else mm$gene_annot[lead_e],
    cs_lead_outcome_rsid = if (is.na(lead_o)) NA_character_ else mm$rsid[lead_o],
    cs_lead_outcome_pos37 = if (is.na(lead_o)) NA_integer_ else mm$p37[lead_o],
    cs_lead_outcome_gene = if (is.na(lead_o)) NA_character_ else mm$gene_annot[lead_o],
    mr_instrument_grch38 = instr38, mr_instrument_grch37 = instr37,
    mr_instrument_id_grch37 = instr37_id, mr_instrument_in_finemap_set = instr_in_set,
    mr_instrument_gene = if (is.na(instr37)) NA_character_ else annot_gene(instr37),
    r2_instr_vs_cs_lead_eqtl = r2_inst_leadE,
    r2_instr_vs_cs_lead_outcome = r2_inst_leadO,
    n_palindromic = sum(mm$palindromic),
    n_cs_outcome_nopal = length(cs_nopal),
    converge_outcome_nopal = conv_nopal,
    max_pip_outcome_nopal = mpip_nopal,
    cs_lead_outcome_nopal_pos37 = lead_nopal_pos37,
    cs_lead_outcome_nopal_gene = lead_nopal_gene,
    argmax_rsid = mm$rsid[argmax1(mm$PIP_eqtl)],
    argmax_pos37 = mm$p37[argmax1(mm$PIP_eqtl)],
    argmax_gene = mm$gene_annot[argmax1(mm$PIP_eqtl)])

  if (length(var_rows) > 0) fwrite(rbindlist(var_rows, fill = TRUE), OUT_V)
  if (length(pair_rows) > 0) fwrite(rbindlist(pair_rows, fill = TRUE), OUT_P)
  if (length(gm_rows) > 0) fwrite(rbindlist(gm_rows, fill = TRUE), OUT_G)
  flush(con)
}

cat("\ntime", format(Sys.time()), "\n")
sink(type = "message"); sink(); close(con)
