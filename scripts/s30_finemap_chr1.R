# -*- coding: utf-8 -*-
# s30_finemap_chr1.R -- Task D: L1 locus fine-mapping (CDC42 vs WNT4 attribution)
#
# Purpose
#   Task D of the 2026-09-24 order: resolve whether the L1 (chr1p36.12) signal that MR
#   attributes to CDC42 actually localises to CDC42 or to WNT4.
#
# Method  (fine-mapping ONLY -- no coloc.susie anywhere; standing prohibition respected)
#   * Variant set + harmonisation: identical to s28_coloc.R (allele-aware key on
#     pos38 + sorted allele pair; outcome keeps its own GRCh38 position; eQTL/LD in GRCh37).
#   * LD: signed r matrix from the SAME OneK1K 980-donor PLINK bfile used for clumping,
#     pruned with --indep-pairwise 200 50 0.9 (primary) / 0.99 (sensitivity).
#   * Fine-mapping: susieR::susie_rss on (a) the outcome and (b) the eQTL, per cell type.
#     * r matches z element order; n supplied; estimate_residual_variance = FALSE is the
#       package-correct choice here because R is NOT in-sample for the GWAS.
#   * Gene annotation from Ensembl GRCh37 REST (queried 2026-09-24).
#
# Outputs (tables/)
#   35_finemap_variant.csv    long: per variant per side per pair -> PIP + CS membership
#   36_finemap_pair.csv       per pair: CS summary, lead gene, PIP mass per gene body
#   37_finemap_outcome.csv    outcome-only fine-mapping on the full L1 window
#
# Usage: Rscript s30_finemap_chr1.R [prune_r2]     (default 0.9)

suppressPackageStartupMessages({
  library(data.table)
  library(susieR)
})

PROJ  <- "D:/endometriosis_project/11_sc_eqtl_mr_project"
PREP  <- file.path(PROJ, "00_data_raw/onek1k/coloc_prep")
TAB   <- file.path(PROJ, "tables")
WORK  <- file.path(PROJ, "00_data_raw/onek1k/finemap_work")
PLINK <- "D:/endometriosis_project/_tools/plink.exe"
BFILE <- "D:/endometriosis_project/_onek1k_plink/plink_merged_980_donors"
dir.create(WORK, showWarnings = FALSE, recursive = TRUE)

LOG <- "D:/endometriosis_project/_s30_finemap.log"
con <- file(LOG, open = "wt", encoding = "UTF-8")
sink(con, split = TRUE)
sink(con, type = "message")

args <- commandArgs(trailingOnly = TRUE)
PRUNE_R2 <- if (length(args) >= 1) as.numeric(args[1]) else 0.9
TAG <- if (PRUNE_R2 == 0.9) "" else sprintf("_r%s", gsub("\\.", "", as.character(PRUNE_R2)))

cat("=== s30 finemap chr1 @", format(Sys.time()), "===\n")
cat("prune r2 threshold =", PRUNE_R2, " tag =", TAG, "\n\n")

N_EQTL <- 980          # OneK1K donors (same convention as s28)
LOCUS  <- "L1"
CHROM  <- "1"
LO37   <- 22000000
HI37   <- 22900000

# --- gene bodies, Ensembl GRCh37 REST, queried 2026-09-24 -----------------------
GENES <- data.table(
  gene  = c("LINC00339", "CDC42", "WNT4"),
  ensid = c("ENSG00000218510", "ENSG00000070831", "ENSG00000162552"),
  start37 = c(22351681L, 22379120L, 22443798L),
  end37   = c(22357716L, 22419437L, 22470462L)
)
cat("gene bodies (GRCh37):\n"); print(GENES)

# Flanking labels are assigned by NEAREST gene and carry the distance in kb. A plain
# "first gene within 25 kb" rule mislabels the CDC42-WNT4 gap: e.g. GRCh37:22,436,446
# sits 17.0 kb past CDC42 but only 7.4 kb before WNT4, so it must read flank5_WNT4_7kb.
annot_gene <- function(p37) {
  out <- rep(NA_character_, length(p37))
  for (i in seq_len(nrow(GENES))) {
    out[p37 >= GENES$start37[i] & p37 <= GENES$end37[i]] <- GENES$gene[i]
  }
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

# --- LD builder: signed r matrix on the OneK1K reference ------------------------
# pval is used ONLY to cap the retained set at 1000 variants when pruning leaves more
# than that; the cap is by outcome p ascending (same convention as s28). Documented
# limitation: fine-mapping PIPs are therefore conditional on the retained set.
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

# --- harmonisation: identical to s28_coloc.R -----------------------------------
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
  m <- m[!duplicated(rsid)]        # rsid is the LD dimname -> must be unique
  m[, p37 := as.integer(sub("^[^:]+:", "", snp_grch37))]
  m[, N_out := n_cases + n_controls]
  setorder(m, pos38)
  m
}

# --- SuSiE on one side ---------------------------------------------------------
# warnings are captured (not just printed) because "IBSS algorithm did not converge"
# is the decisive diagnostic for LD-reference / summary-stat inconsistency and must
# be carried into the output table rather than lost in the log.
run_susie <- function(beta, se, R, n, L = 10) {
  z <- beta / se
  z[!is.finite(z)] <- 0
  wmsg <- character(0)
  r <- tryCatch(
    withCallingHandlers(
      susie_rss(z = z, R = R, n = n, L = L,
                estimate_residual_variance = FALSE, check_prior = TRUE),
      warning = function(w) {
        wmsg <<- c(wmsg, conditionMessage(w))
        invokeRestart("muffleWarning")
      }),
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

pip_of  <- function(r) if (is.null(r$.error)) susie_get_pip(r) else rep(NA_real_, 0)
cs_list <- function(r) if (is.null(r$.error)) r$sets$cs else NULL
# safe max: returns NA_real_ when everything is NA / empty
mx <- function(x) { x <- x[is.finite(x)]; if (length(x) == 0) NA_real_ else max(x) }
argmax1 <- function(x) { x[!is.finite(x)] <- -Inf; if (length(x) == 0) NA_integer_ else which.max(x) }

# --- load inputs ---------------------------------------------------------------
disc_all <- fread(file.path(TAB, "15_discovery_significant_with_locus.csv"))
disc <- disc_all[locus_id == LOCUS]
mr31 <- fread(file.path(TAB, "31_final_MR_with_method.csv"))
setnames(mr31, old = c("p", "p_outcome"), new = c("p_mr", "p_outcome_mr"), skip_absent = TRUE)
disc <- merge(disc, mr31[method == "Wald ratio", .(locus_id, gene, cell_type,
                                                  instrument_grch38 = instrument,
                                                  p_mr, F_stat_mr = F_stat)],
              by = c("locus_id", "gene", "cell_type"), all.x = TRUE)
cat("\nL1 pairs =", nrow(disc), "\n"); print(disc[, .(gene, cell_type, instrument_grch38, F_stat_mr)])

# authoritative GRCh37 instrument coordinates (forward-liftover record, per MEMORY B2)
instr37_map <- unique(fread(file.path(TAB, "01b_discovery_instruments_grch38.csv"))[
  , .(gene, cell_type, pos37, pos38, variant_id, variant_id_grch37)])
cat("\n01b instrument records for L1 genes:\n")
print(instr37_map[gene %in% c("CDC42", "LINC00339")][order(gene, cell_type)])

b  <- fread(file.path(PREP, sprintf("bim_universe_%s.csv", LOCUS)))
gw <- fread(file.path(PREP, sprintf("gwas_%s.csv", LOCUS)))

# explicit GRCh37 window trim: s26 read the GWAS on a loose GRCh38 window; the eQTL
# universe must be the same GRCh37 locus window used for coloc (22.0-22.9 Mb)
b <- b[pos >= LO37 & pos <= HI37]
cat("bim universe in window =", nrow(b), "\n")

# =============================================================================
# PART 1 -- outcome-only fine-mapping on the full L1 window
# =============================================================================
cat("\n========== PART 1: outcome-only fine-mapping ==========\n")
g0 <- copy(gw)
g0[, ea := toupper(ea)]; g0[, oa := toupper(oa)]
g0 <- g0[!is.na(pos38) & !is.na(beta) & !is.na(se) & se > 0]
bb <- data.table(snp37 = b$snp, pos37 = as.integer(b$pos),
                 A1 = toupper(b$a1), A2 = toupper(b$a2))
# map bim GRCh37 positions to GRCh38 via the forward-lifted eQTL table (same source as s26)
eq_map <- unique(fread(file.path(PREP, sprintf("eqtl_%s_CDC42.csv", LOCUS)))[
  , .(snp_grch37, snp_grch38)])
bb <- merge(bb, eq_map, by.x = "snp37", by.y = "snp_grch37", all.x = TRUE)
bb[, pos38_b := as.integer(sub("^[^:]+:", "", snp_grch38))]
bb <- bb[!is.na(pos38_b)]
bb[, allele_key := paste0(pos38_b, "_", pmin(A1, A2), "_", pmax(A1, A2))]
g0[, allele_key := paste0(pos38, "_", pmin(ea, oa), "_", pmax(ea, oa))]
g0 <- g0[!duplicated(allele_key)]
o1 <- merge(bb[, .(snp37, pos37, A1, A2, pos38 = pos38_b, allele_key)], g0,
            by = "allele_key", suffixes = c("_b", "_o"))
# merge() suffixes BOTH non-key pos38 columns -> keep the bim-derived one and rename it
# back to pos38 (same pitfall as s28's pos38_e/pos38_o; both are identical by construction)
if ("pos38_b" %in% names(o1)) setnames(o1, "pos38_b", "pos38")
if ("pos38_o" %in% names(o1)) o1[, pos38_o := NULL]
if (!"pos38" %in% names(o1)) stop("pos38 column lost after merge")
o1 <- o1[(ea == A1 & oa == A2) | (ea == A2 & oa == A1)]
o1[, flip := (ea == A2) & (oa == A1)]
# NOTE: no eQTL MAF available on this path, so palindromic A/T, C/G variants are
# dropped outright rather than resolved by allele frequency (conservative, documented).
o1[, palindromic := paste0(pmin(A1, A2), pmax(A1, A2)) %in% c("AT", "CG")]
n_pal <- nrow(o1[palindromic == TRUE])
o1 <- o1[palindromic == FALSE]
o1[, bo := ifelse(flip, -beta, beta)]
o1[, p37 := pos37]
o1 <- o1[!duplicated(rsid)]
setorder(o1, pos38)
cat("outcome-only intersection (bim x gwas, window-trimmed) =", nrow(o1),
    " (dropped palindromic =", n_pal, ")\n")

# keep the pruned set (plink) then cap by outcome p
LD1 <- ld_signed(o1$snp37, "outcome", o1$p)
out_rows <- list(); out_pair <- list()
if (is.null(LD1)) {
  cat("!! outcome LD unavailable\n")
} else {
  cat("outcome LD =", nrow(LD1), "x", ncol(LD1), "\n")
  ids <- rownames(LD1)
  oo <- o1[match(ids, snp37)]
  nn <- as.numeric(median(oo$n_cases + oo$n_controls))
  cat("outcome n (median) =", nn, "\n")
  # NOTE: dimnames must be set ON the matrix that stays named LD (copy-on-modify would
  # leave LD with the original chr:pos dimnames and break rsid lookups later).
  dimnames(LD1) <- list(oo$rsid, oo$rsid)
  R <- LD1
  r <- run_susie(oo$bo, oo$se, R, nn, L = 10)
  cat("outcome converge flag:", warn_flag(r), "\n")
  if (!is.null(r$.error)) {
    cat("!! outcome susie_rss error:", r$.error, "\n")
  } else {
    pip <- susie_get_pip(r)
    cs  <- r$sets$cs
    cat("outcome: niter =", r$niter, " nCS =", length(cs),
        " max PIP =", round(mx(pip), 4), "\n")
    oo[, PIP := pip]
    oo[, gene_annot := annot_gene(p37)]
    oo[, in_cs := ""]
    if (length(cs) > 0) {
      for (j in seq_along(cs)) {
        idx <- cs[[j]]
        nm <- names(cs)[j]
        oo[idx, in_cs := ifelse(in_cs == "", nm, paste(in_cs, nm, sep = ";"))]
        lead <- idx[argmax1(pip[idx])]
        cat(sprintf("  CS %s: size=%d  lead=%s (GRCh37:%d, %s)  PIP_lead=%.4f\n",
                    nm, length(idx), oo$rsid[lead], oo$p37[lead],
                    oo$gene_annot[lead], pip[lead]))
      }
    }
    oo[, side := "outcome"]
    setcolorder(oo, c("side", "rsid", "snp37", "p37", "pos38", "gene_annot", "PIP", "in_cs"))
    out_rows[[1]] <- oo[, .(side, rsid, snp37 = snp37, pos37 = p37, pos38, gene_annot,
                            PIP, in_cs, beta = bo, se)]
    gm <- oo[, .(PIP_mass = sum(PIP, na.rm = TRUE), n_var = .N), by = gene_of(gene_annot)]
    setnames(gm, "gene_of", "gene_region")
    cat("-- outcome PIP mass by gene region --\n"); print(gm[order(-PIP_mass)])
    out_pair[[1]] <- data.table(
      scope = "outcome_only", nsnp = nrow(oo), niter = r$niter, n_cs = length(cs),
      converge = warn_flag(r),
      max_pip = mx(pip), argmax_snp = oo$rsid[argmax1(pip)],
      argmax_pos37 = oo$p37[argmax1(pip)], argmax_gene = oo$gene_annot[argmax1(pip)],
      pip_mass_CDC42 = oo[gene_annot == "CDC42", sum(PIP)],
      pip_mass_WNT4  = oo[gene_annot == "WNT4", sum(PIP)],
      pip_mass_LINC00339 = oo[gene_annot == "LINC00339", sum(PIP)],
      pip_mass_flank = oo[!gene_annot %in% GENES$gene, sum(PIP)],
      n_var_CDC42 = oo[gene_annot == "CDC42", .N],
      n_var_WNT4  = oo[gene_annot == "WNT4", .N])
  }
}

# =============================================================================
# PART 2 -- per (gene x cell type) eQTL fine-mapping + outcome on same set
# =============================================================================
cat("\n========== PART 2: eQTL fine-mapping per pair ==========\n")
f_e_cache <- list()
var_rows <- if (length(out_rows)) out_rows else list()
pair_rows <- if (length(out_pair)) out_pair else list()

OUT_V <- file.path(TAB, sprintf("35_finemap_variant%s.csv", TAG))
OUT_P <- file.path(TAB, sprintf("36_finemap_pair%s.csv", TAG))

for (i in seq_len(nrow(disc))) {
  g  <- disc$gene[i]; ct <- disc$cell_type[i]
  tag <- sprintf("%s_%s", g, ct)
  cat(sprintf("\n[%02d/%02d] %s / %s\n", i, nrow(disc), g, ct))

  f_e <- file.path(PREP, sprintf("eqtl_%s_%s.csv", LOCUS, g))
  if (!file.exists(f_e)) { cat("  ! missing eqtl file\n"); next }
  if (is.null(f_e_cache[[g]])) f_e_cache[[g]] <- fread(f_e)
  e_all <- f_e_cache[[g]]

  m <- harmonise(e_all, b, gw, cell = ct)
  if (is.null(m)) { cat("  ! harmonise failed\n"); next }
  cat("  harmonised nsnp =", nrow(m), "\n")

  LD <- ld_signed(m$snp_grch37, sprintf("%s_%s", g, ct), m$pval_nominal)
  if (is.null(LD)) { cat("  ! LD unavailable\n"); next }
  ids <- rownames(LD)
  mm <- m[match(ids, snp_grch37)]
  cat("  LD =", nrow(LD), "x", ncol(LD), " ; N_out median =",
      round(median(mm$N_out)), "\n")

  # dimnames must be set on the object kept as LD (see PART 1 note)
  dimnames(LD) <- list(mm$rsid, mm$rsid)
  R <- LD
  N_out_med <- as.numeric(median(mm$N_out))

  # --- eQTL side (quant; n = 980 donors, same convention as s28) ---
  re <- run_susie(mm$slope, mm$slope_se, R, N_EQTL, L = 10)
  # --- outcome side on the SAME variant set ---
  ro <- run_susie(mm$b2, mm$se, R, N_out_med, L = 10)
  cat("  converge: eQTL =", warn_flag(re), " | outcome =", warn_flag(ro), "\n")

  rec <- function(r) {
    if (is.null(r$.error)) return(list(pip = susie_get_pip(r), cs = r$sets$cs, niter = r$niter, err = NA_character_))
    list(pip = rep(NA_real_, nrow(mm)), cs = NULL, niter = NA_integer_, err = r$.error)
  }
  E <- rec(re); O <- rec(ro)
  if (!is.na(E$err)) cat("  eQTL susie error:", E$err, "\n")
  if (!is.na(O$err)) cat("  outcome susie error:", O$err, "\n")

  mm[, PIP_eqtl := E$pip]
  mm[, PIP_outcome := O$pip]
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
  if (!is.na(lead_e)) cat(sprintf("    eQTL CS lead: %s GRCh37:%d  %s  PIP=%.4f   (r2 with MR instrument below)\n",
        mm$rsid[lead_e], mm$p37[lead_e], mm$gene_annot[lead_e], mm$PIP_eqtl[lead_e]))
  if (!is.na(lead_o)) cat(sprintf("    out  CS lead: %s GRCh37:%d  %s  PIP=%.4f\n",
        mm$rsid[lead_o], mm$p37[lead_o], mm$gene_annot[lead_o], mm$PIP_outcome[lead_o]))

  # --- MR instrument on this set + LD to CS leads ---
  # GRCh37 position comes from 01b (authoritative forward-liftover record), NOT from a
  # nearest-neighbour match against the retained set: nearest-neighbour silently returns
  # a variant up to tens of kb away and yields a meaningless r2 ~ 0.
  instr38 <- disc$instrument_grch38[i]
  instr37 <- NA_integer_; instr37_id <- NA_character_; instr_in_set <- FALSE
  r2_inst_leadE <- NA_real_; r2_inst_leadO <- NA_real_
  im <- instr37_map[gene == g & cell_type == ct]
  if (nrow(im) >= 1) {
    instr37 <- as.integer(im$pos37[1])
    instr37_id <- as.character(im$variant_id_grch37[1])
    instr_in_set <- instr37_id %in% mm$snp_grch37
    if (instr_in_set) {
      rid <- mm$rsid[match(instr37_id, mm$snp_grch37)]
      if (!is.na(lead_e)) r2_inst_leadE <- LD[rid, mm$rsid[lead_e]]^2
      if (!is.na(lead_o)) r2_inst_leadO <- LD[rid, mm$rsid[lead_o]]^2
    }
    cat(sprintf("    MR instrument %s -> %s (GRCh37:%d, %s)  in_finemap_set=%s  r2(inst,CS_eqTL)=%s  r2(inst,CS_out)=%s\n",
                instr38, instr37_id, instr37, annot_gene(instr37), instr_in_set,
                ifelse(is.na(r2_inst_leadE), "NA", sprintf("%.3f", r2_inst_leadE)),
                ifelse(is.na(r2_inst_leadO), "NA", sprintf("%.3f", r2_inst_leadO))))
  } else {
    cat("    note: no 01b record for this gene/cell_type\n")
  }

  # --- outcome-side palindromic sensitivity ---
  nopal <- which(!mm$palindromic)
  cs_nopal <- NULL; conv_nopal <- NA_character_; mpip_nopal <- NA_real_
  lead_nopal_pos37 <- NA_integer_; lead_nopal_gene <- NA_character_
  if (length(nopal) >= 20) {
    rn <- run_susie(mm$b2[nopal], mm$se[nopal], LD[nopal, nopal, drop = FALSE], N_out_med, L = 10)
    conv_nopal <- warn_flag(rn)
    if (is.null(rn$.error)) {
      pn <- susie_get_pip(rn); cs_nopal <- rn$sets$cs
      mpip_nopal <- mx(pn)
      if (length(cs_nopal) > 0) {
        i0 <- sort(unique(unlist(cs_nopal)))
        lead_nopal_pos37 <- mm$p37[nopal][i0][argmax1(pn[i0])]
        lead_nopal_gene  <- mm$gene_annot[nopal][i0][argmax1(pn[i0])]
      }
      cat(sprintf("    [sens] outcome without %d palindromic SNPs: nsnp=%d  converge=%s  nCS=%d  maxPIP=%s  CS_lead=%s\n",
                  sum(mm$palindromic), length(nopal), conv_nopal, length(cs_nopal),
                  ifelse(is.na(mpip_nopal), "NA", sprintf("%.4f", mpip_nopal)),
                  ifelse(is.na(lead_nopal_gene), "none",
                         sprintf("%d %s", lead_nopal_pos37, lead_nopal_gene))))
    } else {
      cat("    [sens] outcome without palindromic SNPs: error ", rn$.error, "\n")
    }
  }

  gm_e <- mm[, .(PIP_mass = sum(PIP_eqtl, na.rm = TRUE), n_var = .N), by = gene_of(gene_annot)]
  gm_o <- mm[, .(PIP_mass = sum(PIP_outcome, na.rm = TRUE), n_var = .N), by = gene_of(gene_annot)]
  setnames(gm_e, "gene_of", "gene_region"); setnames(gm_o, "gene_of", "gene_region")
  cat("  eQTL PIP mass:"); print(gm_e[order(-PIP_mass)])
  cat("  out  PIP mass:"); print(gm_o[order(-PIP_mass)])

  var_rows[[length(var_rows) + 1]] <- mm[, .(scope = sprintf("%s|%s", g, ct), side = "eqtl",
      rsid, snp37 = snp_grch37, pos37 = p37, pos38, gene_annot, PIP = PIP_eqtl, in_cs = in_cs_eqtl,
      beta = slope, se = slope_se)]
  var_rows[[length(var_rows) + 1]] <- mm[, .(scope = sprintf("%s|%s", g, ct), side = "outcome",
      rsid, snp37 = snp_grch37, pos37 = p37, pos38, gene_annot, PIP = PIP_outcome, in_cs = in_cs_outcome,
      beta = b2, se = se)]

  pair_rows[[length(pair_rows) + 1]] <- data.table(
    scope = tag, nsnp = nrow(mm), niter_eqtl = E$niter, niter_outcome = O$niter,
    converge_eqtl = warn_flag(re), converge_outcome = warn_flag(ro),
    n_cs_eqtl = length(E$cs), n_cs_outcome = length(O$cs),
    max_pip_eqtl = mx(mm$PIP_eqtl),
    max_pip_outcome = mx(mm$PIP_outcome),
    cs_lead_eqtl_rsid = if (is.na(lead_e)) NA_character_ else mm$rsid[lead_e],
    cs_lead_eqtl_pos37 = if (is.na(lead_e)) NA_integer_ else mm$p37[lead_e],
    cs_lead_eqtl_gene = if (is.na(lead_e)) NA_character_ else mm$gene_annot[lead_e],
    cs_lead_outcome_rsid = if (is.na(lead_o)) NA_character_ else mm$rsid[lead_o],
    cs_lead_outcome_pos37 = if (is.na(lead_o)) NA_integer_ else mm$p37[lead_o],
    cs_lead_outcome_gene = if (is.na(lead_o)) NA_character_ else mm$gene_annot[lead_o],
    mr_instrument_grch38 = instr38,
    mr_instrument_grch37 = instr37, mr_instrument_id_grch37 = instr37_id,
    mr_instrument_in_finemap_set = instr_in_set,
    mr_instrument_gene = if (is.na(instr37)) NA_character_ else annot_gene(instr37),
    r2_instr_vs_cs_lead_eqtl = r2_inst_leadE,
    r2_instr_vs_cs_lead_outcome = r2_inst_leadO,
    pip_mass_eqtl_CDC42 = gm_e[gene_region == "CDC42", sum(PIP_mass)],
    pip_mass_eqtl_WNT4  = gm_e[gene_region == "WNT4", sum(PIP_mass)],
    pip_mass_eqtl_LINC  = gm_e[gene_region == "LINC00339", sum(PIP_mass)],
    pip_mass_eqtl_flank = mm[!gene_annot %in% GENES$gene, sum(PIP_eqtl, na.rm = TRUE)],
    pip_mass_out_CDC42 = gm_o[gene_region == "CDC42", sum(PIP_mass)],
    pip_mass_out_WNT4  = gm_o[gene_region == "WNT4", sum(PIP_mass)],
    pip_mass_out_LINC  = gm_o[gene_region == "LINC00339", sum(PIP_mass)],
    pip_mass_out_flank = mm[!gene_annot %in% GENES$gene, sum(PIP_outcome, na.rm = TRUE)],
    # palindromic sensitivity for the outcome side
    n_palindromic = sum(mm$palindromic),
    converge_outcome_nopal = conv_nopal, n_cs_outcome_nopal = length(cs_nopal),
    max_pip_outcome_nopal = mpip_nopal,
    cs_lead_outcome_nopal_pos37 = lead_nopal_pos37,
    cs_lead_outcome_nopal_gene = lead_nopal_gene)

  # incremental write
  if (length(var_rows) > 0) fwrite(rbindlist(var_rows, fill = TRUE), OUT_V)
  if (length(pair_rows) > 0) fwrite(rbindlist(pair_rows, fill = TRUE), OUT_P)
  flush(con)
}

if (length(var_rows) > 0) {
  V <- rbindlist(var_rows, fill = TRUE)
  fwrite(V, OUT_V)
  cat("\nwrote", OUT_V, "rows =", nrow(V), "\n")
}
if (length(pair_rows) > 0) {
  P <- rbindlist(pair_rows, fill = TRUE)
  fwrite(P, OUT_P)
  cat("wrote", OUT_P, "rows =", nrow(P), "\n")
}

cat("\ntime", format(Sys.time()), "\n")
sink(type = "message"); sink(); close(con)
