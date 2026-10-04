# -*- coding: utf-8 -*-
# s28_coloc.R -- colocalization via coloc.abf + coloc.susie
# Rules:
#   * eQTL/LD in GRCh37 (eQTL side is the reference); outcome keeps its own GRCh38 position
#   * allele-aware merge key (pos38 + sorted allele pair) to kill multi-allelic duplicates
#   * outcome n_cases/n_controls VARY per variant -> N and s must be per-variant vectors
#   * CDC42 and LINC00339 are run and reported SEPARATELY (red line 17-1)
#   * all output written to files; run R through PowerShell; ASCII-only log text
# Usage: Rscript s28_coloc.R [row_index ...]    (no args = all rows)
suppressPackageStartupMessages({
  library(data.table)
  library(coloc)
})

PROJ  <- "D:/endometriosis_project/11_sc_eqtl_mr_project"
PREP  <- file.path(PROJ, "00_data_raw/onek1k/coloc_prep")
TAB   <- file.path(PROJ, "tables")
WORK  <- file.path(PROJ, "00_data_raw/onek1k/coloc_work")
PLINK <- "D:/endometriosis_project/_tools/plink.exe"
BFILE <- "D:/endometriosis_project/_onek1k_plink/plink_merged_980_donors"
dir.create(WORK, showWarnings = FALSE, recursive = TRUE)

LOG <- "D:/endometriosis_project/_s28_coloc.log"
con <- file(LOG, open = "wt", encoding = "UTF-8")
sink(con, split = TRUE)
sink(con, type = "message")

cat("=== s28 coloc ===\ntime ", format(Sys.time()), "\n\n")
N_EQTL <- 980

LOCI <- list(
  L1 = list(chrom = "1",  lo37 = 22000000, hi37 = 22900000),
  L2 = list(chrom = "10", lo37 = 26900000, hi37 = 27900000),
  L3 = list(chrom = "2",  lo37 = 69300000, hi37 = 70200000)
)

disc <- fread(file.path(TAB, "15_discovery_significant_with_locus.csv"))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  disc <- disc[as.integer(args)]
  cat("subset rows: ", paste(args, collapse = ","), "\n")
}
cat("n_pairs = ", nrow(disc), "\n")

# ---- incremental output + resume (R may segfault at random points: rlang/cli heap issue) ----
OUT <- file.path(TAB, "33_coloc_results.csv")
RESUME <- length(args) == 0
res <- list()
done_key <- character(0)
if (file.exists(OUT)) {
  prev <- tryCatch(fread(OUT), error = function(z) NULL)
  if (!is.null(prev) && nrow(prev) > 0) {
    res[[1]] <- prev
    done_key <- unique(paste(prev$locus_id, prev$gene, prev$cell_type))
    cat("resume: ", length(done_key), " finished pair(s) already in 33_coloc_results.csv\n")
  }
}
write_out <- function() {
  R <- rbindlist(res, fill = TRUE)
  if (nrow(R) > 0) {
    R[verdict == "", verdict := fifelse(
      is.na(PP.H4), "undetermined",
      fifelse(PP.H4 > 0.8, "strong_shared",
              fifelse(PP.H4 >= 0.5, "moderate_shared", "no_shared")))]
  }
  fwrite(R, OUT)
  invisible(R)
}

ld_of <- function(locus, snps37, tag, pval) {
  f_in  <- file.path(WORK, sprintf("prune_%s.in", tag))
  f_out <- file.path(WORK, sprintf("prune_%s", tag))
  fwrite(data.table(V1 = snps37), f_in, col.names = FALSE, sep = "\t")
  cmd1 <- sprintf('"%s" --bfile "%s" --allow-no-sex --chr %s --extract "%s" --indep-pairwise 200 50 0.9 --out "%s"',
                  PLINK, BFILE, LOCI[[locus]]$chrom, f_in, f_out)
  system(cmd1, ignore.stdout = TRUE, ignore.stderr = TRUE)
  keep <- file.path(WORK, sprintf("prune_%s.prune.in", tag))
  if (!file.exists(keep)) return(NULL)
  k <- fread(keep, header = FALSE)$V1
  if (length(k) > 1000) {
    dt <- data.table(snp = snps37, p = pval)
    k <- dt[snp %in% k][order(p)][1:1000]$snp
  }
  if (length(k) < 2) return(NULL)
  f_k <- file.path(WORK, sprintf("keep_%s.txt", tag))
  fwrite(data.table(V1 = k), f_k, col.names = FALSE, sep = "\t")
  ld_o <- file.path(WORK, sprintf("ld_%s", tag))
  cmd2 <- sprintf('"%s" --bfile "%s" --allow-no-sex --chr %s --extract "%s" --write-snplist --r square --out "%s"',
                  PLINK, BFILE, LOCI[[locus]]$chrom, f_k, ld_o)
  system(cmd2, ignore.stdout = TRUE, ignore.stderr = TRUE)
  f_ld <- paste0(ld_o, ".ld"); f_sp <- paste0(ld_o, ".snplist")
  if (!file.exists(f_ld) || !file.exists(f_sp)) return(NULL)
  m <- as.matrix(fread(f_ld, header = FALSE))
  ids <- fread(f_sp, header = FALSE)$V1
  if (nrow(m) != length(ids) || nrow(m) < 2) return(NULL)
  dimnames(m) <- list(ids, ids)
  ok <- apply(m, 1, function(z) all(is.finite(z)))
  m <- m[ok, ok, drop = FALSE]
  m
}

for (i in seq_len(nrow(disc))) {
  g  <- disc$gene[i]; ct <- disc$cell_type[i]; lo <- disc$locus_id[i]
  tag <- sprintf("%s_%s_%s", lo, g, ct)
  if (RESUME && paste(lo, g, ct) %in% done_key) { cat(sprintf("\n[%02d/%02d] %s / %s / %s -- skip (done)\n", i, nrow(disc), lo, g, ct)); next }
  cat(sprintf("\n[%02d/%02d] %s / %s / %s\n", i, nrow(disc), lo, g, ct))

  f_e <- file.path(PREP, sprintf("eqtl_%s_%s.csv", lo, g))
  f_g <- file.path(PREP, sprintf("gwas_%s.csv", lo))
  f_b <- file.path(PREP, sprintf("bim_universe_%s.csv", lo))
  if (!file.exists(f_e) || !file.exists(f_g) || !file.exists(f_b)) { cat("  ! missing input, skip\n"); next }

  e  <- fread(f_e); n_e_all <- uniqueN(e$snp_grch37); e <- e[cell_type == ct]
  if (nrow(e) == 0) { cat("  ! no eQTL rows for this cell type, skip\n"); next }
  b  <- fread(f_b); gw <- fread(f_g)

  e <- merge(e, b, by.x = "snp_grch37", by.y = "snp", all.x = TRUE)
  e <- e[!is.na(a1) & !is.na(a2)]
  e[, pos38 := as.integer(sub("^[^:]+:", "", snp_grch38))]
  e <- e[!is.na(pos38)]
  e[, a1 := toupper(a1)]; e[, a2 := toupper(a2)]
  e[, allele_key := paste0(pos38, "_", pmin(a1, a2), "_", pmax(a1, a2))]

  gw[, ea := toupper(ea)]; gw[, oa := toupper(oa)]
  gw <- gw[!is.na(pos38) & !is.na(beta) & !is.na(se) & se > 0]
  gw[, allele_key := paste0(pos38, "_", pmin(ea, oa), "_", pmax(ea, oa))]
  gw <- gw[!duplicated(allele_key)]
  n_gw_win <- nrow(gw)

  m <- merge(e, gw, by = "allele_key", suffixes = c("_e", "_o"))
  if ("pos38_e" %in% names(m)) { setnames(m, "pos38_e", "pos38"); m[, pos38_o := NULL] }
  cat(sprintf("  nsnp eqtl_cis=%d  eqtl_ct=%d  gwas_win=%d  allele-aware intersect=%d\n",
              n_e_all, nrow(e), n_gw_win, nrow(m)))
  if (nrow(m) < 10) { cat("  ! intersect too small, skip\n"); next }

  m[, flip := (ea == a2) & (oa == a1)]
  m <- m[(ea == a1 & oa == a2) | (ea == a2 & oa == a1)]
  m[, palindromic := paste0(pmin(a1, a2), pmax(a1, a2)) %in% c("AT", "CG")]
  m[, af_mism := abs(af - eaf)]
  m[, af_flip := abs(af - (1 - eaf))]
  n_pal_bad <- nrow(m[palindromic == TRUE & pmin(af_mism, af_flip) > 0.05 &
                        abs(af_mism - af_flip) < 0.05])
  m <- m[!(palindromic == TRUE & pmin(af_mism, af_flip) > 0.05 &
             abs(af_mism - af_flip) < 0.05)]
  m[palindromic == TRUE & af_flip < af_mism, flip := !flip]
  m[, b2 := ifelse(flip, -beta, beta)]
  cat(sprintf("  harmonised=%d (flip=%d; palindromic=%d, dropped_af_ambiguous=%d)\n",
              nrow(m), sum(m$flip), sum(m$palindromic), n_pal_bad))
  if (nrow(m) < 10) { cat("  ! harmonised too small, skip\n"); next }

  m[, maf_e := pmin(af, 1 - af)]
  m[, maf_o := pmin(eaf, 1 - eaf)]
  m <- m[is.finite(maf_e) & is.finite(maf_o) & maf_e > 0 & maf_o > 0 &
           is.finite(slope) & is.finite(slope_se) & slope_se > 0]
  m <- m[!is.na(rsid) & rsid != ""]
  setorder(m, pos38)
  m[, N_out := n_cases + n_controls]
  m[, s_out := n_cases / (n_cases + n_controls)]
  cat(sprintf("  final nsnp=%d  rsid_ok=%d  N_out=%d..%d  s_out=%.3f..%.3f\n",
              nrow(m), sum(m$rsid != ""), min(m$N_out), max(m$N_out), min(m$s_out), max(m$s_out)))

  # coloc 5.2.3: coloc.abf requires a SCALAR s for type="cc" (vector s -> coercion error);
  # N stays per-variant (the harmonised file varies n_cases/n_controls by variant).
  # s is therefore taken as the window-aggregate case fraction; this is a documented
  # approximation, NOT a per-variant s (the per-variant values are kept in s_out for audit).
  s_scalar <- sum(m$n_cases) / sum(m$N_out)
  cat(sprintf("  s_scalar=%.4f (window-aggregate case fraction; per-variant s_out %.3f..%.3f)\n",
              s_scalar, min(m$s_out), max(m$s_out)))

  d1 <- list(beta = m$slope, varbeta = m$slope_se^2, snp = m$rsid, position = m$pos38,
             MAF = m$maf_e, N = N_EQTL, type = "quant")
  d2 <- list(beta = m$b2, varbeta = m$se^2, snp = m$rsid, position = m$pos38,
             MAF = m$maf_o, N = m$N_out, type = "cc", s = s_scalar)

  base <- function(method, nsnp_used, nsnp_ld, verdict, note) data.table(
    locus_id = lo, gene = g, cell_type = ct, method = method,
    nsnp_eqtl_cis = n_e_all, nsnp_eqtl_ct = nrow(e), nsnp_gwas_win = n_gw_win,
    nsnp_intersect = nrow(m), nsnp_used = nsnp_used, nsnp_ld = nsnp_ld,
    PP.H0 = NA_real_, PP.H1 = NA_real_, PP.H2 = NA_real_, PP.H3 = NA_real_, PP.H4 = NA_real_,
    n_snps_strong = NA_integer_, verdict = verdict, note = note)

  ab <- tryCatch(coloc.abf(d1, d2), error = function(z) { cat("   abf err: ", conditionMessage(z), "\n"); NULL })
  if (!is.null(ab)) {
    s <- ab$summary
    rr <- base("coloc.abf", as.integer(s["nsnps"]), NA_integer_, "", "ok")
    rr[, `:=`(PP.H0 = s["PP.H0.abf"], PP.H1 = s["PP.H1.abf"], PP.H2 = s["PP.H2.abf"],
              PP.H3 = s["PP.H3.abf"], PP.H4 = s["PP.H4.abf"])]
    res[[length(res) + 1]] <- rr
    cat(sprintf("  abf: nsnps=%d  H0=%.4f H1=%.4f H2=%.4f H3=%.4f H4=%.4f\n",
                as.integer(s["nsnps"]), s["PP.H0.abf"], s["PP.H1.abf"],
                s["PP.H2.abf"], s["PP.H3.abf"], s["PP.H4.abf"]))
  } else {
    res[[length(res) + 1]] <- base("coloc.abf", 0L, NA_integer_, "abf_failed", "coloc.abf error")
  }

  ld <- ld_of(lo, m$snp_grch37, tag, m$pval_nominal)
  if (is.null(ld)) {
    res[[length(res) + 1]] <- base("coloc.susie", 0L, 0L, "susie_skipped", "LD unavailable (PLINK)")
    cat("  ! LD unavailable -> fallback to abf only\n")
  } else {
    cat(sprintf("  LD ok: %d x %d (pruned r2<0.9)\n", nrow(ld), ncol(ld)))
    mp <- setNames(m$rsid, m$snp_grch37)
    ids <- rownames(ld)
    keep <- ids[ids %in% names(mp)]
    if (length(keep) < 2) {
      res[[length(res) + 1]] <- base("coloc.susie", 0L, nrow(ld), "susie_skipped", "LD x intersect too small")
      cat("  ! LD/intersect overlap too small\n")
    } else {
      mm <- m[match(keep, snp_grch37)]
      rn <- unname(mp[keep])
      L <- ld[keep, keep, drop = FALSE]
      dimnames(L) <- list(rn, rn)
      # susie (unlike coloc.abf) needs a SCALAR N -> use median per-variant N of the LD subset
      N2s <- as.numeric(median(mm$N_out))
      d1s <- list(beta = mm$slope, varbeta = mm$slope_se^2, snp = rn, position = mm$pos38,
                  MAF = mm$maf_e, N = N_EQTL, type = "quant", LD = L)
      d2s <- list(beta = mm$b2, varbeta = mm$se^2, snp = rn, position = mm$pos38,
                  MAF = mm$maf_o, N = N2s, type = "cc", s = s_scalar, LD = L)
      cat(sprintf("  susie inputs: nsnp=%d  N2_scalar=%.0f\n", length(rn), N2s))
      # coloc 5.2.3: runsusie() takes a SINGLE dataset -> call it twice, then coloc.susie(S1,S2)
      cs <- tryCatch({
        S1 <- runsusie(d1s)
        S2 <- runsusie(d2s)
        if (is.null(S1) || is.null(S2)) {
          list(ok = FALSE, msg = "runsusie returned NULL (no credible set)")
        } else {
          list(ok = TRUE, cs = coloc.susie(S1, S2))
        }
      }, error = function(z) list(ok = FALSE, msg = conditionMessage(z)))
      if (cs$ok) {
        sm <- cs$cs$summary
        if (is.null(sm) || nrow(sm) == 0) {
          res[[length(res) + 1]] <- base("coloc.susie", nrow(mm), nrow(ld), "susie_nosignal", "no reportable signal pair")
          cat("  susie: no reportable signal pair\n")
        } else {
          for (j in seq_len(nrow(sm))) {
            rr <- base(sprintf("coloc.susie#%d", j), nrow(mm), nrow(ld), "", paste0("hit", j))
            nstrong <- NA_integer_
            # coloc.susie()$results is a data.table: snp, SNP.PP.H4.row1 .. SNP.PP.H4.rowK
            col_nm <- paste0("SNP.PP.H4.row", j)
            spp <- tryCatch(cs$cs$results[[col_nm]], error = function(z) NULL)
            if (!is.null(spp)) nstrong <- as.integer(sum(spp > 0.8, na.rm = TRUE))
            rr[, `:=`(PP.H0 = sm[["PP.H0.abf"]][j], PP.H1 = sm[["PP.H1.abf"]][j],
                      PP.H2 = sm[["PP.H2.abf"]][j], PP.H3 = sm[["PP.H3.abf"]][j],
                      PP.H4 = sm[["PP.H4.abf"]][j], n_snps_strong = nstrong)]
            res[[length(res) + 1]] <- rr
            cat(sprintf("  susie#%d: H0=%.4f H3=%.4f H4=%.4f nSNP.PP.H4>0.8=%s\n", j,
                        sm[["PP.H0.abf"]][j], sm[["PP.H3.abf"]][j], sm[["PP.H4.abf"]][j],
                        ifelse(is.na(nstrong), "NA", as.character(nstrong))))
          }
        }
      } else {
        res[[length(res) + 1]] <- base("coloc.susie", nrow(mm), nrow(ld), "susie_failed", substr(cs$msg, 1, 200))
        cat("  ! susie failed: ", substr(cs$msg, 1, 160), "\n")
      }
    }
  }
  write_out()   # incremental: survive a mid-run segfault
  flush(con)
}

R <- write_out()
cat("\nresult rows = ", nrow(R), "\n")
if (RESUME) {
  chk <- unique(paste(R$locus_id, R$gene, R$cell_type))
  cat("distinct pairs in output = ", length(chk), " (expected ", nrow(disc), ")\n", sep = "")
  miss <- setdiff(paste(disc$locus_id, disc$gene, disc$cell_type), chk)
  if (length(miss) > 0) cat("!! MISSING PAIRS: ", paste(miss, collapse = " | "), "\n")
}
if (nrow(R) > 0) { cat("\n--- verdict counts (per row, incl. every susie pair) ---\n"); print(R[, .N, by = .(method, verdict)][order(method, verdict)]) }
cat("\ntime ", format(Sys.time()), "\n")
sink(type = "message"); sink(); close(con)
