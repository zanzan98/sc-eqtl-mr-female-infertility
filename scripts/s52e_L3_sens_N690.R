# -*- coding: utf-8 -*-
# s52e_L3_sens_N690.R -- L3 九组参数矩阵：N=690 与 N=980 并排（带诊断）
suppressPackageStartupMessages({ library(data.table); library(coloc) })
PROJ <- "D:/endometriosis_project/11_sc_eqtl_mr_project"
DIR3 <- file.path(PROJ, "00_data_raw/onek1k/coloc_sens_prep_recheck")
TAB  <- file.path(PROJ, "tables")
LOG  <- "D:/endometriosis_project/_s52e_L3_sens_N690.log"
con <- file(LOG, open = "wt", encoding = "UTF-8"); sink(con, split = TRUE)

cat("=== s52e L3 九组矩阵 N=690 vs N=980 ===\n"); cat("time", format(Sys.time()), "\n\n")

old <- fread(file.path(TAB, "51b_coloc_sensitivity_L3_recheck.csv"))
cat("old rows=", nrow(old), " cols=", paste(names(old), collapse=","), "\n", sep="")

run_abf_sens <- function(e, g, p12, Nq) {
  e <- copy(e); g <- copy(g)
  e[, a1 := toupper(a1)][, a2 := toupper(a2)]
  e[, allele_key := paste0(pos38, "_", pmin(a1, a2), "_", pmax(a1, a2))]
  g[, ea := toupper(ea)][, oa := toupper(oa)]
  g[, allele_key := paste0(pos38, "_", pmin(ea, oa), "_", pmax(ea, oa))]
  g <- g[!duplicated(allele_key)]
  m <- merge(e, g, by = "allele_key", all = FALSE)
  if (nrow(m) < 10) return(NULL)
  if ("pos38.x" %in% names(m)) { m[, pos38 := pos38.x]; m[, c("pos38.x","pos38.y") := NULL] }
  m[, flip := (ea == a2) & (oa == a1)]
  m <- m[(ea == a1 & oa == a2) | (ea == a2 & oa == a1)]
  m[, palindromic := paste0(pmin(a1,a2), pmax(a1,a2)) %in% c("AT","CG")]
  m[, af_mism := abs(af - eaf)][, af_flip := abs(af - (1 - eaf))]
  m <- m[!(palindromic == TRUE & pmin(af_mism, af_flip) > 0.05 & abs(af_mism - af_flip) < 0.05)]
  m[palindromic == TRUE & af_flip < af_mism, flip := !flip]
  m[, b2 := ifelse(flip, -beta, beta)]
  m[, maf_e := pmin(af, 1-af)][, maf_o := pmin(eaf, 1-eaf)]
  m <- m[is.finite(maf_e) & is.finite(maf_o) & maf_e > 0 & maf_o > 0 &
           is.finite(slope) & is.finite(slope_se) & slope_se > 0]
  m <- m[!is.na(rsid) & rsid != ""]
  setorder(m, pos38)
  if (nrow(m) < 10) return(NULL)
  No <- m$n_cases + m$n_controls
  d1 <- list(beta=m$slope, varbeta=m$slope_se^2, snp=m$rsid, position=m$pos38,
             MAF=m$maf_e, N=Nq, type="quant")
  d2 <- list(beta=m$b2, varbeta=m$se^2, snp=m$rsid, position=m$pos38,
             MAF=m$maf_o, N=No, type="cc", s=sum(m$n_cases)/sum(No))
  ab2 <- tryCatch(suppressWarnings(coloc.abf(d1, d2, p12 = p12)), error=function(z) NULL)
  if (is.null(ab2)) return(NULL)
  s <- ab2$summary
  list(nsnp = as.integer(s["nsnps"]), H0 = s[["PP.H0.abf"]], H1 = s[["PP.H1.abf"]],
       H2 = s[["PP.H2.abf"]], H3 = s[["PP.H3.abf"]], H4 = s[["PP.H4.abf"]])
}

WINDOWS <- c("500Kb", "1Mb", "2Mb"); P12S <- c("1e-4"=1e-4, "1e-5"=1e-5, "1e-6"=1e-6)
rows3 <- list()
for (w in WINDOWS) {
  tag <- sprintf("ANXA4_Mono_NC_%s", w)
  e <- fread(file.path(DIR3, sprintf("eqtl_%s.csv", tag)))
  g <- fread(file.path(DIR3, sprintf("gwas_%s.csv", tag)))
  cat(sprintf("\n[%s] eQTL=%d GWAS=%d\n", w, nrow(e), nrow(g)))
  for (pn in names(P12S)) {
    r <- run_abf_sens(e, g, P12S[[pn]], 690L)
    if (is.null(r)) { cat(sprintf("  %s p12=%s FAILED\n", w, pn)); next }
    sel <- old[window == w & abs(p12 - P12S[[pn]]) < 1e-18]
    o4 <- if (nrow(sel) == 1) sel$PP.H4 else NA_real_
    o3 <- if (nrow(sel) == 1) sel$PP.H3 else NA_real_
    cat(sprintf("  %-6s p12=%-5s nsnp=%-5d H4(690)=%.5f  H3(690)=%.5f | H4(980)=%.5f H3(980)=%.5f  dH4=%+.5f\n",
                w, pn, r$nsnp, r$H4, r$H3, o4, o3, r$H4 - o4))
    rows3[[length(rows3) + 1]] <- data.table(
      locus="L3", gene="ANXA4", cell_type="Mono_NC", window=w, p12=pn, nsnp_N690=r$nsnp,
      PP.H0_N690=r$H0, PP.H1_N690=r$H1, PP.H2_N690=r$H2, PP.H3_N690=r$H3, PP.H4_N690=r$H4,
      PP.H4_delivered_N980=o4, delta_H4 = r$H4 - o4,
      robust_N690 = ifelse(r$H4 >= 0.8, "稳健", "敏感性不稳健"),
      robust_N980 = ifelse(o4 >= 0.8, "稳健", "敏感性不稳健"))
  }
}
cat("\nrows3 收集数 =", length(rows3), "\n")
if (length(rows3)) {
  R3 <- rbindlist(rows3)
  fwrite(R3, file.path(TAB, "55b_coloc_sensitivity_L3_N690.csv"))
  R3[, flip := robust_N690 != robust_N980]
  cat("\n--- 分层翻转统计 ---\n"); print(R3[, .(n=.N, n_flip=sum(flip)), by=window])
  cat("\n--- dH4 汇总 ---\n"); print(R3[, .(max_abs_dH4=max(abs(delta_H4)))])
}
cat("\ntime", format(Sys.time()), "\n"); sink(); close(con)
