# -*- coding: utf-8 -*-
# s51c_coloc_sens_L2L3.R -- 任务A：L2 / L3 的 coloc.abf 参数敏感性矩阵（3 窗口 x 3 p12）
#
# 严格复用 s34_coloc_sensitivity.R 的 run_abf()（逐字符相同口径），仅改动：
#   1) LOCI 列表 -> L2 (YME1L1/CD4_NC) + L3-recheck (ANXA4/Mono_NC)
#   2) 输入目录按位点分别指定：L2 = coloc_sens_prep，L3 = coloc_sens_prep_recheck
#   3) 输出到**新文件**，不覆盖已交付的 tables/40_coloc_sensitivity_matrix.csv
#
# 用法: Rscript s51c_coloc_sens_L2L3.R

suppressPackageStartupMessages({
  library(data.table)
  library(coloc)
})

PROJ <- "D:/endometriosis_project/11_sc_eqtl_mr_project"
DIR_L2 <- file.path(PROJ, "00_data_raw/onek1k/coloc_sens_prep")
DIR_L3 <- file.path(PROJ, "00_data_raw/onek1k/coloc_sens_prep_recheck")
OUT_L2 <- file.path(PROJ, "tables/51_coloc_sensitivity_L2.csv")
OUT_L3 <- file.path(PROJ, "tables/51b_coloc_sensitivity_L3_recheck.csv")

LOG <- "D:/endometriosis_project/_s51c_coloc_sens.log"
con <- file(LOG, open = "wt", encoding = "UTF-8")
sink(con, split = TRUE)

N_EQTL <- 980

LOCI <- list(
  YME1L1_CD4_NC = list(nm="YME1L1_CD4_NC", locus="L2", gene="YME1L1", ct="CD4_NC",
                       chrom="10", lead38=27154694L, dir=DIR_L2),
  ANXA4_Mono_NC = list(nm="ANXA4_Mono_NC", locus="L3", gene="ANXA4",  ct="Mono_NC",
                       chrom="2",  lead38=69729589L, dir=DIR_L3)
)
WINDOWS <- c("500Kb", "1Mb", "2Mb")
P12S <- c("1e-4"=1e-4, "1e-5"=1e-5, "1e-6"=1e-6)

cat("=== s51c L2/L3 coloc sensitivity ===\n")
cat("time", format(Sys.time()), "\n\n")

# ---- run_abf: 逐字符复刻 s34_coloc_sensitivity.R ----
run_abf <- function(e, g, p12) {
  e <- copy(e); g <- copy(g)
  e[, a1 := toupper(a1)]; e[, a2 := toupper(a2)]
  e[, allele_key := paste0(pos38, "_", pmin(a1, a2), "_", pmax(a1, a2))]
  g[, ea := toupper(ea)]; g[, oa := toupper(oa)]
  g[, allele_key := paste0(pos38, "_", pmin(ea, oa), "_", pmax(ea, oa))]
  g <- g[!duplicated(allele_key)]
  m <- merge(e, g, by = "allele_key", all = FALSE)
  if (nrow(m) < 10) return(NULL)
  if ("pos38.x" %in% names(m)) { m[, pos38 := pos38.x]; m[, c("pos38.x","pos38.y") := NULL] }
  m[, flip := (ea == a2) & (oa == a1)]
  m <- m[(ea == a1 & oa == a2) | (ea == a2 & oa == a1)]
  m[, palindromic := paste0(pmin(a1, a2), pmax(a1, a2)) %in% c("AT", "CG")]
  m[, af_mism := abs(af - eaf)]
  m[, af_flip := abs(af - (1 - eaf))]
  m <- m[!(palindromic == TRUE & pmin(af_mism, af_flip) > 0.05 & abs(af_mism - af_flip) < 0.05)]
  m[palindromic == TRUE & af_flip < af_mism, flip := !flip]
  m[, b2 := ifelse(flip, -beta, beta)]
  m[, maf_e := pmin(af, 1 - af)]
  m[, maf_o := pmin(eaf, 1 - eaf)]
  m <- m[is.finite(maf_e) & is.finite(maf_o) & maf_e > 0 & maf_o > 0 &
           is.finite(slope) & is.finite(slope_se) & slope_se > 0]
  m <- m[!is.na(rsid) & rsid != ""]
  setorder(m, pos38)
  if (nrow(m) < 10) return(NULL)
  m[, N_out := n_cases + n_controls]
  s_scalar <- sum(m$n_cases) / sum(m$N_out)
  d1 <- list(beta = m$slope, varbeta = m$slope_se^2, snp = m$rsid, position = m$pos38,
             MAF = m$maf_e, N = N_EQTL, type = "quant")
  d2 <- list(beta = m$b2, varbeta = m$se^2, snp = m$rsid, position = m$pos38,
             MAF = m$maf_o, N = m$N_out, type = "cc", s = s_scalar)
  ab <- tryCatch(coloc.abf(d1, d2, p12 = p12), error = function(z) NULL)
  if (is.null(ab)) return(NULL)
  s <- ab$summary
  list(nsnp = as.integer(s["nsnps"]), H0 = s["PP.H0.abf"], H1 = s["PP.H1.abf"],
       H2 = s["PP.H2.abf"], H3 = s["PP.H3.abf"], H4 = s["PP.H4.abf"])
}

res <- list()
for (nm in names(LOCI)) {
  L <- LOCI[[nm]]
  cat(sprintf("\n--- %s (%s) gene=%s ct=%s dir=%s\n", nm, L$locus, L$gene, L$ct, L$dir))
  for (w in WINDOWS) {
    tag <- sprintf("%s_%s_%s", L$gene, L$ct, w)
    f_e <- file.path(L$dir, sprintf("eqtl_%s.csv", tag))
    f_g <- file.path(L$dir, sprintf("gwas_%s.csv", tag))
    if (!file.exists(f_e) || !file.exists(f_g)) {
      cat(sprintf("  ! missing %s (e=%s g=%s)\n", tag, file.exists(f_e), file.exists(f_g)))
      next
    }
    e <- fread(f_e); g <- fread(f_g)
    cat(sprintf("  %s : eQTL=%d  GWAS=%d\n", tag, nrow(e), nrow(g)))
    for (pname in names(P12S)) {
      p12 <- P12S[[pname]]
      r <- run_abf(e, g, p12)
      if (is.null(r)) {
        res[[length(res)+1]] <- data.table(locus=L$locus, gene=L$gene, cell_type=L$ct,
          window=w, p12=pname, nsnp=NA_integer_, PP.H0=NA_real_, PP.H1=NA_real_,
          PP.H2=NA_real_, PP.H3=NA_real_, PP.H4=NA_real_, robust=NA_character_)
        cat(sprintf("    %s p12=%s -> FAILED\n", w, pname))
      } else {
        robust <- ifelse(r$H4 >= 0.8, "稳健", "敏感性不稳健")
        res[[length(res)+1]] <- data.table(locus=L$locus, gene=L$gene, cell_type=L$ct,
          window=w, p12=pname, nsnp=r$nsnp, PP.H0=r$H0, PP.H1=r$H1,
          PP.H2=r$H2, PP.H3=r$H3, PP.H4=r$H4, robust=robust)
        cat(sprintf("    %s p12=%s nsnp=%d PP.H4=%.4f [%s]\n", w, pname, r$nsnp, r$H4, robust))
      }
    }
  }
  flush(con)
}

R <- rbindlist(res)
R[, locus_row := sprintf("%s_%s", gene, cell_type)]
R2 <- R[locus == "L2"]
R3 <- R[locus == "L3"]
fwrite(R2, OUT_L2)
fwrite(R3, OUT_L3)
cat("\n=== 写出 ===\n", OUT_L2, " rows=", nrow(R2), "\n", OUT_L3, " rows=", nrow(R3), "\n", sep="")

cat("\n--- 每位点×窗口 的 PP.H4 min/max（跨 p12）---\n")
print(R[, .(minH4=min(PP.H4, na.rm=TRUE), maxH4=max(PP.H4, na.rm=TRUE),
            any_below08=any(PP.H4 < 0.8, na.rm=TRUE)), by=.(locus, window)])

cat("\ntime", format(Sys.time()), "\n")
sink(); close(con)
