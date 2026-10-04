# -*- coding: utf-8 -*-
# s34_coloc_sensitivity.R -- coloc.abf 参数敏感性（p12 x 窗口）
# 仅 H4_verdict=strong 的 5 个位点，3 窗口(500Kb/1Mb/2Mb) x 3 p12(1e-4/1e-5/1e-6)
# 输出 45 组合的 PP.H4 矩阵；跌破 0.8 标记"敏感性不稳健"。
# 只跑 coloc.abf（主判据）；susie 单独在诊断环节（LINC00339/CD4_NC 等位异质性）处理。
suppressPackageStartupMessages({
  library(data.table)
  library(coloc)
})

PROJ <- "D:/endometriosis_project/11_sc_eqtl_mr_project"
SENS <- file.path(PROJ, "00_data_raw/onek1k/coloc_sens_prep")
OUT  <- file.path(PROJ, "tables/40_coloc_sensitivity_matrix.csv")

LOG <- "D:/endometriosis_project/_s34_coloc_sens.log"
con <- file(LOG, open = "wt", encoding = "UTF-8")
sink(con, split = TRUE)

N_EQTL <- 980

LOCI <- list(
  CDC42_B_IN   = list(gene="CDC42", ct="B_IN",   chrom="1", lead38=22132301),
  CDC42_B_MEM  = list(gene="CDC42", ct="B_MEM",  chrom="1", lead38=22135618),
  CDC42_Mono_C = list(gene="CDC42", ct="Mono_C", chrom="1", lead38=22096228),
  CDC42_Mono_NC= list(gene="CDC42", ct="Mono_NC",chrom="1", lead38=22096228),
  ANXA4_Mono_NC= list(gene="ANXA4", ct="Mono_NC",chrom="2", lead38=69729589)
)
WINDOWS <- c("500Kb", "1Mb", "2Mb")
P12S <- c("1e-4"=1e-4, "1e-5"=1e-5, "1e-6"=1e-6)

cat("=== s34 coloc sensitivity ===\n")
cat("time", format(Sys.time()), "\n\n")

res <- list()

run_abf <- function(e, g, p12) {
  # e: eQTL (snp_grch37, pos38, af, slope, slope_se, a1, a2)
  # g: GWAS (pos38, ea, oa, beta, se, eaf, rsid, n_cases, n_controls)
  e <- copy(e); g <- copy(g)
  e[, a1 := toupper(a1)]; e[, a2 := toupper(a2)]
  e[, allele_key := paste0(pos38, "_", pmin(a1, a2), "_", pmax(a1, a2))]
  g[, ea := toupper(ea)]; g[, oa := toupper(oa)]
  g[, allele_key := paste0(pos38, "_", pmin(ea, oa), "_", pmax(ea, oa))]
  g <- g[!duplicated(allele_key)]
  m <- merge(e, g, by = "allele_key", all = FALSE)
  if (nrow(m) < 10) return(NULL)
  # merge 后 pos38 可能变成 pos38.x/pos38.y，统一回 pos38（e 侧为准）
  if ("pos38.x" %in% names(m)) { m[, pos38 := pos38.x]; m[, c("pos38.x","pos38.y") := NULL] }
  # harmonise: eQTL a1 作为参考 effect allele
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

for (nm in names(LOCI)) {
  L <- LOCI[[nm]]
  for (w in WINDOWS) {
    tag <- sprintf("%s_%s_%s", L$gene, L$ct, w)
    f_e <- file.path(SENS, sprintf("eqtl_%s.csv", tag))
    f_g <- file.path(SENS, sprintf("gwas_%s.csv", tag))
    if (!file.exists(f_e) || !file.exists(f_g)) { cat("  ! missing", tag, "\n"); next }
    e <- fread(f_e); g <- fread(f_g)
    for (pname in names(P12S)) {
      p12 <- P12S[[pname]]
      r <- run_abf(e, g, p12)
      if (is.null(r)) {
        res[[length(res)+1]] <- data.table(locus=nm, gene=L$gene, cell_type=L$ct,
          window=w, p12=pname, nsnp=NA_integer_, PP.H0=NA_real_, PP.H1=NA_real_,
          PP.H2=NA_real_, PP.H3=NA_real_, PP.H4=NA_real_, robust=NA_character_)
        cat(sprintf("  %s %s p12=%s -> FAILED\n", nm, w, pname))
      } else {
        robust <- ifelse(r$H4 >= 0.8, "稳健", "敏感性不稳健")
        res[[length(res)+1]] <- data.table(locus=nm, gene=L$gene, cell_type=L$ct,
          window=w, p12=pname, nsnp=r$nsnp, PP.H0=r$H0, PP.H1=r$H1,
          PP.H2=r$H2, PP.H3=r$H3, PP.H4=r$H4, robust=robust)
        cat(sprintf("  %s %s p12=%s nsnp=%d PP.H4=%.4f [%s]\n",
                    nm, w, pname, r$nsnp, r$H4, robust))
      }
    }
  }
  flush(con)
}

R <- rbindlist(res)
fwrite(R, OUT)
cat("\n=== 45 组合矩阵完成，写入", OUT, "===\n")

# 汇总：每个 位点×窗口 在 3 种 p12 下的 PP.H4 范围
cat("\n--- 每 位点×窗口 的 PP.H4 min/max（跨 p12）---\n")
R2 <- R[, .(minH4=min(PP.H4, na.rm=TRUE), maxH4=max(PP.H4, na.rm=TRUE),
            any_below08=any(PP.H4 < 0.8, na.rm=TRUE)), by=.(locus, window)]
print(R2)

cat("\ntime", format(Sys.time()), "\n")
sink(); close(con)
