# -*- coding: utf-8 -*-
# s52d_N_sensitivity.R -- N 口径灵敏度：coloc.abf 中 eQTL 的 N 用 980 vs 真实细胞类型 N
#
# 背景（2026-09-28 实测）：OneK1K 各细胞类型 eQTL 有效样本量不同（Mono_NC=690, CD4_NC=980…），
# 而全部既有脚本硬编码 N_EQTL <- 980。coloc.abf(type="quant") 的 N 经 sdY.est()
# （lm(2*n*maf*(1-maf) ~ 1/varbeta)，斜率 = sdY^2）进入 ABF，故 N 用错会使 sdY 偏移 sqrt(980/N)。
#
# 本脚本回答：把 N 改成真实细胞类型 N 后，PP.H0–H4 与「稳健性分层」是否改变。
# 输出为新文件，不覆盖任何已交付产物。
#
# Usage: Rscript s52d_N_sensitivity.R

suppressPackageStartupMessages({
  library(data.table)
  library(coloc)
})

PROJ <- "D:/endometriosis_project/11_sc_eqtl_mr_project"
COND <- file.path(PROJ, "00_data_raw/onek1k/cond_coloc_L2L3")
DIR2 <- list(L2 = file.path(PROJ, "00_data_raw/onek1k/coloc_sens_prep"),
             L3 = file.path(PROJ, "00_data_raw/onek1k/coloc_sens_prep_recheck"))
TAB  <- file.path(PROJ, "tables")
LOG  <- "D:/endometriosis_project/_s52d_N_sens.log"
OUT  <- file.path(TAB, "55_N_sensitivity_L2L3.csv")
OUT3 <- file.path(TAB, "55b_coloc_sensitivity_L3_N690.csv")

con <- file(LOG, open = "wt", encoding = "UTF-8")
sink(con, split = TRUE)

cat("=== s52d N 口径灵敏度（coloc.abf 的 eQTL N）===\n")
cat("time", format(Sys.time()), "\n")
cat("coloc", as.character(packageVersion("coloc")), "\n\n")

NEFF <- c(CD4_NC = 980L, Mono_NC = 690L, B_MEM = 970L, Mono_C = 851L, B_IN = 975L,
          CD8_S100B = 959L, CD4_ET = 980L, CD8_NC = 980L, CD8_ET = 980L, NK = 980L)

PAIRS <- list(
  L2 = list(gene="YME1L1", ct="CD4_NC", m="m_L2_CD4_NC.csv"),
  L3 = list(gene="ANXA4",  ct="Mono_NC", m="m_L3_Mono_NC.csv")
)
ZDET <- list()   # sdY 诊断

# ---------- 1) 逐位点：N=980 vs N=N_true，uncond coloc.abf，两个 keep 规则 ----------
run_abf_m <- function(m, p12, keep_rule, Nq, sdY = NULL) {
  k <- m[1 - r^2 >= keep_rule]
  s_scalar <- sum(k$n_cases) / sum(k$n_cases + k$n_controls)
  d1 <- list(beta = k$slope, varbeta = k$slope_se^2, snp = k$rsid, position = k$pos38,
             MAF = k$maf_e, type = "quant")
  if (is.null(sdY)) d1$N <- Nq else d1$sdY <- sdY
  d2 <- list(beta = k$b2, varbeta = k$se^2, snp = k$rsid, position = k$pos38,
             MAF = k$maf_o, N = k$n_cases + k$n_controls, type = "cc", s = s_scalar)
  ab <- tryCatch(suppressWarnings(coloc.abf(d1, d2, p12 = p12)), error = function(z) NULL)
  if (is.null(ab)) return(NULL)
  s <- ab$summary
  list(nsnp = as.integer(s["nsnps"]), H0 = s["PP.H0.abf"], H1 = s["PP.H1.abf"],
       H2 = s["PP.H2.abf"], H3 = s["PP.H3.abf"], H4 = s["PP.H4.abf"])
}

rows <- list()
for (nm in names(PAIRS)) {
  P <- PAIRS[[nm]]
  Ntrue <- NEFF[[P$ct]]
  m <- fread(file.path(COND, P$m))
  ld <- fread(file.path(COND, sprintf("r_to_%s_%s.csv", nm,
                ifelse(nm == "L2", "rs693965", "rs62133984"))))
  m <- merge(m, ld[, .(snp_grch37, r = r_to_cond, r2 = r2_to_cond)], by = "snp_grch37")
  m <- m[!is.na(r)]
  cat(sprintf("\n--- [%s] %s / %s : 交集 SNP=%d  N_true=%d  N_delivered=980\n",
              nm, P$gene, P$ct, nrow(m), Ntrue))
  # sdY 诊断（无 LD、全 SNP 集）
  for (Nq in c(980L, Ntrue)) {
    sd_est <- tryCatch(coloc:::sdY.est(m$slope_se^2, m$maf_e, Nq), error = function(e) NA_real_)
    cat(sprintf("    sdY.est(N=%-4d) = %.6f\n", Nq, sd_est))
    ZDET[[length(ZDET) + 1]] <- data.table(locus = nm, gene = P$gene, cell_type = P$ct,
                                           N = Nq, sdY_est = sd_est)
  }
  for (kr in c(0.01, 0.1)) {
    for (Nq in c(980L, Ntrue)) {
      r <- run_abf_m(m, 1e-5, kr, Nq)
      if (is.null(r)) next
      rows[[length(rows) + 1]] <- data.table(
        locus = nm, gene = P$gene, cell_type = P$ct, keep_rule = kr,
        N_used = Nq, N_true = Ntrue, nsnp = r$nsnp,
        PP.H0 = r$H0, PP.H1 = r$H1, PP.H2 = r$H2, PP.H3 = r$H3, PP.H4 = r$H4,
        verdict = ifelse(r$H4 >= 0.8, "shared(>=0.8)", ifelse(r$H4 >= 0.5, "moderate(0.5-0.8)", "<0.5")))
      cat(sprintf("    keep 1-r2>=%-4g  N=%-4d  nsnp=%-5d  H0=%.4f H1=%.4f H2=%.4f H3=%.4f H4=%.4f  [%s]\n",
                  kr, Nq, r$nsnp, r$H0, r$H1, r$H2, r$H3, r$H4,
                  ifelse(r$H4 >= 0.8, "shared(>=0.8)", ifelse(r$H4 >= 0.5, "moderate(0.5-0.8)", "<0.5"))))
    }
  }
}
R <- rbindlist(rows)
fwrite(R, OUT)
cat("\n写出 ", OUT, " rows=", nrow(R), "\n", sep = "")

# ---------- 2) L3：完整 9 组矩阵改用 N=690，与交付版（N=980）逐格比较 ----------
cat("\n\n===== L3 9 组参数矩阵：N=690 vs 交付版 N=980 =====\n")
run_abf_sens <- function(e, g, p12, Nq) {
  e <- copy(e); g <- copy(g)
  e[, a1 := toupper(a1)]; e[, a2 := toupper(a2)]
  e[, allele_key := paste0(pos38, "_", pmin(a1, a2), "_", pmax(a1, a2))]
  g[, ea := toupper(ea)]; g[, oa := toupper(oa)]
  g[, allele_key := paste0(pos38, "_", pmin(ea, oa), "_", pmax(ea, oa))]
  g <- g[!duplicated(allele_key)]
  m <- merge(e, g, by = "allele_key", all = FALSE)
  if (nrow(m) < 10) return(NULL)
  if ("pos38.x" %in% names(m)) { m[, pos38 := pos38.x]; m[, c("pos38.x", "pos38.y") := NULL] }
  m[, flip := (ea == a2) & (oa == a1)]
  m <- m[(ea == a1 & oa == a2) | (ea == a2 & oa == a1)]
  m[, palindromic := paste0(pmin(a1, a2), pmax(a1, a2)) %in% c("AT", "CG")]
  m[, af_mism := abs(af - eaf)]; m[, af_flip := abs(af - (1 - eaf))]
  m <- m[!(palindromic == TRUE & pmin(af_mism, af_flip) > 0.05 & abs(af_mism - af_flip) < 0.05)]
  m[palindromic == TRUE & af_flip < af_mism, flip := !flip]
  m[, b2 := ifelse(flip, -beta, beta)]
  m[, maf_e := pmin(af, 1 - af)]; m[, maf_o := pmin(eaf, 1 - eaf)]
  m <- m[is.finite(maf_e) & is.finite(maf_o) & maf_e > 0 & maf_o > 0 &
           is.finite(slope) & is.finite(slope_se) & slope_se > 0]
  m <- m[!is.na(rsid) & rsid != ""]
  setorder(m, pos38)
  if (nrow(m) < 10) return(NULL)
  N_out <- m$n_cases + m$n_controls
  d1 <- list(beta = m$slope, varbeta = m$slope_se^2, snp = m$rsid, position = m$pos38,
             MAF = m$maf_e, N = Nq, type = "quant")
  d2 <- list(beta = m$b2, varbeta = m$se^2, snp = m$rsid, position = m$pos38,
             MAF = m$maf_o, N = N_out, type = "cc", s = sum(m$n_cases) / sum(N_out))
  ab <- tryCatch(suppressWarnings(coloc.abf(d1, d2, p12 = p12)), error = function(z) NULL)
  if (is.null(ab)) return(NULL)
  s <- ab$summary
  list(nsnp = as.integer(s["nsnps"]), H4 = s["PP.H4.abf"], H3 = s["PP.H3.abf"],
       H2 = s["PP.H2.abf"], H1 = s["PP.H1.abf"], H0 = s["PP.H0.abf"])
}
WINDOWS <- c("500Kb", "1Mb", "2Mb")
P12S <- c("1e-4" = 1e-4, "1e-5" = 1e-5, "1e-6" = 1e-6)
old <- fread(file.path(TAB, "51b_coloc_sensitivity_L3_recheck.csv"))
rows3 <- list()
for (w in WINDOWS) {
  tag <- sprintf("ANXA4_Mono_NC_%s", w)
  e <- fread(file.path(DIR2$L3, sprintf("eqtl_%s.csv", tag)))
  g <- fread(file.path(DIR2$L3, sprintf("gwas_%s.csv", tag)))
  for (pn in names(P12S)) {
    r <- run_abf_sens(e, g, P12S[[pn]], 690L)
    o4 <- old[window == w & p12 == pn]$PP.H4
    rows3[[length(rows3) + 1]] <- data.table(
      locus = "L3", gene = "ANXA4", cell_type = "Mono_NC", window = w, p12 = pn,
      nsnp_N690 = r$nsnp, PP.H4_N690 = r$H4, PP.H3_N690 = r$H3,
      PP.H4_delivered_N980 = o4, delta_H4 = r$H4 - o4,
      robust_N690 = ifelse(r$H4 >= 0.8, "稳健", "敏感性不稳健"),
      robust_N980 = ifelse(o4 >= 0.8, "稳健", "敏感性不稳健"))
    cat(sprintf("  %-7s p12=%-5s  H4(N=690)=%.4f  H4(N=980)=%.4f  dH4=%+.4f  [%s -> %s]\n",
                w, pn, r$H4, o4, r$H4 - o4,
                ifelse(o4 >= 0.8, "稳健", "不稳健"),
                ifelse(r$H4 >= 0.8, "稳健", "不稳健")))
  }
}
R3 <- rbindlist(rows3)
fwrite(R3, OUT3)
cat("\n写出 ", OUT3, " rows=", nrow(R3), "\n", sep = "")
cat("\n--- 分层是否翻转 ---\n")
R3[, flip := robust_N690 != robust_N980]
print(R3[, .(n = .N, n_flip = sum(flip)), by = .(window)])
cat("\ntime", format(Sys.time()), "\n")
sink(); close(con)
