# -*- coding: utf-8 -*-
# s40d_cond_coloc_negctl.R -- 任务 3.1 补丁：**配对正确**的程序性阴性对照
#
# 目的：修正 s40c 负对照臂的实现缺陷（LD 列硬绑 rs56318008，与条件变异不配对）。
#       本脚本对每个条件变异使用**其自身**的带符号 LD（r_to_negctl.csv），
#       检验「条件机制不会无差别摧毁共定位」。
#
# 与主结果的对照关系：
#   主口径  cond = rs56318008（WNT4 可信集 lead，自身即 CDC42 强 cis-eQTL）
#   负对照  cond = 两性状 |z| 均最小且 MAF>0.05 的变异（B_MEM=rs2473247, Mono_NC=rs12048511）
#
# 方法同 s40c v2：z_cond = (z_j - r_jk z_k)/sqrt(1-r_jk^2)；(beta,varbeta) 等价形式；
#                 coloc.abf 三口径（uncond / cond_eqtl_only / cond_both）；
#                 keep 阈值 1-r2>=0.01（主）与 0.1（敏感性），均基于**该条件变异自身**的 LD。
#
# 输出：tables/43c_cond_coloc_negctl.csv
# 日志：D:\\endometriosis_project\\_s40d_cond_coloc_negctl.log
# 纪律：只读既有产物；只新增文件；结果写文件（R 退出码可能 139，看产物不看退出码）
suppressPackageStartupMessages({
  library(data.table)
  library(coloc)
})

PROJ <- "D:/endometriosis_project/11_sc_eqtl_mr_project"
COND <- file.path(PROJ, "00_data_raw/onek1k/cond_coloc")
TAB  <- file.path(PROJ, "tables")
LOG  <- "D:/endometriosis_project/_s40d_cond_coloc_negctl.log"
OUT  <- file.path(TAB, "43c_cond_coloc_negctl.csv")

con <- file(LOG, open = "wt", encoding = "UTF-8")
sink(con, split = TRUE); sink(con, type = "message")

COND_SNP37 <- "1:22470407"
N_EQTL <- 980
BONF <- 0.05 / 8612
ZBONF <- qnorm(1 - BONF / 2)
WMAP <- character(0)

cat("=== s40d 配对正确的阴性对照（条件变异 = 其自身 LD） ===\n")
cat(sprintf("Bonferroni（8 612 检验）p = %.3g  双边 z 阈值 = %.4f\n", BONF, ZBONF))
cat("time ", format(Sys.time()), "\n\n")

cwrap <- function(expr) {
  withCallingHandlers(expr, warning = function(w) {
    WMAP <<- c(WMAP, conditionMessage(w)); invokeRestart("muffleWarning")
  })
}

ld <- fread(file.path(COND, "r_to_negctl.csv"))
cat("r_to_negctl.csv rows = ", nrow(ld), "  cols = ", ncol(ld), "\n")

res <- list()

run_one <- function(m, keep_rule, cond_row, scope) {
  m[, `:=`(
    slope_c = slope - r * (slope_se / cond_row$slope_se) * cond_row$slope,
    var_e_c = slope_se^2 * (1 - r^2),
    b2_c    = b2    - r * (se     / cond_row$se)     * cond_row$b2,
    var_o_c = se^2 * (1 - r^2))]
  zk_e <- cond_row$slope / cond_row$slope_se
  zk_o <- cond_row$b2 / cond_row$se
  m[, `:=`(z_e = slope / slope_se, z_o = b2 / se,
           z_e_c = (slope / slope_se - r * zk_e) / sqrt(pmax(1 - r^2, 1e-12)),
           z_o_c = (b2    / se       - r * zk_o) / sqrt(pmax(1 - r^2, 1e-12)))]
  chk <- max(abs(m$slope_c / sqrt(m$var_e_c) - m$z_e_c), na.rm = TRUE)
  s_scalar <- sum(m$n_cases) / sum(m$N_out)
  d1  <- list(beta = m$slope,   varbeta = m$slope_se^2, snp = m$rsid, position = m$pos38,
              MAF = m$maf_e, N = N_EQTL, type = "quant")
  d2  <- list(beta = m$b2,      varbeta = m$se^2,       snp = m$rsid, position = m$pos38,
              MAF = m$maf_o, N = m$N_out, type = "cc", s = s_scalar)
  d1c <- list(beta = m$slope_c, varbeta = m$var_e_c,    snp = m$rsid, position = m$pos38,
              MAF = m$maf_e, N = N_EQTL, type = "quant")
  d2c <- list(beta = m$b2_c,    varbeta = m$var_o_c,    snp = m$rsid, position = m$pos38,
              MAF = m$maf_o, N = m$N_out, type = "cc", s = s_scalar)
  mk <- function(a, side) {
    s <- a$summary
    data.table(scope = scope, tag = side, keep_rule = keep_rule,
               cond_snp = cond_row$rsid, cond_snp37 = cond_row$snp_grch37,
               nsnps = as.integer(s["nsnps"]),
               PP.H0 = s["PP.H0.abf"], PP.H1 = s["PP.H1.abf"], PP.H2 = s["PP.H2.abf"],
               PP.H3 = s["PP.H3.abf"], PP.H4 = s["PP.H4.abf"],
               max_infl = max(1 / sqrt(pmax(1 - m$r^2, 1e-300))),
               max_abs_z_e_pre = max(abs(m$z_e)),  n_e_pre_bonf  = sum(abs(m$z_e)   > ZBONF),
               max_abs_z_o_pre = max(abs(m$z_o)),  n_o_pre_bonf  = sum(abs(m$z_o)   > ZBONF),
               max_abs_z_e_post = max(abs(m$z_e_c)), n_e_post_bonf = sum(abs(m$z_e_c) > ZBONF),
               max_abs_z_o_post = max(abs(m$z_o_c)), n_o_post_bonf = sum(abs(m$z_o_c) > ZBONF),
               zcheck = chk)
  }
  rbindlist(list(
    mk(cwrap(coloc.abf(d1,  d2 )), "uncond"),
    mk(cwrap(coloc.abf(d1c, d2 )), "cond_eqtl_only"),
    mk(cwrap(coloc.abf(d1c, d2c)), "cond_both")))
}

for (CT in c("B_MEM", "Mono_NC")) {
  cat("\n----------------------------------------\n[", CT, "]\n", sep = "")
  m0 <- fread(file.path(COND, sprintf("m_%s.csv", CT)))
  m0 <- merge(m0, ld, by = "snp_grch37", all.x = TRUE)
  m0[, z_e := slope / slope_se][, z_o := b2 / se]
  m0 <- m0[!is.na(pos37)]
  cat("  交集 SNP（有 LD 的）= ", nrow(m0), "\n", sep = "")

  # 复现 s40c 的阴性对照选取规则（两性状 |z| 均最小且 MAF>0.05，排除主条件变异）
  cand <- m0[snp_grch37 != COND_SNP37 & maf_e > 0.05 & maf_o > 0.05]
  cand[, score := pmax(abs(z_e), abs(z_o))]
  nc <- cand[order(score)][1]
  cat(sprintf("  阴性对照条件变异 = %s (%s)  |z_eqtl|=%.4f  |z_outcome|=%.4f\n",
              nc$rsid, nc$snp_grch37, abs(nc$z_e), abs(nc$z_o)))
  ldcol <- sprintf("r_to_%s", nc$rsid)
  if (!ldcol %in% names(m0)) { cat("  !! 缺少配对 LD 列 ", ldcol, "，跳过\n", sep = ""); next }
  m0[, r := get(ldcol)]                                # ★ 配对：r_jk 用条件变异自身的 LD
  cat(sprintf("  配对 LD 列 = %s   其自身 r=%.6f (应=1)  与 rs56318008 的 r=%+.4f\n",
              ldcol, m0[snp_grch37 == nc$snp_grch37]$r,
              m0[snp_grch37 == nc$snp_grch37]$r_to_rs56318008_check))
  ck <- m0[snp_grch37 == nc$snp_grch37]
  if (nrow(ck) != 1L) { cat("  !! 条件变异不唯一，跳过\n"); next }

  for (kr in c(0.01, 0.1)) {
    keep <- m0[1 - r^2 >= kr]
    mx <- max(1 / sqrt(pmax(1 - keep$r^2, 1e-300)))
    cat(sprintf("  keep 1-r2>=%-5g -> nsnp=%d（剔除 %d）  max_infl=%.3f（|r|=%.4f）\n",
                kr, nrow(keep), nrow(m0) - nrow(keep), mx, max(abs(keep$r))))
    rr <- run_one(copy(keep), kr, ck, CT)
    res[[length(res) + 1]] <- rr
    cat(sprintf("  [%s|cond=%s|1-r2>=%g] nsnp=%d  uncond H4=%.4f | cond_eqtl_only H4=%.4f | cond_both H4=%.4f (H2=%.4f,H3=%.4f)\n",
                CT, nc$rsid, kr, nrow(keep),
                rr[tag == "uncond"]$PP.H4[1], rr[tag == "cond_eqtl_only"]$PP.H4[1],
                rr[tag == "cond_both"]$PP.H4[1],
                rr[tag == "cond_both"]$PP.H2[1], rr[tag == "cond_both"]$PP.H3[1]))
    cat(sprintf("     条件后 eQTL max|z|=%.3f（超 Bonf 数=%d）；结局 max|z|=%.3f（超 Bonf 数=%d）\n",
                rr[tag == "cond_both"]$max_abs_z_e_post[1], rr[tag == "cond_both"]$n_e_post_bonf[1],
                rr[tag == "cond_both"]$max_abs_z_o_post[1], rr[tag == "cond_both"]$n_o_post_bonf[1]))
  }
  fwrite(rbindlist(res, fill = TRUE), OUT)
  cat("  [incremental write] ", OUT, "\n", sep = "")
}

R <- rbindlist(res, fill = TRUE)
cat("\n===== 阴性对照汇总（PP.H4）=====\n")
print(dcast(R, tag + keep_rule + cond_snp ~ scope, value.var = "PP.H4"))
cat("\n===== 阴性对照 cond_both 分解 =====\n")
print(R[tag == "cond_both", .(scope, keep_rule, cond_snp, nsnps, max_infl,
                              max_abs_z_e_pre, n_e_pre_bonf, max_abs_z_e_post, n_e_post_bonf,
                              PP.H2, PP.H3, PP.H4)])

cat("\n===== coloc 警告（", length(WMAP), " 条）=====\n", sep = "")
if (length(WMAP)) print(unique(WMAP))
cat("\ntime ", format(Sys.time()), "\n")
sink(type = "message"); sink(); close(con)
