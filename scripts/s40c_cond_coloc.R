# -*- coding: utf-8 -*-
# s40c_cond_coloc.R  (v2)  -- 任务 3.1：以 rs56318008 (GRCh37 1:22470407, WNT4 可信集 lead) 为条件的共定位
#
# 方法（GCTA-COJO 型单 SNP 条件分析 + 外参 LD；coloc 侧用条件汇总统计量重跑 coloc.abf）：
#   等位单位：OneK1K .bim 的 A1（eQTL slope 与结局 b2 均已对齐至 A1；由 s40b 独立校验通过）
#   r       = OneK1K 980 供者中 A1 等位剂量的带符号相关系数
#   z_cond  = (z_j - r_jk * z_k) / sqrt(1 - r_jk^2)                <- 精确形式（同队列同 N 边际 z 的 MVN）
#   等价 (beta, varbeta)：beta_cond_j = beta_j - r_jk*(se_j/se_k)*beta_k ; var_cond_j = se_j^2*(1-r_jk^2)
#
# ★ v2 相对 v1 的两处方法学修正：
#   (1) **条件方差膨胀闸门**：|r|->1 时 var_cond -> 0，条件 z 会被 1/sqrt(1-r2) 放大而失控
#       （v1 保留集中仍有 1 个 |r|>0.99 的变异，导致结局侧条件后 max|z| 达 8.8e4 的数值伪影）。
#       v2 对保留集施加**预先设定的** 1-r2 下限，并逐口径报告最大膨胀因子 1/sqrt(1-r2)。
#   (2) **程序性阴性对照**：改用「在两性状中 |z| 均最小且 MAF>0.05」的变异作条件变异，
#       检验该条件机制不会无差别摧毁共定位（期望 PP.H4 与不条件时一致）。
#
# 输出：tables/43_cond_coloc_CDC42.csv        四类口径 × 两个 keep 阈值 × 两个条件变异
#       tables/43b_cond_coloc_variant_detail.csv  主口径下的逐变异条件前后 z
# 纪律：只读既有产物；增量落盘；全部结果写文件（R 退出码可能 139，看产物不看退出码）
# Usage: Rscript s40c_cond_coloc.R
suppressPackageStartupMessages({
  library(data.table)
  library(coloc)
})

PROJ <- "D:/endometriosis_project/11_sc_eqtl_mr_project"
COND <- file.path(PROJ, "00_data_raw/onek1k/cond_coloc")
TAB  <- file.path(PROJ, "tables")
LOG  <- "D:/endometriosis_project/_s40c_cond_coloc.log"
OUT  <- file.path(TAB, "43_cond_coloc_CDC42.csv")
DET  <- file.path(TAB, "43b_cond_coloc_variant_detail.csv")

con <- file(LOG, open = "wt", encoding = "UTF-8")
sink(con, split = TRUE); sink(con, type = "message")

COND_SNP37 <- "1:22470407"
N_EQTL <- 980
BONF <- 0.05 / 8612
ZBONF <- qnorm(1 - BONF / 2)   # ★ 计数必须用 z 阈值（≈5.4529），BONF 本身是 p 阈值，不可直接比 z
WMAP <- character(0)

cat("=== s40c v2 条件共定位（条件变异 = rs56318008 / ", COND_SNP37, "）===\n", sep = "")
cat(sprintf("Bonferroni（8 612 检验）p = %.3g  对应双边 z 阈值 = %.4f\n", BONF, ZBONF))
cat("time ", format(Sys.time()), "\n\n")

# ---- 收集 coloc 的 warning，不让它静默 ----
cwrap <- function(expr) {
  withCallingHandlers(expr, warning = function(w) {
    WMAP <<- c(WMAP, conditionMessage(w)); invokeRestart("muffleWarning")
  })
}

ld <- fread(file.path(COND, "r_to_rs56318008.csv"))
ld <- ld[, .(snp_grch37, r = r_to_cond, r2 = r2_to_cond)]
cat("LD 表 rows = ", nrow(ld), "\n")

res <- list(); det <- list()

run_one <- function(m, tag, keep_rule, cond_row, scope) {
  m[, `:=`(
    slope_c = slope - r * (slope_se / cond_row$slope_se) * cond_row$slope,
    var_e_c = slope_se^2 * (1 - r^2),
    b2_c    = b2    - r * (se     / cond_row$se)     * cond_row$b2,
    var_o_c = se^2 * (1 - r^2)
  )]
  zk_e <- cond_row$slope / cond_row$slope_se
  zk_o <- cond_row$b2 / cond_row$se
  m[, `:=`(z_e = slope / slope_se, z_o = b2 / se,
           z_e_c = (slope / slope_se - r * zk_e) / sqrt(pmax(1 - r^2, 1e-12)),
           z_o_c = (b2 / se - r * zk_o) / sqrt(pmax(1 - r^2, 1e-12)))]
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
    data.table(scope = scope, tag = side, keep_rule = keep_rule, cond_snp = cond_row$rsid,
               cond_snp37 = cond_row$snp_grch37,
               nsnps = as.integer(s["nsnps"]),
               PP.H0 = s["PP.H0.abf"], PP.H1 = s["PP.H1.abf"], PP.H2 = s["PP.H2.abf"],
               PP.H3 = s["PP.H3.abf"], PP.H4 = s["PP.H4.abf"],
               max_infl = max(1 / sqrt(pmax(1 - m$r^2, 1e-300))),
               max_abs_z_e_pre = max(abs(m$z_e)), n_e_pre_bonf = sum(abs(m$z_e) > ZBONF),
               max_abs_z_o_pre = max(abs(m$z_o)), n_o_pre_bonf = sum(abs(m$z_o) > ZBONF),
               max_abs_z_e_post = max(abs(m$z_e_c)), n_e_post_bonf = sum(abs(m$z_e_c) > ZBONF),
               max_abs_z_o_post = max(abs(m$z_o_c)), n_o_post_bonf = sum(abs(m$z_o_c) > ZBONF),
               zcheck = chk)
  }
  out <- list()
  out[[1]] <- mk(cwrap(coloc.abf(d1,  d2 )), "uncond")          # 两侧均不条件（对照）
  out[[2]] <- mk(cwrap(coloc.abf(d1c, d2 )), "cond_eqtl_only")  # 仅条件暴露侧
  out[[3]] <- mk(cwrap(coloc.abf(d1c, d2c)), "cond_both")       # 两侧均条件
  rbindlist(out)
}

for (CT in c("B_MEM", "Mono_NC")) {
  cat("\n----------------------------------------\n[", CT, "]\n", sep = "")
  m0 <- fread(file.path(COND, sprintf("m_%s.csv", CT)))
  m0 <- merge(m0, ld, by = "snp_grch37", all.x = TRUE)
  n0 <- nrow(m0)
  m0 <- m0[!is.na(r)]
  m0[, z_e := slope / slope_se][, z_o := b2 / se]
  cat("  交集 SNP = ", n0, " ；有 LD 的 = ", nrow(m0), "\n", sep = "")

  ck <- m0[snp_grch37 == COND_SNP37]
  if (nrow(ck) != 1L) { cat("  !! 条件变异不唯一，跳过\n"); next }
  # ---- 程序性阴性对照条件变异：两性状 |z| 均最小且 MAF>0.05（排除条件变异本身） ----
  cand <- m0[snp_grch37 != COND_SNP37 & maf_e > 0.05 & maf_o > 0.05]
  cand[, score := pmax(abs(z_e), abs(z_o))]
  nc <- cand[order(score)][1]
  cat(sprintf("  阴性对照条件变异 = %s (%s)  其 |z_eqtl|=%.4f  |z_outcome|=%.4f  r(与rs56318008)=%+.4f\n",
              nc$rsid, nc$snp_grch37, abs(nc$z_e), abs(nc$z_o), nc$r))

  # ---- LD 结构与 keep 阈值 ----
  cat(sprintf("  |r| 分布：>0.99（不含条件变异）= %d ；>0.9 = %d ；>0.8 = %d\n",
              as.integer(sum(abs(m0$r) > 0.99)) - 1L,
              as.integer(sum(abs(m0$r) > 0.9)), as.integer(sum(abs(m0$r) > 0.8))))
  for (kr in c(1e-9, 0.01, 0.1)) {
    keep <- m0[1 - r^2 >= kr]
    mx <- max(1 / sqrt(pmax(1 - keep$r^2, 1e-300)))
    cat(sprintf("  keep 1-r2 >= %-6g -> nsnp=%d（剔除 %d）  最大膨胀因子=%.3f（|r|=%.4f）\n",
                kr, nrow(keep), nrow(m0) - nrow(keep), mx, max(abs(keep$r))))
  }

  for (kr in c(0.01, 0.1)) {
    keep <- m0[1 - r^2 >= kr]
    for (nm in c("rs56318008", "negctl")) {
      cr <- if (nm == "rs56318008") ck else nc
      tag <- sprintf("%s|cond=%s|1-r2>=%g", CT, cr$rsid, kr)
      KK <- copy(keep)
      rr <- run_one(KK, sprintf("%s|cond=%s|1-r2>=%g", CT, cr$rsid, kr), kr, cr, CT)
      res[[length(res) + 1]] <- rr
      cat(sprintf("  [%s|%s|1-r2>=%g] nsnp=%d  uncond H4=%.4f | cond_eqtl_only H4=%.4f | cond_both H4=%.4f (H2=%.4f,H3=%.4f)\n",
                  CT, cr$rsid, kr, nrow(KK),
                  rr[tag == "uncond"]$PP.H4[1], rr[tag == "cond_eqtl_only"]$PP.H4[1],
                  rr[tag == "cond_both"]$PP.H4[1],
                  rr[tag == "cond_both"]$PP.H2[1], rr[tag == "cond_both"]$PP.H3[1]))
      if (abs(kr - 0.01) < 1e-12 && cr$rsid == "rs56318008") {
        m <- KK
        m[, `:=`(abs_z_e = abs(z_e), abs_z_e_c = abs(z_e_c), abs_z_o = abs(z_o), abs_z_o_c = abs(z_o_c))]
        top <- m[order(-abs_z_e)][1:min(15, .N)]
        cat(sprintf("    主口径逐变异：条件后约简 eQTL max|z|=%.3f（超 Bonf z=%.2f 的变异数=%d）；结局侧条件后 max|z|=%.3f（超 Bonf 的变异数=%d）\n",
                    max(m$abs_z_e_c, na.rm = TRUE), ZBONF, as.integer(sum(m$abs_z_e_c > ZBONF, na.rm = TRUE)),
                    max(m$abs_z_o_c, na.rm = TRUE), as.integer(sum(m$abs_z_o_c > ZBONF, na.rm = TRUE))))
        for (i in seq_len(nrow(top))) {
          cat(sprintf("      %-14s pos37=%-12s r=%+.4f  z_e=%+.3f->%+.3f   z_o=%+.3f->%+.3f\n",
                      top$rsid[i], top$snp_grch37[i], top$r[i],
                      top$z_e[i], top$z_e_c[i], top$z_o[i], top$z_o_c[i]))
          det[[length(det) + 1]] <- data.table(
            scope = CT, keep_rule = kr, rank = i, rsid = top$rsid[i],
            snp_grch37 = top$snp_grch37[i], pos38 = top$pos38[i], r_to_cond = top$r[i],
            z_eqtl_pre = top$z_e[i], z_eqtl_post = top$z_e_c[i],
            z_out_pre = top$z_o[i], z_out_post = top$z_o_c[i],
            max_abs_z_eqtl_post = max(m$abs_z_e_c), n_eqtl_post_bonf = sum(m$abs_z_e_c > ZBONF),
            max_abs_z_out_post = max(m$abs_z_o_c), n_out_post_bonf = sum(m$abs_z_o_c > ZBONF))
        }
      }
    }
  }
  fwrite(rbindlist(res, fill = TRUE), OUT)
  fwrite(rbindlist(det, fill = TRUE), DET)
  cat("  [incremental write] ", OUT, " rows=", length(res), "  ", DET, "\n", sep = "")
}

cat("\n===== 汇总（PP.H4）=====\n")
R <- rbindlist(res, fill = TRUE)
print(dcast(R, tag + keep_rule + cond_snp ~ scope, value.var = "PP.H4"))
cat("\n===== 汇总（cond_both 的 H2/H3/H4 分解）=====\n")
print(R[tag == "cond_both", .(scope, keep_rule, cond_snp, nsnps, max_infl,
                              max_abs_z_e_pre, n_e_pre_bonf, max_abs_z_e_post, n_e_post_bonf,
                              max_abs_z_o_pre, n_o_pre_bonf, max_abs_z_o_post, n_o_post_bonf,
                              PP.H2, PP.H3, PP.H4)])

cat("\n对照：既有 34_coloc_summary.csv\n")
old <- tryCatch(fread(file.path(TAB, "34_coloc_summary.csv")), error = function(z) NULL)
if (!is.null(old)) {
  have <- intersect(c("locus_id", "gene", "cell_type", "abf_verdict", "abf_H3", "abf_H4",
                      "abf_H4_over_H3H4"), names(old))
  print(old[gene == "CDC42" & cell_type %in% c("B_MEM", "Mono_NC"), ..have])
}

cat("\n===== coloc 警告（", length(WMAP), " 条）=====\n", sep = "")
if (length(WMAP)) print(unique(WMAP))
cat("\ntime ", format(Sys.time()), "\n")
sink(type = "message"); sink(); close(con)
