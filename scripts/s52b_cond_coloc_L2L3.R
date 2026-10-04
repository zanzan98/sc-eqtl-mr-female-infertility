# -*- coding: utf-8 -*-
# s52b_cond_coloc_L2L3.R -- 任务A：L2/L3 条件共定位（严格复刻 s40c_cond_coloc.R v2）
#
# 条件变异 = 本座 outcome-only susie 精细定位 argmax（规则 R2，已由 L1 定标验证）
#   L2 -> 由 tables/52b_finemap_pair_L2.csv (scope=outcome_only) 的 argmax_rsid/argmax_pos37 决定
#   L3 -> 同上，来自 52b_finemap_pair_L3.csv
#
# 方法（与 s40c 逐字符一致）：
#   z_cond = (z_j - r_jk z_k)/sqrt(1-r_jk^2)
#   beta_cond_j = beta_j - r_jk (se_j/se_k) beta_k ; var_cond_j = se_j^2 (1-r_jk^2)
#   三口径 uncond / cond_eqtl_only / cond_both；keep 阈值 1-r2 >= 0.01 / 0.1
#   条件方差膨胀闸门 + 程序性阴性对照（两性状 |z| 均最小且 MAF>0.05）
#
# Usage: Rscript s52b_cond_coloc_L2L3.R

suppressPackageStartupMessages({
  library(data.table)
  library(coloc)
})

PROJ <- "D:/endometriosis_project/11_sc_eqtl_mr_project"
COND <- file.path(PROJ, "00_data_raw/onek1k/cond_coloc_L2L3")
TAB  <- file.path(PROJ, "tables")
LOG  <- "D:/endometriosis_project/_s52b_cond_coloc.log"
OUT  <- file.path(TAB, "54_cond_coloc_L2L3.csv")
DET  <- file.path(TAB, "54b_cond_coloc_detail_L2L3.csv")
SUM  <- file.path(TAB, "54c_cond_coloc_summary_L2L3.csv")

con <- file(LOG, open = "wt", encoding = "UTF-8")
sink(con, split = TRUE); sink(con, type = "message")

N_EQTL <- 980   # 兼容旧口径；实际按位点用 NEFF[[ct]]（见下），差异已量化（|dH4| <= 4e-4）
BONF <- 0.05 / 8612
ZBONF <- qnorm(1 - BONF / 2)   # ≈5.4529 —— 与 chr1 链同口径，便于三位点并排
WMAP <- character(0)

cat("=== s52b L2/L3 条件共定位（复刻 s40c v2）===\n")
cat(sprintf("Bonferroni(8612) p = %.3g  对应双边 z 阈值 = %.4f\n", BONF, ZBONF))
cat("time ", format(Sys.time()), "\n\n")

cwrap <- function(expr) {
  withCallingHandlers(expr, warning = function(w) {
    WMAP <<- c(WMAP, conditionMessage(w)); invokeRestart("muffleWarning")
  })
}

pairs <- list(
  L2 = list(locus="L2", gene="YME1L1", ct="CD4_NC", m_file="m_L2_CD4_NC.csv",
            pair_file="52b_finemap_pair_L2.csv", N=980L),
  L3 = list(locus="L3", gene="ANXA4",  ct="Mono_NC", m_file="m_L3_Mono_NC.csv",
            pair_file="52b_finemap_pair_L3.csv", N=690L)
)

run_one <- function(m, keep_rule, cond_row, scope, Nq = N_EQTL) {
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
              MAF = m$maf_e, N = Nq, type = "quant")
  d2  <- list(beta = m$b2,      varbeta = m$se^2,       snp = m$rsid, position = m$pos38,
              MAF = m$maf_o, N = m$N_out, type = "cc", s = s_scalar)
  d1c <- list(beta = m$slope_c, varbeta = m$var_e_c,    snp = m$rsid, position = m$pos38,
              MAF = m$maf_e, N = Nq, type = "quant")
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
  rbindlist(list(mk(cwrap(coloc.abf(d1,  d2 )), "uncond"),
                 mk(cwrap(coloc.abf(d1c, d2 )), "cond_eqtl_only"),
                 mk(cwrap(coloc.abf(d1c, d2c)), "cond_both")), fill = TRUE)
}

res <- list(); det <- list(); summ <- list()
for (nm in names(pairs)) {
  P <- pairs[[nm]]
  cat("\n----------------------------------------\n[", P$locus, "] ", P$gene, "/", P$ct, "\n", sep = "")
  pf <- fread(file.path(TAB, P$pair_file))
  p0 <- pf[scope == "outcome_only"]
  if (nrow(p0) != 1L) { cat("  !! 缺 outcome_only 精细定位行，跳过\n"); next }
  COND_RS  <- as.character(p0$argmax_rsid[1])
  COND_P37 <- as.integer(p0$argmax_pos37[1])
  COND_ID37 <- sprintf("%s:%d", ifelse(P$locus == "L2", 10L, 2L), COND_P37)
  cat(sprintf("  outcome-only argmax = %s (GRCh37 %s, %s)  PIP=%.4f  nCS=%d\n",
              COND_RS, COND_ID37, as.character(p0$argmax_gene[1]), p0$max_pip[1], p0$n_cs[1]))

  rfile <- file.path(COND, sprintf("r_to_%s_%s.csv", P$locus, COND_RS))
  if (!file.exists(rfile)) { cat("  !! 缺 LD 表 ", rfile, "，跳过\n"); next }
  ld <- fread(rfile)[, .(snp_grch37, r = r_to_cond, r2 = r2_to_cond)]
  m0 <- fread(file.path(COND, P$m_file))
  m0 <- merge(m0, ld, by = "snp_grch37", all.x = TRUE)
  n0 <- nrow(m0)
  m0 <- m0[!is.na(r)]
  m0[, z_e := slope / slope_se][, z_o := b2 / se]
  cat(sprintf("  交集 SNP = %d ；有 LD 的 = %d\n", n0, nrow(m0)))

  ck <- m0[snp_grch37 == COND_ID37]
  if (nrow(ck) != 1L) { cat("  !! 条件变异不唯一/不在交集，跳过\n"); next }
  cand <- m0[snp_grch37 != COND_ID37 & maf_e > 0.05 & maf_o > 0.05]
  if (nrow(cand) == 0) { cat("  !! 无阴性对照候选，跳过\n"); next }
  cand[, score := pmax(abs(z_e), abs(z_o))]
  nc <- cand[order(score)][1]
  cat(sprintf("  阴性对照条件变异 = %s (%s) 其 |z_e|=%.4f |z_o|=%.4f r(与条件)=%+.4f\n",
              nc$rsid, nc$snp_grch37, abs(nc$z_e), abs(nc$z_o), nc$r))
  cat(sprintf("  |r| 分布：>0.99=%d ; >0.9=%d ; >0.8=%d\n",
              as.integer(sum(abs(m0$r) > 0.99)) - 1L,
              as.integer(sum(abs(m0$r) > 0.9)), as.integer(sum(abs(m0$r) > 0.8))))
  for (kr in c(1e-9, 0.01, 0.1)) {
    keep <- m0[1 - r^2 >= kr]
    mx <- max(1 / sqrt(pmax(1 - keep$r^2, 1e-300)))
    cat(sprintf("  keep 1-r2 >= %-6g -> nsnp=%d（剔除 %d） 最大膨胀因子=%.3f（|r|=%.4f）\n",
                kr, nrow(keep), nrow(m0) - nrow(keep), mx, max(abs(keep$r))))
  }
  for (kr in c(0.01, 0.1)) {
    keep <- m0[1 - r^2 >= kr]
    for (knd in c("lead", "negctl")) {
      cr <- if (knd == "lead") ck else nc
      KK <- copy(keep)
      rr <- run_one(KK, kr, cr, P$locus, P$N)
      rr[, `:=`(tag_pair = sprintf("%s|%s", P$locus, knd))]
      res[[length(res) + 1]] <- rr
      cat(sprintf("  [%s|cond=%s|1-r2>=%g] nsnp=%d  uncond H4=%.4f | cond_eqtl_only H4=%.4f | cond_both H4=%.4f (H2=%.4f,H3=%.4f)\n",
                  P$locus, cr$rsid, kr, nrow(keep),
                  rr[tag == "uncond"]$PP.H4[1], rr[tag == "cond_eqtl_only"]$PP.H4[1],
                  rr[tag == "cond_both"]$PP.H4[1],
                  rr[tag == "cond_both"]$PP.H2[1], rr[tag == "cond_both"]$PP.H3[1]))
      if (abs(kr - 0.01) < 1e-12 && knd == "lead") {
        m <- KK
        m[, `:=`(abs_z_e = abs(z_e), abs_z_e_c = abs(z_e_c), abs_z_o = abs(z_o), abs_z_o_c = abs(z_o_c))]
        top <- m[order(-abs_z_e)][1:min(15, .N)]
        cat(sprintf("    主口径：条件后 eQTL max|z|=%.3f（超 Bonf z=%.2f 个数=%d）；结局侧条件后 max|z|=%.3f（超 Bonf 个数=%d）\n",
                    max(m$abs_z_e_c, na.rm = TRUE), ZBONF, as.integer(sum(m$abs_z_e_c > ZBONF, na.rm = TRUE)),
                    max(m$abs_z_o_c, na.rm = TRUE), as.integer(sum(m$abs_z_o_c > ZBONF, na.rm = TRUE))))
        for (i in seq_len(nrow(top))) {
          det[[length(det) + 1]] <- data.table(
            locus = P$locus, keep_rule = kr, rank = i, rsid = top$rsid[i],
            snp_grch37 = top$snp_grch37[i], pos38 = top$pos38[i], r_to_cond = top$r[i],
            z_eqtl_pre = top$z_e[i], z_eqtl_post = top$z_e_c[i],
            z_out_pre = top$z_o[i], z_out_post = top$z_o_c[i],
            max_abs_z_eqtl_post = max(m$abs_z_e_c), n_eqtl_post_bonf = sum(m$abs_z_e_c > ZBONF),
            max_abs_z_out_post = max(m$abs_z_o_c), n_out_post_bonf = sum(m$abs_z_o_c > ZBONF))
        }
      }
    }
  }
  # 汇总行
  for (kr in c(0.01, 0.1)) {
    a <- rbindlist(res)[tag_pair == sprintf("%s|lead", P$locus) & keep_rule == kr]
    b <- rbindlist(res)[tag_pair == sprintf("%s|negctl", P$locus) & keep_rule == kr]
    summ[[length(summ) + 1]] <- data.table(
      locus = P$locus, gene = P$gene, cell_type = P$ct,
      cond_snp = COND_RS, cond_snp37 = COND_ID37, keep_rule = kr,
      nsnps = a[tag == "uncond"]$nsnps,
      H4_uncond = a[tag == "uncond"]$PP.H4,
      H4_cond_eqtl_only = a[tag == "cond_eqtl_only"]$PP.H4,
      H4_cond_both = a[tag == "cond_both"]$PP.H4,
      H2_cond_both = a[tag == "cond_both"]$PP.H2,
      H3_cond_both = a[tag == "cond_both"]$PP.H3,
      max_infl = a[tag == "cond_both"]$max_infl,
      negctl_H4_uncond = b[tag == "uncond"]$PP.H4,
      negctl_H4_cond_both = b[tag == "cond_both"]$PP.H4,
      negctl_cond_snp = nc$rsid)
  }
  fwrite(rbindlist(res, fill = TRUE), OUT)
  fwrite(rbindlist(det, fill = TRUE), DET)
  fwrite(rbindlist(summ, fill = TRUE), SUM)
  cat("  [incremental write] ", OUT, "\n", sep = "")
}

cat("\n===== 汇总（条件前后 PP.H4）=====\n")
R <- rbindlist(res, fill = TRUE)
print(dcast(R, tag_pair + keep_rule + cond_snp ~ tag, value.var = "PP.H4"))
cat("\n===== 主口径 cond_both 的 H2/H3/H4 与膨胀诊断 =====\n")
print(R[tag == "cond_both", .(tag_pair, keep_rule, nsnps, max_infl, PP.H2, PP.H3, PP.H4)])
cat("\n对照：本座无条件 coloc（34_coloc_summary.csv）\n")
old <- tryCatch(fread(file.path(TAB, "34_coloc_summary.csv")), error = function(z) NULL)
if (!is.null(old)) print(old[gene %in% c("YME1L1", "ANXA4"),
                            .(locus_id, gene, cell_type, abf_H3, abf_H4, abf_verdict, susie_verdict)])

cat("\n===== coloc 警告（", length(WMAP), " 条）=====\n", sep = "")
if (length(WMAP)) print(unique(WMAP))
cat("\ntime ", format(Sys.time()), "\n")
sink(type = "message"); sink(); close(con)
