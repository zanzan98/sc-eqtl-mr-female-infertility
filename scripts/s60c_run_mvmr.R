# -*- coding: utf-8 -*-
# s60c_run_mvmr.R -- 任务B：CDC42 + LINC00339 双暴露 MVMR（逐 OneK1K 细胞类型）
#
# 暴露：CDC42 顺式 eQTL 与 LINC00339 顺式 eQTL（同一细胞类型）
# 结局：GCST90483463 女性不孕症
# 工具：两基因 cis-eQTL P<5e-8 之并集，再按 LD 逐阈值剪枝为独立集
#
# 逐 (细胞类型 × 剪枝阈值) 报告：
#   ① 条件 F 统计量（MVMR::strength_mvmr，Sanderson 2019）
#   ② MVMR-IVW 效应估计与 95% CI（MVMR::ivw_mvmr）
#   ③ LD 条件数（工具子集相关阵 λmax/λmin）与最大两两 r²
#   ④ 可识别性判定：k>=2 且 两暴露条件 F>=10 且 LD 条件数<=30 → 可识别
#
# 读：00_data_raw/mvmr_L1/<cell>.mvmr_input.csv 与 ld/<cell>.ld
# 写：tables/61_mvmr_results.csv、tables/61b_mvmr_ld_diag.csv、scripts/_s60c_run_mvmr.log
suppressMessages({
  library(MVMR)
  library(data.table)
})

PROJ <- "D:/endometriosis_project/11_sc_eqtl_mr_project"
IN <- file.path(PROJ, "00_data_raw/mvmr_L1")
LDD <- file.path(IN, "ld")
TAB <- file.path(PROJ, "tables")
LOG <- file.path(PROJ, "scripts/_s60c_run_mvmr.log")

CELLS <- c("B_IN", "B_MEM", "CD4_ET", "CD4_NC", "CD8_ET", "CD8_NC", "CD8_S100B", "NK")
THRS <- c(0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.5)
F_COND_MIN <- 10
LD_COND_MAX <- 30

L <- character(0)
p <- function(...) {
  s <- paste0(...)
  L <<- c(L, s)
  cat(s, "\n")
}

p("====================================================================================================")
p("s60c — MVMR：CDC42 + LINC00339（逐细胞类型 × LD 剪枝阈值）")
p("====================================================================================================")
p(sprintf("MVMR 包版本 = %s", as.character(packageVersion("MVMR"))))
p(sprintf("阈值集 = %s ；条件 F 门槛 = %g ；LD 条件数门槛 = %g",
          paste(THRS, collapse = ", "), F_COND_MIN, LD_COND_MAX))

res <- list()
diag_all <- list()
k <- 0

# 稳健行构造：强制每列为长度 1（防 strength_mvmr / ivw_mvmr 返回非标量）
mkrow <- function(...) {
  v <- list(...)
  v <- lapply(v, function(z) {
    if (is.null(z)) return(NA)
    z <- unlist(z, use.names = FALSE)
    if (length(z) == 0) return(NA)
    z[1]
  })
  as.data.frame(v, stringsAsFactors = FALSE, check.names = FALSE)
}

for (cell in CELLS) {
  fin <- file.path(IN, sprintf("%s.mvmr_input.csv", cell))
  fld <- file.path(LDD, sprintf("%s.ld", cell))
  if (!file.exists(fin) || !file.exists(fld)) {
    p(sprintf("[%s] 输入缺失，跳过", cell)); next
  }
  D <- as.data.frame(fread(fin))
  names(D)[1] <- "snp"          # 防 utf-8-sig BOM 污染首列名
  LD <- as.matrix(read.table(fld, header = FALSE))
  if (nrow(LD) != nrow(D)) {
    p(sprintf("[%s] **维度不符** nrow(LD)=%d nrow(D)=%d，跳过", cell, nrow(LD), nrow(D))); next
  }
  p("")
  p("----------------------------------------------------------------------------------------------------")
  p(sprintf("[%s] 工具并集（pre-clump P<5e-8）= %d 个变异", cell, nrow(D)))
  # 全局共线性（全并集）
  ev_all <- eigen(LD, symmetric = TRUE, only.values = TRUE)$values
  r2off <- LD^2; diag(r2off) <- NA
  p(sprintf("  全并集 LD：λmax=%.3f λmin=%.3e 条件数=%.3e ；最大两两 r²=%.4f ；平均 r²=%.4f",
            max(ev_all), min(ev_all), max(ev_all)/max(min(ev_all), 1e-300),
            max(r2off, na.rm = TRUE), mean(r2off, na.rm = TRUE)))
  diag_all[[length(diag_all) + 1]] <- data.frame(
    cell = cell, level = "all_preclump", thr = NA, k = nrow(D),
    ld_cond = max(ev_all)/max(min(ev_all), 1e-300),
    max_r2 = max(r2off, na.rm = TRUE), mean_r2 = mean(r2off, na.rm = TRUE))

  ord <- order(pmin(D$p_cdc42, D$p_linc00339))
  for (thr in THRS) {
    keep <- integer(0)
    for (i in ord) {
      if (length(keep) == 0) { keep <- i; next }
      if (all((LD[i, keep]^2) < thr)) keep <- c(keep, i)
    }
    kk <- length(keep)
    if (kk < 2) {
      res[[length(res) + 1]] <- mkrow(
        cell = cell, thr = thr, k = kk, identifiable = "否（工具数 < 2，模型不可识别）",
        F_cond_CDC42 = NA, F_cond_LINC00339 = NA, F_ind_CDC42 = NA, F_ind_LINC00339 = NA,
        b_CDC42 = NA, se_CDC42 = NA, lo_CDC42 = NA, hi_CDC42 = NA, p_CDC42 = NA,
        b_LINC00339 = NA, se_LINC00339 = NA, lo_LINC00339 = NA, hi_LINC00339 = NA, p_LINC00339 = NA,
        ld_cond = NA, max_r2 = NA)
      next
    }
    BX <- as.matrix(D[keep, c("b_cdc42", "b_linc00339")])
    seBX <- as.matrix(D[keep, c("se_cdc42", "se_linc00339")])
    BY <- D$b_out[keep]; seBY <- D$se_out[keep]
    fmt <- format_mvmr(BXGs = BX, BYG = BY, seBXGs = seBX, seBYG = seBY, RSID = D$snp[keep])
    st <- tryCatch(suppressWarnings(strength_mvmr(fmt, gencov = 0)), error = function(e) NULL)
    # ivw_mvmr 返回 summary(lm)$coef 矩阵（非 lm 对象）；capture.output 抑制其打印
    # ★临时变量名不得用 res（res 为累加用列表）
    .co <- NULL
    .ivtxt <- tryCatch(capture.output(.co <- suppressWarnings(ivw_mvmr(fmt, gencov = 0))),
                       error = function(e) NULL)
    co <- .co
    F1 <- if (!is.null(st)) st[[1]] else NA
    F2 <- if (!is.null(st)) st[[2]] else NA
    # 单暴露（边际）F：mean (b/se)^2
    Fi1 <- mean((BX[, 1]/seBX[, 1])^2); Fi2 <- mean((BX[, 2]/seBX[, 2])^2)
    b1 <- se1 <- b2 <- se2 <- pv1 <- pv2 <- NA_real_
    if (!is.null(co) && is.matrix(co) && nrow(co) >= 2) {
      b1 <- co[1, 1]; se1 <- co[1, 2]; pv1 <- co[1, 4]
      b2 <- co[2, 1]; se2 <- co[2, 2]; pv2 <- co[2, 4]
    }
    M <- LD[keep, keep]
    ev <- eigen(M, symmetric = TRUE, only.values = TRUE)$values
    cnd <- max(ev)/max(min(ev), 1e-300)
    r2 <- M^2; diag(r2) <- NA
    ident <- (kk >= 2 && !is.na(F1) && !is.na(F2) && F1 >= F_COND_MIN && F2 >= F_COND_MIN && cnd <= LD_COND_MAX)
    verdict <- if (ident) "可识别" else {
      why <- c()
      if (!is.na(F1) && F1 < F_COND_MIN) why <- c(why, sprintf("CDC42 条件F=%.2f<%.0f", F1, F_COND_MIN))
      if (!is.na(F2) && F2 < F_COND_MIN) why <- c(why, sprintf("LINC00339 条件F=%.2f<%.0f", F2, F_COND_MIN))
      if (cnd > LD_COND_MAX) why <- c(why, sprintf("LD条件数=%.3e>%g", cnd, LD_COND_MAX))
      paste0("否（", paste(why, collapse = "；"), "）")
    }
    p(sprintf("  thr=%-6g k=%-4d 条件F: CDC42=%s LINC00339=%s ｜ 单暴露F: %.1f / %.1f ｜ LD条件数=%.3e max r²=%.4f ｜ %s",
              thr, kk,
              ifelse(is.na(F1), "NA", sprintf("%.3f", F1)),
              ifelse(is.na(F2), "NA", sprintf("%.3f", F2)),
              Fi1, Fi2, cnd, max(r2, na.rm = TRUE), verdict))
    if (!is.null(co) && is.matrix(co) && nrow(co) >= 2) {
      fmtci <- function(b, se) if (is.na(se) || se == 0) sprintf("b=%+.5f (se=不可估计：恰好识别)", b)
                              else sprintf("b=%+.5f (se %.5f, 95%%CI %+.5f..%+.5f, p=%.3g)", b, se, b-1.96*se, b+1.96*se, pv1)
      p(sprintf("           MVMR-IVW：CDC42 %s", fmtci(b1, se1)))
      p(sprintf("                     LINC00339 %s",
                if (is.na(se2) || se2 == 0) sprintf("b=%+.5f (se=不可估计：恰好识别)", b2)
                else sprintf("b=%+.5f (se %.5f, 95%%CI %+.5f..%+.5f, p=%.3g)", b2, se2, b2-1.96*se2, b2+1.96*se2, pv2)))
    }
    res[[length(res) + 1]] <- mkrow(
      cell = cell, thr = thr, k = kk, identifiable = verdict,
      F_cond_CDC42 = F1, F_cond_LINC00339 = F2,
      F_ind_CDC42 = Fi1, F_ind_LINC00339 = Fi2,
      b_CDC42 = b1, se_CDC42 = se1, lo_CDC42 = b1 - 1.96*se1, hi_CDC42 = b1 + 1.96*se1, p_CDC42 = pv1,
      b_LINC00339 = b2, se_LINC00339 = se2, lo_LINC00339 = b2 - 1.96*se2,
      hi_LINC00339 = b2 + 1.96*se2, p_LINC00339 = pv2,
      ld_cond = cnd, max_r2 = max(r2, na.rm = TRUE))
    diag_all[[length(diag_all) + 1]] <- data.frame(
      cell = cell, level = sprintf("pruned_r2_%g", thr), thr = thr, k = kk,
      ld_cond = cnd, max_r2 = max(r2, na.rm = TRUE), mean_r2 = mean(r2, na.rm = TRUE))
  }
}

R <- do.call(rbind, res)
DD <- do.call(rbind, diag_all)
write.csv(R, file.path(TAB, "61_mvmr_results.csv"), row.names = FALSE)
write.csv(DD, file.path(TAB, "61b_mvmr_ld_diag.csv"), row.names = FALSE)

p("")
p("====================================================================================================")
p("结果总览")
p("====================================================================================================")
p(sprintf("%-9s %-8s %4s %10s %10s %12s %9s  %s",
          "cell", "thr", "k", "F_CDC42", "F_LINC", "LD_cond", "max_r2", "可识别"))
for (i in seq_len(nrow(R))) {
  r <- R[i, ]
  p(sprintf("%-9s %-8g %4d %10s %10s %12s %9s  %s",
            r$cell, r$thr, r$k,
            ifelse(is.na(r$F_cond_CDC42), "NA", sprintf("%.3f", r$F_cond_CDC42)),
            ifelse(is.na(r$F_cond_LINC00339), "NA", sprintf("%.3f", r$F_cond_LINC00339)),
            ifelse(is.na(r$ld_cond), "NA", sprintf("%.3e", r$ld_cond)),
            ifelse(is.na(r$max_r2), "NA", sprintf("%.4f", r$max_r2)),
            r$identifiable))
}
p("")
p("★ 结论判定（k>=2 且两暴露条件 F>=10 且 LD 条件数<=30 方可识别）：")
n_ident <- sum(grepl("^可识别", R$identifiable))
p(sprintf("  可识别的 (细胞类型 × 阈值) 组合 = %d / %d", n_ident, nrow(R)))
writeLines(L, LOG)
cat("WROTE ", LOG, "\n")
