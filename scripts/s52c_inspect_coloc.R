# -*- coding: utf-8 -*-
# s52c_inspect_coloc.R -- 检查 coloc.abf 中 N 的进入方式（只读，不改任何产物）
suppressPackageStartupMessages(library(coloc))
out <- "D:/endometriosis_project/_s52c_coloc_source.txt"
con <- file(out, open = "wt", encoding = "UTF-8")
sink(con, split = FALSE)
cat("=== coloc 版本 ===\n"); print(packageVersion("coloc"))
for (fn in c("sdY.est", "abf", "coloc.abf", "p12.est", "sd.prior", "approx.bf.quant")) {
  cat("\n\n===== ", fn, " =====\n", sep = "")
  b <- tryCatch(getFromNamespace(fn, "coloc"), error = function(e) NULL)
  if (is.null(b)) { cat("(不存在)\n"); next }
  print(b)
}
cat("\n\n===== 判断 coloc.abf 是否把 N 传给 sdY.est =====\n")
b <- tryCatch(getFromNamespace("coloc.abf", "coloc"), error = function(e) NULL)
if (!is.null(b)) {
  txt <- paste(deparse(b), collapse = "\n")
  cat("包含 'sdY' 的行：\n")
  for (l in strsplit(txt, "\n")[[1]]) if (grepl("sdY|N \\*|N\\*|beta", l)) cat("   ", l, "\n")
}
sink(); close(con)
cat("WROTE ", out, "\n")
