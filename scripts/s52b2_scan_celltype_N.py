# -*- coding: utf-8 -*-
"""s52b2_scan_celltype_N.py -- 扫描全部细胞类型 parquet 的有效样本量

原理：parquet 有 ma_count（minor allele 计数）与 af（minor allele freq）。
      N_eff = ma_count / af / 2 （每个 variant 应给出同一个整数 N_eff）
若 N_eff != 980，则该项目里用 980 供者 PLINK 面板做 LD 参考的步骤
（条件共定位 signed r、susie_rss 的 R、SMR 的 LD）存在 LD 参考与 eQTL 样本不同源的问题。
"""
import io, os, glob
import numpy as np
import pandas as pd

PROJ = r"D:\endometriosis_project\11_sc_eqtl_mr_project"
PQ = os.path.join(PROJ, "00_data_raw", "onek1k", "parquet_subset")
LOG = r"D:\endometriosis_project\_s52b2_scan_N.log"
L = []

def p(s=""):
    s = str(s); L.append(s)
    try: print(s)
    except UnicodeEncodeError: print(s.encode("ascii", "backslashreplace").decode("ascii"))

p("=== 扫描 parquet 有效样本量 N_eff = ma_count/af/2 ===")
p("（每个 (细胞类型, 染色体) 文件独立统计）")
p("")

files = sorted(glob.glob(os.path.join(PQ, "*.parquet")))
p("parquet 文件总数 = %d" % len(files))
p("")

rows = []
for f in files:
    b = os.path.basename(f)
    df = pd.read_parquet(f, columns=["af", "ma_count"])
    n = len(df)
    a = df["af"].astype(np.float64).values
    m = df["ma_count"].astype(np.float64).values
    del df
    ok = (a > 1e-6) & (a < 1 - 1e-6) & (m > 0)
    ne = m[ok] / a[ok] / 2.0
    ne_m = m[ok] / np.minimum(a[ok], 1 - a[ok]) / 2.0
    # 模型选择：只在 {2N} 网格上取值？取中位数与四分位
    med = float(np.median(ne)); q10 = float(np.quantile(ne, .1)); q90 = float(np.quantile(ne, .9))
    med_ma = float(np.median(ne_m))
    # 检验 N_eff 是否落在整数网格（用中位数的四舍五入值）
    n_round = int(round(med))
    on_grid = float(np.median(np.abs(ne - np.round(ne))))
    rows.append(dict(file=b, n_var=n, n_eff_af_med=med, q10=q10, q90=q90,
                     n_eff_ma_med=med_ma, n_round=n_round, grid_dev=on_grid))
    p("  %-58s n_var=%-8d N_af_med=%8.2f [%.1f, %.1f]  N_ma_med=%8.2f  grid_dev=%.4f"
      % (b, n, med, q10, q90, med_ma, on_grid))

out = pd.DataFrame(rows).sort_values(["n_eff_af_med", "file"])
out.to_csv(os.path.join(PROJ, "tables", "_s52b2_celltype_Neff.csv"), index=False, encoding="utf-8-sig")

p("")
p("=" * 78)
p("汇总：N_eff 取值分布（按 N_af_med 归并）")
p(out.groupby("n_eff_af_med")["file"].apply(lambda s: " ; ".join(s)).to_string())
p("")
p("★ 结论：")
uniq = sorted(out["n_eff_af_med"].unique())
p("  出现的 N_eff 值 = %s" % [round(x, 2) for x in uniq])
p("  LD 参考面板 plink_merged_980_donors 的 N = 980")
diff = out[np.abs(out["n_eff_af_med"] - 980) > 0.5]
p("  ★ N_eff != 980 的文件数 = %d" % len(diff))
for _, r in diff.iterrows():
    p("     %-58s N_eff=%.2f" % (r["file"], r["n_eff_af_med"]))

io.open(LOG, "w", encoding="utf-8", newline="\n").write("\n".join(L) + "\n")
print("WROTE " + LOG)
