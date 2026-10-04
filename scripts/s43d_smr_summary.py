# -*- coding: utf-8 -*-
"""s43d_smr_summary.py — 任务三 3.3 步骤 4：SMR/HEIDI 结果汇总与判读

输入：tables/45_smr_3p3_results.csv（主分析 5e-8）、tables/45b_smr_3p3_results_sens.csv（敏感性 1e-5）
输出：
  tables/46_smr_3p3_gene_summary.csv   基因层面汇总（HEIDI 通过/否决计数、P_SMR 范围、方向一致性）
  tables/45c_smr_3p3_combined.csv      重写：加 heidi_verdict / fdr_smr 列
  日志 scripts/_s43d_smr_summary.log
"""
import os, io, sys, csv
import pandas as pd
import numpy as np

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
PROJ = r"D:\endometriosis_project\11_sc_eqtl_mr_project"
TAB = os.path.join(PROJ, "tables")
LOG = os.path.join(PROJ, "scripts", "_s43d_smr_summary.log")
buf = []
def w(s=""):
    buf.append(str(s)); print(s)

def bh(p):
    p = np.asarray(p, float)
    n = len(p)
    o = np.argsort(p)
    q = np.empty(n)
    q[o] = p[o] * n / (np.arange(n) + 1.0)
    q = np.minimum.accumulate(q[o][::-1])[::-1]
    out = np.empty(n)
    out[o] = np.minimum(q, 1.0)
    return out

def load(fn, tag):
    d = pd.read_csv(os.path.join(TAB, fn))
    for c in ("b_SMR", "se_SMR", "p_SMR", "p_HEIDI", "b_eQTL", "se_eQTL", "p_eQTL",
              "b_GWAS", "se_GWAS", "p_GWAS", "Freq"):
        d[c] = pd.to_numeric(d[c], errors="coerce")
    d["nsnp_HEIDI"] = pd.to_numeric(d["nsnp_HEIDI"], errors="coerce")
    d["heidi_verdict"] = np.where(d["p_HEIDI"].isna(), "NA",
                                  np.where(d["p_HEIDI"] > 0.05, "consistent", "heterogeneity"))
    d["fdr_smr"] = bh(d["p_SMR"].values)
    d["direction"] = np.where(d["b_SMR"] > 0, "positive", "negative")
    d["analysis"] = tag
    return d

w("=" * 92)
w("s43d_smr_summary.py — SMR/HEIDI 汇总（任务三 3.3）")
w("=" * 92)

m = load("45_smr_3p3_results.csv", "main(5e-8)")
s = load("45b_smr_3p3_results_sens.csv", "sens(1e-5)")
w()
w("主分析行数 = %d   敏感性行数 = %d" % (len(m), len(s)))

# ---- 汇总
rows = []
for tag, d in (("main(5e-8)", m), ("sens(1e-5)", s)):
    for g in ("CDC42", "LINC00339"):
        dd = d[d["Gene"] == g]
        if len(dd) == 0:
            continue
        rows.append(dict(
            analysis=tag, gene=g,
            n_cells_tested=len(dd),
            n_smr_p05=int((dd["p_SMR"] < 0.05).sum()),
            n_smr_fdr05=int((dd["fdr_smr"] < 0.05).sum()),
            min_p_SMR=dd["p_SMR"].min(),
            argmin_pSMR_cell=dd.loc[dd["p_SMR"].idxmin(), "cell_type"],
            n_heidi_consistent=int((dd["heidi_verdict"] == "consistent").sum()),
            n_heidi_heterogeneity=int((dd["heidi_verdict"] == "heterogeneity").sum()),
            n_heidi_NA=int((dd["heidi_verdict"] == "NA").sum()),
            b_SMR_sign_positive=int((dd["direction"] == "positive").sum()),
            b_SMR_sign_negative=int((dd["direction"] == "negative").sum()),
            b_SMR_min=dd["b_SMR"].min(), b_SMR_max=dd["b_SMR"].max(),
            min_p_HEIDI=dd["p_HEIDI"].min(), median_p_HEIDI=dd["p_HEIDI"].median(),
            heidi_fail_cells=";".join(sorted(dd.loc[dd["heidi_verdict"] == "heterogeneity", "cell_type"])),
        ))
S = pd.DataFrame(rows)
S.to_csv(os.path.join(TAB, "46_smr_3p3_gene_summary.csv"), index=False, encoding="utf-8")
w()
w("基因层面汇总（tables/46_smr_3p3_gene_summary.csv）：")
w("%-11s %-10s %5s %5s %7s %-12s %6s %6s %5s %-6s %-6s"
  % ("analysis", "gene", "nCell", "SMR<.05", "FDR<.05", "argminP_cell", "HEpass", "HEfail", "HEna", "b>0", "b<0"))
for _, r in S.iterrows():
    w("%-11s %-10s %5d %5d %7d %-12s %6d %6d %5d %-6d %-6d"
      % (r["analysis"], r["gene"], r["n_cells_tested"], r["n_smr_p05"], r["n_smr_fdr05"],
         r["argmin_pSMR_cell"], r["n_heidi_consistent"], r["n_heidi_heterogeneity"],
         r["n_heidi_NA"], r["b_SMR_sign_positive"], r["b_SMR_sign_negative"]))
w()
for _, r in S.iterrows():
    w("  %s / %s: b_SMR ∈ [%.4f, %.4f] ; P_SMR min=%.3e ; P_HEIDI min=%.3e median=%.4g"
      % (r["analysis"], r["gene"], r["b_SMR_min"], r["b_SMR_max"],
         r["min_p_SMR"], r["min_p_HEIDI"], r["median_p_HEIDI"]))
    if r["heidi_fail_cells"]:
        w("      HEIDI 否决(P<=0.05) 细胞: %s" % r["heidi_fail_cells"])

# ---- 重写合并表（含 heidi_verdict / fdr_smr）
comb = pd.concat([m, s], ignore_index=True)
cols = ["analysis", "cell_type", "Gene", "probeID", "Probe_bp", "topSNP", "topSNP_bp",
        "A1", "A2", "Freq", "b_GWAS", "se_GWAS", "p_GWAS", "b_eQTL", "se_eQTL", "p_eQTL",
        "b_SMR", "se_SMR", "p_SMR", "fdr_smr", "direction", "p_HEIDI", "nsnp_HEIDI", "heidi_verdict"]
comb[cols].to_csv(os.path.join(TAB, "45c_smr_3p3_combined.csv"), index=False, encoding="utf-8")
w()
w("重写 tables/45c_smr_3p3_combined.csv（%d 行，含 heidi_verdict / fdr_smr / direction）" % len(comb))

# ---- 一致性检查：主分析与敏感性在同 (cell,gene) 上是否逐位一致
key = ["cell_type", "Gene"]
mm = m.set_index(key); ss = s.set_index(key)
common = mm.index.intersection(ss.index)
w()
w("主分析 ∩ 敏感性 的交集 = %d 对；在交集上 max|Δb_SMR| = %.3e  max|Δp_SMR| = %.3e  max|Δp_HEIDI| = %.3e"
  % (len(common),
     np.nanmax(np.abs(mm.loc[common, "b_SMR"].values - ss.loc[common, "b_SMR"].values)),
     np.nanmax(np.abs(mm.loc[common, "p_SMR"].values - ss.loc[common, "p_SMR"].values)),
     np.nanmax(np.abs(mm.loc[common, "p_HEIDI"].values - ss.loc[common, "p_HEIDI"].values))))
w("   → 交集内逐位一致即证明放宽 --peqtl-smr 不改变 target SNP 选择（仅在无 5e-8 命中时新增对）")

# ---- 逐行明细（主分析）
w()
w("主分析逐行（按 P_SMR 升序）")
w("%-12s %-10s %-14s %9s %11s %11s %5s %-14s" %
  ("cell", "Gene", "topSNP", "b_SMR", "p_SMR", "p_HEIDI", "m", "heidi_verdict"))
for _, r in m.sort_values("p_SMR").iterrows():
    w("%-12s %-10s %-14s %9.4f %11.3e %11.3e %5d %-14s"
      % (r["cell_type"], r["Gene"], r["topSNP"], r["b_SMR"], r["p_SMR"], r["p_HEIDI"],
         int(r["nsnp_HEIDI"]), r["heidi_verdict"]))

with io.open(LOG, "w", encoding="utf-8") as f:
    f.write("\n".join(buf) + "\n")
print("[wrote] %s" % LOG)
