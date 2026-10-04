# -*- coding: utf-8 -*-
"""s25_prep_ivw_inputs.py —— 生成 s06 / s10 可直接消费的工具表

产出三张表（列格式一致，超集）：
  tables/29a_ivw_instruments_main.csv      主口径 r²<0.001（26 行，全 nsnp=1）
  tables/29b_ivw_instruments_sens001.csv   敏感性 r²<0.01（LINC00339/CD4_NC 有 2 个）
  tables/29c_ivw_instruments_unclumped.csv 未 clump（负对照，2,842 行）
列：gene, cell_type, chrom, variant_id(GRCh38), variant_id_grch37, beta, se, af, F,
    pval_nominal, qval(NaN), tss_distance, num_var
"""
import os
import numpy as np
import pandas as pd

PROJ = r"D:\endometriosis_project\11_sc_eqtl_mr_project"
T = os.path.join(PROJ, "tables")
TOT = pd.read_csv(os.path.join(T, "27b_raw_cis_totals.csv"), encoding="utf-8-sig")
PRE = pd.read_csv(os.path.join(T, "27_raw_instruments_preclump.csv"), encoding="utf-8-sig")

JOBS = [("29a_ivw_instruments_main.csv", "28_raw_instruments_clumped_main.csv"),
        ("29b_ivw_instruments_sens001.csv", "28_raw_instruments_clumped_sens001.csv"),
        ("29c_ivw_instruments_unclumped.csv", None)]

for out_name, clump_name in JOBS:
    if clump_name is None:
        d = PRE.copy()
        d["clump_index"] = 1
        d["n_index"] = d.groupby(["gene", "cell_type"])["variant_id_grch37"].transform("size")
    else:
        d = pd.read_csv(os.path.join(T, clump_name), encoding="utf-8-sig")
    d = d.merge(TOT[["gene", "cell_type", "n_cis_variants"]], on=["gene", "cell_type"], how="left")
    d = d.rename(columns={"F_approx": "F", "n_cis_variants": "num_var",
                          "variant_id_grch38": "variant_id"})
    d["qval"] = np.nan
    cols = ["gene", "cell_type", "chrom", "variant_id", "variant_id_grch37",
            "beta", "se", "af", "F", "pval_nominal", "qval", "tss_distance", "num_var",
            "clump_index", "n_index"]
    d = d[[c for c in cols if c in d.columns]]
    d = d.sort_values(["chrom", "gene", "cell_type", "pval_nominal"]).reset_index(drop=True)
    d.to_csv(os.path.join(T, out_name), index=False, encoding="utf-8-sig")
    ns = d.groupby(["gene", "cell_type"]).size()
    print("%-34s 行=%5d  组合=%2d  nsnp=1 的组合=%d  nsnp>=2 的组合=%d"
          % (out_name, len(d), len(ns), int((ns == 1).sum()), int((ns >= 2).sum())))
    if (ns >= 2).any():
        print("   多 SNP 组合：%s" % ns[ns >= 2].to_dict())

print()
print("=== 未 clump 的 nsnp 分布（前 10）===")
d3 = pd.read_csv(os.path.join(T, "29c_ivw_instruments_unclumped.csv"), encoding="utf-8-sig")
print(d3.groupby(["gene", "cell_type"]).size().sort_values(ascending=False).head(10).to_string())
