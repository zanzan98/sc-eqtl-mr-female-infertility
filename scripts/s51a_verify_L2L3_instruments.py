# -*- coding: utf-8 -*-
"""s51a_verify_L2L3_instruments.py  --  任务A 步骤1/2 前置核验（修正版）

★坐标体系实证（由 s51a2_coord_forensics.py 标定）：
    parquet_subset 的 variant_id 与 _onek1k_plink bim SNP ID 一致，均为
        "<chr>:<pos GRCh37>"
    而 15_discovery_significant_with_locus.csv / 31_final_MR_with_method.csv 的
    `instrument` / `variant_id_grch38` 列为 "<chr>:<pos GRCh38>"。
    标定锚点（L1）：CDC42/CD4_NC  GRCh38 1:22088292 ↔ GRCh37 1:22414785（parquet 内命中后者）。

本脚本按 **GRCh37** 在 parquet 内定位 L2/L3 工具变量，核对
    slope / slope_se / pval_nominal / F=(slope/se)^2
与已交付 31_final_MR_with_method.csv 的 F_stat 及由 beta_outcome 反推的效应量是否自洽。
只读；不下载；不写任何已交付产物。
"""
import io
import os

import numpy as np
import pandas as pd

PROJ = r"D:\endometriosis_project\11_sc_eqtl_mr_project"
PQ = os.path.join(PROJ, "00_data_raw", "onek1k", "parquet_subset")
TAB = os.path.join(PROJ, "tables")
LOG = r"D:\endometriosis_project\_s51a_verify_log.txt"

L = []


def p(s=""):
    L.append(str(s))
    print(s)


# 期望值逐字符抄自 tables/31_final_MR_with_method.csv 第 15、16 行
TARGETS = [
    dict(locus="L2", gene="YME1L1", ct="CD4_NC", chrom="10",
         lead37="10:27443623", lead38="10:27154694",
         b_mr=0.1116160716615767, se_mr=0.026844557680235,
         F_exp=204.4232973160344, b_out=-0.0352, se_out=0.0081),
    dict(locus="L3", gene="ANXA4", ct="Mono_NC", chrom="2",
         lead37="2:69956721", lead38="2:69729589",
         b_mr=-0.0754702275441082, se_mr=0.0188839307962625,
         F_exp=96.49130239435812, b_out=0.0385, se_out=0.0088),
]

p("=== s51a L2/L3 工具变量核验（GRCh37 对齐版）===")
p("")
mr31 = pd.read_csv(os.path.join(TAB, "31_final_MR_with_method.csv"))
p("[0] 31_final_MR_with_method.csv 原始行（instrument 列为 GRCh38 字符串）")
p(mr31[mr31["locus_id"].isin(["L2", "L3"])].to_string(index=False))

results = []
for t in TARGETS:
    p("")
    p("=" * 78)
    p("[%s] %s / %s / chr%s" % (t["locus"], t["gene"], t["ct"], t["chrom"]))
    p("   GRCh37 工具 = %s   (31_ 表 instrument = %s, GRCh38)" % (t["lead37"], t["lead38"]))
    f = os.path.join(PQ, "OneK1K_%s.cis_qtl_pairs.chr%s.parquet" % (t["ct"], t["chrom"]))
    df = pd.read_parquet(f,
                         columns=["phenotype_id", "variant_id", "af", "slope",
                                  "slope_se", "pval_nominal"])
    e = df[df["phenotype_id"].astype(str) == t["gene"]].copy()
    e["vid"] = e["variant_id"].astype(str)
    e["pos37"] = e["vid"].str.split(":").str[1].astype(int)
    p("  parquet %s 变异数 = %d   pos37 范围 [%d, %d]"
      % (t["gene"], len(e), e["pos37"].min(), e["pos37"].max()))

    imin = int(e["pval_nominal"].idxmin())
    rm = e.loc[imin]
    p("  该基因最强 cis 变异 = %s  p=%.4e  slope=%+.6f se=%.6f  F=%.4f"
      % (rm["vid"], rm["pval_nominal"], rm["slope"], rm["slope_se"],
         (rm["slope"] / rm["slope_se"]) ** 2))

    q = e[e["vid"] == t["lead37"]]
    if len(q) == 0:
        p("  !! %s 不在 parquet 中" % t["lead37"])
        results.append(dict(locus=t["locus"], gene=t["gene"], cell_type=t["ct"],
                            instrument_grch37=t["lead37"], in_parquet=False))
        del df, e
        continue
    q = q.iloc[0]
    slope = float(q["slope"]); se = float(q["slope_se"]); pval = float(q["pval_nominal"])
    F = (slope / se) ** 2
    # 由 15_ 表的 beta_outcome/se_outcome 反推 MR 效应与（二阶 delta）SE
    b_mr_calc = t["b_out"] / slope
    se_mr_calc = np.sqrt(t["se_out"] ** 2 / slope ** 2 + t["b_out"] ** 2 * se ** 2 / slope ** 4)
    dF = abs(F - t["F_exp"])
    d_mr = abs(b_mr_calc - t["b_mr"])
    d_se = abs(se_mr_calc - t["se_mr"])
    # ★parquet 的 slope/slope_se 是 float32 → F 的绝对误差可达 ~1e-5，必须用相对容差。
    relF = dF / t["F_exp"]
    ok = (relF < 1e-5) and (d_mr < 5e-5) and (d_se < 5e-5)
    p("  --- 比对 ---")
    p("    slope(暴露)   = %+.15f   se = %.15f   p = %.6e" % (slope, se, pval))
    p("    F=(slope/se)^2 = %.6f       31_ 表 F_stat = %.6f    |diff|=%.3e  相对=%.3e"
      % (F, t["F_exp"], dF, relF))
    p("    MR b  = beta_out/slope  = %+.15f   31_ 表 b = %+.15f   |diff|=%.3e"
      % (b_mr_calc, t["b_mr"], d_mr))
    p("    MR se = 二阶 delta      = %.15f   31_ 表 se = %.15f   |diff|=%.3e"
      % (se_mr_calc, t["se_mr"], d_se))
    p("    af = %.8f   →  ==> %s" % (q["af"], "一致 (PASS)" if ok else "!! 不一致 (FAIL)"))
    results.append(dict(locus=t["locus"], gene=t["gene"], cell_type=t["ct"],
                        instrument_grch37=t["lead37"], instrument_grch38=t["lead38"],
                        in_parquet=True, slope=slope, slope_se=se, pval=pval, af=float(q["af"]),
                        F_calc=F, F_31=t["F_exp"], rel_F=relF, b_MR_calc=b_mr_calc, b_MR_31=t["b_mr"],
                        se_MR_calc=se_mr_calc, se_MR_31=t["se_mr"], d_F=dF, d_b=d_mr, d_se=d_se,
                        pass_=ok))
    del df, e

p("")
p("=" * 78)
R = pd.DataFrame(results)
R.to_csv(os.path.join(TAB, "_s51a_verify_L2L3.csv"), index=False, encoding="utf-8-sig")
p("汇总：")
p(R.to_string(index=False))
p("")
p("全部 PASS = %s" % bool(R["pass_"].all() if "pass_" in R.columns else False))
io.open(LOG, "w", encoding="utf-8", newline="\n").write("\n".join(L) + "\n")
print("WROTE " + LOG)
