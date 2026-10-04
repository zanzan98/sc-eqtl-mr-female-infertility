# -*- coding: utf-8 -*-
"""s40a_recon_3p1.py -- 任务 3.1 条件共定位：输入盘点

只读，不改动任何已闭环产物。
输出 -> D:\\endometriosis_project\\_s40a_recon_3p1.log (UTF-8)
"""
import io
import os

import numpy as np
import pandas as pd

PROJ = r"D:\endometriosis_project\11_sc_eqtl_mr_project"
PREP = os.path.join(PROJ, "00_data_raw", "onek1k", "coloc_prep")
TAB = os.path.join(PROJ, "tables")
BIM_REF = r"D:\endometriosis_project\_onek1k_plink\plink_merged_980_donors.bim"
LOG = r"D:\endometriosis_project\_s40a_recon_3p1.log"

TARGET_RS = "rs56318008"
TARGET_POS37 = 22470407

L = []


def p(s=""):
    L.append(str(s))
    print(s)


p("=== s40a 任务3.1 条件共定位：输入盘点 ===")

# ---------- 1) 精细定位表里的 rs56318008 ----------
fp = os.path.join(TAB, "35_finemap_variant.csv")
fm = pd.read_csv(fp)
p("")
p("[1] tables/35_finemap_variant.csv")
p("   shape = %s" % (fm.shape,))
p("   cols  = %s" % list(fm.columns))
sub = fm[fm.apply(lambda r: TARGET_RS in r.astype(str).values, axis=1)]
p("   含 %s 的行 = %d" % (TARGET_RS, len(sub)))
if len(sub):
    for _, r in sub.head(4).iterrows():
        p("      " + " | ".join("%s=%s" % (c, r[c]) for c in fm.columns))

# ---------- 2) 窗口定义 ----------
p("")
p("[2] coloc_prep/window.json")
p("   " + io.open(os.path.join(PREP, "window.json"), encoding="utf-8").read().replace("\n", " "))

# ---------- 3) eQTL 侧 ----------
p("")
p("[3] eQTL 侧：eqtl_L1_CDC42.csv")
e = pd.read_csv(os.path.join(PREP, "eqtl_L1_CDC42.csv"))
p("   shape = %s" % (e.shape,))
p("   cols  = %s" % list(e.columns))
p("   细胞类型 = %s" % sorted(e["cell_type"].unique().tolist()))
for ct in ["B_MEM", "Mono_NC", "B_IN", "Mono_C"]:
    g = e[e["cell_type"] == ct]
    hit = g[g["snp_grch37"] == "1:%d" % TARGET_POS37]
    p("   %-8s rows=%5d  cis唯一位置=%5d  %s@1:%d -> %s"
      % (ct, len(g), g["snp_grch37"].nunique(), TARGET_RS, TARGET_POS37,
         ("存在" if len(hit) else "不存在")))

# ---------- 4) 结局侧 ----------
p("")
p("[4] 结局侧：gwas_L1.csv")
g = pd.read_csv(os.path.join(PREP, "gwas_L1.csv"))
p("   shape = %s" % (g.shape,))
p("   cols  = %s" % list(g.columns))
h = g[g["rsid"].astype(str) == TARGET_RS]
p("   %s -> %s" % (TARGET_RS, "存在" if len(h) else "不存在"))
if len(h):
    p("      " + " | ".join("%s=%s" % (c, h.iloc[0][c]) for c in g.columns))
p("   pos38 范围 = %s .. %s" % (g["pos38"].min(), g["pos38"].max()))

# ---------- 5) LD 参考 bim ----------
p("")
p("[5] OneK1K LD 参考 bim")
bim = pd.read_csv(BIM_REF, sep="\t", header=None, dtype={0: str},
                  names=["chr", "snp", "cm", "pos", "a1", "a2"])
p("   总变异 = %d" % len(bim))
b1 = bim[(bim["chr"].astype(str) == "1") & (bim["pos"] >= 22000000) & (bim["pos"] <= 22900000)]
p("   chr1:22.0-22.9Mb(GRCh37) 变异 = %d" % len(b1))
t = bim[bim["pos"] == TARGET_POS37]
p("   pos=%d -> %s" % (TARGET_POS37, ("存在 " + str(t.iloc[0]["snp"]) +
                                     " A1=%s A2=%s" % (t.iloc[0]["a1"], t.iloc[0]["a2"])) if len(t) else "不存在"))

# ---------- 6) 复现 s28 的 harmonise，检查 rs56318008 是否在交集内 ----------
p("")
p("[6] 复现 s28 harmonise（allele-aware key）并查 rs56318008")
b = pd.read_csv(os.path.join(PREP, "bim_universe_L1.csv"))
p("   bim_universe_L1 rows = %d" % len(b))
gw = g.copy()
gw["ea"] = gw["ea"].astype(str).str.upper()
gw["oa"] = gw["oa"].astype(str).str.upper()
gw = gw[gw["pos38"].notna() & gw["beta"].notna() & gw["se"].notna() & (gw["se"] > 0)]
gw["allele_key"] = (gw["pos38"].astype(int).astype(str) + "_" +
                    np.minimum(gw["ea"], gw["oa"]) + "_" + np.maximum(gw["ea"], gw["oa"]))
gw = gw[~gw["allele_key"].duplicated()]
for ct in ["B_MEM", "Mono_NC"]:
    ec = e[e["cell_type"] == ct].copy()
    ec = ec.merge(b, left_on="snp_grch37", right_on="snp", how="left")
    ec = ec[ec["a1"].notna() & ec["a2"].notna()]
    ec["pos38"] = ec["snp_grch38"].astype(str).str.split(":").str[1]
    ec = ec[ec["pos38"].notna()]
    ec["pos38"] = ec["pos38"].astype(int)
    ec["a1"] = ec["a1"].astype(str).str.upper()
    ec["a2"] = ec["a2"].astype(str).str.upper()
    ec["allele_key"] = (ec["pos38"].astype(str) + "_" +
                        np.minimum(ec["a1"], ec["a2"]) + "_" + np.maximum(ec["a1"], ec["a2"]))
    m = ec.merge(gw, on="allele_key", suffixes=("_e", "_o"))
    p("   %-8s eqtl_ct=%5d  gwas_win=%5d  交集=%5d" % (ct, len(ec), len(gw), len(m)))
    hh = m[m["rsid"].astype(str) == TARGET_RS]
    p("        %s 在交集内 -> %s" % (TARGET_RS, "是" if len(hh) else "否"))
    if len(hh):
        r0 = hh.iloc[0]
        p("        eQTL: slope=%.6f slope_se=%.6f af=%.4f pval=%.3g"
          % (r0["slope"], r0["slope_se"], r0["af"], r0["pval_nominal"]))
        p("        GWAS: beta=%.6f se=%.6f eaf=%.4f p=%.3g n_cases=%s n_controls=%s"
          % (r0["beta"], r0["se"], r0["eaf"], r0["p"],
             r0["n_cases"], r0["n_controls"]))
        p("        pos38=%s  a1=%s a2=%s  ea=%s oa=%s"
          % (r0["pos38_o"], r0["a1"], r0["a2"], r0["ea"], r0["oa"]))

# ---------- 7) 已有条件/精细定位产物 ----------
p("")
p("[7] 已有 36_finemap_pair.csv（chr1）")
fr = pd.read_csv(os.path.join(TAB, "36_finemap_pair.csv"))
p("   cols = %s" % list(fr.columns))
c1 = fr[fr["locus_id"].astype(str).str.startswith("L1")] if "locus_id" in fr.columns else fr
p("   L1 行数 = %d" % len(c1))
p(c1.head(15).to_string())

io.open(LOG, "w", encoding="utf-8", newline="\n").write("\n".join(L) + "\n")
print("WROTE " + LOG)
