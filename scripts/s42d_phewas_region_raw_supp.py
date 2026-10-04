# -*- coding: utf-8 -*-
"""s42d_phewas_region_raw_supp.py -- 任务 3.2 (b) 的**假设自由**佐证

目的
----
(b) 区域臂用 rs56318008 作工具、以 CDC42 的 eQTL 效应标定（因 WNT4 在血液免疫细胞
**无 cis-eQTL**，无法用 WNT4 转录本标定）。为免「标定本身引入 CDC42 假设」，
这里再做一层**完全不依赖任何 eQTL** 的检查：

  把 rs56318008 当作一个普通变异，直接对 2,469 个 FinnGen R12 端点做变异级 PheWAS
  （BH，m = 2,469），看同样的生殖/产科/免疫血液信号是否仍成立。

★ 这是**附加**的辅助检查（不在 (a)/(b)/(c) 原始规格内），仅在汇报中作为 (b) 的
  稳健性支撑，不单独构成结论。

产物：tables/44f_phewas_3p2_region_variant_raw.csv
"""
import io
import os
import sys
import time

import numpy as np
import pandas as pd

PROJ = r"D:\endometriosis_project\11_sc_eqtl_mr_project"
TABX = os.path.join(PROJ, "00_data_raw", "finngen_r12_tabix")
TAB = os.path.join(PROJ, "tables")
MANIFEST = r"D:\endometriosis_project\01_data_raw\finngen_R12_manifest.tsv"
LOG = r"D:\endometriosis_project\_s42d_phewas_region_raw.log"
OUT = os.path.join(TAB, "44f_phewas_3p2_region_variant_raw.csv")

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass
L = []


def p(s=""):
    L.append(str(s))
    print(s)


def bh(pv):
    pv = np.asarray(pv, float)
    ok = np.isfinite(pv)
    q = np.full(pv.shape, np.nan)
    v = pv[ok]
    n = v.size
    o = np.argsort(v)
    r = v[o]
    qq = np.minimum.accumulate((r * n / np.arange(1, n + 1))[::-1])[::-1]
    q[ok] = np.clip(qq, 0, 1)[np.argsort(o)]
    return q


DOM = {
    "免疫/血液/肿瘤": [
        "I Certain infectious and parasitic diseases (AB1_)",
        "III Diseases of the blood and blood-forming organs and certain disorders involving the immune mechanism (D3_)",
        "II Neoplasms from hospital discharges (CD2_)",
        "II Neoplasms, from cancer register (ICD-O-3)",
        "Diseases marked as autimmune origin"],
    "女性生殖/不孕": ["XIV Diseases of the genitourinary system (N14_)"],
    "妊娠/分娩/产褥": ["XV Pregnancy, childbirth and the puerperium (O15_)",
                       "XVI Certain conditions originating in the perinatal period (P16_)"],
    "肌肉骨骼/结缔组织": ["XIII Diseases of the musculoskeletal system and connective tissue (M13_)",
                          "Rheuma endpoints"],
}
C2D = {c: d for d, cs in DOM.items() for c in cs}

p("=== s42d 区域变异（rs56318008）变异级原始 PheWAS %s ===" % time.strftime("%Y-%m-%d %H:%M:%S"))

lg = pd.read_csv(os.path.join(TABX, "phewas_3snp_all_endpoints.csv"))
lg = lg[lg["pos38"].astype(int) == 22143914].copy()
lg = lg.rename(columns={"endpoint": "phenocode", "pval": "pval_raw", "beta": "beta_raw",
                        "sebeta": "se_raw", "af_alt": "af_raw"})
p("端点 = %d（rs56318008 每端点 1 行）" % len(lg))

man = pd.read_csv(MANIFEST, sep="\t")[["phenocode", "phenotype", "category", "num_cases", "num_controls"]]
d = man.merge(lg[["phenocode", "pval_raw", "beta_raw", "se_raw", "af_raw"]], on="phenocode", how="inner")
d["domain"] = np.where(d["category"].isin(C2D), d["category"].map(C2D), "其他")
d["p_FDR_raw"] = bh(d["pval_raw"].values)
d["risk_if_cdc42_inhibited"] = np.where(d["beta_raw"] < 0, "increased", "decreased")
d.sort_values("p_FDR_raw").to_csv(OUT, index=False, encoding="utf-8-sig")

sig = d[d["p_FDR_raw"] < 0.05]
p("")
p("[总体] 检验数 m = %d ；FDR<0.05 = %d（方向↑ %d / ↓ %d）；min p = %.3e"
  % (len(d), len(sig), int((sig["beta_raw"] < 0).sum()), int((sig["beta_raw"] > 0).sum()),
     float(d["pval_raw"].min())))
p("")
p("[按域]")
for dm, g in d.groupby("domain"):
    s = g[g["p_FDR_raw"] < 0.05]
    p("   %-14s n=%4d  FDR<0.05=%3d（↑ %3d / ↓ %3d）  min_p=%.3e"
      % (dm, len(g), len(s), int((s["beta_raw"] < 0).sum()), int((s["beta_raw"] > 0).sum()),
         float(g["pval_raw"].min())))

# 与 (a)/(c) 的 12 个风险升高端点对照
old = pd.read_csv(os.path.join(TAB, "42_phewas_safety_assessment.csv"))
inc = old[old["risk_if_cdc42_inhibited"] == "increased"][["endpoint", "cell_type"]]
eps = sorted(set(inc["endpoint"]))
p("")
p("[对照] 42_ 的 12 个「风险升高」行涉及 %d 个端点" % len(eps))
p("   %-46s %12s %12s %12s %6s" % ("endpoint", "beta_raw", "se_raw", "p_raw", "q_raw"))
nq = 0
for e in eps:
    r = d[d["phenocode"] == e]
    if not len(r):
        p("   %-46s  !!! 缺失" % e)
        continue
    r = r.iloc[0]
    if r["p_FDR_raw"] < 0.05:
        nq += 1
    p("   %-46s %+12.4f %12.5f %12.3e %6s"
      % (e, r["beta_raw"], r["se_raw"], r["pval_raw"],
         "%.2e" % r["p_FDR_raw"]))
p("   → 12 端点中 %d/%d 在**变异级** BH<0.05（Mono_NC 与 B_MEM 共用同一边异，故按端点计 %d 个）"
  % (nq, len(eps), len(eps)))

p("")
p("产物：%s" % OUT)
p("VERDICT = DONE")
io.open(LOG, "w", encoding="utf-8", newline="\n").write("\n".join(L) + "\n")
print("WROTE " + LOG)
