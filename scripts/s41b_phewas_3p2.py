# -*- coding: utf-8 -*-
"""s41b_phewas_3p2.py -- 任务 3.2：条件 PheWAS（a）与 区域工具 PheWAS（b）

设计
----
暴露侧效应（OneK1K，980 供者，GRCh37；已对齐至 bim A1）：
  CDC42 / B_MEM   工具 rs12037376 (1:22462111)  β=0.3595215 se=0.04347525  F=31.13
  CDC42 / Mono_NC 工具 rs10917151 (1:22422721)  β=0.3296235 se=0.05788484  F=26.12
  条件变异 rs56318008 (1:22470407，WNT4 可信集 lead)
       B_MEM  slope=0.344998 se=0.044752 (z=+7.7091)
       Mono_NC slope=0.317485 se=0.060542 (z=+5.2440)
结局侧效应：FinnGen R12 PheWeb `/api/variant/...`（2,469 端点；已与 42_…csv 的 beta_out 逐位核对）

手臂（arm）
----------
(a) 条件 PheWAS：
    uncond          工具边际 β_exp × 结局边际 β_out（= 复现原始 4,938 检验的回归锚点）
    cond_eqtl_only  工具条件 β_exp|k × 结局边际 β_out
    cond_both       工具条件 β_exp|k × 结局条件 β_out|k   ← 教科书式（两侧同 r 条件）
(b) 区域工具 PheWAS（不含任何条件）：
    region          条件变异 rs56318008 自身作为区域工具，分别用 B_MEM / Mono_NC 的 eQTL 效应

条件公式（GCTA-COJO 型；外参 LD = OneK1K 980 欧洲血统，与主结局 EUR 血统一致）
  β_cond = β_j − r_jk·(se_j/se_k)·β_k ;  var_cond = se_j²·(1−r_jk²)
  z_cond = (z_j − r_jk·z_k)/sqrt(1−r_jk²)                    （两条路径互相核验）
  Wald: β_MR = β_o/β_e ; se_MR = sqrt(var_o/β_e² + β_o²·var_e/β_e⁴)   ← 二阶 delta（硬纪律）
  F = (β_e_c / se_e_c)²

FDR：按「家族（family）」把两个细胞类型合并后做 BH（与原始 42_ 口径一致：4,938 检验）
     同时给出按 arm 单独的 FDR 供参考。

域计数：**必须用 manifest 的 `category` 字段**（CSV 的 category_prefix 是 phenocode 混合编码）
       并断言复现 42b 的 484/144/111/273/1457（×2 = 968/288/222/546/2914）。

红线：不改动任何已闭环产物；只新增 tables/44*；结果写文件再读。
"""
import io
import json
import os
import subprocess
import sys
import time

import numpy as np
import pandas as pd
from scipy import stats

PROJ = r"D:\endometriosis_project\11_sc_eqtl_mr_project"
COND = os.path.join(PROJ, "00_data_raw", "onek1k", "cond_coloc")
PHW = os.path.join(PROJ, "00_data_raw", "finngen_r12_pheweb")
TAB = os.path.join(PROJ, "tables")
MANIFEST = r"D:\endometriosis_project\01_data_raw\finngen_R12_manifest.tsv"
LOG = r"D:\endometriosis_project\_s41b_phewas_3p2.log"

OUT_ALL = os.path.join(TAB, "44_phewas_3p2_arms.csv")
OUT_SUM = os.path.join(TAB, "44b_phewas_3p2_summary.csv")
OUT_DOM = os.path.join(TAB, "44c_phewas_3p2_domains.csv")
OUT_R12 = os.path.join(TAB, "44d_phewas_3p2_risk12_tracking.csv")

BONF_NOTE = "BH-FDR"
L = []


def p(s=""):
    L.append(str(s))
    print(s)


def free_gb():
    try:
        out = subprocess.check_output(
            ["powershell", "-NoProfile", "-Command",
             "(Get-CimInstance Win32_OperatingSystem).FreePhysicalMemory"],
            text=True, stderr=subprocess.DEVNULL, timeout=60).strip()
        return int(out) / 1024.0 / 1024.0
    except Exception:
        return float("nan")


def bh(pv):
    pv = np.asarray(pv, dtype=float)
    ok = np.isfinite(pv)
    q = np.full(pv.shape, np.nan)
    if ok.sum() == 0:
        return q
    v = pv[ok]
    n = v.size
    order = np.argsort(v)
    ranked = v[order]
    qq = ranked * n / (np.arange(1, n + 1))
    qq = np.minimum.accumulate(qq[::-1])[::-1]
    qq = np.clip(qq, 0, 1)
    tmp = np.empty_like(qq)
    tmp[order] = qq
    q[ok] = tmp
    return q


def wald(be, vare, bo, varo):
    """二阶 delta 的 Wald ratio。返回 (beta, se, p, F)"""
    be = np.asarray(be, float); vare = np.asarray(vare, float)
    bo = np.asarray(bo, float); varo = np.asarray(varo, float)
    with np.errstate(divide="ignore", invalid="ignore"):
        beta = bo / be
        se = np.sqrt(varo / be ** 2 + bo ** 2 * vare / be ** 4)
        z = beta / se
        pv = 2 * stats.norm.sf(np.abs(z))
        F = be ** 2 / vare
    bad = ~(np.isfinite(beta) & np.isfinite(se) & (se > 0) & np.isfinite(pv))
    beta[bad] = np.nan; se[bad] = np.nan; pv[bad] = np.nan; F[bad] = np.nan
    return beta, se, pv, F


T0 = time.time()
p("=== s41b 任务 3.2 条件 PheWAS + 区域工具 PheWAS ===")
p("time %s   可用内存 %.2f GB" % (time.strftime("%Y-%m-%d %H:%M:%S"), free_gb()))

# ------------------------------------------------------------------ 0) 输入
manifest = pd.read_csv(MANIFEST, sep="\t")
manifest = manifest[["phenocode", "phenotype", "category", "num_cases", "num_controls"]].copy()
p("")
p("[0] manifest 端点 = %d   唯一 phenocode = %d" % (len(manifest), manifest["phenocode"].nunique()))

var = {}
for rs in ["rs12037376", "rs10917151", "rs56318008"]:
    d = pd.read_csv(os.path.join(PHW, "var_%s.tsv" % rs), sep="\t")
    d = d[["phenocode", "phenostring", "beta", "sebeta", "pval", "maf", "n_case", "n_control"]].copy()
    d = d.rename(columns={"beta": "b_" + rs, "sebeta": "se_" + rs, "pval": "p_" + rs,
                          "maf": "maf_" + rs})
    # 去掉与 manifest 不匹配的行（API 多返回 1 条）
    d = d[d["phenocode"].isin(set(manifest["phenocode"]))]
    d = d[~d["phenocode"].duplicated()]
    var[rs] = d
    p("   %-12s 行=%d  maf=%.4f  n_case 中位=%s" % (rs, len(d), d["maf_" + rs].iloc[0],
                                                    int(d["n_case"].median())))

df = manifest.merge(var["rs12037376"], on="phenocode", how="inner")
df = df.merge(var["rs10917151"], on="phenocode", how="inner")
df = df.merge(var["rs56318008"], on="phenocode", how="inner")
df = df.drop(columns=[c for c in ["phenostring"] if c in df.columns])
p("   合并后端点 = %d" % len(df))

# ------------------------------------------------------------------ 1) 等位对齐自校验
p("")
p("[1] 等位对齐自校验（eQTL af(效应等位) vs FinnGen maf）")
mB = pd.read_csv(os.path.join(COND, "m_B_MEM.csv"))
ldtab = pd.read_csv(os.path.join(COND, "r_to_rs56318008.csv"))
ldtab["snp_grch37"] = ldtab["snp_grch37"].astype(str)
POS37 = {"rs12037376": "1:22462111", "rs10917151": "1:22422721", "rs56318008": "1:22470407"}
exp = {}
for rs, sid in POS37.items():
    qm = mB[mB["snp_grch37"] == sid]
    qld = ldtab[ldtab["snp_grch37"] == sid]
    exp[rs] = dict(sid=sid, a1=qm.iloc[0]["a1"], a2=qm.iloc[0]["a2"], af=float(qm.iloc[0]["af"]))
    maf_api = float(df["maf_" + rs].iloc[0])
    d_align = abs(exp[rs]["af"] - maf_api)
    d_flip = abs(exp[rs]["af"] - (1 - maf_api))
    ok = d_align < 0.05 and d_align < d_flip
    p("   %-12s bim A1=%s  eQTL af=%.4f  FinnGen maf=%.4f  |对齐|=%.4f  |若翻转|=%.4f  -> %s"
      % (rs, exp[rs]["a1"], exp[rs]["af"], maf_api, d_align, d_flip, "ALIGNED" if ok else "!! 需复核"))
    if not ok:
        p("   !! 等位对齐判定失败，终止")
        io.open(LOG, "w", encoding="utf-8", newline="\n").write("\n".join(L) + "\n")
        sys.exit(5)

# ------------------------------------------------------------------ 2) 暴露参数与 r
p("")
p("[2] 暴露参数与 LD（r vs rs56318008）")
mM = pd.read_csv(os.path.join(COND, "m_Mono_NC.csv"))
mM["snp_grch37"] = mM["snp_grch37"].astype(str)
SID = {"rs12037376": "1:22462111", "rs10917151": "1:22422721", "rs56318008": "1:22470407"}


def etab(m):
    o = {}
    for rs, sid in SID.items():
        q = m[m["snp_grch37"] == sid]
        o[rs] = dict(slope=float(q.iloc[0]["slope"]), se=float(q.iloc[0]["slope_se"]),
                     z=float(q.iloc[0]["slope"]) / float(q.iloc[0]["slope_se"]))
    return o


E = {"B_MEM": etab(mB), "Mono_NC": etab(mM)}
R = {}
for rs in ["rs12037376", "rs10917151"]:
    q = ldtab[ldtab["snp_grch37"] == SID[rs]]
    R[rs] = dict(r=float(q.iloc[0]["r_to_cond"]), r2=float(q.iloc[0]["r2_to_cond"]))
    p("   %-12s r=%.6f  r2=%.6f" % (rs, R[rs]["r"], R[rs]["r2"]))
for ct in ["B_MEM", "Mono_NC"]:
    for rs in SID:
        p("   %-8s %-12s slope=%.6f se=%.6f z=%+.4f" % (ct, rs, E[ct][rs]["slope"],
                                                         E[ct][rs]["se"], E[ct][rs]["z"]))

# 工具的条件 F（弱工具预警）
p("")
p("[3] 条件后工具强度（F=(β_cond/se_cond)²）")
weak = {}
for ct, ins in [("B_MEM", "rs12037376"), ("Mono_NC", "rs10917151")]:
    k = "rs56318008"
    r = R[ins]["r"]
    be, se = E[ct][ins]["slope"], E[ct][ins]["se"]
    bk, sek = E[ct][k]["slope"], E[ct][k]["se"]
    be_c = be - r * (se / sek) * bk
    se_c = se * np.sqrt(1 - r ** 2)
    F_c = (be_c / se_c) ** 2
    weak[ct] = dict(be_c=be_c, se_c=se_c, F_c=F_c)
    p("   %-8s 工具 %-12s 边际 F=%.2f  ->  条件 F=%.2f  (β_c=%.6f se_c=%.6f, 条件 z=%+.4f)"
      % (ct, ins, (be / se) ** 2, F_c, be_c, se_c, be_c / se_c))
p("   ★ F_cond < 10 提示条件工具偏弱，条件臂的阴性需与功效一并陈述（正式功效分析属 3.5）")

# ------------------------------------------------------------------ 4) 建 arm
p("")
p("[4] 计算各 arm")
ARMS = []
for ct, ins in [("B_MEM", "rs12037376"), ("Mono_NC", "rs10917151")]:
    k = "rs56318008"
    r = R[ins]["r"]
    be, se_e = E[ct][ins]["slope"], E[ct][ins]["se"]
    bk, se_k = E[ct][k]["slope"], E[ct][k]["se"]
    be_c = be - r * (se_e / se_k) * bk
    var_e_c = se_e ** 2 * (1 - r ** 2)
    bo = df["b_" + ins].values.astype(float)
    vo = df["se_" + ins].values.astype(float) ** 2
    bok = df["b_" + k].values.astype(float)
    vok = df["se_" + k].values.astype(float) ** 2
    bo_c = bo - r * (np.sqrt(vo) / np.sqrt(vok)) * bok
    vo_c = vo * (1 - r ** 2)

    b1, s1, p1, F1 = wald(np.full(len(df), be), np.full(len(df), se_e ** 2), bo, vo)
    b2, s2, p2, F2 = wald(np.full(len(df), be_c), np.full(len(df), var_e_c), bo, vo)
    b3, s3, p3, F3 = wald(np.full(len(df), be_c), np.full(len(df), var_e_c), bo_c, vo_c)
    for tag, bb, ss, pp in [("uncond", b1, s1, p1), ("cond_eqtl_only", b2, s2, p2),
                            ("cond_both", b3, s3, p3)]:
        ARMS.append(dict(family=tag, cell_type=ct, instrument=ins, beta_MR=bb, se_MR=ss,
                         p_MR=pp, F=F1 if tag == "uncond" else F2))
    p("   %-8s %-14s 工具=%s  条件F=%.2f" % (ct, "cond_*", ins, weak[ct]["F_c"]))

# (b) 区域工具：rs56318008 自身
for ct in ["B_MEM", "Mono_NC"]:
    be = E[ct]["rs56318008"]["slope"]; se_e = E[ct]["rs56318008"]["se"]
    bo = df["b_rs56318008"].values.astype(float)
    vo = df["se_rs56318008"].values.astype(float) ** 2
    b, s, pv, F = wald(np.full(len(df), be), np.full(len(df), se_e ** 2), bo, vo)
    ARMS.append(dict(family="region", cell_type=ct, instrument="rs56318008",
                     beta_MR=b, se_MR=s, p_MR=pv, F=F))
    p("   %-8s %-14s 工具=rs56318008  F=%.2f" % (ct, "region", (be / se_e) ** 2))

# ------------------------------------------------------------------ 5) 组装 + FDR
p("")
p("[5] BH-FDR（按 family 合并两细胞类型，m = 该 family 有效检验数）")
res = []
for a in ARMS:
    d = df[["phenocode", "phenotype", "category", "num_cases", "num_controls"]].copy()
    d["family"] = a["family"]; d["cell_type"] = a["cell_type"]; d["instrument"] = a["instrument"]
    d["beta_MR"] = a["beta_MR"]; d["se_MR"] = a["se_MR"]; d["p_MR"] = a["p_MR"]; d["F"] = a["F"]
    res.append(d)
R_all = pd.concat(res, ignore_index=True)

fam_sum = []
for fam, g in R_all.groupby("family", sort=False):
    n_fam = int(g["p_MR"].notna().sum())
    q = bh(g["p_MR"].values)
    R_all.loc[g.index, "p_FDR_family"] = q
    # 单 arm FDR（参考）
    R_all.loc[g.index, "p_FDR_arm"] = np.nan
    for (fam2, ct), g2 in g.groupby(["family", "cell_type"], sort=False):
        R_all.loc[g2.index, "p_FDR_arm"] = bh(g2["p_MR"].values)
    sig = (q < 0.05) & np.isfinite(q)
    up = sig & (R_all.loc[g.index, "beta_MR"].values < 0)
    dn = sig & (R_all.loc[g.index, "beta_MR"].values > 0)
    fam_sum.append(dict(family=fam, n_tests=n_fam, m_family=int(np.isfinite(q).sum()),
                        n_FDR05=int(sig.sum()), n_risk_up=int(up.sum()), n_risk_down=int(dn.sum()),
                        min_p_MR=float(np.nanmin(g["p_MR"])),
                        min_p_MR_risk_up=float(np.nanmin(g.loc[g["beta_MR"] < 0, "p_MR"]))
                        if (g["beta_MR"] < 0).any() else np.nan))
    p("   %-15s m=%d  FDR<0.05=%d（↑风险 %d / ↓保护 %d）  min_p=%.3e  min_p(↑)=%s"
      % (fam, n_fam, int(sig.sum()), int(up.sum()), int(dn.sum()), float(np.nanmin(g["p_MR"])),
         ("%.3e" % np.nanmin(g.loc[g["beta_MR"] < 0, "p_MR"])) if (g["beta_MR"] < 0).any() else "—"))
SUM = pd.DataFrame(fam_sum)
R_all.rename(columns={"num_cases": "n_cases_manifest"}).to_csv(OUT_ALL, index=False, encoding="utf-8-sig")
SUM.to_csv(OUT_SUM, index=False, encoding="utf-8-sig")

# ------------------------------------------------------------------ 6) 域计数（manifest category）
p("")
p("[6] 域计数（manifest `category`）")
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
    "肌肉骨骼/结缔组织": [
        "XIII Diseases of the musculoskeletal system and connective tissue (M13_)",
        "Rheuma endpoints"],
}
CAT2DOM = {c: d for d, cs in DOM.items() for c in cs}
df["domain"] = np.where(df["category"].isin(CAT2DOM), df["category"].map(CAT2DOM), "其他")
chk = df.groupby("domain").size().to_dict()
p("   端点域计数 = %s" % chk)
want = {"免疫/血液/肿瘤": 484, "女性生殖/不孕": 144, "妊娠/分娩/产褥": 111,
        "肌肉骨骼/结缔组织": 273, "其他": 1457}
ok_dom = all(chk.get(k, 0) == v for k, v in want.items())
p("   与 42b 口径核对（×2 应=968/288/222/546/2914）：%s" % ("一致 ✔" if ok_dom else "!! 不一致"))
for k, v in want.items():
    p("     %-14s 期望 %4d  实测 %4d" % (k, v, chk.get(k, 0)))

R_all["domain"] = R_all["phenocode"].map(dict(zip(df["phenocode"], df["domain"])))
rows = []
for (fam, ctd), g in R_all.groupby(["family", "domain"], sort=False):
    sig = g["p_FDR_family"] < 0.05
    rows.append(dict(family=fam, domain=ctd, n_tests=int(g["p_MR"].notna().sum()),
                     n_sig=int(sig.sum()),
                     n_risk_up=int((sig & (g["beta_MR"] < 0)).sum()),
                     n_risk_down=int((sig & (g["beta_MR"] > 0)).sum()),
                     min_p_MR=float(np.nanmin(g["p_MR"])) if g["p_MR"].notna().any() else np.nan))
DOMS = pd.DataFrame(rows).sort_values(["family", "domain"])
DOMS.to_csv(OUT_DOM, index=False, encoding="utf-8-sig")

# ------------------------------------------------------------------ 7) 回归锚点
p("")
p("[7] 回归锚点：uncond arm 应复现既有 4,938 检验结果")
old = pd.read_csv(os.path.join(TAB, "42_phewas_safety_assessment.csv"))
old_sig = old[["rsid", "endpoint", "cell_type", "beta_MR", "se_MR", "p_MR", "p_FDR"]].copy()
sig_unc = R_all[(R_all["family"] == "uncond") & (R_all["p_FDR_family"] < 0.05)]
p("   既有 FDR<0.05 行 = %d ；本次 uncond 复算 FDR<0.05 行 = %d" % (len(old_sig), len(sig_unc)))
mg = old_sig.merge(R_all[R_all["family"] == "uncond"],
                   left_on=["rsid", "endpoint", "cell_type"],
                   right_on=["instrument", "phenocode", "cell_type"], how="left", suffixes=("", "_n"))
p("   键匹配 = %d / %d" % (int(mg["beta_MR_n"].notna().sum()), len(old_sig)))
for c in ["beta_MR", "se_MR", "p_MR", "p_FDR"]:
    d = (mg[c] - mg[c + "_n"]).abs()
    p("     |Δ %-5s| max = %.3e  (>1e-6:%d)" % (c, np.nanmax(d), int((d > 1e-6).sum())))

# ------------------------------------------------------------------ 8) 12 个风险升高端点跟踪
p("")
p("[8] 12 个「风险升高」端点在各 arm 的表现")
r12 = pd.read_csv(os.path.join(PROJ, "supplementary_data",
                              "Fig5bc_phewas_risk_increasing_12.csv"))
keys = r12[["phenotype", "cell_type"]].drop_duplicates()
tr = []
for _, k in keys.iterrows():
    sub = R_all[(R_all["phenotype"] == k["phenotype"]) & (R_all["cell_type"] == k["cell_type"])]
    row = dict(phenotype=k["phenotype"], cell_type=k["cell_type"])
    for fam in ["uncond", "cond_eqtl_only", "cond_both"]:
        q = sub[sub["family"] == fam]
        if len(q):
            row["b_" + fam] = float(q.iloc[0]["beta_MR"])
            row["p_" + fam] = float(q.iloc[0]["p_MR"])
            row["fdr_" + fam] = float(q.iloc[0]["p_FDR_family"])
    q = sub[sub["family"] == "region"]
    if len(q):
        row["b_region"] = float(q.iloc[0]["beta_MR"])
        row["p_region"] = float(q.iloc[0]["p_MR"])
        row["fdr_region"] = float(q.iloc[0]["p_FDR_family"])
    row["domain"] = sub["domain"].iloc[0] if len(sub) else ""
    tr.append(row)
TR = pd.DataFrame(tr)
TR.to_csv(OUT_R12, index=False, encoding="utf-8-sig")
for _, r in TR.iterrows():
    p("   %-52s %-8s | uncond b=%+.3f fdr=%.2e | cond_both b=%+.3f fdr=%.2e | region b=%+.3f fdr=%.2e"
      % (str(r["phenotype"])[:52], r["cell_type"], r.get("b_uncond", np.nan),
         r.get("fdr_uncond", np.nan), r.get("b_cond_both", np.nan), r.get("fdr_cond_both", np.nan),
         r.get("b_region", np.nan), r.get("fdr_region", np.nan)))

p("")
p("产物：%s\n      %s\n      %s\n      %s" % (OUT_ALL, OUT_SUM, OUT_DOM, OUT_R12))
p("总耗时 %.1f s   末次可用内存 %.2f GB" % (time.time() - T0, free_gb()))
p("VERDICT = %s" % ("PASS" if ok_dom else "WARN(域计数不一致)"))
io.open(LOG, "w", encoding="utf-8", newline="\n").write("\n".join(L) + "\n")
print("WROTE " + LOG)
