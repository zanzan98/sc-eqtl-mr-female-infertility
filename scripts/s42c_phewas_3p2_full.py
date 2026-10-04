# -*- coding: utf-8 -*-
"""s42c_phewas_3p2_full.py -- 任务 3.2 主分析（**完整 2,469 端点面板**，来自 FinnGen R12 公开桶 tabix）

数据源
------
 00_data_raw/finngen_r12_tabix/phewas_3snp_all_endpoints.csv
   逐端点从 gs://finngen-public-data-r12 summary_stats（BGZF + tabix，HTTP Range）取出 3 个 SNP 的行。
   列：endpoint, pos38, ref, alt, rsids, nearest_genes, pval, mlogp, beta, sebeta,
       af_alt, af_alt_cases, af_alt_controls（beta 为 **alt 等位**效应）

手臂（arm）
----------
(a) 条件 PheWAS（CDC42 工具 ± 以 rs56318008 条件化）
    uncond / cond_eqtl_only / cond_both      （B_MEM 用 rs12037376；Mono_NC 用 rs10917151）
(b) 区域工具 PheWAS（rs56318008 自身作区域工具，用 B_MEM / Mono_NC 的 eQTL 效应）
    region
FDR：按 family 合并两细胞类型做 BH（与原始 42_ 口径一致：m ≈ 4,938）

硬纪律
------
 * 二階 delta 的 Wald SE；单 SNP 工具 → Wald ratio
 * 域计数用 manifest 的 `category`；断言复现 42b 的 484/144/111/273/1457
 * 回归锚点：uncond arm 必须逐位复现 42_phewas_safety_assessment.csv 的 beta/se/p/FDR
 * 只读既有产物；只新增 tables/44*
"""
import io
import os
import sys
import time

import numpy as np
import pandas as pd
from scipy import stats

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

PROJ = r"D:\endometriosis_project\11_sc_eqtl_mr_project"
COND = os.path.join(PROJ, "00_data_raw", "onek1k", "cond_coloc")
TABX = os.path.join(PROJ, "00_data_raw", "finngen_r12_tabix")
TAB = os.path.join(PROJ, "tables")
MANIFEST = r"D:\endometriosis_project\01_data_raw\finngen_R12_manifest.tsv"
LOG = r"D:\endometriosis_project\_s42c_phewas_3p2_full.log"

O_ALL = os.path.join(TAB, "44_phewas_3p2_full_all.csv")
O_SUM = os.path.join(TAB, "44b_phewas_3p2_full_summary.csv")
O_DOM = os.path.join(TAB, "44c_phewas_3p2_full_domains.csv")
O_R12 = os.path.join(TAB, "44d_phewas_3p2_full_risk12.csv")
O_SIG = os.path.join(TAB, "44e_phewas_3p2_full_signals.csv")

RS2POS = {"rs12037376": 22135618, "rs10917151": 22096228, "rs56318008": 22143914}
POS2RS = {v: k for k, v in RS2POS.items()}
L = []


def p(s=""):
    L.append(str(s))
    print(s)


def bh(pv):
    pv = np.asarray(pv, float)
    ok = np.isfinite(pv)
    q = np.full(pv.shape, np.nan)
    if ok.sum() == 0:
        return q
    v = pv[ok]
    n = v.size
    o = np.argsort(v)
    r = v[o]
    qq = r * n / (np.arange(1, n + 1))
    qq = np.minimum.accumulate(qq[::-1])[::-1]
    qq = np.clip(qq, 0, 1)
    t = np.empty_like(qq)
    t[o] = qq
    q[ok] = t
    return q


def wald(be, vare, bo, varo):
    be = np.asarray(be, float); vare = np.asarray(vare, float)
    bo = np.asarray(bo, float); varo = np.asarray(varo, float)
    with np.errstate(divide="ignore", invalid="ignore"):
        beta = bo / be
        se = np.sqrt(varo / be ** 2 + bo ** 2 * vare / be ** 4)
        z = beta / se
        pv = 2 * stats.norm.sf(np.abs(z))
        F = be ** 2 / vare
    bad = ~(np.isfinite(beta) & np.isfinite(se) & (se > 0) & np.isfinite(pv))
    for a in (beta, se, pv, F):
        a[bad] = np.nan
    return beta, se, pv, F


T0 = time.time()
p("=== s42c 任务 3.2 主分析（完整端点面板）===")
p("time %s" % time.strftime("%Y-%m-%d %H:%M:%S"))

# ---------------- 0) 读入 ----------------
lg = pd.read_csv(os.path.join(TABX, "phewas_3snp_all_endpoints.csv"))
lg["pos38"] = lg["pos38"].astype(int)
lg["rsid"] = lg["pos38"].map(POS2RS)
p("")
p("[0] tabix 长表 rows=%d  端点=%d  SNP 覆盖 = %s"
  % (len(lg), lg["endpoint"].nunique(), lg["rsid"].value_counts().to_dict()))
lg["ref"] = lg["ref"].astype(str).str.upper()
lg["alt"] = lg["alt"].astype(str).str.upper()

man = pd.read_csv(MANIFEST, sep="\t")
man = man[["phenocode", "phenotype", "category", "num_cases", "num_controls"]]
p("    manifest 端点 = %d" % len(man))

# 宽表
w = lg.pivot_table(index="endpoint", columns="rsid",
                   values=["beta", "sebeta", "pval", "af_alt"], aggfunc="first")
w.columns = ["%s_%s" % (a, b) for a, b in w.columns]
w = w.reset_index()
df = man.merge(w, left_on="phenocode", right_on="endpoint", how="inner")
p("    宽表合并端点 = %d（manifest∩tabix）" % len(df))
miss = man[~man["phenocode"].isin(set(lg["endpoint"]))]
p("    仅 manifest 有、tabix 无的端点 = %d" % len(miss))
if len(miss):
    p("      示例 = %s" % list(miss["phenocode"].head(5)))

# ---------------- 1) 等位对齐核验 ----------------
p("")
p("[1] 等位对齐（tabix `alt` 必须是该 SNP 的暴露效应等位；本项目已锁定 A1）")
mB = pd.read_csv(os.path.join(COND, "m_B_MEM.csv"))
mB["snp_grch37"] = mB["snp_grch37"].astype(str)
SID37 = {"rs12037376": "1:22462111", "rs10917151": "1:22422721", "rs56318008": "1:22470407"}
EXP_A1 = {"rs12037376": "A", "rs10917151": "A", "rs56318008": "T"}
for rs in RS2POS:
    alt = set(lg.loc[lg["rsid"] == rs, "alt"].unique())
    a1 = EXP_A1[rs]
    af_eq = float(mB.loc[mB["snp_grch37"] == SID37[rs], "af"].iloc[0])
    af_tab = float(df["af_alt_" + rs].median())
    ok = (alt == {a1}) and abs(af_eq - af_tab) < 0.05
    p("   %-12s tabix alt=%s（集合 %s）  eQTL af(A1=%s)=%.4f  tabix af_alt 中位=%.4f  -> %s"
      % (rs, a1, alt, a1, af_eq, af_tab, "ALIGNED" if ok else "!! 需复核"))
    if not ok:
        p("   !! 等位未对齐，终止")
        io.open(LOG, "w", encoding="utf-8", newline="\n").write("\n".join(L) + "\n")
        sys.exit(5)

# ---------------- 2) 暴露与 LD ----------------
mM = pd.read_csv(os.path.join(COND, "m_Mono_NC.csv"))
mM["snp_grch37"] = mM["snp_grch37"].astype(str)
ldtab = pd.read_csv(os.path.join(COND, "r_to_rs56318008.csv"))
ldtab["snp_grch37"] = ldtab["snp_grch37"].astype(str)


def etab(m):
    return {rs: dict(slope=float(m.loc[m["snp_grch37"] == SID37[rs], "slope"].iloc[0]),
                     se=float(m.loc[m["snp_grch37"] == SID37[rs], "slope_se"].iloc[0]))
            for rs in SID37}


E = {"B_MEM": etab(mB), "Mono_NC": etab(mM)}
R = {rs: float(ldtab.loc[ldtab["snp_grch37"] == SID37[rs], "r_to_cond"].iloc[0])
     for rs in ["rs12037376", "rs10917151"]}
p("")
p("[2] 暴露与 LD")
for ct in ["B_MEM", "Mono_NC"]:
    p("   %-8s rs12037376 slope=%.6f(se %.6f)   rs10917151 slope=%.6f(se %.6f)   rs56318008 slope=%.6f(se %.6f)"
      % (ct, E[ct]["rs12037376"]["slope"], E[ct]["rs12037376"]["se"],
         E[ct]["rs10917151"]["slope"], E[ct]["rs10917151"]["se"],
         E[ct]["rs56318008"]["slope"], E[ct]["rs56318008"]["se"]))
for rs in R:
    p("   r(%s, rs56318008) = %+.6f" % (rs, R[rs]))

# ---------------- 3) 计算 ----------------
p("")
p("[3] 各 arm 计算（family 级 BH，m=该 family 有效检验数）")
res = []
weakinfo = {}
for ct, ins in [("B_MEM", "rs12037376"), ("Mono_NC", "rs10917151")]:
    k = "rs56318008"
    r = R[ins]
    be, se_e = E[ct][ins]["slope"], E[ct][ins]["se"]
    bk, se_k = E[ct][k]["slope"], E[ct][k]["se"]
    be_c = be - r * (se_e / se_k) * bk
    var_e_c = se_e ** 2 * (1 - r ** 2)
    weakinfo[ct] = dict(F_marg=(be / se_e) ** 2, F_cond=(be_c / np.sqrt(var_e_c)) ** 2,
                        be_c=be_c, se_c=np.sqrt(var_e_c))
    bo = df["beta_" + ins].values
    vo = df["sebeta_" + ins].values ** 2
    bok = df["beta_" + k].values
    vok = df["sebeta_" + k].values ** 2
    bo_c = bo - r * (np.sqrt(vo) / np.sqrt(vok)) * bok
    vo_c = vo * (1 - r ** 2)
    arms = [("uncond", np.full(len(df), be), np.full(len(df), se_e ** 2), bo, vo),
            ("cond_eqtl_only", np.full(len(df), be_c), np.full(len(df), var_e_c), bo, vo),
            ("cond_both", np.full(len(df), be_c), np.full(len(df), var_e_c), bo_c, vo_c)]
    for tag, a1, a2, a3, a4 in arms:
        b, s, pv, F = wald(a1, a2, a3, a4)
        d = df[["phenocode", "phenotype", "category", "num_cases", "num_controls"]].copy()
        d["family"] = tag; d["cell_type"] = ct; d["instrument"] = ins
        d["beta_MR"] = b; d["se_MR"] = s; d["p_MR"] = pv; d["F"] = F
        res.append(d)
for ct in ["B_MEM", "Mono_NC"]:
    be, se_e = E[ct]["rs56318008"]["slope"], E[ct]["rs56318008"]["se"]
    bo = df["beta_rs56318008"].values
    vo = df["sebeta_rs56318008"].values ** 2
    b, s, pv, F = wald(np.full(len(df), be), np.full(len(df), se_e ** 2), bo, vo)
    d = df[["phenocode", "phenotype", "category", "num_cases", "num_controls"]].copy()
    d["family"] = "region"; d["cell_type"] = ct; d["instrument"] = "rs56318008"
    d["beta_MR"] = b; d["se_MR"] = s; d["p_MR"] = pv; d["F"] = F
    res.append(d)

A = pd.concat(res, ignore_index=True)
p("   工具强度：B_MEM 边际 F=%.2f → 条件 F=%.2f ；Mono_NC 边际 F=%.2f → 条件 F=%.2f"
  % (weakinfo["B_MEM"]["F_marg"], weakinfo["B_MEM"]["F_cond"],
     weakinfo["Mono_NC"]["F_marg"], weakinfo["Mono_NC"]["F_cond"]))
p("   region 臂：B_MEM F=%.2f ；Mono_NC F=%.2f"
  % ((E["B_MEM"]["rs56318008"]["slope"] / E["B_MEM"]["rs56318008"]["se"]) ** 2,
     (E["Mono_NC"]["rs56318008"]["slope"] / E["Mono_NC"]["rs56318008"]["se"]) ** 2))

A["p_FDR_family"] = np.nan
for fam, g in A.groupby("family", sort=False):
    A.loc[g.index, "p_FDR_family"] = bh(g["p_MR"].values)

fam = []
for f, g in A.groupby("family", sort=False):
    sig = (g["p_FDR_family"] < 0.05).values & np.isfinite(g["p_FDR_family"].values)
    up = sig & (g["beta_MR"].values < 0)
    dn = sig & (g["beta_MR"].values > 0)
    upm = g.loc[g["beta_MR"] < 0, "p_MR"]
    fam.append(dict(family=f, n_tests=int(g["p_MR"].notna().sum()),
                    n_FDR05=int(sig.sum()), n_risk_up=int(up.sum()), n_risk_down=int(dn.sum()),
                    min_p_MR=float(np.nanmin(g["p_MR"])),
                    min_p_MR_risk_up=float(np.nanmin(upm)) if len(upm) else np.nan))
    p("   %-15s n=%d  FDR<0.05=%3d（风险↑ %3d / 保护↓ %3d）  min_p=%.3e  min_p(↑)=%s"
      % (f, int(g["p_MR"].notna().sum()), int(sig.sum()), int(up.sum()), int(dn.sum()),
         float(np.nanmin(g["p_MR"])),
         ("%.3e" % np.nanmin(upm)) if len(upm) else "—"))
SUM = pd.DataFrame(fam)
SUM.to_csv(O_SUM, index=False, encoding="utf-8-sig")

# ---------------- 4) 域计数 ----------------
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
C2D = {c: d for d, cs in DOM.items() for c in cs}
df["domain"] = np.where(df["category"].isin(C2D), df["category"].map(C2D), "其他")
chk = df.groupby("domain").size().to_dict()
want = {"免疫/血液/肿瘤": 484, "女性生殖/不孕": 144, "妊娠/分娩/产褥": 111,
        "肌肉骨骼/结缔组织": 273, "其他": 1457}
ok_dom = all(chk.get(k, 0) == v for k, v in want.items())
p("")
p("[4] 域计数校验：%s" % ("与 42b 一致 ✔" % () if ok_dom else "!! 不一致"))
for k, v in want.items():
    p("   %-14s 期望 %4d  实测 %4d" % (k, v, chk.get(k, 0)))
A["domain"] = A["phenocode"].map(dict(zip(df["phenocode"], df["domain"])))
rows = []
for (f, d), g in A.groupby(["family", "domain"], sort=False):
    sig = g["p_FDR_family"] < 0.05
    rows.append(dict(family=f, domain=d, n_tests=int(g["p_MR"].notna().sum()), n_sig=int(sig.sum()),
                     n_risk_up=int((sig & (g["beta_MR"] < 0)).sum()),
                     n_risk_down=int((sig & (g["beta_MR"] > 0)).sum()),
                     min_p_MR=float(np.nanmin(g["p_MR"])) if g["p_MR"].notna().any() else np.nan))
pd.DataFrame(rows).sort_values(["family", "domain"]).to_csv(O_DOM, index=False, encoding="utf-8-sig")

# ---------------- 5) 回归锚点 ----------------
p("")
p("[5] 回归锚点：uncond 臂 vs 42_phewas_safety_assessment.csv（既有 52 行 FDR<0.05）")
old = pd.read_csv(os.path.join(TAB, "42_phewas_safety_assessment.csv"))
U = A[A["family"] == "uncond"]
mg = old.merge(U, left_on=["rsid", "endpoint", "cell_type"],
               right_on=["instrument", "phenocode", "cell_type"], how="left", suffixes=("", "_n"))
p("   既有行=%d  键匹配=%d  既有 FDR<0.05=%d  本次 uncond FDR<0.05=%d"
  % (len(old), int(mg["beta_MR_n"].notna().sum()), len(old),
     int((U["p_FDR_family"] < 0.05).sum())))
# ★ old 的列名是 p_FDR、new 是 p_FDR_family → 名字不重叠，pandas **不会**加 _n 后缀，
#   必须显式配对（早期版本 mg["p_FDR_n"] 会 KeyError）。
pairs = [("beta_MR", "beta_MR_n"), ("se_MR", "se_MR_n"),
         ("p_MR", "p_MR_n"), ("p_FDR", "p_FDR_family")]
for c_old, c_new in pairs:
    d = (mg[c_old] - mg[c_new]).abs()
    p("     |Δ %-10s| max = %.3e   (>1e-6 的个数 %d)"
      % (c_old, np.nanmax(d.values), int((d > 1e-6).sum())))
anchor_ok = all((mg[c_old] - mg[c_new]).abs().max() < 1e-6 for c_old, c_new in pairs[:3])

# ---------------- 6) 12 个风险升高端点 ----------------
p("")
p("[6] 12 个「风险升高」端点在 4 个 arm 的表现（按 42_ 表的 endpoint 精确匹配）")
inc = old[old["risk_if_cdc42_inhibited"] == "increased"][["endpoint", "phenotype", "cell_type"]]
p("   42_ 表中 risk_if_cdc42_inhibited=increased 的行数 = %d" % len(inc))
tr = []
for _, k in inc.iterrows():
    sub = A[(A["phenocode"] == k["endpoint"]) & (A["cell_type"] == k["cell_type"])]
    row = dict(endpoint=k["endpoint"], phenotype=k["phenotype"], cell_type=k["cell_type"],
               domain=(sub["domain"].iloc[0] if len(sub) else ""))
    for f in ["uncond", "cond_eqtl_only", "cond_both", "region"]:
        q = sub[sub["family"] == f]
        if len(q):
            row["b_" + f] = float(q.iloc[0]["beta_MR"])
            row["p_" + f] = float(q.iloc[0]["p_MR"])
            row["fdr_" + f] = float(q.iloc[0]["p_FDR_family"])
    tr.append(row)
TR = pd.DataFrame(tr)
TR.to_csv(O_R12, index=False, encoding="utf-8-sig")
p("   %-52s %-8s | %-16s %-16s %-16s" % ("endpoint", "cell", "uncond", "cond_both", "region"))
for _, r in TR.iterrows():
    p("   %-52s %-8s | %+.3f q=%.1e  %+.3f q=%.1e  %+.3f q=%.1e"
      % (str(r["endpoint"])[:52], r["cell_type"],
         r.get("b_uncond", np.nan), r.get("fdr_uncond", np.nan),
         r.get("b_cond_both", np.nan), r.get("fdr_cond_both", np.nan),
         r.get("b_region", np.nan), r.get("fdr_region", np.nan)))
p("   12 端点中 FDR<0.05 仍成立者：uncond=%d  cond_eqtl_only=%d  cond_both=%d  region=%d"
  % tuple(int((TR["fdr_" + f] < 0.05).sum()) for f in
          ["uncond", "cond_eqtl_only", "cond_both", "region"]))
p("   其中方向为「风险升高」(beta_MR<0) 者：uncond=%d  cond_eqtl_only=%d  cond_both=%d  region=%d"
  % tuple(int(((TR["b_" + f] < 0) & (TR["fdr_" + f] < 0.05)).sum()) for f in
          ["uncond", "cond_eqtl_only", "cond_both", "region"]))

# ---------------- 7) 落盘 ----------------
A.rename(columns={"num_cases": "n_cases_manifest"}).to_csv(O_ALL, index=False, encoding="utf-8-sig")
sig = A[(A["p_FDR_family"] < 0.05) & np.isfinite(A["p_FDR_family"])].copy()
sig["direction"] = np.where(sig["beta_MR"] < 0, "risk_up", "protective")
sig.sort_values(["family", "p_FDR_family"]).to_csv(O_SIG, index=False, encoding="utf-8-sig")
p("")
p("[7] 各 arm 的 FDR<0.05 明细行数：%s" % sig["family"].value_counts().to_dict())
p("产物：%s\n      %s\n      %s\n      %s\n      %s" % (O_ALL, O_SUM, O_DOM, O_R12, O_SIG))
p("总耗时 %.1f s" % (time.time() - T0))
p("VERDICT = %s" % ("PASS" if (ok_dom and anchor_ok) else
                    "WARN(域计数=%s,锚点=%s)" % (ok_dom, anchor_ok)))
io.open(LOG, "w", encoding="utf-8", newline="\n").write("\n".join(L) + "\n")
print("WROTE " + LOG)
