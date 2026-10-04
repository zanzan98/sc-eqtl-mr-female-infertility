# -*- coding: utf-8 -*-
"""s60d_mvmr_decomposition.py -- 任务B 第3步：效应分解对照（同一工具集下 单暴露 IVW vs 双暴露 MVMR）

对 s60c 判定为「可识别」的每个 (细胞类型 × 剪枝阈值) 组合，在**完全相同的工具集**上比较：
  ① 单暴露 IVW（CDC42 单独 / LINC00339 单独）—— 未调整
  ② 双暴露 MVMR-IVW（CDC42 + LINC00339）—— 互为调整（取自 tables/61）
从而量化「条件后效应如何分解」。

单暴露 IVW（与 MVMR 同权重口径）：β = Σ(w_k β_Xk β_Yk)/Σ(w_k β_Xk²)，w_k = 1/se_Yk²
                                  se = sqrt(1/Σ(w_k β_Xk²))

读：00_data_raw/mvmr_L1/<cell>.mvmr_input.csv、ld/<cell>.ld、tables/61_mvmr_results.csv
写：tables/62_mvmr_decomposition.csv
"""
import io
import math
import os
import sys

import numpy as np
import pandas as pd

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
PROJ = r"D:\endometriosis_project\11_sc_eqtl_mr_project"
IN = os.path.join(PROJ, "00_data_raw", "mvmr_L1")
LDD = os.path.join(IN, "ld")
TAB = os.path.join(PROJ, "tables")
LOG = os.path.join(PROJ, "scripts", "_s60d_mvmr_decomposition.log")
THRS = [0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.5]

buf = []
def w(s=""):
    buf.append(str(s)); print(s)

w("=" * 110)
w("s60d — MVMR 效应分解（同一工具集：单暴露 IVW vs 双暴露 MVMR）")
w("=" * 110)

R = pd.read_csv(os.path.join(TAB, "61_mvmr_results.csv"))
R["ident"] = R["identifiable"].astype(str).str.startswith("可识别")

rows = []
for _, r in R.iterrows():
    cell, thr = r["cell"], r["thr"]
    fin = os.path.join(IN, "%s.mvmr_input.csv" % cell)
    fld = os.path.join(LDD, "%s.ld" % cell)
    if not (os.path.exists(fin) and os.path.exists(fld)):
        continue
    D = pd.read_csv(fin, encoding="utf-8-sig")
    D.columns = ["snp"] + list(D.columns[1:])
    LD = pd.read_csv(fld, sep=r"\s+", header=None).values
    if LD.shape[0] != len(D):
        continue
    ordr = np.argsort(np.minimum(D["p_cdc42"].values, D["p_linc00339"].values), kind="stable")
    keep = []
    for i in ordr:
        if not keep:
            keep.append(i); continue
        if np.all((LD[i, keep] ** 2) < thr):
            keep.append(i)
    kk = len(keep)
    if kk < 2:
        continue
    S = D.iloc[keep].reset_index(drop=True)
    wy = 1.0 / S["se_out"].values ** 2

    def ivw(bx, sebx):
        num = np.sum(wy * bx * S["b_out"].values)
        den = np.sum(wy * bx ** 2)
        b = num / den
        se = np.sqrt(1.0 / den)
        z = abs(b / se)
        return b, se, math.erfc(z / math.sqrt(2.0))   # 双侧正态 p

    b1, se1, p1 = ivw(S["b_cdc42"].values, S["se_cdc42"].values)
    b2, se2, p2 = ivw(S["b_linc00339"].values, S["se_linc00339"].values)

    # ---- 双暴露 MVMR-IVW 独立复算 + 稳健（sandwich）SE ----
    B = np.column_stack([S["b_cdc42"].values, S["b_linc00339"].values])
    G = S["b_out"].values
    Wm = np.diag(wy)
    BtWB = B.T @ Wm @ B
    BtWB_inv = np.linalg.pinv(BtWB)
    bm = BtWB_inv @ B.T @ Wm @ G
    resid = G - B @ bm
    meat = B.T @ Wm @ np.diag(resid ** 2) @ Wm @ B
    cov_sand = BtWB_inv @ meat @ BtWB_inv          # 异质性稳健
    cov_fix = BtWB_inv                             # 无过识别残差时退化
    se_sand = np.sqrt(np.diag(cov_sand))
    se_fix = np.sqrt(np.diag(cov_fix))

    rows.append(dict(
        cell=cell, thr=thr, k=kk, identifiable=r["identifiable"],
        uni_CDC42_b=b1, uni_CDC42_se=se1, uni_CDC42_p=p1,
        uni_LINC00339_b=b2, uni_LINC00339_se=se2, uni_LINC00339_p=p2,
        mvmr_CDC42_b=r["b_CDC42"], mvmr_CDC42_se=r["se_CDC42"], mvmr_CDC42_p=r["p_CDC42"],
        mvmr_LINC00339_b=r["b_LINC00339"], mvmr_LINC00339_se=r["se_LINC00339"],
        mvmr_LINC00339_p=r["p_LINC00339"],
        re_CDC42_b=bm[0], re_CDC42_se_fix=se_fix[0], re_CDC42_se_sand=se_sand[0],
        re_LINC00339_b=bm[1], re_LINC00339_se_fix=se_fix[1], re_LINC00339_se_sand=se_sand[1],
        F_cond_CDC42=r["F_cond_CDC42"], F_cond_LINC00339=r["F_cond_LINC00339"],
        ld_cond=r["ld_cond"]))

T = pd.DataFrame(rows)
T.to_csv(os.path.join(TAB, "62_mvmr_decomposition.csv"), index=False, encoding="utf-8-sig")

idn = T[T["identifiable"].astype(str).str.startswith("可识别")].copy()
w("可识别组合数 = %d（含恰好识别 k=2 时 se 不可估计者）" % len(idn))
w("")
w("%-9s %-7s %4s | %-28s | %-28s" % ("cell", "thr", "k", "CDC42  单→双", "LINC00339 单→双"))
w("-" * 110)
for _, r in idn.iterrows():
    def f(b, se):
        return ("%+.4f" % b) + (" (se %.4f)" % se if pd.notna(se) else " (se NA)")
    w("%-9s %-7g %4d | %-28s | %-28s"
      % (r["cell"], r["thr"], r["k"],
         "%s → %s" % (f(r["uni_CDC42_b"], r["uni_CDC42_se"]), f(r["mvmr_CDC42_b"], r["mvmr_CDC42_se"])),
         "%s → %s" % (f(r["uni_LINC00339_b"], r["uni_LINC00339_se"]),
                      f(r["mvmr_LINC00339_b"], r["mvmr_LINC00339_se"]))))
w("")
w("已写 tables/62_mvmr_decomposition.csv（%d 行）" % len(T))
w("")
w("★ 独立复算 vs MVMR::ivw_mvmr 一致性与稳健 SE（k>=3 才有过识别残差）")
w("%-9s %-6s %4s | %-26s | %-30s" % ("cell", "thr", "k", "CDC42 (包 / 复算)", "LINC00339 (包 / 复算)"))
w("-" * 110)
for _, r in idn.iterrows():
    w("%-9s %-6g %4d | %-26s | %-30s"
      % (r["cell"], r["thr"], r["k"],
         "%+.4f / %+.4f" % (r["mvmr_CDC42_b"], r["re_CDC42_b"]),
         "%+.4f / %+.4f" % (r["mvmr_LINC00339_b"], r["re_LINC00339_b"])))
w("")
w("（复算点估计与包内一致即为通过；稳健 SE 见 CSV 的 re_*_se_sand 列）")
io.open(LOG, "w", encoding="utf-8", newline="\n").write("\n".join(buf) + "\n")
print("WROTE " + LOG)
