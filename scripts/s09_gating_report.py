# -*- coding: utf-8 -*-
"""s09_gating_report.py —— 阳性对照门控的形式化判定

两层判据（不可混用 BH 口径）：
 ① 预设面板（13 个 a-priori 检验）：小样本定向假设检验 → 报告名义 P + 符号检验（方向一致性），
    **不**用全量 8,596 检验的 BH 值判定（那会把"是否全基因组最强"错当成"基因是否因果"）。
 ② 已知 RA 位点回收率：以全量检验的 FDR<0.05 率为背景，用 Fisher 精确检验比较
    "已知 RA 基因命中率 vs 背景命中率" → 这是管道有效性的主判据。
"""
import os, traceback
import numpy as np
import pandas as pd
from scipy import stats

OUT1 = r"D:\endometriosis_project\11_sc_eqtl_mr_project\tables\03_门控判定_阳性对照.csv"
OUT2 = r"D:\endometriosis_project\_gating_stats.txt"
ROOT = r"D:\endometriosis_project\11_sc_eqtl_mr_project"
L = []
try:
    d = pd.read_csv(os.path.join(ROOT, r"tables\02_阳性对照_RA全量MR.csv"))
    v = d[d['p_MR'].notna()].copy()
    N = len(v)
    n_bg = int((v['p_FDR'] < 0.05).sum())
    L.append("=== 基础计数 ===")
    L.append("工具变量行 = %d ; 有效检验 = %d ; FDR<0.05 = %d (%.2f%%)" % (len(d), N, n_bg, 100.0 * n_bg / N))
    L.append("")

    # ---------- ① 预设面板 ----------
    preset = ["GLIPR1", "XBP1", "IKZF1", "IL2RA"]
    g = v[v['gene'].isin(preset)].copy().sort_values('p_MR')
    L.append("=== ① 预设阳性对照面板（a-priori 定向检验）===")
    L.append("检验数 = %d（该面板仅作方法学质控，不进入主分析筛选）" % len(g))
    L.append("  %-7s %-10s %-15s %9s %8s %11s %10s" % ("gene", "cell", "variant", "beta_MR", "z_MR", "p_MR", "p_FDR"))
    for _, r in g.iterrows():
        L.append("  %-7s %-10s %-15s %9.4f %8.2f %11.2e %10.2e" % (
            r['gene'], r['cell_type'], str(r['variant_id']), r['beta_MR'], r['z_MR'], r['p_MR'], r['p_FDR']))
    n_nom = int((g['p_MR'] < 0.05).sum())
    n_neg = int((g['beta_MR'] < 0).sum())
    L.append("")
    L.append("  名义 P<0.05 = %d / %d (%.1f%%)；期望背景 5%%" % (n_nom, len(g), 100.0 * n_nom / len(g)))
    L.append("  方向一致（beta_MR<0）= %d / %d" % (n_neg, len(g)))
    # 二项检验：命中数 vs p=0.05
    bt = stats.binomtest(n_nom, len(g), 0.05, alternative='greater')
    L.append("  二项检验（H0: 命中率=5%%）p = %.4f" % bt.pvalue)
    # 符号检验：方向一致性
    st = stats.binomtest(n_neg, len(g), 0.5, alternative='two-sided')
    L.append("  符号检验（H0: 正负各半）p = %.4f" % st.pvalue)
    L.append("  → 判定：%s" % ("通过（名义富集 + 方向一致）" if (bt.pvalue < 0.05 and st.pvalue < 0.05)
                              else ("部分通过（仅方向一致）" if st.pvalue < 0.05 else "不通过")))
    g.to_csv(OUT1, index=False, encoding='utf-8-sig')
    L.append("  已写出 %s" % os.path.basename(OUT1))
    L.append("")

    # ---------- ② 已知 RA 位点回收率 ----------
    known = {
        'HLA': ["HLA-DRB1", "HLA-DRB5", "HLA-DQA1", "HLA-DQA2", "HLA-DQB1", "HLA-DPB1", "HLA-B", "HLA-A"],
        'nonHLA': ["PADI4", "PADI2", "FCRL3", "FCRL5", "TNFRSF14", "FCGR2A", "FCGR2B", "RGS1", "MMEL1",
                   "IL12RB2", "ORMDL3", "TNFSF13B", "ERAP1", "IL2RA", "IKZF1", "CTLA4", "TRAF1",
                   "ANKRD55", "UBE2L3", "BLK", "ETS1", "SLC22A5"],
    }
    L.append("=== ② 已知 RA 易感位点回收率（管道有效性主判据）===")
    L.append("注：清单为**人工整理的经典 RA GWAS 位点（非穷尽）**，仅用于门控评价")

    def rate(genes, label):
        avail = [x for x in genes if x in set(v['gene'])]
        absent = [x for x in genes if x not in set(v['gene'])]
        hit = []
        for x in avail:
            s = v[v['gene'] == x]
            if s['p_FDR'].min() < 0.05:
                hit.append(x)
        return avail, absent, hit

    res_rows = []
    for lbl, genes in known.items():
        avail, absent, hit = rate(genes, lbl)
        L.append("")
        L.append("  [%s] 在工具变量集中可用 %d / 不在集合 %d" % (lbl, len(avail), len(absent)))
        L.append("       回收（best FDR<0.05）= %d / %d = %.1f%%" % (len(hit), len(avail), 100.0 * len(hit) / max(len(avail), 1)))
        L.append("       命中: %s" % ", ".join(sorted(hit)))
        miss = [x for x in avail if x not in hit]
        L.append("       未命中: %s" % (", ".join(sorted(miss)) if miss else "无"))
        res_rows.append((lbl, len(hit), len(avail), 100.0 * len(hit) / max(len(avail), 1)))

    tot_hit = sum(r[1] for r in res_rows)
    tot_avail = sum(r[2] for r in res_rows)
    L.append("")
    L.append("  合计回收 = %d / %d = %.1f%%（背景 FDR<0.05 率 = %.2f%%）" % (
        tot_hit, tot_avail, 100.0 * tot_hit / max(tot_avail, 1), 100.0 * n_bg / N))
    # Fisher 精确检验：已知基因命中率 vs 背景命中率
    table = [[tot_hit, tot_avail - tot_hit], [n_bg, N - n_bg]]
    _, fisher_p = stats.fisher_exact(table, alternative='greater')
    L.append("  Fisher 精确检验（已知基因命中率 > 背景命中率）: p = %.3e" % fisher_p)
    L.append("  富集倍数 = %.1f x" % ((100.0 * tot_hit / max(tot_avail, 1)) / (100.0 * n_bg / N)))
    # 仅非 HLA 基因（更严格，排除 MHC 长 LD 影响）
    nh = [r for r in res_rows if r[0] == 'nonHLA'][0]
    table2 = [[nh[1], nh[2] - nh[1]], [n_bg, N - n_bg]]
    _, fisher_p2 = stats.fisher_exact(table2, alternative='greater')
    L.append("  ★ 仅非 HLA 基因（排除 MHC 长 LD 影响）: %d/%d = %.1f%%, Fisher p = %.3e, 富集 %.1f x" % (
        nh[1], nh[2], nh[3], fisher_p2, nh[3] / (100.0 * n_bg / N)))
    L.append("")

    # ---------- ③ 稳健性 ----------
    L.append("=== ③ 稳健性检查 ===")
    v['chr'] = v['variant_id'].astype(str).str.split(':').str[0]
    v['pos'] = v['variant_id'].astype(str).str.split(':').str[1].astype(float)
    ismhc = (v['chr'] == '6') & (v['pos'] >= 25e6) & (v['pos'] <= 35e6)
    L.append("  FDR<0.05 中 MHC 占 %d，非 MHC 占 %d" % (int(((v['p_FDR'] < 0.05) & ismhc).sum()),
                                                     int(((v['p_FDR'] < 0.05) & ~ismhc).sum())))
    L.append("  剔除 MHC 后重算 BH：")
    sub = v[~ismhc].copy()
    pv = sub['p_MR'].values
    ns = len(pv); order = np.argsort(pv); ranked = pv[order]
    q = np.minimum.accumulate((ranked * ns / np.arange(1, ns + 1))[::-1])[::-1]
    tmp = np.empty_like(q); tmp[order] = np.clip(q, 0, 1)
    sub['p_FDR_nomhc'] = tmp
    L.append("    非 MHC 检验 = %d ; 剔除 MHC 后 FDR<0.05 = %d" % (ns, int((sub['p_FDR_nomhc'] < 0.05).sum())))
    nhit = []
    for x in known['nonHLA']:
        s = sub[sub['gene'] == x]
        if len(s) and s['p_FDR_nomhc'].min() < 0.05:
            nhit.append(x)
    L.append("    非 HLA 已知基因回收（剔除 MHC 后）= %d / %d = %.1f%%" % (len(nhit), nh[2], 100.0 * len(nhit) / max(nh[2], 1)))
    L.append("    命中: %s" % ", ".join(sorted(nhit)))
except Exception:
    L.append("!!! EXCEPTION !!!")
    L.append(traceback.format_exc())
open(OUT2, "w", encoding="utf-8").write("\n".join(L))
