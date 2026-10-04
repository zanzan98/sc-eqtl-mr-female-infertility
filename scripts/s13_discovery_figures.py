#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""s13_discovery_figures.py —— Discovery 阶段诊断图

产出（figures/，PDF + PNG @600 dpi，全部文字 ASCII）：
  FigD1_discovery_volcano     火山图：MR 效应量 vs -log10(P)，标注 MHC 与 FDR 显著
  FigD2_discovery_manhattan   14 细胞类型小倍数图（chromosome x within-chromosome position）
  FigD3_discovery_forest      15 个 FDR 显著对的森林图（按基因座分组）
  FigD4_discovery_celltype    细胞类型级负荷：名义显著数 / FDR 显著数 / 最小 qval
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import _fig_style as S

T = r'D:\endometriosis_project\11_sc_eqtl_mr_project\tables'
LOG = r'D:\endometriosis_project\_fig_disc_out.txt'
L = []


def rd(n):
    return pd.read_csv(os.path.join(T, n))


def main():
    S.setup()
    d = rd('10_discovery_MR_main.csv')
    d['chr_i'] = pd.to_numeric(d['chr'], errors='coerce')
    d = d[d['chr_i'].notna()].copy()
    d['chr_i'] = d['chr_i'].astype(int)
    d['neglogp'] = -np.log10(d['p'].clip(lower=1e-320))
    d['is_MHC'] = d['is_MHC'].astype(bool)
    sig = d[d['qval'] < 0.05].copy()
    p_thr = float(sig['p'].max())
    L.append('全表 %d 行 ; FDR<0.05 = %d ; BH 阈值对应 p = %.4g (-log10=%.3f)'
             % (len(d), len(sig), p_thr, -np.log10(p_thr)))

    # ---------------------------------------------------------------- FigD1
    fig, ax = plt.subplots(figsize=(4.6, 3.5))
    ns = d[~d['is_MHC'] & (d['p'] >= 0.05) & (d['qval'] >= 0.05)]
    nm = d[~d['is_MHC'] & (d['p'] < 0.05) & (d['qval'] >= 0.05)]
    mh = d[d['is_MHC']]
    ax.scatter(ns['b'], ns['neglogp'], s=3, c=S.C['ns'], lw=0, alpha=0.5, label='not significant', rasterized=True)
    ax.scatter(nm['b'], nm['neglogp'], s=3.5, c=S.C['nominal'], lw=0, alpha=0.75, label='nominal P < 0.05', rasterized=True)
    ax.scatter(mh['b'], mh['neglogp'], s=4, c=S.C['mhc'], lw=0, alpha=0.8, label='MHC region (chr6:25-35 Mb)')
    ax.scatter(sig['b'], sig['neglogp'], s=16, c=S.C['sig'], lw=0, marker='D',
               label='FDR < 0.05 (BH P threshold = %.2g)' % p_thr, zorder=5)
    ax.axhline(-np.log10(p_thr), ls='--', lw=0.6, c=S.C['grey'])
    ax.text(ax.get_xlim()[1] - 0.004, -np.log10(p_thr) - 0.45, 'FDR = 0.05',
            fontsize=6, color=S.C['grey'], va='top', ha='right')
    ax.axvline(0, ls='-', lw=0.5, c=S.C['light'])
    _lab_off = {'ANXA4': (5, 3), 'YME1L1': (5, 3)}
    for _, r in sig.sort_values('p').groupby('gene').head(2).iterrows():
        ax.annotate('%s (%s)' % (r['gene'], r['cell_type']), (r['b'], r['neglogp']),
                    textcoords='offset points',
                    xytext=_lab_off.get(r['gene'], (4, 3)),
                    fontsize=5.8, color=S.C['dark'])
    ax.set_xlabel('MR effect size (beta per SD of expression, log-odds of female infertility)')
    ax.set_ylabel('-log10(P)')
    ax.set_title('Discovery MR: blood immune-cell cis-eQTL vs female infertility', fontsize=8.5)
    ax.legend(loc='upper left', markerscale=2.2)
    S.finalize(ax)
    S.save(fig, 'FigD1_discovery_volcano')

    # ---------------------------------------------------------------- FigD2
    cells = sorted(d['cell_type'].unique())
    maxpos = d.groupby(['chr_i'])['tss_pos_grch38'].max()
    d['frac'] = d.apply(lambda r: (r['tss_pos_grch38'] / maxpos.get(r['chr_i'], np.nan))
                        if np.isfinite(r['tss_pos_grch38']) else np.nan, axis=1)
    d['xpos'] = d['chr_i'] + d['frac'].fillna(0.5) * 0.86

    fig, axes = plt.subplots(4, 4, figsize=(7.2, 5.4), sharey=True)
    axes = axes.ravel()
    for i, cell in enumerate(cells):
        ax = axes[i]
        g = d[d['cell_type'] == cell]
        ax.scatter(g.loc[~g['is_MHC'], 'xpos'], g.loc[~g['is_MHC'], 'neglogp'],
                   s=2.2, c=S.C['ns'], lw=0, alpha=0.65, rasterized=True)
        ax.scatter(g.loc[g['is_MHC'], 'xpos'], g.loc[g['is_MHC'], 'neglogp'],
                   s=2.6, c=S.C['mhc'], lw=0, alpha=0.85, rasterized=True)
        gs = g[g['qval'] < 0.05]
        if len(gs):
            ax.scatter(gs['xpos'], gs['neglogp'], s=12, c=S.C['sig'], lw=0, marker='D', zorder=5)
            thr = -np.log10(gs['p'].max())
        else:
            thr = -np.log10(0.05 / max(len(g), 1))
        ax.axhline(thr, ls='--', lw=0.5, c=S.C['grey'])
        ax.set_title('%s (n = %d%s)' % (cell, len(g), ', sig = %d' % len(gs) if len(gs) else ''),
                     fontsize=6.4, pad=2)
        ax.set_xticks([1, 6, 12, 18, 22])
        ax.set_xticklabels(['1', '6', '12', '18', '22'], fontsize=5.5)
        ax.tick_params(labelsize=5.5)
        S.finalize(ax)
    for j in range(len(cells), 16):
        axes[j].axis('off')
    for k in (12, 13, 14, 15):
        axes[k].set_xlabel('chromosome', fontsize=6)
    for k in range(0, 16, 4):
        axes[k].set_ylabel('-log10(P)', fontsize=6.5)
    fig.suptitle('Per-cell-type Manhattan (x = chromosome; position normalized within chromosome)',
                 fontsize=8.5, y=0.995)
    fig.tight_layout(rect=[0, 0, 1, 0.975])
    S.save(fig, 'FigD2_discovery_manhattan')

    # ---------------------------------------------------------------- FigD3
    s15 = rd('15_discovery_significant_with_locus.csv')
    loc_color = {'L1': S.C['locus1'], 'L2': S.C['locus2'], 'L3': S.C['locus3']}
    s15 = s15.sort_values(['locus_id', 'p'], ascending=[True, True]).reset_index(drop=True)
    n = len(s15)
    fig, ax = plt.subplots(figsize=(5.6, 3.6))
    for i, r in s15.iterrows():
        y = n - 1 - i
        c = loc_color.get(r['locus_id'], S.C['grey'])
        ax.errorbar(r['b'], y, xerr=1.96 * r['se'], fmt='D', ms=3.2, color=c,
                    ecolor=c, elinewidth=0.7, capsize=1.6, lw=0.7)
    ax.axvline(0, ls='-', lw=0.6, c=S.C['dark'])
    ax.set_yticks(range(n))
    ax.set_yticklabels(['%s | %s' % (r['gene'], r['cell_type']) for _, r in s15.iloc[::-1].iterrows()],
                       fontsize=6)
    ax.set_ylim(-0.8, n - 0.2)
    ax.set_xlabel('MR effect size (beta, 95% CI)')
    handles = [plt.Line2D([], [], marker='D', ls='', color=loc_color[k],
                          label='%s  %s' % (k, lab)) for k, lab in
               [('L1', 'chr1:22.03 Mb (CDC42/LINC00339)'),
                ('L2', 'chr10:27.13 Mb (YME1L1)'),
                ('L3', 'chr2:69.74 Mb (ANXA4)')]]
    ax.legend(handles=handles, loc='lower right', fontsize=5.8)
    ax.set_title('Significant gene x cell type pairs (FDR < 0.05), grouped by locus', fontsize=8.5)
    S.finalize(ax)
    S.save(fig, 'FigD3_discovery_forest')

    # ---------------------------------------------------------------- FigD4
    cs = rd('17_discovery_celltype_summary.csv').sort_values('qval_min')
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 2.9),
                             gridspec_kw={'width_ratios': [1.35, 1]})
    ax = axes[0]
    y = np.arange(len(cs))
    ax.barh(y, cs['n_nominal_main'], color=S.C['nominal'], height=0.62, label='nominal P < 0.05')
    ax.barh(y, cs['n_FDR_main'], color=S.C['sig'], height=0.62, label='FDR < 0.05')
    ax.set_yticks(y)
    ax.set_yticklabels(cs['cell_type'], fontsize=6.5)
    ax.invert_yaxis()
    ax.set_xlabel('number of gene x cell type pairs')
    ax.legend(loc='lower right')
    S.finalize(ax)

    ax = axes[1]
    ax.barh(y, -np.log10(cs['qval_min'].clip(lower=1e-320)), color=S.C['mhc'], height=0.62)
    ax.axvline(-np.log10(0.05), ls='--', lw=0.6, c=S.C['sig'])
    ax.text(-np.log10(0.05) + 0.15, len(cs) - 1.2, 'FDR = 0.05', fontsize=6, color=S.C['sig'])
    ax.set_yticks(y)
    ax.set_yticklabels([])
    ax.invert_yaxis()
    ax.set_xlabel('-log10(minimum q value)')
    S.finalize(ax)
    fig.suptitle('Cell-type burden of discovery signal', fontsize=8.5, y=0.99)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    S.save(fig, 'FigD4_discovery_celltype')

    L.append('图已输出到 %s' % S.FIGDIR)
    open(LOG, 'w', encoding='utf-8').write('\n'.join(L))


if __name__ == '__main__':
    main()
