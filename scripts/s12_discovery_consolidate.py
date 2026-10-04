#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""s12_discovery_consolidate.py —— Discovery 阶段结果整合

产出四张表：
  15_discovery_significant_with_locus.csv  —— 显著对 + 基因座归属 + 敏感性口径是否保留
  16_discovery_locus_summary.csv           —— 基因座级汇总（独立信号数，避免用 pair 数夸大）
  17_discovery_celltype_summary.csv        —— 细胞类型级汇总
  18_replication_candidate_pool.csv        —— 进入 Replication 的候选池（含 GRCh37/38 双坐标与等位）
"""
import os
import numpy as np
import pandas as pd

ROOT = r'D:\endometriosis_project\11_sc_eqtl_mr_project'
T = os.path.join(ROOT, 'tables')
LOCUS_WIN = 1_000_000  # 同染色体 1 Mb 内视为同一基因座


def assign_locus(d):
    d = d.copy()
    d['_chr'] = d['chr'].astype(str)
    d['_pos'] = pd.to_numeric(d['tss_pos_grch38'], errors='coerce')
    d = d.sort_values(['_chr', '_pos']).reset_index(drop=True)
    lid = []
    cur = -1
    prev_c, prev_p = None, None
    for c, p in zip(d['_chr'], d['_pos']):
        if prev_c is None or c != prev_c or not np.isfinite(p) or not np.isfinite(prev_p) or (p - prev_p) > LOCUS_WIN:
            cur += 1
        lid.append('L%d' % (cur + 1))
        prev_c, prev_p = c, p
    d['locus_id'] = lid
    return d


def main():
    main_all = pd.read_csv(os.path.join(T, '10_discovery_MR_main.csv'))
    qval_all = pd.read_csv(os.path.join(T, '12_discovery_sensitivity_qval05.csv'))
    main_sig = main_all[main_all['qval'] < 0.05].copy()
    qval_sig = qval_all[qval_all['qval'] < 0.05].copy()

    key_main = set(zip(main_sig.gene, main_sig.cell_type))
    key_qval = set(zip(qval_sig.gene, qval_sig.cell_type))
    print('主口径显著 = %d ; qval口径显著 = %d ; 交集 = %d' % (len(key_main), len(key_qval), len(key_main & key_qval)))
    print('仅主口径有 =', sorted(key_main - key_qval))
    print('仅 qval 口径有 =', sorted(key_qval - key_main))

    s = assign_locus(main_sig)
    s['in_qval_sensitivity'] = [ (g, c) in key_qval for g, c in zip(s.gene, s.cell_type) ]
    s['locus_label'] = s.apply(lambda r: '%s:%.2fMb' % (r['chr'], r['tss_pos_grch38'] / 1e6)
                               if np.isfinite(r['tss_pos_grch38']) else 'NA', axis=1)
    # ★ 结局自身在该工具变量上的 p 值：区分「MR 显著且结局该位点达 GWS」vs「仅 MR 显著、结局为提示性」
    _vl = pd.read_csv(os.path.join(T, '_disc_main_variantlevel.csv'), dtype={'variant_id': str})
    _vl = _vl[['gene', 'cell_type', 'p_out', 'beta_out', 'se_out']].rename(
        columns={'p_out': 'p_outcome', 'beta_out': 'beta_outcome', 'se_out': 'se_outcome'})
    s = s.merge(_vl, on=['gene', 'cell_type'], how='left')
    s['outcome_GWS'] = s['p_outcome'] < 5e-8
    cols = ['locus_id', 'gene', 'cell_type', 'method', 'nsnp', 'b', 'se', 'z', 'p', 'qval',
            'qval_noMHC', 'is_MHC', 'F_stat', 'F_min', 'eQTL_pval_nominal', 'eQTL_qval',
            'p_outcome', 'outcome_GWS', 'beta_outcome', 'se_outcome',
            'in_qval_sensitivity', 'chr', 'locus_label', 'tss_pos_grch38', 'tss_distance',
            'variant_id_grch38', 'n_var_tested', 'allele_src']
    s = s[[c for c in cols if c in s.columns]]
    s.to_csv(os.path.join(T, '15_discovery_significant_with_locus.csv'), index=False, encoding='utf-8-sig')

    # 基因座级汇总
    ls = (s.groupby('locus_id')
            .agg(chr=('chr', 'first'), locus_label=('locus_label', 'first'),
                 genes=('gene', lambda x: '|'.join(sorted(set(x)))),
                 n_cell_types=('cell_type', 'nunique'),
                 n_pairs=('gene', 'size'),
                 p_min=('p', 'min'), qval_min=('qval', 'min'),
                 b_range=('b', lambda x: '%.3f~%.3f' % (x.min(), x.max())),
                 F_stat_min=('F_stat', 'min'))
            .reset_index().sort_values('p_min'))
    ls.to_csv(os.path.join(T, '16_discovery_locus_summary.csv'), index=False, encoding='utf-8-sig')

    # 细胞类型级汇总
    cs = []
    for cell, g in main_all.groupby('cell_type'):
        gq = qval_all[qval_all['cell_type'] == cell]
        cs.append(dict(cell_type=cell,
                       n_tests_main=len(g), n_nominal_main=int((g['p'] < 0.05).sum()),
                       n_FDR_main=int((g['qval'] < 0.05).sum()),
                       n_tests_qval=len(gq), n_FDR_qval=int((gq['qval'] < 0.05).sum()),
                       p_min=float(g['p'].min()), qval_min=float(g['qval'].min()),
                       F_stat_median=float(g['F_stat'].median()),
                       n_MHC_pairs=int(g['is_MHC'].sum())))
    csd = pd.DataFrame(cs).sort_values('qval_min')
    csd.to_csv(os.path.join(T, '17_discovery_celltype_summary.csv'), index=False, encoding='utf-8-sig')

    # 复现候选池：带 GRCh37 坐标 + 等位（等位取自 s06 逐变异表的 harmonise 结果）
    iv = pd.read_csv(os.path.join(T, '01b_discovery_instruments_grch38.csv'),
                     dtype={'variant_id': str, 'variant_id_grch37': str})
    vl = pd.read_csv(os.path.join(T, '_disc_main_variantlevel.csv'), dtype={'variant_id': str})
    vl = vl[['gene', 'cell_type', 'variant_id', 'effect_allele', 'other_allele',
             'beta_out', 'se_out', 'p_out', 'status', 'allele_source']]
    pool = s.merge(iv[['gene', 'cell_type', 'variant_id_grch37', 'af', 'beta', 'se',
                       'pval_nominal', 'qval', 'ma_samples', 'ma_count']],
                   on=['gene', 'cell_type'], how='left', suffixes=('', '_iv'))
    pool = pool.merge(vl, on=['gene', 'cell_type'], how='left', suffixes=('', '_vl'))
    pool.to_csv(os.path.join(T, '18_replication_candidate_pool.csv'), index=False, encoding='utf-8-sig')
    print('复现候选池 = %d 行 ; 缺失等位 = %d' % (len(pool), int(pool['effect_allele'].isna().sum())))

    print('独立基因座数 = %d' % s['locus_id'].nunique())
    print(ls.to_string(index=False))
    print()
    print(csd.to_string(index=False))


if __name__ == '__main__':
    main()
