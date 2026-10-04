#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""s10_discovery_aggregate.py —— 把 s06 的逐变异 MR 结果聚合为 gene x cell_type 结果表

方法学：
  * nsnp == 1  → Wald ratio（单工具变量）
  * nsnp >= 2  → 逆方差加权 IVW（固定效应，主口径）+ Cochran Q / I²；
                 若 Q 显著（p<0.05），额外给出乘法随机效应 IVW-MRE 的 p 供稳健性判断
  * F_stat = 工具变量 F 均值（另报 F_min；单 SNP 时二者相等）
  * is_MHC = 基因 TSS（GRCh38）落在 chr6:25,000,000-35,000,000
             TSS = 变异 GRCh38 位置 - tss_distance（TensorQTL 口径）
  * qval   = BH-FDR（对全部有效 gene x cell_type 检验）
  * qval_noMHC = 剔除 MHC 后重算的 BH（MHC 行置 NaN，n 仍按全表记账 → 更保守）

输出列严格包含用户指定字段：gene, cell_type, method, nsnp, b, se, p, qval, is_MHC, F_stat
"""
import argparse
import os
import numpy as np
import pandas as pd
from scipy import stats

MHC_CHR = '6'
MHC_LO, MHC_HI = 25_000_000, 35_000_000


def bh(p):
    """Benjamini-Hochberg：m 取实际参与检验数（NaN 剔除），在 p 升序空间单调化后映回原序"""
    p = np.asarray(p, dtype=float)
    out = np.full(p.size, np.nan)
    ok = ~np.isnan(p)
    m = int(ok.sum())
    if m == 0:
        return out
    pv = p[ok]
    order = np.argsort(pv)
    ranked = pv[order]
    q = np.minimum.accumulate((ranked * m / np.arange(1, m + 1))[::-1])[::-1]
    tmp = np.empty_like(q)
    tmp[order] = np.clip(q, 0, 1)
    out[ok] = tmp
    return out


def agg_one(b, se):
    nsnp = len(b)
    if nsnp == 1:
        return dict(method='Wald ratio', nsnp=1, b=float(b[0]), se=float(se[0]),
                    p=float(2 * stats.norm.sf(abs(b[0] / se[0]))),
                    p_IVW_MRE=np.nan, Q=np.nan, Q_df=0, Q_p=np.nan, I2=np.nan)
    w = 1.0 / se ** 2
    sw = w.sum()
    bi = float((w * b).sum() / sw)
    sei = float(np.sqrt(1.0 / sw))
    p_fe = float(2 * stats.norm.sf(abs(bi / sei)))
    Q = float((w * (b - bi) ** 2).sum())
    df = nsnp - 1
    Qp = float(stats.chi2.sf(Q, df)) if df > 0 else np.nan
    I2 = float(max(0.0, (Q - df) / Q) * 100) if Q > 0 else 0.0
    if df > 0 and Q > df:
        scale = np.sqrt(Q / df)
        p_mre = float(2 * stats.norm.sf(abs(bi / (sei * scale))))
    else:
        p_mre = p_fe
    return dict(method='IVW', nsnp=nsnp, b=bi, se=sei, p=p_fe,
                p_IVW_MRE=p_mre, Q=Q, Q_df=df, Q_p=Qp, I2=I2)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--in', dest='inp', required=True, help='s06 逐变异输出')
    ap.add_argument('--instruments', required=True, help='配对用的工具变量表（GRCh38 版）')
    ap.add_argument('--out', dest='outp', required=True)
    ap.add_argument('--sig-out', dest='sig', default=None)
    ap.add_argument('--label', default='DISC_MAIN')
    a = ap.parse_args()

    d = pd.read_csv(a.inp, dtype={'variant_id': str, 'variant_id_grch37': str})
    print('逐变异表: %d 行' % len(d))
    print(d['status'].value_counts().to_string())
    if 'allele_source' in d:
        print('等位来源:', d['allele_source'].value_counts().to_dict())

    ok = d[d['p_MR'].notna()].copy()
    print('有效逐变异检验 = %d ; 唯一位点 = %d' % (len(ok), ok['variant_id'].nunique()))

    rows = []
    for (gene, cell), g in ok.groupby(['gene', 'cell_type'], sort=False):
        r = agg_one(g['beta_MR'].values, g['se_MR'].values)
        chr38 = g['variant_id'].str.split(':').str[0]
        r.update(dict(
            gene=gene, cell_type=cell,
            chr=str(chr38.mode().iloc[0]) if len(chr38.mode()) else None,
            variant_id_grch38=g['variant_id'].iloc[0],
            F_stat=float(g['F'].mean()), F_min=float(g['F'].min()),
            allele_src=str(g['allele_source'].mode().iloc[0]) if 'allele_source' in g and len(g['allele_source'].mode()) else None,
        ))
        rows.append(r)

    res = pd.DataFrame(rows)

    # s06 未透传 tss_distance / num_var → 从工具变量表按 (gene, cell_type) 回填
    iv = pd.read_csv(a.instruments, dtype={'variant_id': str, 'variant_id_grch37': str})
    ev = (iv.groupby(['gene', 'cell_type'])
            .agg(eQTL_pval_nominal=('pval_nominal', 'min'),
                 eQTL_qval=('qval', 'min'),
                 eQTL_af=('af', 'first'),
                 tss_distance=('tss_distance', 'first'),
                 n_var_tested=('num_var', 'first'))
            .reset_index())
    res = res.merge(ev, on=['gene', 'cell_type'], how='left')
    # TSS(GRCh38) = 变异 GRCh38 位置 - tss_distance（TensorQTL 口径）
    _p38 = pd.to_numeric(res['variant_id_grch38'].str.split(':').str[1], errors='coerce')
    res['tss_pos_grch38'] = _p38 - pd.to_numeric(res['tss_distance'], errors='coerce')

    res['z'] = res['b'] / res['se']
    res['p_bonf'] = np.minimum(1.0, res['p'] * len(res))
    res['qval'] = bh(res['p'].values)
    res['is_MHC'] = (res['chr'].astype(str) == MHC_CHR) & \
                    (res['tss_pos_grch38'] >= MHC_LO) & (res['tss_pos_grch38'] <= MHC_HI)
    res['qval_noMHC'] = bh(np.where(res['is_MHC'], np.nan, res['p'].values))
    res['sig_P_FDR'] = res['qval'] < 0.05
    res['sig_P_FDR_noMHC'] = res['qval_noMHC'] < 0.05
    res['sig_qval_eQTL'] = res['eQTL_qval'] < 0.05

    res = res.sort_values('p').reset_index(drop=True)
    order = ['gene', 'cell_type', 'method', 'nsnp', 'b', 'se', 'p', 'qval', 'is_MHC', 'F_stat',
             'sig_P_FDR', 'qval_noMHC', 'sig_P_FDR_noMHC', 'z', 'p_bonf',
             'eQTL_pval_nominal', 'eQTL_qval', 'sig_qval_eQTL', 'eQTL_af',
             'F_min', 'p_IVW_MRE', 'Q', 'Q_df', 'Q_p', 'I2',
             'chr', 'tss_pos_grch38', 'tss_distance', 'variant_id_grch38', 'n_var_tested', 'allele_src']
    res = res[[c for c in order if c in res.columns]]

    os.makedirs(os.path.dirname(a.outp), exist_ok=True)
    res.to_csv(a.outp, index=False, encoding='utf-8-sig')
    print('已写出 %s (%d 行)' % (a.outp, len(res)))

    sig = res[res['qval'] < 0.05]
    if a.sig:
        sig.to_csv(a.sig, index=False, encoding='utf-8-sig')
        print('P_FDR<0.05 = %d 已写出 %s' % (len(sig), a.sig))

    print('=' * 70)
    print('method:', res['method'].value_counts().to_dict())
    print('is_MHC:', res['is_MHC'].value_counts().to_dict())
    print('名义 p<0.05 = %d ; BH<0.05 = %d ; BH(剔MHC)<0.05 = %d ; Bonferroni<0.05 = %d'
          % (int((res['p'] < 0.05).sum()), int((res['qval'] < 0.05).sum()),
             int((res['qval_noMHC'] < 0.05).sum()), int((res['p_bonf'] < 0.05).sum())))
    print()
    print('== Top 25 ==')
    print(res.head(25).to_string(index=False))


if __name__ == '__main__':
    main()
