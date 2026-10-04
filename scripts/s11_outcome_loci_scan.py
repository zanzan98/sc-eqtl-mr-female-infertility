#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""s11_outcome_loci_scan.py —— 结局文件自身的全基因组显著位点扫描（QC + 位点归属）

目的（两条，都很硬）：
  1) QC：确认 GCST90483463 确为真实的全基因组关联结果（存在一批 GWS 峰，且位点数目量级合理）
  2) 位点归属：把 Discovery MR 显著基因的 TSS 与结局自身 GWS 峰做重叠判定，
     从而区分「MR 命中恰在结局 GWS 峰内」vs「MR 命中不在任何 GWS 峰内」

算法：流式扫描 *.h.tsv.gz，保留 p<5e-8 的变异，按 1 Mb 窗口 (pos//1e6) 取窗口内最小 p 者。
输出：tables/14_outcome_gws_windows.csv
"""
import gzip, sys, os
import numpy as np
import pandas as pd

OUTCOME = r'D:\endometriosis_project\11_sc_eqtl_mr_project\00_data_raw\gwas\GCST90483463.h.tsv.gz'
OUT = r'D:\endometriosis_project\11_sc_eqtl_mr_project\tables\14_outcome_gws_windows.csv'
GWS = 5e-8
WIN = 1_000_000

ALIAS = {
    'chromosome': ['chromosome', 'chr', 'chrom'],
    'pos': ['base_pair_location', 'pos'],
    'ea': ['effect_allele', 'ea'],
    'oa': ['other_allele', 'oa'],
    'beta': ['beta', 'b'],
    'se': ['standard_error', 'se'],
    'eaf': ['effect_allele_frequency', 'eaf'],
    'p': ['p_value', 'pval', 'p'],
}


def fnum(x):
    try:
        return float(x)
    except Exception:
        return float('nan')


def main():
    best = {}
    n = 0
    ngws = 0
    with gzip.open(OUTCOME, 'rt', errors='replace') as fh:
        hdr = fh.readline().lstrip('#').rstrip('\n').split('\t')
        h = [c.strip().lower() for c in hdr]
        ix = {}
        for k, alts in ALIAS.items():
            for a in alts:
                if a in h:
                    ix[k] = h.index(a); break
        print('表头映射 =', ix, flush=True)
        for ln in fh:
            n += 1
            f = ln.rstrip('\n').split('\t')
            if len(f) <= max(ix.values()):
                continue
            p = fnum(f[ix['p']])
            if not (p < GWS):
                continue
            ngws += 1
            pos = int(float(f[ix['pos']]))
            chrom = f[ix['chromosome']].replace('chr', '')
            key = (chrom, pos // WIN)
            cur = best.get(key)
            if cur is None or p < cur['p']:
                best[key] = dict(chr=chrom, win_mb=pos // WIN,
                                 top_pos=pos, p=p,
                                 beta=fnum(f[ix['beta']]) if 'beta' in ix else np.nan,
                                 se=fnum(f[ix['se']]) if 'se' in ix else np.nan,
                                 ea=f[ix['ea']] if 'ea' in ix else '',
                                 oa=f[ix['oa']] if 'oa' in ix else '',
                                 eaf=fnum(f[ix['eaf']]) if 'eaf' in ix else np.nan)
            if n % 5_000_000 == 0:
                print('扫描 %d 行, GWS=%d, 窗口=%d' % (n, ngws, len(best)), flush=True)

    d = pd.DataFrame(list(best.values()))
    d = d.sort_values('p')
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    d.to_csv(OUT, index=False, encoding='utf-8-sig')
    print('总行数 = %d' % n)
    print('GWS 变异 = %d ; 1Mb 窗口 = %d' % (ngws, len(d)))
    print('已写出 %s' % OUT)
    print(d.head(40).to_string(index=False))


if __name__ == '__main__':
    main()
