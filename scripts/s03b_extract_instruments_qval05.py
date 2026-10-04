#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""s03b_extract_instruments_qval05.py —— 敏感性分析工具变量集：TensorQTL 置换 qval < 0.05

与主口径 s03（pval_nominal < 5e-8）的区别：
  * 主口径 = 名义 P<5e-8（8,958 个 pair，8,958 个变异）
  * 本口径 = OneK1K 官方置换 FDR qval<0.05（~14,485 个 pair）
两者都取自同一个 top 表（每基因-细胞类型仅 1 个 top 变异），故均为单 SNP 工具变量。
"""
import zipfile, io, csv, os

ROOT = r'D:\endometriosis_project\11_sc_eqtl_mr_project'
Z = os.path.join(ROOT, r'00_data_raw\onek1k\top_eQTL_summary.zip')
OUT = os.path.join(ROOT, r'tables\01c_discovery_instruments_qval05.csv')
REP = r'D:\endometriosis_project\_qval_instruments_report.txt'
QTHR = 0.05
FMIN = 10.0

recs = []
with zipfile.ZipFile(Z) as zf:
    for nm in zf.namelist():
        if nm.startswith('__MACOSX') or nm.endswith('/'):
            continue
        base = nm.split('/')[-1]
        cell = base.split('.sig_cis_qtl_pairs')[0].replace('OneK1K_', '')
        with zf.open(nm) as fh:
            for r in csv.DictReader(io.TextIOWrapper(fh, errors='replace'), delimiter='\t'):
                try:
                    q = float(r['qval']); p = float(r['pval_nominal'])
                    b = float(r['slope']); se = float(r['slope_se'])
                except Exception:
                    continue
                if q >= QTHR or not (se > 0):
                    continue
                F = (b / se) ** 2
                if F < FMIN:
                    continue
                recs.append(dict(
                    cell_type=cell, gene=r['phenotype_id'], variant_id=r['variant_id'],
                    tss_distance=int(float(r['tss_distance'])),
                    af=float(r['af']), beta=b, se=se, z=b / se, F=F,
                    pval_nominal=p, qval=q,
                    pval_perm=float(r['pval_perm']), pval_beta=float(r['pval_beta']),
                    ma_samples=int(float(r['ma_samples'])), ma_count=int(float(r['ma_count'])),
                    num_var=int(float(r['num_var'])),
                ))

recs.sort(key=lambda x: (x['gene'], x['cell_type']))
os.makedirs(os.path.dirname(OUT), exist_ok=True)
with open(OUT, 'w', newline='', encoding='utf-8-sig') as f:
    w = csv.DictWriter(f, fieldnames=list(recs[0].keys())); w.writeheader()
    for r in recs:
        w.writerow(r)

import collections
rep = []
rep.append('输出: %s' % OUT)
rep.append('qval<0.05 且 F>10 工具变量数 = %d' % len(recs))
rep.append('涉及基因 = %d' % len(set(r['gene'] for r in recs)))
rep.append('细胞类型 = %d' % len(set(r['cell_type'] for r in recs)))
rep.append('F min/median = %.2f / %.1f' % (min(r['F'] for r in recs), sorted(r['F'] for r in recs)[len(recs) // 2]))
pps = sorted(r['pval_nominal'] for r in recs)
rep.append('pval_nominal min/median/max = %.3e / %.3e / %.3e' % (pps[0], pps[len(pps) // 2], pps[-1]))
n_lt5e8 = sum(1 for r in recs if r['pval_nominal'] < 5e-8)
rep.append('其中 pval_nominal<5e-8 的 = %d (%.1f%%)' % (n_lt5e8, 100.0 * n_lt5e8 / len(recs)))
rep.append('新增（仅在 qval 口径、不在 P<5e-8 口径）= %d' % (len(recs) - n_lt5e8))
rep.append('')
rep.append('== 每细胞类型工具变量数 ==')
cc = collections.Counter(r['cell_type'] for r in recs)
for k, v in sorted(cc.items(), key=lambda x: -x[1]):
    rep.append('  %-12s %5d' % (k, v))

open(REP, 'w', encoding='utf-8').write('\n'.join(rep))
print('ok')
