"""从 OneK1K top cis-eQTL 表抽取 P<5e-8 工具变量，输出 MR 输入表。
字段说明（源文件列）：
  phenotype_id  = gene symbol
  variant_id    = chr:pos (OneK1K)  ★ 坐标系为 GRCh37 (hg19)，非 GRCh38 —— 见 00b_执行期更正_2026-09-22.md §1
  slope         = beta（每等位基因表达效应，标准化后）
  slope_se      = SE
  pval_nominal  = 名义 p
  qval          = 置换 FDR
  af            = 等位频率（注意：等位定义须用 plink bim 的 A1 核对）
"""
import zipfile, io, csv, os

Z = r'D:\endometriosis_project\11_sc_eqtl_mr_project\00_data_raw\onek1k\top_eQTL_summary.zip'
OUT = r'D:\endometriosis_project\11_sc_eqtl_mr_project\tables\01_discovery_instruments.csv'
PC_GENES = ['ESR1','FSHB','SYNE1','KDR','CDC42','GREB1','VEZT','WNT4','PGR','CDKN2B',
            'TRPA1','CALCA','CALCB','CALCRL','RAMP1','TACR1','SCN9A','TRPV1']
THR = 5e-8

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
                    p = float(r['pval_nominal'])
                except Exception:
                    continue
                if p >= THR:
                    continue
                b = float(r['slope']); se = float(r['slope_se'])
                recs.append(dict(
                    cell_type=cell, gene=r['phenotype_id'], variant_id=r['variant_id'],
                    tss_distance=int(float(r['tss_distance'])),
                    af=float(r['af']), beta=b, se=se, z=b/se, F=(b/se)**2,
                    pval_nominal=p, qval=float(r['qval']),
                    pval_perm=float(r['pval_perm']), pval_beta=float(r['pval_beta']),
                    ma_samples=int(float(r['ma_samples'])), ma_count=int(float(r['ma_count'])),
                    num_var=int(float(r['num_var'])),
                ))

recs.sort(key=lambda x: (x['gene'], x['cell_type']))
cols = list(recs[0].keys())
os.makedirs(os.path.dirname(OUT), exist_ok=True)
with open(OUT, 'w', newline='', encoding='utf-8-sig') as f:
    w = csv.DictWriter(f, fieldnames=cols); w.writeheader()
    for r in recs: w.writerow(r)

rep = []
rep.append(f"输出: {OUT}")
rep.append(f"P<5e-8 工具变量数 = {len(recs):,}")
rep.append(f"涉及基因 = {len(set(r['gene'] for r in recs)):,}")
rep.append(f"细胞类型 = {len(set(r['cell_type'] for r in recs)):,}")
rep.append(f"F 最小值 = {min(r['F'] for r in recs):.2f}  (全部 >10? {all(r['F']>10 for r in recs)})")
rep.append(f"F 中位数 = {sorted(r['F'] for r in recs)[len(recs)//2]:.1f}")
rep.append(f"|tss_distance| 最大值 = {max(abs(r['tss_distance']) for r in recs):,} bp")
rep.append("")
rep.append("== 阳性对照 / 本项目关注基因是否在 P<5e-8 集合内 ==")
gset = set(r['gene'] for r in recs)
for g in PC_GENES:
    hit = [r for r in recs if r['gene'] == g]
    if hit:
        rep.append(f"  {g:8s} 命中 {len(hit):3d} 个细胞类型: " +
                   ", ".join(f"{h['cell_type']}(p={h['pval_nominal']:.2e})" for h in hit))
    else:
        rep.append(f"  {g:8s} 无 P<5e-8 cis-eQTL")
rep.append("")
rep.append("== 每细胞类型工具变量数 ==")
import collections
cc = collections.Counter(r['cell_type'] for r in recs)
for k, v in sorted(cc.items(), key=lambda x: -x[1]):
    rep.append(f"  {k:10s} {v:5d}")

open(r'D:\endometriosis_project\_instruments_report.txt','w',encoding='utf-8').write("\n".join(rep))
print("ok")
