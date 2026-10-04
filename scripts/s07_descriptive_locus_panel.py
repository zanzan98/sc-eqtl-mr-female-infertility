"""描述性位点面板：对已知内异症/不孕症 GWAS 基因 + 本项目关注基因，
报告其在 OneK1K 各免疫细胞类型中的最小 cis-eQTL P、qval、F、变异位点。
★ 用途：描述性位点重叠分析 —— 不作阳性对照门控判定（见方案 §2 问题 5 修正）。
"""
import zipfile, io, csv, collections, os

Z = r'D:\endometriosis_project\11_sc_eqtl_mr_project\00_data_raw\onek1k\top_eQTL_summary.zip'
OUT = r'D:\endometriosis_project\11_sc_eqtl_mr_project\tables\00d_descriptive_locus_panel.csv'
REP = r'D:\endometriosis_project\_locus_panel.txt'

# 已知内异症/不孕症 GWAS 位点基因（Open Targets FinnGen R12 可信集 L2G 首选基因）
GWAS_GENES = ['SYNE1','FSHB','KDR','CCDC170','CDC42','ESR1','GREB1','VEZT','ETAA1',
              'NFE2L3','PGR','CD109','CDKN2B','WNT4']
# 本项目主线关注基因
PROJECT_GENES = ['TRPA1','TRPV1','SCN9A','CALCA','CALCB','CALCRL','RAMP1','TACR1']
PANEL = GWAS_GENES + PROJECT_GENES

rows = []
with zipfile.ZipFile(Z) as zf:
    for nm in zf.namelist():
        if nm.startswith('__MACOSX') or nm.endswith('/'):
            continue
        cell = nm.split('/')[-1].split('.sig_cis_qtl_pairs')[0].replace('OneK1K_', '')
        with zf.open(nm) as fh:
            for r in csv.DictReader(io.TextIOWrapper(fh, errors='replace'), delimiter='\t'):
                g = r['phenotype_id']
                if g not in PANEL:
                    continue
                rows.append(dict(
                    gene=g, cell_type=cell, variant_id=r['variant_id'],
                    tss_distance=r['tss_distance'], af=r['af'],
                    beta=r['slope'], se=r['slope_se'],
                    F='%.1f' % ((float(r['slope'])/float(r['slope_se']))**2),
                    pval_nominal=r['pval_nominal'], qval=r['qval'],
                    pval_perm=r['pval_perm'], pval_beta=r['pval_beta'],
                    p_lt_5e8=str(float(r['pval_nominal']) < 5e-8),
                    q_lt_0p05=str(float(r['qval']) < 0.05),
                ))

# 每基因取 p 最小的那一行作为"最佳"
best = {}
for r in rows:
    g = r['gene']
    if g not in best or float(r['pval_nominal']) < float(best[g]['pval_nominal']):
        best[g] = r

rows.sort(key=lambda x: (PANEL.index(x['gene']), float(x['pval_nominal'])))
os.makedirs(os.path.dirname(OUT), exist_ok=True)
cols = list(rows[0].keys()) if rows else ['gene','cell_type']
with open(OUT, 'w', newline='', encoding='utf-8-sig') as f:
    w = csv.DictWriter(f, fieldnames=cols); w.writeheader()
    for r in rows: w.writerow(r)

L = []
L.append("=" * 100)
L.append("描述性位点面板 —— OneK1K top cis-eQTL（每基因每细胞类型 1 行 = 该基因最强 cis 变异）")
L.append("★ 本表为描述性分析，ESR1/FSHB/KDR 等未达 P<5e-8 属'暴露数据缺失'，不判定门控失败")
L.append("=" * 100)
L.append("")
L.append(f"{'gene':9s} {'出现的细胞类型数':>8s} {'最佳细胞类型':>12s} {'最佳p':>11s} {'qval':>10s} {'F':>8s} {'P<5e-8':>7s} {'q<0.05':>7s}  {'最佳变异':>14s}")
L.append("-" * 110)
for g in PANEL:
    sub = [r for r in rows if r['gene'] == g]
    if not sub:
        L.append(f"{g:9s} {'0 (不在表中)':>8s}")
        continue
    b = best[g]
    L.append(f"{g:9s} {len(sub):>8d} {b['cell_type']:>12s} {float(b['pval_nominal']):>11.3e} "
             f"{float(b['qval']):>10.3e} {b['F']:>8s} {b['p_lt_5e8']:>7s} {b['q_lt_0p05']:>7s}  {b['variant_id']:>14s}")

L.append("")
L.append("== 关键判断 ==")
for g in PANEL:
    sub = [r for r in rows if r['gene'] == g]
    n5 = sum(1 for r in sub if r['p_lt_5e8'] == 'True')
    nq = sum(1 for r in sub if r['q_lt_0p05'] == 'True')
    if not sub:
        L.append(f"  {g:9s} 未出现在 OneK1K top cis-eQTL 表中 → 该细胞类型集合内无可报告的 cis-eQTL")
    else:
        L.append(f"  {g:9s} 细胞类型={len(sub):2d}  P<5e-8 达标={n5:2d}  qval<0.05 达标={nq:2d}")

open(REP, 'w', encoding='utf-8').write("\n".join(L))
print("ok rows=", len(rows))
