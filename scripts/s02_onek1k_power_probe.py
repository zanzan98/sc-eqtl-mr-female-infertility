import zipfile, io, csv, collections, statistics, json

z = r'D:\endometriosis_project\11_sc_eqtl_mr_project\00_data_raw\onek1k\top_eQTL_summary.zip'
out = r'D:\endometriosis_project\11_sc_eqtl_mr_project\tables\00b_onek1k_top_power_probe.csv'
rep = r'D:\endometriosis_project\_onek1k_power.txt'

rows = []
with zipfile.ZipFile(z) as zf:
    for nm in zf.namelist():
        if nm.startswith('__MACOSX') or nm.endswith('/'):
            continue
        # filename: OneK1K_<CELLTYPE>.sig_cis_qtl_pairs.chr<N>.csv
        base = nm.split('/')[-1]
        cell = base.split('.sig_cis_qtl_pairs')[0].replace('OneK1K_', '')
        with zf.open(nm) as fh:
            txt = io.TextIOWrapper(fh, errors='replace')
            rd = csv.DictReader(txt, delimiter='\t')
            if rd.fieldnames is None or 'pval_nominal' not in rd.fieldnames:
                # try comma
                fh.seek(0)
                txt = io.TextIOWrapper(fh, errors='replace')
                rd = csv.DictReader(txt)
            for r in rd:
                try:
                    p = float(r['pval_nominal'])
                except Exception:
                    continue
                rows.append((cell, r.get('phenotype_id',''), r.get('variant_id',''),
                             p, float(r['slope']), float(r['slope_se']),
                             float(r.get('qval') or 'nan'), float(r.get('af') or 'nan'),
                             float(r.get('tss_distance') or 'nan')))

print("total rows:", len(rows), file=open(rep,'w',encoding='utf-8'))

by = collections.defaultdict(list)
for c, g, v, p, b, se, q, af, d in rows:
    by[c].append(p)

lines = []
lines.append(f"total rows = {len(rows):,}")
lines.append("")
hdr = f"{'celltype':12s} {'n_pairs':>8s} {'p<5e-8':>7s} {'p<1e-6':>7s} {'p<1e-5':>7s} {'p<0.05':>7s} {'median_p':>11s} {'min_p':>11s}"
lines.append(hdr)
lines.append("-"*len(hdr))
summary = []
for c in sorted(by):
    ps = sorted(by[c])
    n = len(ps)
    n8  = sum(1 for p in ps if p < 5e-8)
    n6  = sum(1 for p in ps if p < 1e-6)
    n5  = sum(1 for p in ps if p < 1e-5)
    n05 = sum(1 for p in ps if p < 5e-2)
    lines.append(f"{c:12s} {n:8d} {n8:7d} {n6:7d} {n5:7d} {n05:7d} {statistics.median(ps):11.3e} {ps[0]:11.3e}")
    summary.append(dict(cell_type=c, n_top_pairs=n, n_p_lt_5e8=n8, n_p_lt_1e6=n6,
                        n_p_lt_1e5=n5, n_p_lt_0p05=n05,
                        median_pval=statistics.median(ps), min_pval=ps[0]))

tot8 = sum(1 for _,_,_,p,_,_,_,_,_ in rows if p < 5e-8)
lines.append("")
lines.append(f"ALL cell types combined: pairs={len(rows):,}  p<5e-8 = {tot8:,}")
lines.append(f"distinct genes (any cell type) = {len(set(g for _,g,_,_,_,_,_,_,_ in rows)):,}")
lines.append(f"distinct genes with p<5e-8 = {len(set(g for _,g,_,p,_,_,_,_,_ in rows if p<5e-8)):,}")

with open(rep,'w',encoding='utf-8') as f:
    f.write("\n".join(lines))

import csv as _csv
with open(out,'w',newline='',encoding='utf-8-sig') as f:
    w = _csv.DictWriter(f, fieldnames=list(summary[0].keys()))
    w.writeheader()
    for s in summary: w.writerow(s)
print("ok")
