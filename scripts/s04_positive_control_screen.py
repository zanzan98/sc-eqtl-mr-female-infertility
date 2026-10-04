"""在 OneK1K P<5e-8 工具变量表中筛选候选阳性对照基因，并统计其跨细胞类型强度。"""
import csv, collections, os

IN = r'D:\endometriosis_project\11_sc_eqtl_mr_project\tables\01_discovery_instruments.csv'
OUT = r'D:\endometriosis_project\11_sc_eqtl_mr_project\tables\00c_positive_control_screen.csv'
REP = r'D:\endometriosis_project\_pc_screen.txt'

# 候选池：用户指定 + 免疫/PCOS 经典基因 + 本项目关注基因
USER_PC = ['GLIPR1','XBP1','IL2RA','PTPN22','IKZF1']
IMMUNE_CLASSIC = ['ORMDL3','GSDMB','SH2B3','ATXN2','IL18RAP','CD226','IRF5','STAT4',
                  'TNFSF13B','ERAP1','TYK2','JAZF1','IL23R','CCR6','BACH2','PRDM1',
                  'CD40','TNFRSF14','FCGR2A','SUOX','DENND1A','THADA','FSHR','LHCGR',
                  'INSR','YAP1','RAB5B','HMGA2','TOX3','AMH','AMHR2','KISS1','LEPR']
PROJECT = ['TRPA1','CALCA','CALCB','CALCRL','RAMP1','TACR1','SCN9A','TRPV1',
           'ESR1','FSHB','KDR','GREB1','VEZT','WNT4','PGR','CDKN2B','SYNE1','CDC42']

CAND = USER_PC + IMMUNE_CLASSIC + PROJECT

per_gene = collections.defaultdict(list)   # gene -> list of rows
nsnp = collections.defaultdict(set)

with open(IN, encoding='utf-8-sig') as f:
    for r in csv.DictReader(f):
        r['pval_nominal'] = float(r['pval_nominal']); r['F'] = float(r['F'])
        r['beta'] = float(r['beta']); r['se'] = float(r['se'])
        per_gene[r['gene']].append(r)
        nsnp[r['gene']].add(r['variant_id'])

rows_out = []
for g in CAND:
    hits = per_gene.get(g, [])
    if not hits:
        rows_out.append(dict(gene=g, group=('user' if g in USER_PC else ('immune' if g in IMMUNE_CLASSIC else 'project')),
                             n_celltypes_p5e8=0, best_celltype='', best_snp='', best_beta='', best_se='',
                             best_F='', best_p='', n_distinct_snp=0, n_celltypes_total_tested=None,
                             qualifies_gating=False))
        continue
    b = min(hits, key=lambda x: x['pval_nominal'])
    rows_out.append(dict(gene=g, group=('user' if g in USER_PC else ('immune' if g in IMMUNE_CLASSIC else 'project')),
                         n_celltypes_p5e8=len(hits), best_celltype=b['cell_type'], best_snp=b['variant_id'],
                         best_beta=b['beta'], best_se=b['se'], best_F=b['F'], best_p=b['pval_nominal'],
                         n_distinct_snp=len(nsnp[g]), n_celltypes_total_tested=None,
                         qualifies_gating=True))

cols = ['gene','group','n_celltypes_p5e8','best_celltype','best_snp','best_beta','best_se','best_F','best_p',
        'n_distinct_snp','n_celltypes_total_tested','qualifies_gating']
os.makedirs(os.path.dirname(OUT), exist_ok=True)
with open(OUT,'w',newline='',encoding='utf-8-sig') as f:
    w = csv.DictWriter(f, fieldnames=cols); w.writeheader()
    for r in rows_out: w.writerow(r)

L = []
L.append("== 用户指定阳性对照候选（GLIPR1/XBP1/IL2RA/PTPN22/IKZF1）==")
for g in USER_PC:
    h = per_gene.get(g, [])
    if h:
        b = min(h, key=lambda x: x['pval_nominal'])
        L.append(f"  {g:9s} 合格  达标细胞类型数={len(h):2d}  最强: {b['cell_type']:9s} {b['variant_id']:14s} "
                 f"beta={b['beta']:+.3f} se={b['se']:.4f} F={b['F']:.1f} p={b['pval_nominal']:.2e}")
    else:
        L.append(f"  {g:9s} **不合格**（无 P<5e-8 cis-eQTL）")

L.append("")
L.append("== 免疫经典基因（可供补充阳性对照）==")
for g in IMMUNE_CLASSIC:
    h = per_gene.get(g, [])
    if h:
        b = min(h, key=lambda x: x['pval_nominal'])
        L.append(f"  {g:9s} 合格  达标细胞类型数={len(h):2d}  最强: {b['cell_type']:9s} {b['variant_id']:14s} "
                 f"beta={b['beta']:+.3f} se={b['se']:.4f} F={b['F']:.1f} p={b['pval_nominal']:.2e}")

L.append("")
L.append("== 本项目关注基因（描述性位点重叠用，不作门控）==")
for g in PROJECT:
    h = per_gene.get(g, [])
    if h:
        b = min(h, key=lambda x: x['pval_nominal'])
        L.append(f"  {g:9s} 达标细胞类型数={len(h):2d}  最强: {b['cell_type']:9s} p={b['pval_nominal']:.2e}")
    else:
        L.append(f"  {g:9s} 无 P<5e-8 cis-eQTL")

L.append("")
L.append("== 参考：工具变量表中最强/覆盖最广的基因（按达标细胞类型数排序，前 25）==")
rank = sorted(per_gene.items(), key=lambda kv: -len(kv[1]))[:25]
for g, h in rank:
    b = min(h, key=lambda x: x['pval_nominal'])
    L.append(f"  {g:12s} 细胞类型数={len(h):2d}  最强 {b['cell_type']:9s} p={b['pval_nominal']:.2e} F={b['F']:.0f}")

open(REP,'w',encoding='utf-8').write("\n".join(L))
print("ok")
