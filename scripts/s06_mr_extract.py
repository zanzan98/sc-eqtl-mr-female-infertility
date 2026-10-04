#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
s06_mr_extract.py —— 通用 cis-eQTL → 结局 的 Wald ratio MR（流式抽取 + harmonise + FDR）

同一脚本用于：
  * 阳性对照门控（--genes GLIPR1,XBP1,IKZF1,IL2RA,...  --outcome GCST90132222 类风湿关节炎）
  * 主分析（不传 --genes，跑全部 8,958 个工具变量）

等位判定优先级：
  1) --bim-zip 给定 → 用 PLINK bim 的 A1 作为效应等位（OneK1K 官方口径：effect allele = A1）
  2) 未给定 → 用 AF 匹配启发式：比较 eQTL 的 af 与结局的 effect_allele_frequency，
     取 |af - af_out| 与 |af - (1-af_out)| 更小者；此路径会在输出中标注 allele_source=af_heuristic
"""
import argparse, csv, gzip, math, os, sys
import numpy as np
import pandas as pd
from scipy import stats

PALINDROMIC = {('A', 'T'), ('T', 'A'), ('C', 'G'), ('G', 'C')}
COMP = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

# GWAS-SSF / 常见表头的别名映射
ALIAS = {
    'chromosome': ['chromosome', 'chr', 'chrom', '#chrom', 'hm_chrom'],
    'pos': ['base_pair_location', 'pos', 'bp', 'position', 'hm_pos'],
    'ea': ['effect_allele', 'ea', 'a1', 'hm_effect_allele', 'alt', 'allele1'],
    'oa': ['other_allele', 'oa', 'a2', 'hm_other_allele', 'ref', 'allele0'],
    'beta': ['beta', 'b', 'effect', 'hm_beta', 'logOR'],
    'se': ['standard_error', 'se', 'stderr', 'sebeta', 'hm_se'],
    'eaf': ['effect_allele_frequency', 'eaf', 'freq', 'af', 'hm_effect_allele_frequency', 'af_alt'],
    'p': ['p_value', 'pval', 'p', 'p_value_neg_log10', 'hm_pval'],
    'n': ['n', 'samplesize', 'n_total'],
    'nc': ['n_case', 'n_cases', 'num_cases'],
    'nk': ['n_control', 'n_controls', 'num_controls'],
}


def map_header(hdr):
    h = [c.strip().lower().lstrip('#') for c in hdr]
    m = {}
    for canon, alts in ALIAS.items():
        for a in alts:
            if a in h:
                m[canon] = h.index(a)
                break
    return m


# GWAS Catalog harmonised 文件用 'NA' 表示缺失；直接 float() 会抛 ValueError。
NA_TOKENS = {'', 'na', 'n/a', 'nan', 'none', 'null', '.', '-'}


def sf(x):
    """容错数值解析：缺失/非数值一律返回 NaN，绝不抛异常。"""
    if x is None:
        return float('nan')
    s = str(x).strip()
    if s.lower() in NA_TOKENS:
        return float('nan')
    try:
        return float(s)
    except ValueError:
        return float('nan')


def load_bim(zip_path):
    import zipfile
    out = {}
    with zipfile.ZipFile(zip_path) as zf:
        for nm in zf.namelist():
            # 排除 macOS 资源叉（__MACOSX/._*.bim）：其内容为二进制元数据，非 PLINK bim
            if not nm.lower().endswith('.bim') or '__MACOSX' in nm:
                continue
            for ln in zf.read(nm).decode('utf-8', errors='replace').splitlines():
                f = ln.split()
                if len(f) < 6:
                    continue
                try:
                    bp = int(f[3])
                except ValueError:
                    continue
                out[(f[0], bp)] = (f[4].upper(), f[5].upper())
    return out


def harmonise(be, ee, oe, eafe, bo, eo, oo, eafo, band=0.08):
    ee, oe, eo, oo = ee.upper(), oe.upper(), eo.upper(), oo.upper()
    if not ee or not eo:
        return None, None, None, None, 'missing_allele'
    if (ee, oe) in PALINDROMIC and abs(eafe - 0.5) < band:
        return None, None, None, None, 'palindromic_ambiguous'
    if ee == eo and oe == oo:
        return be, ee, oe, bo, 'ok'
    if ee == oo and oe == eo:
        return be, ee, oe, -bo, 'flipped'
    if COMP.get(ee) == eo and COMP.get(oe) == oo:
        return be, ee, oe, bo, 'complement'
    if COMP.get(ee) == oo and COMP.get(oe) == eo:
        return be, ee, oe, -bo, 'complement_flipped'
    return None, None, None, None, 'allele_mismatch'


def wald(be, see, bo, seo):
    beta = bo / be
    se = math.sqrt(seo**2 / be**2 + (bo**2) * (see**2) / (be**4))
    z = beta / se if se > 0 else float('nan')
    return beta, se, (2 * stats.norm.sf(abs(z)) if se > 0 else float('nan'))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--instruments', required=True)
    ap.add_argument('--outcome-gz', required=True)
    ap.add_argument('--bim-zip', default=None)
    ap.add_argument('--genes', default=None, help='逗号分隔，只跑这些基因')
    ap.add_argument('--out', required=True)
    ap.add_argument('--label', default='outcome')
    args = ap.parse_args()

    iv = pd.read_csv(args.instruments, dtype={'variant_id': str, 'chr': str})
    if args.genes:
        gs = set(x.strip() for x in args.genes.split(','))
        iv = iv[iv['gene'].isin(gs)]
    # 丢弃位置缺失（如 liftover 无映射）的工具变量——否则 astype(int) 会因 NaN 失败
    n0 = len(iv)
    iv = iv[iv['variant_id'].astype(str).str.contains(':', na=False)].copy()
    if len(iv) < n0:
        print(f'[iv] 丢弃位置缺失的工具变量 {n0 - len(iv)} 行')
    iv['pos_i'] = iv['variant_id'].str.split(':').str[1].astype(int)
    iv['chr_s'] = iv['variant_id'].str.split(':').str[0]
    # ★ 双坐标系支持：结局按 variant_id(GRCh38) 匹配；PLINK bim 的 A1/A2 以原始（GRCh37）位置为键。
    #   若提供 variant_id_grch37 列，则 bim 查表用它，避免坐标系错配导致查表全落空。
    if 'variant_id_grch37' in iv.columns:
        _b37 = iv['variant_id_grch37'].astype(str)
        iv = iv[_b37.str.contains(':', na=False)].copy()
        _b37 = iv['variant_id_grch37'].astype(str)
        iv['bim_chr'] = _b37.str.split(':').str[0]
        iv['bim_pos'] = _b37.str.split(':').str[1].astype(int)
        print('[iv] 启用 variant_id_grch37 作为 bim 查表键（双坐标系模式）')
    else:
        iv['bim_chr'] = iv['chr_s']
        iv['bim_pos'] = iv['pos_i']
    want = {}
    for r in iv.itertuples(index=False):
        want.setdefault((r.chr_s, r.pos_i), []).append(r)
    print(f'[iv] 行={len(iv):,} 基因={iv["gene"].nunique():,} 唯一位点={len(want):,}')

    bim = load_bim(args.bim_zip) if args.bim_zip and os.path.exists(args.bim_zip) else {}

    print(f'[outcome] 扫描 {args.outcome_gz}')
    hits, hdr_seen = [], None
    with gzip.open(args.outcome_gz, 'rt', errors='replace') as fh:
        first = fh.readline()
        hdr = first.lstrip('#').rstrip('\n').split('\t')
        hdr_seen = '|'.join(hdr)
        ix = map_header(hdr)
        print(f'[outcome] 表头映射 = {ix}')
        need = {'chromosome', 'pos', 'ea', 'oa', 'beta', 'se'}
        if not need.issubset(ix):
            sys.exit(f'FATAL: 表头缺字段 {need - set(ix)}；实际表头：{hdr_seen}')
        offset = 1 if not first.startswith('#') else 1
        for ln in fh:
            f = ln.rstrip('\n').split('\t')
            if len(f) <= max(ix.values()):
                continue
            c = f[ix['chromosome']].replace('chr', '')
            try:
                p = int(float(f[ix['pos']]))
            except ValueError:
                continue
            k = (c, p)
            if k in want:
                hits.append(dict(chr=c, pos=p,
                                 ea_out=f[ix['ea']], oa_out=f[ix['oa']],
                                 beta_out=sf(f[ix['beta']]),
                                 se_out=sf(f[ix['se']]),
                                 eaf_out=sf(f[ix['eaf']]) if 'eaf' in ix else float('nan'),
                                 p_out=sf(f[ix['p']]) if 'p' in ix else float('nan')))
    print(f'[outcome] 命中位点行数 = {len(hits):,}')

    hdf = pd.DataFrame(hits)
    recs = []
    for k, rows in want.items():
        sub = hdf[(hdf['chr'] == k[0]) & (hdf['pos'] == k[1])] if not hdf.empty else hdf
        if sub.empty:
            for r in rows:
                recs.append(dict(gene=r.gene, cell_type=r.cell_type, variant_id=r.variant_id,
                                 status='outcome_missing'))
            continue
        o = sub.iloc[0]
        if o.beta_out != o.beta_out or o.se_out != o.se_out:
            # 结局该位点的 beta/se 为 NA（GWAS Catalog 常见）→ 单列标记，不混入 outcome_missing
            for r in rows:
                recs.append(dict(gene=r.gene, cell_type=r.cell_type, variant_id=r.variant_id,
                                 status='outcome_na'))
            continue
        for r in rows:
            bk = (r.bim_chr, r.bim_pos)
            if bk in bim:
                ea, oa = bim[bk]; src = 'bim_a1'
            else:
                # AF 启发式：选与 eQTL af 更一致的取向
                if o.eaf_out != o.eaf_out:
                    # 结局缺 EAF（NA）→ 无法判定取向，按结局效应等位对齐并在输出中显式标注
                    ea, oa = o.ea_out, o.oa_out
                    src = 'af_unavailable_ambiguous'
                else:
                    d1 = abs(r.af - o.eaf_out)
                    d2 = abs(r.af - (1 - o.eaf_out))
                    if d1 <= d2:
                        ea, oa = o.ea_out, o.oa_out
                    else:
                        ea, oa = o.oa_out, o.ea_out
                    src = 'af_heuristic'
            be, ee, oe, bo, st = harmonise(r.beta, ea, oa, r.af,
                                           o.beta_out, o.ea_out, o.oa_out, o.eaf_out)
            if be is None:
                recs.append(dict(gene=r.gene, cell_type=r.cell_type, variant_id=r.variant_id,
                                 status=st, allele_source=src))
                continue
            b, s, p = wald(be, r.se, bo, o.se_out)
            recs.append(dict(gene=r.gene, cell_type=r.cell_type, variant_id=r.variant_id,
                             effect_allele=ee, other_allele=oe,
                             beta_exp=r.beta, se_exp=r.se, F=r.F,
                             beta_out=bo, se_out=o.se_out, p_out=o.p_out,
                             beta_MR=b, se_MR=s, p_MR=p, z_MR=b/s if s > 0 else np.nan,
                             status=st, allele_source=src,
                             outcome=args.label))
    res = pd.DataFrame(recs)
    ok = res['p_MR'].notna() if 'p_MR' in res else pd.Series(dtype=bool)
    n = int(ok.sum())
    res['p_FDR'] = np.nan
    if n:
        pv = res.loc[ok, 'p_MR'].values
        order = np.argsort(pv); ranked = pv[order]
        q = np.minimum.accumulate((ranked * n / np.arange(1, n + 1))[::-1])[::-1]
        tmp = np.empty_like(q); tmp[order] = np.clip(q, 0, 1)
        res.loc[ok, 'p_FDR'] = tmp
    res = res.sort_values('p_MR') if 'p_MR' in res else res
    res.to_csv(args.out, index=False, encoding='utf-8-sig')
    print('=' * 60)
    print(f'[out] {args.out}')
    print(f'有效检验 = {n:,} / {len(res):,}')
    if 'status' in res:
        print(res['status'].value_counts().to_string())
    if n:
        sig = res[res['p_FDR'] < 0.05]
        print(f'FDR<0.05 = {len(sig):,}')
        print(sig.head(20).to_string(index=False))


if __name__ == '__main__':
    main()
