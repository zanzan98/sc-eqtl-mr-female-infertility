#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""s22_p2_link_discovery.py —— raw(tar.gz) × Discovery 工具表 连接性实测

对每个 parquet 成员：读 7 列 → 筛出该 (cell_type, chr) 下的 Discovery 工具基因 →
与 01b/01d 逐条比对（键 = variant_id_grch37 chr:pos；校验 af / slope / se / tss_distance /
该基因在 raw 中的 argmin(pval_nominal) 是否即 Discovery 的 top 变异）。

★ 增量落盘 + --resume：任何中断均可继续（state 只记已完成成员，不存大对象）。
"""
import argparse, os, re, shutil, sys, time, json
import tarfile

import pandas as pd
import pyarrow.parquet as pq

PAT = re.compile(r'^\./OneK1K_eQTL_summary_Xue_et_al_bioRxiv_2024/OneK1K_(?P<ct>[A-Za-z0-9_]+)\.cis_qtl_pairs\.(?P<chr>chr[0-9]+)\.parquet$')

HDR = ('set,cell_type,chr37,gene,variant_id_grch37,variant_id_grch38,found,'
       'af_raw,af_top,d_af,slope_raw,beta_top,slope_se_raw,se_top,'
       'tss_dist_raw,tss_dist_top,pval_raw,pval_top,qval_top,'
       'is_raw_argmin_pval,raw_rows_for_gene,min_pval_raw_for_gene\n')


def _num(v):
    try:
        if v is None:
            return None
        return float(v)
    except Exception:
        return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--tar', required=True)
    ap.add_argument('--inst-b', required=True)
    ap.add_argument('--inst-d', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--state', required=True)
    ap.add_argument('--tmp', required=True)
    ap.add_argument('--time-cap', type=float, default=500.0)
    a = ap.parse_args()

    b = pd.read_csv(a.inst_b, dtype={'variant_id': str, 'variant_id_grch37': str})
    b['set'] = 'main_P5e8'
    d = pd.read_csv(a.inst_d, dtype={'variant_id': str, 'variant_id_grch37': str})
    d['set'] = 'sens_qval05'
    ins = pd.concat([b, d], ignore_index=True)
    ins['v37'] = ins['variant_id_grch37'].astype(str)
    ins['chr37'] = ins['v37'].str.split(':').str[0].astype(str)
    ins['pos37'] = ins['v37'].str.split(':').str[1].astype(str)
    ins['v37k'] = ins['chr37'] + ':' + ins['pos37']
    ins = ins[['set', 'cell_type', 'chr37', 'v37k', 'gene', 'variant_id', 'af', 'beta', 'se',
               'tss_distance', 'pval_nominal', 'qval']]
    print('instrument rows=%d (main=%d sens=%d)' % (
        len(ins), int((ins['set'] == 'main_P5e8').sum()), int((ins['set'] == 'sens_qval05').sum())))

    done = set()
    if os.path.isfile(a.state):
        with open(a.state, 'r', encoding='utf-8') as f:
            st = json.load(f)
        done = set(tuple(x) for x in st.get('done', []))
    else:
        with open(a.out, 'w', encoding='utf-8') as f:
            f.write(HDR)

    os.makedirs(a.tmp, exist_ok=True)
    tmpf = os.path.join(a.tmp, 'cur.parquet')
    t0 = time.time()
    n_new = 0
    n_skip = 0
    remaining = []
    try:
        tf = tarfile.open(a.tar, mode='r|gz')
        for m in tf:
            if not m.isfile():
                continue
            mm = PAT.match(m.name)
            if not mm:
                continue
            ct = mm.group('ct')
            ch = mm.group('chr')
            if (ct, ch) in done:
                n_skip += 1
                continue
            if time.time() - t0 > a.time_cap:
                remaining.append((ct, ch))
                continue
            fh = tf.extractfile(m)
            with open(tmpf, 'wb') as g:
                shutil.copyfileobj(fh, g, 1 << 22)
            try:
                fh.close()
            except Exception:
                pass
            cnum = ch.replace('chr', '')
            sub_ins = ins[(ins['cell_type'] == ct) & (ins['chr37'] == cnum)]
            out_rows = []
            if len(sub_ins):
                t = pq.read_table(tmpf, columns=['phenotype_id', 'variant_id', 'tss_distance', 'af',
                                                 'pval_nominal', 'slope', 'slope_se']).to_pandas()
                genes = set(sub_ins['gene'].astype(str))
                gsub = t[t['phenotype_id'].isin(genes)]
                if len(gsub):
                    idx = gsub.groupby('phenotype_id')['pval_nominal'].idxmin()
                    topm = {r['phenotype_id']: r for _, r in gsub.loc[idx].iterrows()}
                    cnt = gsub.groupby('phenotype_id').size().to_dict()
                    mn = gsub.groupby('phenotype_id')['pval_nominal'].min().to_dict()
                    mraw = gsub.set_index(['phenotype_id', 'variant_id'])
                else:
                    topm, cnt, mn, mraw = {}, {}, {}, None
                for _, r in sub_ins.iterrows():
                    g = str(r['gene'])
                    v = str(r['v37k'])
                    rr = None
                    if mraw is not None and (g, v) in mraw.index:
                        rr = mraw.loc[(g, v)]
                        if hasattr(rr, 'iloc') and getattr(rr, 'ndim', 1) == 2:
                            rr = rr.iloc[0]
                    afr, slr, ser = (None, None, None) if rr is None else (
                        _num(rr['af']), _num(rr['slope']), _num(rr['slope_se']))
                    tsr, pvr = (None, None) if rr is None else (
                        _num(rr['tss_distance']), _num(rr['pval_nominal']))
                    trow = topm.get(g)
                    out_rows.append([
                        r['set'], ct, cnum, g, v, str(r['variant_id']), int(rr is not None),
                        afr, _num(r['af']), (abs(afr - _num(r['af'])) if afr is not None and _num(r['af']) is not None else None),
                        slr, _num(r['beta']), ser, _num(r['se']),
                        tsr, _num(r['tss_distance']), pvr, _num(r['pval_nominal']), _num(r['qval']),
                        int(rr is not None and trow is not None and str(trow['variant_id']) == v),
                        int(cnt.get(g, 0)), (float(mn[g]) if g in mn else None),
                    ])
            with open(a.out, 'a', encoding='utf-8') as f:
                for row in out_rows:
                    f.write(','.join('' if x is None else str(x) for x in row) + '\n')
            done.add((ct, ch))
            n_new += 1
            st = {'done': [list(x) for x in done]}
            with open(a.state, 'w', encoding='utf-8') as f:
                json.dump(st, f, ensure_ascii=False)
            try:
                os.remove(tmpf)
            except Exception:
                pass
        try:
            tf.close()
        except Exception:
            pass
    except Exception as e:
        sys.stderr.write('TAR_ERR ' + repr(e) + '\n')
    finally:
        if os.path.exists(tmpf):
            try:
                os.remove(tmpf)
            except Exception:
                pass
    print('new=%d skipped=%d remaining_after_cap=%d elapsed=%.1f' %
          (n_new, n_skip, len(remaining), time.time() - t0))
    if remaining[:5]:
        print('next:', remaining[:5])


if __name__ == '__main__':
    main()
