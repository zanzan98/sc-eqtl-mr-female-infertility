#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""s21_p2_member_profile.py —— 逐成员（cell_type × chr）画像，带断点续跑

对 308 个 parquet 成员逐个：流式抽取 → 读元数据与 3 列 → 追加写 CSV → 删除临时文件。
★ 增量落盘 + --resume，任何中断都可继续，不重跑已完成成员。
"""
import argparse, os, re, shutil, sys, time
import tarfile

import pandas as pd
import pyarrow.parquet as pq

PAT = re.compile(r'^\./OneK1K_eQTL_summary_Xue_et_al_bioRxiv_2024/OneK1K_(?P<ct>[A-Za-z0-9_]+)\.cis_qtl_pairs\.(?P<chr>chr[0-9]+)\.parquet$')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--tar', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--tmp', required=True)
    ap.add_argument('--time-cap', type=float, default=540.0)
    ap.add_argument('--resume', action='store_true')
    a = ap.parse_args()

    done = set()
    if a.resume and os.path.isfile(a.out):
        try:
            d = pd.read_csv(a.out)
            done = set(zip(d['cell_type'].astype(str), d['chr'].astype(str)))
        except Exception:
            done = set()
    else:
        with open(a.out, 'w', encoding='utf-8') as f:
            f.write('cell_type,chr,file_size,rows,n_genes,n_variants,max_vars_per_gene,'
                    'median_vars_per_gene,n_rows_pval_lt_5e8,n_genes_pval_lt_5e8,'
                    'min_af,max_af,sec\n')

    os.makedirs(a.tmp, exist_ok=True)
    tmpf = os.path.join(a.tmp, 'cur.parquet')
    t0 = time.time()
    n_done_new = 0
    n_skip = 0
    remaining = []
    t_start_all = time.time()
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
            if time.time() - t_start_all > a.time_cap:
                remaining.append((ct, ch))
                continue
            t1 = time.time()
            fh = tf.extractfile(m)
            with open(tmpf, 'wb') as g:
                shutil.copyfileobj(fh, g, 1 << 22)
            try:
                fh.close()
            except Exception:
                pass
            rec = [ct, ch, os.path.getsize(tmpf)]
            try:
                pf = pq.ParquetFile(tmpf)
                rows = int(pf.metadata.num_rows)
                t = pf.read(columns=['phenotype_id', 'variant_id', 'af', 'pval_nominal']).to_pandas()
                ng = int(t['phenotype_id'].nunique())
                nv = int(t['variant_id'].nunique())
                vc = t.groupby('phenotype_id')['variant_id'].nunique()
                hit = t['pval_nominal'] < 5e-8
                nh = int(hit.sum())
                ngh = int(t.loc[hit, 'phenotype_id'].nunique())
                rec += [rows, ng, nv,
                        int(vc.max()) if len(vc) else 0,
                        float(vc.median()) if len(vc) else 0.0,
                        nh, ngh,
                        float(t['af'].min()), float(t['af'].max()),
                        round(time.time() - t1, 2)]
            except Exception as e:
                rec += ['ERR:' + repr(e)] + [''] * 8
            with open(a.out, 'a', encoding='utf-8') as f:
                f.write(','.join(str(x) for x in rec) + '\n')
            n_done_new += 1
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
          (n_done_new, n_skip, len(remaining), time.time() - t0))
    if remaining[:5]:
        print('next:', remaining[:5])


if __name__ == '__main__':
    main()
