#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""s20_p2_extract_and_profile.py —— 暂停点 2：从 tar.gz 流式抽取样本成员并画像

不做全量解压：仅抽取 README / plink_merged.bim / 4 个 parquet 样本到临时目录，
画像后删除。产出 JSON 证据文件。
"""
import argparse, json, os, shutil, sys, time
import tarfile

import pandas as pd
import pyarrow.parquet as pq

GRCH38 = {1:248956422,2:242193529,3:198295559,4:190214555,5:181538259,6:170805979,
          7:159345973,8:145138636,9:138394717,10:133797422,11:135086622,12:133275309,
          13:114364328,14:107043718,15:101991189,16:90338345,17:83257441,18:80373285,
          19:58617616,20:64444167,21:46709983,22:50818468}
GRCH37 = {1:249250621,2:243199373,3:198022430,4:191154276,5:180915260,6:171115067,
          7:159138663,8:146364022,9:141213431,10:135534747,11:135006516,12:133851895,
          13:115169878,14:107349540,15:102531392,16:90354753,17:81195210,18:78077248,
          19:59128983,20:63025520,21:48129895,22:51304566}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--tar', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--tmp', required=True)
    ap.add_argument('--keep-tmp', action='store_true')
    ap.add_argument('--reuse', action='store_true', help='临时目录已有样本时跳过解压')
    a = ap.parse_args()

    base = './OneK1K_eQTL_summary_Xue_et_al_bioRxiv_2024/'
    targets = {
        base + 'README': 'README',
        base + 'plink_merged.bim': 'plink_merged.bim',
        base + 'OneK1K_Mono_NC.cis_qtl_pairs.chr15.parquet': 'p1_Mono_NC_chr15.parquet',
        base + 'OneK1K_CD4_NC.cis_qtl_pairs.chr18.parquet': 'p2_CD4_NC_chr18.parquet',
        base + 'OneK1K_DC.cis_qtl_pairs.chr21.parquet': 'p3_DC_chr21.parquet',
        base + 'OneK1K_CD4_SOX4.cis_qtl_pairs.chr2.parquet': 'p4_CD4_SOX4_chr2.parquet',
    }
    if os.path.isdir(a.tmp):
        shutil.rmtree(a.tmp, ignore_errors=True)
    os.makedirs(a.tmp, exist_ok=True)

    R = {'tar': a.tar}
    done = {}
    t0 = time.time()
    n_seen = 0
    decomp = 0
    need = [v for v in targets.values() if v != 'README']
    if a.reuse and all(os.path.isfile(os.path.join(a.tmp, v)) for v in need):
        R['reused_tmp'] = True
        for k, v in targets.items():
            p = os.path.join(a.tmp, v)
            if os.path.isfile(p):
                done[k] = {'path': p, 'size': os.path.getsize(p), 'tar_index': None}
        R['extract_seconds'] = 0.0
        R['members_traversed'] = None
        R['decompressed_bytes_traversed'] = None
    else:
        if os.path.isdir(a.tmp):
            shutil.rmtree(a.tmp, ignore_errors=True)
        os.makedirs(a.tmp, exist_ok=True)
        try:
            tf = tarfile.open(a.tar, mode='r|gz')
            for m in tf:
                n_seen += 1
                if m.isfile():
                    decomp += max(m.size, 0)
                if not m.isfile() or m.name not in targets:
                    if len(done) == len(targets):
                        break
                    continue
                dst = os.path.join(a.tmp, targets[m.name])
                fh = tf.extractfile(m)
                with open(dst, 'wb') as g:
                    shutil.copyfileobj(fh, g, 1 << 22)
                try:
                    fh.close()
                except Exception:
                    pass
                done[m.name] = {'path': dst, 'size': os.path.getsize(dst), 'tar_index': n_seen}
                if len(done) == len(targets):
                    break
            try:
                tf.close()
            except Exception:
                pass
        except Exception as e:
            R['extract_error'] = repr(e)
        R['extract_seconds'] = round(time.time() - t0, 2)
        R['members_traversed'] = n_seen
        R['decompressed_bytes_traversed'] = decomp
    R['extracted'] = done

    # ---------- README ----------
    try:
        with open(os.path.join(a.tmp, 'README'), 'rb') as f:
            raw = f.read()
        R['README_text'] = raw.decode('utf-8', 'replace')
    except Exception as e:
        R['README_text'] = 'ERR ' + repr(e)

    # ---------- bim ----------
    try:
        bim = pd.read_csv(os.path.join(a.tmp, 'plink_merged.bim'), sep=r'\s+', header=None,
                          names=['chr', 'variant_id', 'cm', 'pos', 'a1', 'a2'],
                          dtype={'chr': str, 'variant_id': str}, engine='c')
        bim['pos'] = pd.to_numeric(bim['pos'], errors='coerce')
        g = bim.groupby('chr')['pos']
        per = {k: {'n': int(v), 'min': int(g.min()[k]), 'max': int(g.max()[k])} for k, v in g.size().items()}
        R['bim'] = {
            'n_rows': int(len(bim)),
            'n_chr': int(bim['chr'].nunique()),
            'n_unique_variant_id': int(bim['variant_id'].nunique()),
            'n_rsid': int(bim['variant_id'].str.startswith('rs').sum()),
            'per_chr': per,
            'sample_rows': bim.head(3).to_dict('records'),
            'alleles_a1_sample': {str(k): int(v) for k, v in bim['a1'].value_counts().head(6).items()},
        }
        over38 = []
        over37 = []
        for k, d in per.items():
            try:
                ci = int(k)
            except Exception:
                continue
            if ci in GRCH38 and d['max'] > GRCH38[ci]:
                over38.append(k)
            if ci in GRCH37 and d['max'] > GRCH37[ci]:
                over37.append(k)
        R['bim']['chr_exceeding_GRCh38_length'] = sorted(over38, key=lambda x: int(x))
        R['bim']['n_chr_exceeding_GRCh38'] = len(over38)
        R['bim']['chr_exceeding_GRCh37_length'] = sorted(over37, key=lambda x: int(x))
        R['bim']['n_chr_exceeding_GRCh37'] = len(over37)
    except Exception as e:
        R['bim'] = {'error': repr(e)}

    # ---------- parquet ----------
    P = {}
    for key, fn in [('p1_Mono_NC_chr15', 'p1_Mono_NC_chr15.parquet'),
                    ('p2_CD4_NC_chr18', 'p2_CD4_NC_chr18.parquet'),
                    ('p3_DC_chr21', 'p3_DC_chr21.parquet'),
                    ('p4_CD4_SOX4_chr2', 'p4_CD4_SOX4_chr2.parquet')]:
        fp = os.path.join(a.tmp, fn)
        d = {}
        try:
            pf = pq.ParquetFile(fp)
            d['file'] = fn
            d['file_size'] = os.path.getsize(fp)
            d['n_rows'] = int(pf.metadata.num_rows)
            d['n_row_groups'] = int(pf.metadata.num_row_groups)
            d['n_cols'] = int(pf.metadata.num_columns)
            d['schema'] = [{'name': f.name, 'type': str(f.type)} for f in pf.schema_arrow]
            d['col_names'] = pf.schema_arrow.names
            d['created_by'] = str(pf.metadata.created_by)
            d['kv_metadata'] = {(k.decode('utf-8', 'replace') if isinstance(k, bytes) else str(k)):
                                (v.decode('utf-8', 'replace') if isinstance(v, bytes) else str(v))
                                for k, v in (pf.metadata.metadata or {}).items()}
            t = pf.read_row_group(0).slice(0, 3).to_pandas()
            d['head3'] = json.loads(t.to_json(orient='records', date_format='iso'))
        except Exception as e:
            d['error'] = repr(e)
        P[key] = d
    R['parquet'] = P

    # ---------- 抽样成员全列画像（p3 最小） ----------
    try:
        fp = os.path.join(a.tmp, 'p3_DC_chr21.parquet')
        cols = pq.ParquetFile(fp).schema_arrow.names
        want = [c for c in ['gene_id', 'gene_name', 'variant_id', 'phenotype_id', 'af', 'maf',
                            'slope', 'slope_se', 'pval_nominal', 'tss_distance', 'chr', 'pos'] if c in cols]
        d = pq.read_table(fp, columns=want).to_pandas()
        prof = {'file': 'p3_DC_chr21.parquet', 'n_rows': int(len(d)), 'cols_used': want}
        for c in want:
            try:
                prof['nunique_' + c] = int(d[c].nunique())
                prof['head_' + c] = [str(x) for x in d[c].head(5).tolist()]
                if c in ('pos', 'pval_nominal', 'slope', 'af', 'maf', 'tss_distance', 'slope_se'):
                    prof['range_' + c] = [str(d[c].min()), str(d[c].max())]
            except Exception as e:
                prof['err_' + c] = repr(e)
        R['profile_p3'] = prof
    except Exception as e:
        R['profile_p3'] = {'error': repr(e)}

    with open(a.out, 'w', encoding='utf-8') as f:
        json.dump(R, f, ensure_ascii=False, indent=2, default=str)
    if not a.keep_tmp:
        shutil.rmtree(a.tmp, ignore_errors=True)
    print('WROTE ' + a.out)
    print('extracted=%d members_traversed=%d secs=%.1f' % (len(done), n_seen, R['extract_seconds']))


if __name__ == '__main__':
    main()
