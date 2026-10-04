"""对照实验：Range（断点续传）请求是否导致 EBI / molgenis 挂起。"""
import time, os, traceback

OUT = r'D:\endometriosis_project\_resume_ab.txt'
lines = []

try:
    import requests

    TARGETS = [
        ('EBI_463', 'https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90483001-GCST90484000/GCST90483463/harmonised/GCST90483463.h.tsv.gz'),
        ('molgenis', 'https://molgenis26.gcc.rug.nl/downloads/1m-scbloodnl/eqtls_20201106_genome_wide.tar.gz'),
        ('zenodo_plink', 'https://zenodo.org/records/18870747/files/plink_genotype_merged_980_donors.zip?download=1'),
    ]
    WINDOW, CAP = 20, 20 * 1024 * 1024

    for name, url in TARGETS:
        for label, hdrs in (('no_range', {}), ('range_0_', {'Range': 'bytes=0-'})):
            got, status, note = 0, '', ''
            t0 = time.time()
            try:
                with requests.get(url, headers=hdrs, stream=True, timeout=(15, 20)) as r:
                    status = r.status_code
                    note = r.headers.get('Content-Range') or r.headers.get('Content-Length') or ''
                    for chunk in r.iter_content(1 << 16):
                        got += len(chunk)
                        if got >= CAP or time.time() - t0 > WINDOW:
                            break
                err = ''
            except Exception as e:
                err = f'{type(e).__name__}:{str(e)[:60]}'
            dt = max(time.time() - t0, 0.001)
            lines.append(f"{name:14s} {label:10s} HTTP={str(status):>4s} hdr={str(note)[:24]:24s} "
                         f"bytes={got:>11,}  KBps={got/1024/dt:8.1f}  {err}")
except Exception:
    lines.append('FATAL:\n' + traceback.format_exc())

with open(OUT, 'w', encoding='utf-8') as f:
    f.write("\n".join(lines))
