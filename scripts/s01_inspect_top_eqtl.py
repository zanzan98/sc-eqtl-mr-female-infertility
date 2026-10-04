import zipfile, os, sys, io, gzip, csv

z = r'D:\endometriosis_project\11_sc_eqtl_mr_project\00_data_raw\onek1k\top_eQTL_summary.zip'
out = r'D:\endometriosis_project\_top_eqtl_inspect.txt'
lines = []
with zipfile.ZipFile(z) as zf:
    lines.append("=== members ===")
    for i in zf.infolist():
        lines.append(f"{i.file_size/1e6:10.2f} MB  {i.filename}")
    names = zf.namelist()
    for nm in names[:6]:
        lines.append(f"\n=== preview: {nm} ===")
        with zf.open(nm) as fh:
            if nm.endswith('.gz'):
                fh = gzip.open(fh, 'rt', errors='replace')
            else:
                fh = io.TextIOWrapper(fh, errors='replace')
            for k, ln in enumerate(fh):
                if k >= 4: break
                lines.append(ln.rstrip()[:400])
            # count rows
        with zf.open(nm) as fh2:
            if nm.endswith('.gz'):
                fh2 = gzip.open(fh2, 'rt', errors='replace')
            else:
                fh2 = io.TextIOWrapper(fh2, errors='replace')
            n = sum(1 for _ in fh2)
        lines.append(f"    rows = {n:,}")
with open(out, 'w', encoding='utf-8') as f:
    f.write("\n".join(lines))
print("ok")
