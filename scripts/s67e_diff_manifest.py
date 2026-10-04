# -*- coding: utf-8 -*-
"""s67e_diff_manifest.py — 对比 s38 重跑前后的 manifest，确认旧行逐位不变、新行已加入"""
import csv
import io
import os

OLD = r"D:/_transfer_logs/_supp_backup_20260929_s1s2/_manifest.csv"
NEW = r"D:/endometriosis_project/11_sc_eqtl_mr_project/supplementary_data/_manifest.csv"
OUT = r"D:/endometriosis_project/11_sc_eqtl_mr_project/logs/_s67e_manifest_diff.log"


def rd(p):
    with io.open(p, encoding="utf-8-sig", newline="") as fh:
        return {r["file"]: r for r in csv.DictReader(fh)}


o, n = rd(OLD), rd(NEW)
L = []
L.append("old rows = %d ; new rows = %d" % (len(o), len(n)))
added = sorted(set(n) - set(o))
removed = sorted(set(o) - set(n))
L.append("added   = %s" % (added or "(none)"))
L.append("removed = %s" % (removed or "(none)"))
changed = [k for k in sorted(set(o) & set(n)) if o[k]["sha256"] != n[k]["sha256"]]
L.append("changed = %d %s" % (len(changed), changed or ""))
for k in added:
    L.append("  + %-52s %9s B  %s" % (k, n[k]["bytes"], n[k]["label"]))
L.append("VERDICT: %s" % ("PASS" if not removed and not changed else "CHECK"))
with io.open(OUT, "w", encoding="utf-8", newline="\n") as fh:
    fh.write("\n".join(L) + "\n")
print("\n".join(L))
