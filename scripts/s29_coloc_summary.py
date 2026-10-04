# -*- coding: utf-8 -*-
"""
s29_coloc_summary.py -- collapse tables/33_coloc_results.csv (per-row: 1 abf + K susie pairs)
into one row per locus x gene x cell_type, and join the MR result (31_final_MR_with_method.csv).

Output: tables/34_coloc_summary.csv
"""
import csv
import os
from collections import defaultdict

PROJ = r"D:\endometriosis_project\11_sc_eqtl_mr_project"
TAB = os.path.join(PROJ, "tables")

rows = list(csv.DictReader(open(os.path.join(TAB, "33_coloc_results.csv"), encoding="utf-8-sig")))

F = lambda x: float(x) if x not in ("", None) else None


def verdict(h4):
    if h4 is None:
        return "undetermined"
    if h4 > 0.8:
        return "strong_shared"
    if h4 >= 0.5:
        return "moderate_shared"
    return "no_shared"


groups = defaultdict(lambda: {"abf": None, "susie": []})
for r in rows:
    key = (r["locus_id"], r["gene"], r["cell_type"])
    if r["method"] == "coloc.abf":
        groups[key]["abf"] = r
    else:
        groups[key]["susie"].append(r)

# MR results (Wald / IVW) for the same pair
mr = {}
for r in csv.DictReader(open(os.path.join(TAB, "31_final_MR_with_method.csv"), encoding="utf-8-sig")):
    mr.setdefault((r["locus_id"], r["gene"], r["cell_type"]), []).append(r)

out = []
for key in sorted(groups, key=lambda k: (k[0], k[1], k[2])):
    g = groups[key]
    ab = g["abf"]
    su = g["susie"]
    h4s = [F(x["PP.H4"]) for x in su if F(x["PP.H4"]) is not None]
    h3s = [F(x["PP.H3"]) for x in su if F(x["PP.H3"]) is not None]
    best = max(su, key=lambda x: F(x["PP.H4"]) or -1) if su else None
    abf_h4 = F(ab["PP.H4"]) if ab else None
    abf_h3 = F(ab["PP.H3"]) if ab else None
    m = mr.get(key, [])
    mw = next((x for x in m if x["method"] == "Wald ratio"), None)
    mi = next((x for x in m if x["method"].startswith("IVW")), None)
    out.append({
        "locus_id": key[0], "gene": key[1], "cell_type": key[2],
        "nsnp_intersect": ab["nsnp_intersect"] if ab else "",
        "nsnp_used_abf": ab["nsnp_used"] if ab else "",
        "nsnp_used_susie": best["nsnp_used"] if best else "",
        "abf_H0": ab["PP.H0"] if ab else "",
        "abf_H1": ab["PP.H1"] if ab else "",
        "abf_H2": ab["PP.H2"] if ab else "",
        "abf_H3": ab["PP.H3"] if ab else "",
        "abf_H4": ab["PP.H4"] if ab else "",
        "abf_H4_over_H3H4": (round(abf_h4 / (abf_h3 + abf_h4), 4)
                             if abf_h3 is not None and abf_h4 is not None and (abf_h3 + abf_h4) > 0 else ""),
        "abf_verdict": ab["verdict"] if ab else "",
        "susie_n_pairs": len(su),
        "susie_maxH4": f"{max(h4s):.6g}" if h4s else "",
        "susie_n_pairs_H4_gt_0.8": sum(1 for x in h4s if x > 0.8),
        "susie_minH3": f"{min(h3s):.6g}" if h3s else "",
        "susie_best_note": best["note"] if best else "",
        "susie_verdict": verdict(max(h4s)) if h4s else "",
        "MR_method": (mw or mi or {}).get("method", ""),
        "MR_b": (mw or mi or {}).get("b", ""),
        "MR_se": (mw or mi or {}).get("se", ""),
        "MR_p": (mw or mi or {}).get("p", ""),
        "MR_instrument": (mw or mi or {}).get("instrument", ""),
        "IVW_tested": "yes" if mi else "no",
        "note": "",
    })

# --- notes ---
for r in out:
    n = []
    if r["abf_verdict"] == "no_shared" and r["susie_n_pairs_H4_gt_0.8"] > 0:
        n.append("abf(H3) 与 susie(部分对 H4>0.8) 不一致，以 abf 为主判据、susie 列全部 signal 对")
    if r["gene"] in ("CDC42", "LINC00339"):
        n.append("红线 17-①: 同一 chr1 座内 CDC42 与 LINC00339 分别报告，禁止合并")
    r["note"] = "; ".join(n)

cols = list(out[0].keys())
with open(os.path.join(TAB, "34_coloc_summary.csv"), "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=cols)
    w.writeheader()
    w.writerows(out)

print("wrote 34_coloc_summary.csv  rows =", len(out))
print()
hdr = f"{'locus':<3} {'gene':<10} {'cell_type':<11} {'H3_abf':>7} {'H4_abf':>7} {'verdict_abf':<16} {'nsusie':>6} {'maxH4_susie':>11} {'n>0.8':>5}"
print(hdr)
print("-" * len(hdr))
for r in out:
    print(f"{r['locus_id']:<3} {r['gene']:<10} {r['cell_type']:<11} "
          f"{(r['abf_H3'] if r['abf_H3']=='' else format(F(r['abf_H3']),'.4f')):>7} "
          f"{(r['abf_H4'] if r['abf_H4']=='' else format(F(r['abf_H4']),'.4f')):>7} "
          f"{r['abf_verdict']:<16} {r['susie_n_pairs']:>6} {r['susie_maxH4']:>11} {r['susie_n_pairs_H4_gt_0.8']:>5}")
