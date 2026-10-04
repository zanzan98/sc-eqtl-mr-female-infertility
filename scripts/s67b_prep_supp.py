# -*- coding: utf-8 -*-
"""Prep + verify data for supplementary Fig S1 / S2 (read-only)."""
import csv
import io
import os

P = "D:/endometriosis_project/11_sc_eqtl_mr_project/"
OUT = P + "logs/_s67b_prep.txt"
buf = []
def W(s=""):
    buf.append(str(s))

def rd(rel):
    with io.open(P + rel, encoding="utf-8-sig", newline="") as fh:
        return list(csv.DictReader(fh))

W("=" * 78)
W("A. tables/20_celltype_power.csv  (FULL)")
W("=" * 78)
rows = rd("tables/20_celltype_power.csv")
W("n_rows = %d ; fields = %s" % (len(rows), list(rows[0].keys())))
for r in rows:
    W("  %-10s instr=%-6s genes=%-6s tested=%-6s nominal=%-5s FDR=%-3s pct=%-6s flag=%s"
      % (r["cell_type"], r["n_instruments"], r["n_genes"], r["n_tested"],
         r["n_nominal"], r["n_FDR05"], r["pct_nominal"], r["power_flag"]))
W("  flags present: %s" % sorted(set(r["power_flag"] for r in rows)))
W("  n cell types = %d" % len(rows))
W("")

W("=" * 78)
W("B. tables/17_discovery_celltype_summary.csv  (FULL)")
W("=" * 78)
r17 = rd("tables/17_discovery_celltype_summary.csv")
W("n_rows = %d ; fields = %s" % (len(r17), list(r17[0].keys())))
for r in r17:
    W("  %-10s tests_m=%-6s nom=%-5s FDR=%-3s tests_q=%-6s FDRq=%-3s Fmed=%.2f MHC=%-4s"
      % (r["cell_type"], r["n_tests_main"], r["n_nominal_main"], r["n_FDR_main"],
         r["n_tests_qval"], r["n_FDR_qval"], float(r["F_stat_median"]), r["n_MHC_pairs"]))
W("")

W("=" * 78)
W("C. tables/39_sensitivity_vs_main_GCST90483469.csv  -> the 15/15 concordance")
W("=" * 78)
r39 = rd("tables/39_sensitivity_vs_main_GCST90483469.csv")
W("total rows = %d ; fields = %s" % (len(r39), list(r39[0].keys())))
both = [r for r in r39 if r["both_FDR05"].strip().lower() == "true"]
m_only = [r for r in r39 if r["main_FDR05"].strip().lower() == "true" and r["sens_FDR05"].strip().lower() != "true"]
s_only = [r for r in r39 if r["sens_FDR05"].strip().lower() == "true" and r["main_FDR05"].strip().lower() != "true"]
W("main_FDR05 & sens_FDR05 (both) = %d" % len(both))
W("main only = %d ; sens only = %d" % (len(m_only), len(s_only)))
W("sign_consistent among 'both' = %d / %d"
  % (sum(1 for r in both if r["sign_consistent"].strip().lower() == "true"), len(both)))
W("")
W("  %-11s %-9s %10s %10s %11s %11s %12s %10s %5s %5s"
  % ("gene", "cell_type", "b_main", "b_sens", "p_main", "p_sens", "qval", "F_stat", "sign", "both"))
for r in sorted(both, key=lambda x: -(abs(float(x["b_main"])))):
    W("  %-11s %-9s %10.6f %10.6f %11.4e %11.4e %12.4e %10.2f %5s %5s"
      % (r["gene"], r["cell_type"], float(r["b_main"]), float(r["b_sens"]),
         float(r["p_main"]), float(r["p_sens"]), float(r["qval"]), float(r["F_stat"]),
         r["sign_consistent"], r["both_FDR05"]))
W("")

# gene / cell_type composition of the 15
import collections
W("  gene counts: %s" % dict(collections.Counter(r["gene"] for r in both)))
W("  cell_type counts: %s" % dict(collections.Counter(r["cell_type"] for r in both)))
W("")

# join to locus
r15 = rd("tables/15_discovery_significant_with_locus.csv")
W("=" * 78)
W("D. tables/15_discovery_significant_with_locus.csv (15 pairs) joined to 39_")
W("=" * 78)
W("n_rows = %d" % len(r15))
idx = {(r["gene"], r["cell_type"]): r for r in both}
ok = 0
for r in r15:
    k = (r["gene"], r["cell_type"])
    m = idx.get(k)
    if m:
        ok += 1
    W("  %-4s %-10s %-8s b_main=%9.6f p=%10.3e | sens: %s"
      % (r.get("locus_id"), r["gene"], r["cell_type"], float(r["b"]), float(r["p"]),
         ("b=%9.6f p=%10.3e q=%10.3e" % (float(m["b_sens"]), float(m["p_sens"]), float(m["p_FDR"]))) if m else "MISSING"))
W("joined = %d / %d" % (ok, len(r15)))

# pearson r between b_main and b_sens over the 15
xs = [float(r["b_main"]) for r in both]
ys = [float(r["b_sens"]) for r in both]
n = len(xs)
mx, my = sum(xs) / n, sum(ys) / n
sxy = sum((a - mx) * (b - my) for a, b in zip(xs, ys))
sxx = sum((a - mx) ** 2 for a in xs)
syy = sum((b - my) ** 2 for b in ys)
r = sxy / (sxx * syy) ** 0.5
slope = sxy / sxx
W("")
W("Pearson r(b_main, b_sens) = %.6f ; slope = %.6f ; n = %d" % (r, slope, n))
W("sign consistency: all same sign = %s" % all(a * b > 0 for a, b in zip(xs, ys)))
W("range b_main = [%.6f, %.6f] ; b_sens = [%.6f, %.6f]" % (min(xs), max(xs), min(ys), max(ys)))

with io.open(OUT, "w", encoding="utf-8", newline="\n") as fh:
    fh.write("\n".join(buf) + "\n")
print("DONE", os.path.getsize(OUT))
