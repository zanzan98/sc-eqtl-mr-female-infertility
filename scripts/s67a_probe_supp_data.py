# -*- coding: utf-8 -*-
"""Probe candidate data sources for supplementary Fig S1 / S2."""
import io
import os

P = "D:/endometriosis_project/11_sc_eqtl_mr_project/"
OUT = P + "logs/_s67a_probe.txt"

CAND = [
    ("tables/17_discovery_celltype_summary.csv", 8),
    ("tables/20_celltype_power.csv", 8),
    ("tables/_s60a_cis_tool_counts.csv", 20),
    ("tables/01_discovery_instruments.csv", 4),
    ("tables/01c_discovery_instruments_qval05.csv", 4),
    ("tables/10_discovery_MR_main.csv", 3),
    ("tables/11_discovery_significant_P_FDR05.csv", 3),
    ("tables/12_discovery_sensitivity_qval05.csv", 3),
    ("tables/13_discovery_sensitivity_qval05_significant.csv", 3),
    ("tables/15_discovery_significant_with_locus.csv", 3),
    ("tables/16_discovery_locus_summary.csv", 8),
    ("tables/24b_replication_unique_tests.csv", 3),
    ("tables/31_final_MR_with_method.csv", 3),
    ("tables/38_discovery_MR_sensitivity_GCST90483469.csv", 3),
    ("tables/39_sensitivity_vs_main_GCST90483469.csv", 3),
    ("supplementary_data/Fig2c_replication_comparison.csv", 3),
    ("supplementary_data/Fig5e_power_class_968.csv", 3),
    ("supplementary_data/_manifest.csv", 60),
]

lines = []
def P_(s=""):
    lines.append(str(s))

for rel, nshow in CAND:
    fp = P + rel
    P_("=" * 78)
    P_("FILE %s" % rel)
    if not os.path.exists(fp):
        P_("  MISSING")
        continue
    P_("  bytes = %d" % os.path.getsize(fp))
    with io.open(fp, encoding="utf-8", errors="replace") as fh:
        rows = fh.read().splitlines()
    P_("  lines = %d" % len(rows))
    for i, r in enumerate(rows[:nshow]):
        P_("   [%02d] %s" % (i, r[:600]))
    if len(rows) > nshow:
        P_("   ... (%d more)" % (len(rows) - nshow))

with io.open(OUT, "w", encoding="utf-8", newline="\n") as fh:
    fh.write("\n".join(lines) + "\n")
print("DONE", os.path.getsize(OUT))
