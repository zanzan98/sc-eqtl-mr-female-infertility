# -*- coding: utf-8 -*-
"""s69c_hash_report.py — 汇总当前交付物 sha256 与登记信息（供 65_ 交付报告引用）"""
import hashlib, io, os, zipfile
ROOT = r"D:/endometriosis_project/11_sc_eqtl_mr_project"
FIG = os.path.join(ROOT, "figures")
TAB = os.path.join(ROOT, "tables")
SUP = os.path.join(ROOT, "supplementary_data")

TARGETS = []
for n in ["FigS1_instrument_power_gating.pdf", "FigS1_instrument_power_gating.png",
          "FigS2_sensitivity_outcome_concordance.pdf", "FigS2_sensitivity_outcome_concordance.png",
          "FigS3_manhattan_by_celltype.pdf", "FigS3_manhattan_by_celltype.png",
          "FigS4_discovery_volcano.pdf", "FigS4_discovery_volcano.png",
          "FigOverview_all_figures.pdf", "FigOverview_all_figures.png"]:
    TARGETS.append(os.path.join(FIG, n))
for n in ["63_figure_inventory.csv", "64a_FigS1_celltype_power_gating.csv",
          "64b_FigS2_sensitivity_vs_main_15pairs.csv",
          "64d_FigS3S4_discovery_thresholds.csv"]:
    TARGETS.append(os.path.join(TAB, n))
TARGETS.append(os.path.join(ROOT, "Supplementary_Data.zip"))
TARGETS.append(os.path.join(SUP, "_manifest.csv"))
for n in ["38_英文正文_论文体重写_v5.md", "39_英文正文_论文体v6.md",
          "38_中文版_论文体重写_v5_正文与图注.md", "39_中文版_论文体v7_正文与图注.md",
          "43_论文初稿_带图_v5.md", "39b_技术文档_SM与SN_v8.md"]:
    TARGETS.append(os.path.join(ROOT, n))

out = []
for p in TARGETS:
    if not os.path.exists(p):
        out.append("%-52s MISSING" % os.path.relpath(p, ROOT).replace("\\", "/"))
        continue
    h = hashlib.sha256(open(p, "rb").read()).hexdigest()
    out.append("%-52s %10d  %s" % (os.path.relpath(p, ROOT).replace("\\", "/"), os.path.getsize(p), h))

# manifest 概况
mp = os.path.join(SUP, "_manifest.csv")
lines = io.open(mp, encoding="utf-8-sig").read().splitlines()
rows = [l for l in lines[1:] if l.strip()]
from collections import Counter
kinds = Counter(r.split(",")[1] for r in rows if len(r.split(",")) > 1)
out.append("")
out.append("manifest: total_lines=%d  data_rows=%d  kinds=%s" % (len(lines), len(rows), dict(kinds)))

# zip 内条目数
z = zipfile.ZipFile(os.path.join(ROOT, "Supplementary_Data.zip"))
out.append("zip entries = %d ; namelist head = %s" % (len(z.namelist()), z.namelist()[:4]))

# inventory 行
inv = io.open(os.path.join(TAB, "63_figure_inventory.csv"), encoding="utf-8-sig").read()
out.append("")
out.append("--- 63_figure_inventory.csv (raw) ---")
out.append(inv)

p = os.path.join(ROOT, "logs", "_s69c_hash_report.txt")
with io.open(p, "w", encoding="utf-8", newline="\n") as fh:
    fh.write("\n".join(out) + "\n")
print("wrote", p)
