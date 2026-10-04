#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""s65_figure_overview.py —— 一页式「图件总览」拼版图（Nature 体例）

内容：主图 Fig1-5（5） + 补充图 Fig. S1-S2（2） + 探索性诊断图 Fig. D1-D5（5）
  = 12 格，三区（Main / Supplementary / Exploratory diagnostics）。

规格（沿用 scripts/_fig_style.py）：
  NPG 十色 ; Arial ; 600 dpi ; 刻度朝内 ; 轴线 0.5 pt ; **全部文字 ASCII** ; 双重字形核验
  版心 240 mm（总览用宽版，非投稿单图宽），含三区标题 + 逐图标签 + 元信息行。

输出：figures/FigOverview_all_figures.{pdf,png}
日志：logs/s65_figure_overview.log
"""
import os
import io
import sys
import textwrap

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from PIL import Image

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _fig_style as FS

ROOT = FS.ROOT
FIGDIR = FS.FIGDIR
LOG = os.path.join(ROOT, "logs", "s65_figure_overview.log")
MM = 1.0 / 25.4  # mm -> inch

L = []
def p(*a):
    s = " ".join(str(x) for x in a); L.append(s)
    try: print(s)
    except UnicodeEncodeError: print(s.encode("ascii", "backslashreplace").decode("ascii"))

# ------------------------------------------------------------------ 内容表
MAIN = [
    ("Fig1_study_design",                 "Fig. 1  Study design and workflow",       "schematic; acronym list",  "p4"),
    ("Fig2_discovery_replication",        "Fig. 2  Discovery and replication",       "panels a-c",               "p8"),
    ("Fig3_coloc_sensitivity",            "Fig. 3  Colocalization and robustness",   "panels a-d",               "p10"),
    ("Fig4_chr1p36_12_finemap",           "Fig. 4  LD structure and attribution",    "panels a-f",               "p13"),
    ("Fig5_opentargets_phewas_safety",    "Fig. 5  Open Targets and PheWAS safety",  "panels a-f",               "p15"),
]
SUPP = [
    ("FigS1_instrument_power_gating",
     "Fig. S1  Instrument power and gating",
     "panels a-d | 182.4 x 70.6 mm"),
    ("FigS2_sensitivity_outcome_concordance",
     "Fig. S2  Sensitivity outcome GCST90483469",
     "panels a-b | 182.5 x 77.7 mm"),
    ("FigS3_manhattan_by_celltype",
     "Fig. S3  Manhattan by cell type",
     "panels a-n | 182.4 x 135.5 mm"),
    ("FigS4_discovery_volcano",
     "Fig. S4  Discovery volcano",
     "single panel | single col 89.0 x 68.6 mm"),
]
DIAG = [
    ("FigD1_discovery_volcano",    "Fig. D1  Discovery volcano",        "superseded by Fig. S4"),
    ("FigD2_discovery_manhattan",  "Fig. D2  Manhattan (14 cell types)", "superseded by Fig. S3"),
    ("FigD3_discovery_forest",     "Fig. D3  Forest (15 pairs)",        "superseded by Fig. 2b"),
    ("FigD4_discovery_celltype",   "Fig. D4  Cell-type burden",         "superseded by Fig. S1d"),
    ("FigD5_chr1_LD_structure",    "Fig. D5  chr1 LD structure",        "superseded by Fig. 4a"),
]

W_MM, H_MM = 240.0, 205.0
MG = 6.0            # 页边距
GAP = 4.0           # 格间距
LBL_H = 9.0         # 每格底部标签区高度


def load_thumb(name, cache={}):
    """读 PNG 并等比缩到长边 <= 1400 px（44 mm @600 dpi 约需 1040 px）。"""
    if name in cache:
        return cache[name]
    path = os.path.join(FIGDIR, name + ".png")
    assert os.path.exists(path), "缺图件：%s" % path
    im = Image.open(path).convert("RGB")
    if max(im.size) > 1400:
        r = 1400.0 / max(im.size)
        im = im.resize((max(1, int(im.size[0] * r)), max(1, int(im.size[1] * r))),
                       Image.LANCZOS)
    arr = np.asarray(im)
    cache[name] = arr
    return arr


def fit_rect(x, y, cw, ch, aspect, valign="bottom"):
    """在格子 (x,y,cw,ch) [mm, 左下原点] 内等比放入宽高比 aspect=w/h 的图，返回 mm 矩形。
    valign='bottom' -> 底部对齐（图注紧贴图下沿，网格更整齐）。"""
    if aspect >= cw / ch:
        w = cw; h = cw / aspect
    else:
        h = ch; w = ch * aspect
    ry = y if valign == "bottom" else y + (ch - h) / 2.0
    return x + (cw - w) / 2.0, ry, w, h


def add_thumb(fig, arr, x, y, cw, ch, label, meta, accent):
    h, w = arr.shape[:2]
    rx, ry, rw, rh = fit_rect(x, y, cw, ch, w / float(h), valign="bottom")
    ax = fig.add_axes([rx * MM / (W_MM * MM), ry * MM / (H_MM * MM),
                       rw * MM / (W_MM * MM), rh * MM / (H_MM * MM)])
    ax.imshow(arr, interpolation="lanczos", aspect="auto")
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_linewidth(0.5); sp.set_color(accent)
    ax.set_facecolor("white")
    # 标签（图区下方）
    fig.text((x + cw / 2.0) * MM / (W_MM * MM), (y - 3.0) * MM / (H_MM * MM),
             label, ha="center", va="top", fontsize=7.2, fontweight="bold", color="#2B2B2B")
    fig.text((x + cw / 2.0) * MM / (W_MM * MM), (y - 6.2) * MM / (H_MM * MM),
             meta, ha="center", va="top", fontsize=6.8, color="#7F7F7F")


def add_placeholder(fig, x, y, cw, ch, label, status, note, accent):
    ax = fig.add_axes([x * MM / (W_MM * MM), y * MM / (H_MM * MM),
                       cw * MM / (W_MM * MM), ch * MM / (H_MM * MM)])
    ax.set_xticks([]); ax.set_yticks([]); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    for sp in ax.spines.values():
        sp.set_linewidth(0.7); sp.set_color(accent); sp.set_linestyle((0, (3, 2)))
    ax.add_patch(Rectangle((0.008, 0.008), 0.984, 0.984, transform=ax.transAxes,
                           fill=False, lw=0.7, ls=(0, (3, 2)), ec=accent))
    ax.text(0.5, 0.70, status, ha="center", va="center", fontsize=7.6,
            fontweight="bold", color=accent)
    ax.text(0.5, 0.42, "\n".join(textwrap.wrap(note, 62)),
            ha="center", va="center", fontsize=7.0, color="#2B2B2B", linespacing=1.5)
    fig.text((x + cw / 2.0) * MM / (W_MM * MM), (y - 3.0) * MM / (H_MM * MM),
             label, ha="center", va="top", fontsize=7.2, fontweight="bold", color="#2B2B2B")


def main():
    import atexit, traceback
    atexit.register(lambda: io.open(LOG, "w", encoding="utf-8", newline="\n").write("\n".join(L) + "\n"))
    fam = FS.setup()
    p("== 图件总览拼版图 ==")
    p("字体族 = %s" % fam)
    p("版心 %.0f x %.0f mm" % (W_MM, H_MM))

    fig = plt.figure(figsize=(W_MM * MM, H_MM * MM), facecolor="white")

    def FX(mm): return mm * MM / (W_MM * MM)
    def FY(mm): return mm * MM / (H_MM * MM)

    # ---------------- 页眉 ----------------
    fig.text(FX(MG), FY(H_MM - 5.0),
             "Figure inventory - single-cell eQTL Mendelian randomization for female infertility",
             ha="left", va="top", fontsize=9.2, fontweight="bold", color="#2B2B2B")
    fig.text(FX(MG), FY(H_MM - 10.0),
             "chr1p36.12 | 5 main figures | 4 supplementary figures | "
             "5 exploratory diagnostic figures | generated 2026-09-29",
             ha="left", va="top", fontsize=7.0, color="#7F7F7F")
    fig.add_artist(plt.Line2D([FX(MG), FX(W_MM - MG)], [FY(H_MM - 12.5)] * 2,
                              color="#B8B8B8", lw=0.5, transform=fig.transFigure))

    acc_main = FS.C["mhc"]      # #3C5488
    acc_supp = FS.C["locus2"]   # #00A087（已生成）
    acc_diag = FS.C["grey"]     # #7F7F7F

    # ---------------- A. 主图 ----------------
    yA_top = H_MM - 16.0
    fig.text(FX(MG), FY(yA_top), "A   MAIN FIGURES (Fig. 1-5)  -  submitted set",
             ha="left", va="top", fontsize=8.0, fontweight="bold", color=acc_main)
    rhA = 70.0
    yA = yA_top - 4.0 - rhA
    cwA = (W_MM - 2 * MG - 4 * GAP) / 5.0
    for i, (stem, label, meta, pg) in enumerate(MAIN):
        x = MG + i * (cwA + GAP)
        arr = load_thumb(stem)
        add_thumb(fig, arr, x, yA + LBL_H, cwA, rhA - LBL_H, label,
                  "%s | %s" % (meta, pg), acc_main)
        p("  A%d %-32s %dx%d px" % (i + 1, stem, arr.shape[1], arr.shape[0]))

    # ---------------- B. 补充图 ----------------
    yB_top = yA - 4.0
    fig.text(FX(MG), FY(yB_top),
             "B   SUPPLEMENTARY FIGURES (Fig. S1-S4)  -  generated 2026-09-29",
             ha="left", va="top", fontsize=8.0, fontweight="bold", color=acc_supp)
    rhB = 46.0
    yB = yB_top - 4.0 - rhB
    cwB = (W_MM - 2 * MG - 3 * GAP) / 4.0
    for i, (stem, label, meta) in enumerate(SUPP):
        x = MG + i * (cwB + GAP)
        arr = load_thumb(stem)
        add_thumb(fig, arr, x, yB + LBL_H, cwB, rhB - LBL_H, label, meta, acc_supp)
        p("  B%d %-38s %dx%d px" % (i + 1, stem, arr.shape[1], arr.shape[0]))

    # ---------------- C. 诊断图 ----------------
    yC_top = yB - 4.0
    fig.text(FX(MG), FY(yC_top),
             "C   SUPERSEDED / EXPLORATORY DIAGNOSTIC FIGURES (Fig. D1-D5)  -  "
             "content absorbed into the main and supplementary figures; not for submission",
             ha="left", va="top", fontsize=8.0, fontweight="bold", color=acc_diag)
    rhC = 38.0
    yC = yC_top - 4.0 - rhC
    cwC = (W_MM - 2 * MG - 4 * GAP) / 5.0
    for i, (stem, label, note) in enumerate(DIAG):
        x = MG + i * (cwC + GAP)
        arr = load_thumb(stem)
        add_thumb(fig, arr, x, yC + LBL_H, cwC, rhC - LBL_H, label, note, acc_diag)
        p("  C%d %-32s %dx%d px" % (i + 1, stem, arr.shape[1], arr.shape[0]))

    # ---------------- 页脚 ----------------
    fig.text(FX(MG), FY(yC - 3.0),
             "Source: figures/*.png (600 dpi) | machine-readable inventory: tables/63_figure_inventory.csv | "
             "main figures plotted 2026-09-28, supplementary figures 2026-09-29, diagnostics 2026-09-22",
             ha="left", va="top", fontsize=6.8, color="#7F7F7F")

    paths, wmm, hmm = FS.save(fig, "FigOverview_all_figures")
    p("输出 %s" % " / ".join(paths))
    p("实测 %.1f x %.1f mm" % (wmm, hmm))
    p("RESULT=PASS")


if __name__ == "__main__":
    import traceback
    try:
        main()
    except BaseException:
        L.append(""); L.append("==== TRACEBACK ===="); L.append(traceback.format_exc())
        io.open(LOG, "w", encoding="utf-8", newline="\n").write("\n".join(L) + "\n")
        raise
