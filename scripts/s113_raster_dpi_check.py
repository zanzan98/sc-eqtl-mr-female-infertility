# -*- coding: utf-8 -*-
"""
s113_raster_dpi_check.py —— 只读：figures/ 各 PDF 内嵌位图的"有效分辨率"。
Nature：图像 ≥300 dpi（导出建议 ≥450 dpi）。有效 dpi = 原生像素 / 显示物理尺寸(inch)。
★ 只读；结果写文件再读。
"""
import io
import os
import pymupdf

ROOT = r"D:/endometriosis_project/11_sc_eqtl_mr_project"
FIGD = os.path.join(ROOT, "figures")
LOG = os.path.join(ROOT, "logs", "_s113_raster_dpi.txt")

NAMES = [
    "Fig1_study_design", "Fig2_discovery_replication", "Fig3_coloc_sensitivity",
    "Fig4_chr1p36_12_finemap", "Fig5_opentargets_phewas_safety",
    "FigS1_instrument_power_gating", "FigS2_sensitivity_outcome_concordance",
    "FigS3_manhattan_by_celltype", "FigS4_discovery_volcano",
    "FigD1_discovery_volcano", "FigD2_discovery_manhattan",
    "FigD3_discovery_forest", "FigD4_discovery_celltype",
    "FigD5_chr1_LD_structure", "FigOverview_all_figures",
]
L = []


def w(s=""):
    L.append(str(s))


def main():
    w("=" * 96)
    w("figures/ PDF 内嵌位图有效分辨率（有效 dpi = 原生像素 / 显示尺寸(in)）")
    w("=" * 96)
    w("%-38s %-4s %-6s %-12s %-10s %s" %
      ("图", "页", "图序", "原生px", "显示mm", "有效dpi"))
    w("-" * 96)
    worst = []
    for n in NAMES:
        p = os.path.join(FIGD, n + ".pdf")
        if not os.path.exists(p):
            continue
        d = pymupdf.open(p)
        for pno in range(d.page_count):
            pg = d[pno]
            for i, im in enumerate(pg.get_images(full=True)):
                xref = im[0]
                px_w, px_h = im[2], im[3]
                rects = pg.get_image_rects(xref)
                if not rects:
                    continue
                r = rects[0]
                wmm, hmm = r.width / 72.0 * 25.4, r.height / 72.0 * 25.4
                dpi = px_w / (r.width / 72.0) if r.width else 0
                w("%-38s %-4d %-6d %-12s %-10s %.0f" %
                  (n, pno + 1, i + 1, "%dx%d" % (px_w, px_h),
                   "%.1fx%.1f" % (wmm, hmm), dpi))
                worst.append((dpi, n, pno + 1, i + 1, px_w, px_h, wmm, hmm))
        d.close()
    w("")
    w("=" * 96)
    w("判定")
    w("=" * 96)
    if not worst:
        w("· 无内嵌位图")
    else:
        worst.sort()
        lo = worst[0]
        w("· 最低有效 dpi = %.0f  @ %s p%d img%d  (%dx%d px 显示 %.1fx%.1f mm)" %
          (lo[0], lo[1], lo[2], lo[3], lo[4], lo[5], lo[6], lo[7]))
        for t in worst:
            flag = "  <-- 低于 450" if t[0] < 450 else ""
            if t[0] < 450:
                w("    %.0f dpi  %s p%d img%d" % (t[0], t[1], t[2], t[3]) + flag)
        w("· 内嵌位图总数 = %d（去重按出现次数）" % len(worst))
    with io.open(LOG, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(L) + "\n")
    print("done -> %s" % LOG)


if __name__ == "__main__":
    main()
