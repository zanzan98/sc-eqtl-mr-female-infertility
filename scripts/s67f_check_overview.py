# -*- coding: utf-8 -*-
"""s67f_check_overview.py — 总览图 PDF 字形/几何核验"""
import io
import os

ROOT = "D:/endometriosis_project/11_sc_eqtl_mr_project"
pdf = os.path.join(ROOT, "figures", "FigOverview_all_figures.pdf")
png = os.path.join(ROOT, "figures", "FigOverview_all_figures.png")
L = []
b = open(pdf, "rb").read()
L.append("FigOverview PDF bytes = %d ; LastResortHE = %d ; LastResort* = %d"
         % (len(b), b.count(b"LastResortHE"), b.count(b"LastResort")))
from PIL import Image
im = Image.open(png)
dpi = im.info.get("dpi", (600, 600))
L.append("PNG %dx%d px @ dpi=%.0f -> %.1f x %.1f mm"
         % (im.size[0], im.size[1], dpi[0],
            im.size[0] / dpi[0] * 25.4, im.size[1] / dpi[1] * 25.4))
L.append("VERDICT: %s" % ("PASS" if b.count(b"LastResortHE") == 0 else "FAIL"))
p = os.path.join(ROOT, "logs", "_s67f_overview_check.log")
with io.open(p, "w", encoding="utf-8", newline="\n") as fh:
    fh.write("\n".join(L) + "\n")
print("\n".join(L))
