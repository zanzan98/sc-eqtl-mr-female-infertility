# -*- coding: utf-8 -*-
"""s67g_crop_s4.py — 原分辨率裁切 Fig S4 局部，用于「文字压点」目视判据"""
import os
from PIL import Image
FIG = r"D:/endometriosis_project/11_sc_eqtl_mr_project/figures"
PREV = os.path.join(FIG, "_preview_FigS4_crops.png")

src = Image.open(os.path.join(FIG, "FigS4_discovery_volcano.png")).convert("RGB")
W, H = src.size
dpi = 600.0
mmw, mmh = W / dpi * 25.4, H / dpi * 25.4
print("FigS4 PNG %dx%d px -> %.1f x %.1f mm" % (W, H, mmw, mmh))

# 归一化裁切框（左, 上, 右, 下），top-origin
BOXES = [
    ("A left-labels", (0.00, 0.25, 0.55, 0.62)),
    ("B legend",      (0.02, 0.02, 0.55, 0.22)),
    ("C right-ann",   (0.45, 0.02, 1.00, 0.55)),
]
crops = []
for name, (l, t, r, b) in BOXES:
    box = (int(l * W), int(t * H), int(r * W), int(b * H))
    c = src.crop(box)
    sc = 2 if max(c.size) < 1400 else 1
    if sc > 1:
        c = c.resize((c.size[0] * sc, c.size[1] * sc), Image.LANCZOS)
    crops.append((name, c))
    print("  %-14s box=%s -> %dx%d" % (name, box, c.size[0], c.size[1]))

pad = 14
totw = max(c.size[0] for _, c in crops)
toth = sum(c.size[1] for _, c in crops) + pad * (len(crops) + 1)
canvas = Image.new("RGB", (totw + 2 * pad, toth), "white")
y = pad
for _, c in crops:
    canvas.paste(c, (pad, y))
    y += c.size[1] + pad
canvas.save(PREV)
print("wrote", PREV, canvas.size)
