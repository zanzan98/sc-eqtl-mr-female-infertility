# -*- coding: utf-8 -*-
"""s67d_crop_check.py — 裁剪 FigS2 中缝区域，用于目视排查「文字压文字」"""
import io
import os
from PIL import Image

FIGD = "D:/endometriosis_project/11_sc_eqtl_mr_project/figures"
p = os.path.join(FIGD, "FigS2_sensitivity_outcome_concordance.png")
im = Image.open(p)
W, H = im.size
# 中缝：a 面板右边界 到 b 面板标签区
box = (int(W * 0.26), 0, int(W * 0.52), H)
crop = im.crop(box)
crop = crop.resize((crop.size[0] // 2, crop.size[1] // 2), Image.LANCZOS)
o = os.path.join(FIGD, "_preview_FigS2_gap.png")
crop.save(o)
print("full=%dx%d  crop=%s  saved=%s" % (W, H, box, o))
