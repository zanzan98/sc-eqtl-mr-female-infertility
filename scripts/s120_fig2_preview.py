# -*- coding: utf-8 -*-
"""s120_fig2_preview.py —— 生成 Fig2 重绘后的预览图（只读 figures/，不重绘）。

产出（到 D:/_transfer_logs/_fig2_preview_20261001/）：
  Fig2_full_preview.png      整图，缩到宽 2000 px（用 Image.BOX，项目纪律）
  Fig2_zoom_bc_gap.png       (b) 右部 + 列间空白带 + (c) 左部，1:1 像素，看间隙
  Fig2_zoom_panel_a.png      (a) 面板，1:1 像素，看背景浅色系与阳性点对比

坐标换算：用 matplotlib 的 get_tightbbox() 精确反推 PNG 像素坐标
（fig.savefig 用 bbox_inches='tight', pad_inches=0.02, dpi=600）。
"""
import io
import os
import sys

os.environ["FIGPAL"] = "gb12"
os.environ["FIG_NOTES"] = "0"
os.environ.pop("FIG_OUTDIR", None)

HERE = r"D:/endometriosis_project/11_sc_eqtl_mr_project/scripts"
sys.path.insert(0, HERE)

import matplotlib                       # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt         # noqa: E402
from PIL import Image                   # noqa: E402
import s37_final_figures as S           # noqa: E402

OUT = r"D:/_transfer_logs/_fig2_preview_20261001"
os.makedirs(OUT, exist_ok=True)
LOG = os.path.join(OUT, "_preview_log.txt")
DPI = 600
PAD_IN = 0.02
L = []


def w(s=""):
    L.append(str(s))


SCALE = S.CAL.get("Fig2_discovery_replication", 1.0)
fig = S.fig2(SCALE)
fig.canvas.draw()
r = fig.canvas.get_renderer()
bb = fig.get_tightbbox(r)          # inches, 未含 pad

X0 = bb.x0 - PAD_IN
Y0 = bb.y0 - PAD_IN
W_IN = bb.width + 2 * PAD_IN
H_IN = bb.height + 2 * PAD_IN


def to_px(x_in, y_in):
    """figure 坐标(inch, 原点左下) -> PNG 像素(原点左上)"""
    return ((x_in - X0) * DPI, (bb.y1 + PAD_IN - y_in) * DPI)


def ax_rect(ax):
    p = ax.get_position()
    w_in, h_in = fig.get_size_inches()
    x0, y0 = to_px(p.x0 * w_in, p.y1 * h_in)
    x1, y1 = to_px(p.x1 * w_in, p.y0 * h_in)
    return (x0, y0, x1, y1)


axm, axf, axr = fig.axes
w("tight bbox = %.3f x %.3f in  (pad %.2f)  scale=%.5f" % (W_IN, H_IN, PAD_IN, SCALE))
w("PNG 预期尺寸 = %d x %d px" % (round(W_IN * DPI), round(H_IN * DPI)))
for nm, ax in (("a", axm), ("b", axf), ("c", axr)):
    w("  panel %s  PNG px rect = (%.0f, %.0f) - (%.0f, %.0f)" % ((nm,) + ax_rect(ax)))

src = os.path.join(S.FS.FIGDIR, "Fig2_discovery_replication.png")
im = Image.open(src).convert("RGB")
w("实际 PNG 尺寸 = %d x %d px" % im.size)
plt.close(fig)

# ---- 1. 整图预览 ----
tw = 2000
th = round(im.size[1] * tw / im.size[0])
im.resize((tw, th), Image.BOX).save(os.path.join(OUT, "Fig2_full_preview.png"))
w("-> Fig2_full_preview.png  %d x %d" % (tw, th))

# ---- 2. (b) 右部 + 空白带 + (c) 左部 ----
rb = ax_rect(axf)
rc = ax_rect(axr)
x0 = int(rb[2] - 0.30 * (rb[2] - rb[0]))       # (b) 轴右 30%
x1 = int(rc[0] + 0.30 * (rc[2] - rc[0]))       # (c) 轴左 30%
y0 = int(min(rb[1], rc[1]) - 40)               # 含面板字母余量
y1 = int(max(rb[3], rc[3]) + 10)
x0, x1 = max(0, x0), min(im.size[0], x1)
y0, y1 = max(0, y0), min(im.size[1], y1)
im.crop((x0, y0, x1, y1)).save(os.path.join(OUT, "Fig2_zoom_bc_gap.png"))
w("-> Fig2_zoom_bc_gap.png  crop=(%d, %d, %d, %d)  %d x %d"
  % (x0, y0, x1, y1, x1 - x0, y1 - y0))

# ---- 3. (a) 面板 ----
ra = ax_rect(axm)
x0 = max(0, int(ra[0] - 0.03 * im.size[0]))
x1 = min(im.size[0], int(ra[2] + 0.01 * im.size[0]))
y0 = max(0, int(ra[1] - 30))
y1 = min(im.size[1], int(ra[3] + 10))
im.crop((x0, y0, x1, y1)).save(os.path.join(OUT, "Fig2_zoom_panel_a.png"))
w("-> Fig2_zoom_panel_a.png  crop=(%d, %d, %d, %d)  %d x %d"
  % (x0, y0, x1, y1, x1 - x0, y1 - y0))

w("")
w("PDF = %s" % os.path.join(S.FS.FIGDIR, "Fig2_discovery_replication.pdf"))

with io.open(LOG, "w", encoding="utf-8", newline="\n") as fh:
    fh.write("\n".join(L) + "\n")
print("done -> %s" % OUT)
