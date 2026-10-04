# -*- coding: utf-8 -*-
"""s130_letter_candidates.py —— 生成 Fig2 面板字母 B/C 的**位置候选对比图**（只读正本）。

三个候选（dy 全部保持 1.44 mm，即字母底距绘图区上缘 1.44 mm）：
  C1  dx = -0.0500  ->  字母左缘距绘图区左缘约 2.5 mm （紧贴绘图区左上角）
  C2  dx = -0.1000  ->  约 5.0 mm
  C3  dx = -0.1595  ->  约 8.0 mm （= 当前值，供对照）

★ 纪律：**必须设 FIG_OUTDIR 指向临时目录**，否则 `_fig_style.save` 会静默覆写 figures/ 正本。
★ 每个候选都量「字母 -> 最近文字」的间隙（文本对象口径），保证候选本身不造成重叠。
"""
import io
import os
import sys

CAND = [("C1 dx=0.0000 -> 0.0mm (letter left = plot left)", 0.0000),
        ("C2 dx=-0.0500 -> 2.5mm", -0.0500),
        ("C3 dx=-0.1000 -> 5.0mm", -0.1000),
        ("C4 dx=-0.1595 -> 8.0mm (current)", -0.1595),
        ("C5 dx=-0.6400 -> 32.0mm (align leftmost text)", -0.6400)]

TMP = r"D:/_transfer_logs/_fig2_letter_cand"
os.makedirs(TMP, exist_ok=True)

os.environ["FIGPAL"] = "gb12"
os.environ["FIG_NOTES"] = "0"
os.environ["FIG_OUTDIR"] = TMP          # ★ 防覆写正本
HERE = r"D:/endometriosis_project/11_sc_eqtl_mr_project/scripts"
sys.path.insert(0, HERE)

import matplotlib                       # noqa: E402
matplotlib.use("Agg")
import numpy as np                      # noqa: E402
from PIL import Image, ImageDraw, ImageFont   # noqa: E402
import s37_final_figures as S           # noqa: E402

SCALE = 1.02351
K = 600 / 25.4
LOG = os.path.join(TMP, "_candidates_log.txt")
L = []


def w(s=""):
    L.append(str(s))


w("候选对比生成   FIGDIR=%s（应为临时目录）" % S.FS.FIGDIR)
w("")

SHOTS = []
for name, dxv in CAND:
    S.FIG2_BC_LABEL_DX = dxv
    fig = S.fig2(SCALE)
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    fd = float(fig.dpi)
    bb = fig.get_tightbbox(r)

    # --- 间隙测量（画布内相对量，差值口径） ---
    def gap_of(letter_ax, lab):
        t = [x for x in letter_ax.texts
             if x.get_text().strip().upper() == lab][0]
        lb = t.get_window_extent(r)
        best = None
        for ax2 in fig.axes:
            pool = list(ax2.texts) + list(ax2.get_xticklabels()) + list(ax2.get_yticklabels())
            for t2 in pool:
                if t2 is t or not t2.get_text().strip():
                    continue
                b2 = t2.get_window_extent(r)
                if b2.width <= 0 or b2.height <= 0:
                    continue
                g = max(max(lb.x0 - b2.x1, b2.x0 - lb.x1),
                        max(lb.y0 - b2.y1, b2.y0 - lb.y1)) * 25.4 / fd
                if best is None or g < best[0]:
                    best = (g, t2.get_text().strip())
        return best

    axm, axf, axr = fig.axes
    for lab, ax in (("B", axf), ("C", axr)):
        b = gap_of(ax, lab)
        w("%-32s 面板 %s 最近文字间隙 %.2f mm  <- %r" % (name, lab, b[0], b[1][:40]))
    fig.savefig(os.path.join(TMP, "cand_%s.png" % name.split()[0]),
                dpi=600, bbox_inches="tight", pad_inches=0.02)
    SHOTS.append((name, os.path.join(TMP, "cand_%s.png" % name.split()[0])))
    matplotlib.pyplot.close(fig)

# --- 拼图：以 B 面板为代表（C 与 B 几何同构），裁到含数据区，便于判断"离图多远" ---
im0 = Image.open(SHOTS[0][1]).convert("RGB")
W, H = im0.size


def crop(x0mm, y0mm, x1mm, y1mm):
    return (max(0, int(x0mm * K)), max(0, int(y0mm * K)),
            min(W, int(x1mm * K)), min(H, int(y1mm * K)))


BX = crop(0, 64, 96, 90)          # B 面板：y 标签区 + 到绘图区数据区
BW, BH = BX[2] - BX[0], BX[3] - BX[1]
CAP = 42
PAD = 16
tot_h = 12 + len(SHOTS) * (CAP + BH + PAD)
canvas = Image.new("RGB", (BW + 24, tot_h), "white")
dr = ImageDraw.Draw(canvas)
try:
    font = ImageFont.truetype("C:/Windows/Fonts/arialbd.ttf", 28)
except Exception:
    font = ImageFont.load_default()

y = 12
for name, p in SHOTS:
    dr.text((16, y), "%s    [panel B]" % name, fill="black", font=font)
    y += CAP
    canvas.paste(Image.open(p).convert("RGB").crop(BX), (12, y))
    y += BH + PAD

OUT = os.path.join(TMP, "Fig2_letter_candidates.png")
canvas.save(OUT)
w("")
w("-> %s   %d x %d   (裁剪 x 0-96 mm, y 64-90 mm，1:1 像素)" % (OUT, canvas.size[0], canvas.size[1]))
for name, p in SHOTS:
    w("   %s -> %s" % (name, os.path.basename(p)))

with io.open(LOG, "w", encoding="utf-8", newline="\n") as fh:
    fh.write("\n".join(L) + "\n")
print("done -> %s" % OUT)
