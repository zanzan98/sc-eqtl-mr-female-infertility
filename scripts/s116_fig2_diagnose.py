# -*- coding: utf-8 -*-
"""s116_fig2_diagnose.py —— Fig2 布局诊断（只读，不 save）。

回答三问：
  1. (b)/(c) 面板的 y 轴标签区实际占用多少 mm？两者之间的**真实间隙**是多少？
  2. 间隙不足的根因是 wspace 太小，还是 (c) 标签本身太长？
  3. 背景点/阳性点各有多少个？点径现状多少？
"""
import collections
import io
import os
import sys

os.environ.setdefault("FIGPAL", "gb12")
os.environ.setdefault("FIG_NOTES", "0")
os.environ.setdefault("FIG_OUTDIR", r"D:/_transfer_logs/_p0_scratch")

HERE = r"D:/endometriosis_project/11_sc_eqtl_mr_project/scripts"
sys.path.insert(0, HERE)

import matplotlib                      # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt        # noqa: E402
import s37_final_figures as S          # noqa: E402

LOG = r"D:/_transfer_logs/_p0_scratch/_s116_fig2_layout.txt"
L = []


def w(s=""):
    L.append(str(s))


fig = S.fig2(1.0)
fig.canvas.draw()
r = fig.canvas.get_renderer()
dpi = float(fig.dpi)


def mm(px):
    return px / dpi * 25.4


def box(t):
    b = t.get_window_extent(r)
    return b.x0, b.x1, b.y0, b.y1


def span(txts):
    xs0 = []
    xs1 = []
    for t in txts:
        if not t.get_text().strip():
            continue
        b = box(t)
        xs0.append(b[0])
        xs1.append(b[1])
    if not xs0:
        return None
    return min(xs0), max(xs1)


axm, axf, axr = fig.axes
pw, ph = fig.get_size_inches()
w("figsize = %.3f x %.3f in   dpi = %.0f" % (float(pw), float(ph), dpi))
w("")
for nm, ax in (("axm(a)", axm), ("axf(b)", axf), ("axr(c)", axr)):
    p = ax.get_position()
    w("%-7s  axes 区 [%.2f, %.2f] mm  (宽 %.2f mm)   纵向 [%.2f, %.2f] mm"
      % (nm, mm(p.x0 * pw * dpi), mm(p.x1 * pw * dpi), mm(p.width * pw * dpi),
         mm(p.y0 * ph * dpi), mm(p.y1 * ph * dpi)))
w("")
w("gridspec: 2 行 x 2 列, height_ratios=[0.92, 1.25], hspace=0.42, wspace=0.54")
w("")

# ---- (b) 与 (c) 的 y 刻度标签跨度 ----
sb = span(axf.get_yticklabels())
sc = span(axr.get_yticklabels())
if sb is None or sc is None:
    raise SystemExit("no yticklabels")
bx0, bx1 = mm(sb[0]), mm(sb[1])
cx0, cx1 = mm(sc[0]), mm(sc[1])
axf_x0 = mm(axf.get_position().x0 * pw * dpi)
axf_x1 = mm(axf.get_position().x1 * pw * dpi)
axr_x0 = mm(axr.get_position().x0 * pw * dpi)
axr_x1 = mm(axr.get_position().x1 * pw * dpi)

w("(b) 面板 y 刻度标签区：%.2f -> %.2f mm   跨度 %.2f mm" % (bx0, bx1, bx1 - bx0))
w("(c) 面板 y 刻度标签区：%.2f -> %.2f mm   跨度 %.2f mm  ★" % (cx0, cx1, cx1 - cx0))
w("(b) axes 区：%.2f -> %.2f mm" % (axf_x0, axf_x1))
w("(c) axes 区：%.2f -> %.2f mm" % (axr_x0, axr_x1))
w("")
w("★ 关键间隙 = (c) 标签左端 - (b) axes 右缘 = %.2f - %.2f = **%.2f mm**"
  % (cx0, axf_x1, cx0 - axf_x1))
w("★ (c) 标签左端 - (b) 标签左端 = %.2f mm" % (cx0 - bx0))
w("★ 列间空白带宽 = (c) axes 左缘 - (b) axes 右缘 = %.2f mm" % (axr_x0 - axf_x1))
w("")
if cx0 < axf_x1:
    w(">>> 判定：**(c) 的 y 刻度标签已经压进 (b) 面板的绘图区 %.2f mm**"
      % (axf_x1 - cx0))
else:
    w(">>> 判定：(c) 标签未进入 (b) 绘图区；剩余视觉间隙仅 %.2f mm" % (cx0 - axf_x1))
w("")

# ---- 最长 / 最宽的 y 标签 ----
def widest(txts, n=6):
    rows = []
    for t in txts:
        s = t.get_text()
        if not s.strip():
            continue
        b = box(t)
        rows.append((mm(b[1] - b[0]), s))
    rows.sort(key=lambda x: -x[0])
    return rows[:n]


w("(b) 最宽的 6 个 y 标签：")
for wd, s in widest(axf.get_yticklabels()):
    w("   %6.2f mm  %r" % (wd, s))
w("")
w("(c) 最宽的 6 个 y 标签：")
for wd, s in widest(axr.get_yticklabels()):
    w("   %6.2f mm  %r" % (wd, s))
w("")

# ---- (a) 面板与 (b)/(c) 的纵向关系 ----
axm_y0 = mm(axm.get_position().y0 * ph * dpi)
axf_y1 = mm(axf.get_position().y1 * ph * dpi)
w("(a) 底 %.2f mm  vs  (b) 顶 %.2f mm  ->  纵向带隙 %.2f mm"
  % (axm_y0, axf_y1, axm_y0 - axf_y1))
w("")

# ---- 点数与点径 ----
disc = S.load("10_discovery_MR_main.csv")
allf = S.load("31_final_MR_with_method.csv")
fin = [x for x in allf if x["thr"].startswith("P<5e-8")]
ivw = [x for x in allf if not x["thr"].startswith("P<5e-8")][:1]
sig = [x for x in disc if x["sig_P_FDR"] == "True"]
w("背景点（全部 discovery 行）= %d 个 ；FDR 显著（上色）= %d 个" % (len(disc), len(sig)))
w("(b) 森林行 = %d（fin %d + ivw %d）" % (len(fin) + len(ivw), len(fin), len(ivw)))
w("(a) 背景点 ms = 0.85 pt（直径）；阳性点 ms = 3.6 pt")
w("")
w("阳性点色（GENE_C）：")
for g in ("CDC42", "LINC00339", "YME1L1", "ANXA4", "WNT4"):
    w("   %-10s %s" % (g, S.GENE_C.get(g, "(未定义)")))
w("   %-10s %s  (region / nr 共用)" % ("region/nr", S.FS.C["region"]))
plt.close(fig)

with io.open(LOG, "w", encoding="utf-8", newline="\n") as fh:
    fh.write("\n".join(L) + "\n")
print("done -> %s" % LOG)
