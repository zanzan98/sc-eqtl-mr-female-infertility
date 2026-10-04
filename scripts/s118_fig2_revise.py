# -*- coding: utf-8 -*-
"""s118_fig2_revise.py —— Fig2 排版与配色专改（图件细改 · Fig2 · 问题1+2+3）。

重出 `figures/Fig2_discovery_replication.{pdf,png}`，并做**三重客观验收**：
  A. 幅宽/幅高（fit_save 标定到 182.4 mm 硬窗）
  B. (b)/(c) 面板间隙复测：**(c) 的 y 刻度标签左端 与 (b) 面板绘图区右缘**的距离
     （改前实测 -8.66 mm = 压进 (b) 8.66 mm）
  C. 像素级配色验证：BG6 六色与 5 个基因色的**实际成像像素数**必须都 > 0
     （防止"改了常量但没上图像"）

★ 环境：FIGPAL=gb12、FIG_NOTES=0，且**不设 FIG_OUTDIR** ⇒ 落盘到项目 figures/。
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
import s37_final_figures as S           # noqa: E402

LOG = r"D:/endometriosis_project/11_sc_eqtl_mr_project/logs/_s118_fig2_revise.txt"
NAME = "Fig2_discovery_replication"
L = []


def w(s=""):
    L.append(str(s))


w("=" * 100)
w("Fig2 排版与配色专改（问题1 布局 / 问题2 背景配色 / 问题3 面板字母大写）")
w("FIGDIR = %s" % S.FS.FIGDIR)
w("=" * 100)

# ---------------- A. 重出 ----------------
paths, wmm, hmm = S.fit_save(S.fig2, NAME)
w("")
w("A. 重出结果")
w("   %s" % NAME)
w("   实测 %.2f x %.2f mm   （硬窗 181.5-183.3 mm 宽；高 <= 170 mm）" % (wmm, hmm))
w("   scale = %.5f" % S.CAL.get(NAME, float("nan")))
w("   宽度合格 = %s ；高度合格 = %s" % (181.5 <= wmm <= 183.3, hmm <= 170.0))

# ---------------- B. (b)/(c) 间隙复测 ----------------
scale = S.CAL.get(NAME, 1.0)
fig = S.fig2(scale)
fig.canvas.draw()
r = fig.canvas.get_renderer()
dpi = float(fig.dpi)


def mm(px):
    return px / dpi * 25.4


def span(txts):
    xs = [t.get_window_extent(r) for t in txts if t.get_text().strip()]
    return (min(b.x0 for b in xs), max(b.x1 for b in xs)) if xs else None


axm, axf, axr = fig.axes
pw, ph = fig.get_size_inches()
sb = span(axf.get_yticklabels())
sc = span(axr.get_yticklabels())
bx0, bx1 = mm(sb[0]), mm(sb[1])
cx0, cx1 = mm(sc[0]), mm(sc[1])
axf_x0, axf_x1 = mm(axf.get_position().x0 * pw * dpi), mm(axf.get_position().x1 * pw * dpi)
axr_x0, axr_x1 = mm(axr.get_position().x0 * pw * dpi), mm(axr.get_position().x1 * pw * dpi)

w("")
w("B. (b)/(c) 面板间隙复测   scale = %.5f" % scale)
w("   (b) axes 区        : %.2f -> %.2f mm   (宽 %.2f)" % (axf_x0, axf_x1, axf_x1 - axf_x0))
w("   (c) axes 区        : %.2f -> %.2f mm   (宽 %.2f)" % (axr_x0, axr_x1, axr_x1 - axr_x0))
w("   (b) y 刻度标签区   : %.2f -> %.2f mm   (最宽 %.2f)" % (bx0, bx1, bx1 - bx0))
w("   (c) y 刻度标签区   : %.2f -> %.2f mm   (最宽 %.2f)" % (cx0, cx1, cx1 - cx0))
w("   列间绝对空白带     : %.2f mm  ( = (c) axes 左缘 - (b) axes 右缘 )" % (axr_x0 - axf_x1))
gap = cx0 - axf_x1
w("   ★ 关键间隙 = (c) 标签左端 - (b) axes 右缘 = %.2f - %.2f = **%.2f mm**"
  % (cx0, axf_x1, gap))
w("   ★ 改前 = **-8.66 mm**（压进 (b) 8.66 mm）")
w("   >>> %s" % ("通过：(c) 标签不再进入 (b) 绘图区"
                if gap >= 0 else "未通过：(c) 标签仍压进 (b) %.2f mm" % (-gap)))
w("   纵向带隙 (a)底-(b)顶 = %.2f mm"
  % (mm(axm.get_position().y0 * ph * dpi) - mm(axf.get_position().y1 * ph * dpi)))

# ---------------- D. (a) 左移 / 纸面居中复测（2026-10-01 请求 J） ----------------
# 口径说明：
#   (a) 的**绘图区**左右缘本就与 (b) 左缘 / (c) 右缘等宽对齐（改动前实测同为 px 850 -> 4198）。
#   它「看着偏右」是因为左侧 36 mm 全被 (b) 的长 y 标签占掉，右侧仅 5.5 mm ⇒ 纸面「左重右轻」。
#   本段同时报告两个视角：
#     ① 绘图区对齐差 = (a).x0 - (b).x0   （改后应 ≈ -17.75 mm，因为主动左移）
#     ② 纸面条带留白 = 用**整张 tight bbox 的有效墨迹边界**算左右留白
#        —— 左边界取 min(a 左缘, b 标签左缘)，右边界取 max(a 右缘, c 右缘)。
tb = fig.get_tightbbox(r)
a_x0 = mm(axm.get_position().x0 * pw * dpi)
a_x1 = mm(axm.get_position().x1 * pw * dpi)
w("")
w("D. (a) 左移 / 纸面居中复测（请求 J）")
w("   (a) 绘图区         : %.2f -> %.2f mm   (宽 %.2f)" % (a_x0, a_x1, a_x1 - a_x0))
w("   绘图区对齐差 (a).x0-(b).x0 = %.2f mm   （改后应为负 = (a) 已左移）" % (a_x0 - axf_x0))
w("   理论错位 = -%.2f mm（FIG2_A_SHIFT=%.5f * %.2f mm 画布宽）"
  % (S.FIG2_A_SHIFT * pw * 25.4, S.FIG2_A_SHIFT, pw * 25.4))
# 「全图墨迹缘」= 所有面板墨迹的最外缘（(a) 左缘 vs (b) 标签左端 取左； (a) 右缘 vs (c) axes 右缘 取右）
left_ink = min(a_x0, bx0)
right_ink = max(a_x1, axr_x1)
mL = a_x0 - left_ink
mR = right_ink - a_x1
w("   全图墨迹缘         : %.2f -> %.2f mm" % (left_ink, right_ink))
w("   ★ (a) 左侧留白 = (a).x0 - 墨迹左缘 = %.2f - (%.2f) = **%.2f mm**" % (a_x0, left_ink, mL))
w("   ★ (a) 右侧留白 = 墨迹右缘 - (a).x1 = %.2f - %.2f = **%.2f mm**" % (right_ink, a_x1, mR))
w("   >>> 左右留白差 = %.2f mm   %s" % (abs(mL - mR),
                                       "通过：(a) 在纸面内基本居中（差 <= 3 mm）"
                                       if abs(mL - mR) <= 3.0 else "未通过：仍左右失衡"))
w("   改前对照：左侧留白 0.0 mm（(b) 标签即最左墨迹）／右侧留白 %.2f mm"
  % (right_ink - 163.78))
w("   注：(b) 的 y 标签左缘 = %.2f mm，(a).x0 = %.2f mm ⇒ 上下错位 %.2f mm（用户已选纸面居中）"
  % (bx0, a_x0, a_x0 - bx0))
w("   tight bbox(in) = %.3f x %.3f   (= %.2f x %.2f mm 含 pad 0.02 in)"
  % (tb.width, tb.height, tb.width * 25.4, tb.height * 25.4))

# ---------------- E. PNG 像素级核验（独立第二口径：真实成像，不含布局假设） ----------------
# ★ bbox_inches='tight' 会裁剪并补 0.02 in pad ⇒ PNG 像素坐标 ≠ 画布坐标，必须做偏移换算。
_bb = tb
_PAD = 0.02
_X0 = _bb.x0 - _PAD
_W_IN = _bb.width + 2 * _PAD
_H_IN = _bb.height + 2 * _PAD
_DPI = float(fig.dpi)


def to_px(x_in, y_in):
    return ((x_in - _X0) * _DPI, (_bb.y1 + _PAD - y_in) * _DPI)


def band_rows(ax):
    p = ax.get_position()
    x0, y0 = to_px(p.x0 * pw, p.y1 * ph)
    x1, y1 = to_px(p.x1 * pw, p.y0 * ph)
    return (int(round(min(y0, y1))), int(round(max(y0, y1))))


w("")
w("E. PNG 像素级 (a) 居中核验（独立口径：直接扫非白像素）")
try:
    import numpy as np
    from PIL import Image
    png = os.path.join(S.FS.FIGDIR, NAME + ".png")
    im = Image.open(png).convert("RGB")
    arr = np.asarray(im).astype(int)
    H, W = arr.shape[:2]
    k = 25.4 / 600.0
    nonwhite = (arr.sum(axis=2) < 720)      # 墨迹掩膜（白底 3*255=765）
    colany = nonwhite.any(axis=0)
    x0i = int(np.argmax(colany))
    x1i = int(W - 1 - np.argmax(colany[::-1]))

    def band_x(r0, r1):
        sub = nonwhite[max(0, r0):min(H, r1)]
        c = sub.any(axis=0)
        if not c.any():
            return (float("nan"), float("nan"))
        return (int(np.argmax(c)) * k, int(W - 1 - np.argmax(c[::-1])) * k)

    w("   PNG 预期 %d x %d px（实测 %d x %d）" % (round(_W_IN * 600), round(_H_IN * 600), W, H))
    w("   全图墨迹列范围     : %.2f -> %.2f mm" % (x0i * k, x1i * k))
    ax0p, ax1p = band_x(*band_rows(axm))
    b0p, b1p = band_x(*band_rows(axf))
    c0p, c1p = band_x(*band_rows(axr))
    w("   (a) 行带墨迹列范围 : %.2f -> %.2f mm" % (ax0p, ax1p))
    w("   (b) 行带墨迹列范围 : %.2f -> %.2f mm" % (b0p, b1p))
    w("   (c) 行带墨迹列范围 : %.2f -> %.2f mm" % (c0p, c1p))
    mL2, mR2 = ax0p - x0i * k, x1i * k - ax1p
    w("   ★ PNG 口径 (a) 左留白 = %.2f mm ／右留白 = %.2f mm  ⇒ 差 %.2f mm  %s"
      % (mL2, mR2, abs(mL2 - mR2),
         "通过（<=5 mm）" if abs(mL2 - mR2) <= 5.0 else "未通过"))
    w("   ★ (a) 行带墨迹宽 = %.2f mm vs (b)+(c) 横跨 %.2f ~ %.2f mm"
      % (ax1p - ax0p, b0p, c1p))
    w("   注：(a) 行带含 x 轴刻度数字，(b) 行带含其长 y 标签 ⇒ 两者边界不同源，"
      "纸面居中以上方 Text/axes 口径为准。")
except Exception as e:
    w("   [!! PNG 核验失败] %r" % (e,))
plt.close(fig)

# ---------------- C. 像素级配色验证 ----------------
w("")
w("C. 像素级配色验证（容差 <=4/255，读 600 dpi PNG）")
try:
    import numpy as np
    from PIL import Image
    png = os.path.join(S.FS.FIGDIR, NAME + ".png")
    a = np.asarray(Image.open(png).convert("RGB")).astype(int)

    def near(hexv, tol=4):
        c = np.array([int(hexv[i:i + 2], 16) for i in (1, 3, 5)])
        return int((np.abs(a - c).max(axis=2) <= tol).sum())

    w("   PNG = %s   尺寸 %d x %d px" % (os.path.basename(png), a.shape[1], a.shape[0]))
    w("   -- 背景浅色系 BG6（应全部 > 0）--")
    allbg = True
    for c in S.BG6:
        n = near(c)
        allbg &= n > 0
        w("      %s  %8d px" % (c, n))
    w("   -- 前景基因色（应全部 > 0）--")
    for g, c in S.GENE_C.items():
        w("      %-10s %s  %8d px" % (g, c, near(c)))
    w("   -- 区域归因 / 阈值线色 --")
    w("      region/nr  %s  %8d px" % (S.FS.C["region"], near(S.FS.C["region"])))
    w("   >>> 背景 6 色全部上图像 = %s" % allbg)
except Exception as e:
    w("   [!! 像素验证失败] %r" % (e,))

with io.open(LOG, "w", encoding="utf-8", newline="\n") as fh:
    fh.write("\n".join(L) + "\n")
print("done -> %s" % LOG)
