# -*- coding: utf-8 -*-
"""s108_fig1_overview_v2.py —— Fig1 **V2「图形化研究设计总览」**（graphical-abstract 版）

★ 背景（用户 2026-09-30 指令）：
  V1（`s107_fig1_overview.py`）保留不动；另出一版 ——
  「文字太多、信息量很大、排版拥挤；有的文字是不是可以用图形或一些图标表示？
    既简洁又能让人看懂我们做了什么、怎么做的」。

★ V2 相对 V1 的改动（**体例不变**：3 横带 x 3 列组 + 真实数据小图）：
  ① 每组标题行 = 「短语(左) + **大字号数字**(右，10 pt 粗体，取所在带强调色)」；
  ② V1 的 2–5 行 5.6 pt 说明句 -> 压成 **1 行 5.2 pt**（必要时 1 行 5.0 pt 灰字注）；
  ③ 关键量一律改由**大号数字**承担：
     980 / 14 / 40,024 / 26 / 5-5-5 / 12-14 / 8:1 / 52->0；
  ④ 阈值与方法名改为**图形化小胶囊 chip**（P < 5e-8 | r2 < 0.001 | 10 Mb）；
  ⑤ 14 个细胞类型：由「14 个带名字标签」改为 **14 个细胞点阵**（圆 + 核，按谱系着色）
     + 1 行谱系说明 —— 用图形换掉 14 个文字对象；
  ⑥ 列组之间加**极浅 chevron**（V1 没有）暗示带内左->右流程。
  ★ **完全复用 V1 的图标 / 小图函数**（`import s107_fig1_overview as V1`），
    V1 文件零改动，两者共用同一套真实数据与同一套幅宽反解。

★ 版面自检：直接调用 `V1.audit()`（同一套 BANDS / COLS 网格，四道判据：
  文字越出画布 / 小图足迹越出画布 / 越出列组 / 越出所属横带）。
  `set FIG1_AUDIT=1` 触发打印。

★ 幅宽：与 V1 同一反解公式（`_W_TARGET_MM=182.4`）⇒ `scale=1.0` 直接落 182.37 mm。

用法（★ **必须显式设 FIG_OUTDIR，否则会静默覆写 figures/ 正本**）：
  set FIGPAL=gb12 & set FIG_NOTES=0 & set FIG_OUTDIR=D:\\_transfer_logs\\_fig1_draft_v2
  python s108_fig1_overview_v2.py
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import matplotlib                                          # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt                            # noqa: E402
from matplotlib.patches import Circle, FancyBboxPatch, Rectangle   # noqa: E402

import s107_fig1_overview as V1                            # noqa: E402  ← V1 零改动
FS = V1.FS
CW, CH = V1.CW, V1.CH
FIGW, FIGH = V1.FIGW, V1.FIGH
_W_TARGET_MM, _W_LO_MM, _W_HI_MM = V1._W_TARGET_MM, V1._W_LO_MM, V1._W_HI_MM

BAND_X0, BAND_X1 = V1.BAND_X0, V1.BAND_X1
COLS = V1.COLS
WASH = V1.WASH
ACC = V1.ACC
EDGE, BODY = V1.EDGE, V1.BODY
hx, blend = V1.hx, V1.blend

# ★ 横带几何沿用 V1（这样 V1.audit 的 BANDS / COLS 判据可直接复用）
BANDS = V1.BANDS                                    # a(42.2,62.0) b(22.1,40.7) c(2.0,20.6)
# 2026-10-01（图件细改 · 问题3）：面板字母统一改为**大写**（用户明确要求，
#   注意与 Nature 官方的小写规范相反；回退只需把这里的三处前缀改回小写）。
#   key 仍用小写，因为 ACC[key] 的颜色索引依赖它。
BAND_LABEL = {"a": "A. Data", "b": "B. Discovery", "c": "C. Validation"}

# 每带的标题行 / 副题行 y（自上而下）
HEAD_Y = {"a": 57.4, "b": 36.9, "c": 16.3}
SUB_Y = {"a": 55.9, "b": 35.6, "c": 14.6}

# 谱系着色（14 个细胞点阵按 7 个谱系着色）
LIN_COLOR = {
    "B": FS.C["region"], "CD4T": FS.C["nominal"], "CD8T": FS.C["sig"],
    "NK": FS.C["mhc"], "Mono": FS.C["locus2"], "DC": FS.C["risk_down"],
    "Plasma": FS.C["green_dark"],
}
CELLS_R1 = ["B", "B", "CD4T", "CD4T", "CD4T", "CD8T", "CD8T"]
CELLS_R2 = ["CD8T", "NK", "NK", "Mono", "Mono", "DC", "Plasma"]


# ======================================================================
# 版面元件
# ======================================================================
def head(ax, x0, x1, y, name, big, col):
    """组标题行：短语左对齐 + 大字号数字右对齐（同一基线，省一行高度）。"""
    ax.text(x0, y, name, ha="left", va="top", fontsize=7.0, fontweight="bold",
            color=BODY, zorder=4)
    ax.text(x1, y, big, ha="right", va="top", fontsize=10.0, fontweight="bold",
            color=col, zorder=4)


def sub(ax, cx, y, text, color=None):
    """副题：一行 5.2 pt（V2 每组的说明句只剩这一行）。"""
    ax.text(cx, y, text, ha="center", va="top", fontsize=5.2,
            color=color or "#555555", zorder=4)


def chip(ax, cx, cy, text, fc, ec, tc="#202020", fs=5.0):
    """图形化小胶囊（承载阈值 / 方法名，替代整句文字）。"""
    w = 0.50 * (fs / 5.0) * len(text) + 1.9
    h = 1.9 * (fs / 5.0)
    ax.add_patch(FancyBboxPatch((cx - w / 2, cy - h / 2), w, h,
                                boxstyle="round,pad=0,rounding_size=0.55",
                                fc=fc, ec=ec, lw=0.5, zorder=3))
    ax.text(cx, cy, text, ha="center", va="center", fontsize=fs, color=tc, zorder=4)
    return w


def chip_row(ax, x0, x1, cy, items, fc, ec, fs=5.0, gap=1.3):
    """一行胶囊，按文本宽度自适应居中排布。"""
    ws = [0.50 * (fs / 5.0) * len(t) + 1.9 for t in items]
    total = sum(ws) + gap * (len(items) - 1)
    x = x0 + ((x1 - x0) - total) / 2.0
    for t, w in zip(items, ws):
        chip(ax, x + w / 2, cy, t, fc, ec, fs=fs)
        x += w + gap


def dot_cell(ax, cx, cy, r, col, ec):
    """单个细胞图标：胞体 + 细胞核。"""
    ax.add_patch(Circle((cx, cy), r, fc=hx(blend(col, 0.62)), ec=ec, lw=0.55, zorder=3))
    ax.add_patch(Circle((cx, cy), r * 0.36, fc=ec, ec="none", alpha=0.85, zorder=4))


def chevron(ax, cx, cy, s=0.9, col="#C2C2C2"):
    """极浅 chevron，暗示带内左->右的推进顺序（不承载数据）。"""
    ax.plot([cx - s / 2, cx + s / 2, cx - s / 2], [cy + s, cy, cy - s],
            color=col, lw=0.9, solid_capstyle="round", zorder=2)


def glyph_female(ax, cx, cy, r=1.0, fc=None, ec=None):
    """女性符号 icon（圆 + 竖杆 + 横杆）：承担 band a 的结局节点与 b2 的结局节点。

    ★ 为什么**放弃画子宫**（三轮实测，均为小尺寸 2–3 逻辑单位下的失败）：
      ① V1 `glyph_uterus`（输卵管弧 + 卵巢椭圆在宫体**顶端**）-> 读成一张**脸**；
      ② 附件下移到宫体两侧 -> 读成**钥匙 / 蘑菇**；
      ③ 宫底凹口压深 + 短弧连宫角 -> 读成**心形 / 气球**，且短弧肉眼仍与宫体断开、
         读作悬空的「眉毛」。
      结论：这个尺度下**不存在**可靠的器官剪影 => 改用几何符号，语义交给标签行承担。
    """
    fc = fc or hx(blend(FS.C["sig"], 0.30))
    ec = ec or FS.C["sig"]
    ax.add_patch(Circle((cx, cy + 0.52 * r), r, fc=fc, ec=ec, lw=0.7, zorder=3))
    ax.plot([cx, cx], [cy - 0.48 * r, cy - 1.05 * r], color=ec, lw=0.8,
            solid_capstyle="round", zorder=3)
    ax.plot([cx - 0.58 * r, cx + 0.58 * r], [cy - 0.74 * r, cy - 0.74 * r],
            color=ec, lw=0.8, solid_capstyle="round", zorder=3)


# ======================================================================
# 主图
# ======================================================================
def fig1(scale=1.0):
    FS.setup()
    fig = plt.figure(figsize=(FIGW * scale, FIGH * scale), facecolor="white")
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, CW)
    ax.set_ylim(0, CH)
    ax.axis("off")
    # 背景矩形：保证 tight bbox = 整幅（否则 inset axes 的 figure 分数坐标会错位）
    ax.add_patch(Rectangle((0, 0), CW, CH, fc="white", ec="none", zorder=-20))

    for key, _title, y0, y1 in BANDS:
        ax.add_patch(FancyBboxPatch((BAND_X0, y0), BAND_X1 - BAND_X0, y1 - y0,
                                    boxstyle="round,pad=0,rounding_size=0.9",
                                    fc=WASH[key], ec=EDGE, lw=0.7, zorder=-10))
        # P1（2026-10-01）：面板字母/带标题 8.6 -> **8.0 pt**（Nature：面板字母 8 pt 粗体）。
        #   注：本处为「字母 + 带名」合并写法（"a. Data"），字号按面板字母口径取 8 pt。
        ax.text(BAND_X0 + 1.8, y1 - 1.0, BAND_LABEL[key], ha="left", va="top",
                fontsize=8.0, fontweight="bold", color=ACC[key], zorder=4)
        ax.plot([BAND_X0 + 1.2, BAND_X1 - 1.2], [y1 - 3.35, y1 - 3.35],
                color=ACC[key], lw=0.7, alpha=0.55, zorder=1)

    # 列组间的极浅 chevron（x 间隙 33.0–34.4 / 65.0–66.4）
    # ★ a 带不加：它是 2 组、间隔 ~16 单位，chevron 会孤立悬在留白中（改用下面的长箭头）
    for cx in ((COLS[0][1] + COLS[1][0]) / 2, (COLS[1][1] + COLS[2][0]) / 2):
        for cy in (29.5, 9.0):
            chevron(ax, cx, cy, s=0.85)

    # ================================================== a1  sc-eQTL 暴露
    c1m = (COLS[0][0] + COLS[0][1]) / 2          # 17.9
    head(ax, COLS[0][0] + 2.0, COLS[0][1] - 2.0, HEAD_Y["a"],
         "sc-eQTL exposure", "980", ACC["a"])
    sub(ax, c1m, SUB_Y["a"], "PBMC donors; 14 immune cell types; 1 cohort")
    for row, cy in ((CELLS_R1, 51.6), (CELLS_R2, 48.8)):
        for i, lin in enumerate(row):
            cx = 7.2 + i * 3.3
            dot_cell(ax, cx, cy, 1.05, LIN_COLOR[lin], LIN_COLOR[lin])
    sub(ax, c1m, 46.4, "7 lineages: B / CD4T / CD8T / NK / Mono / DC / Plasma")
    chip(ax, c1m, 44.4, "replication: 119 donors (nominal consistency only)",
         "#FFFFFF", EDGE, tc="#6A6A6A", fs=4.8)

    # ================================================== a2  结局 GWAS
    c2m = (COLS[1][0] + COLS[1][1]) / 2          # 49.7
    c2h = (COLS[2][0] + COLS[2][1]) / 2          # 82.0
    # 由 a1 指向 a2 的极浅长箭头：把「暴露 -> 结局」的设计关系画出来，
    # 同时填补两组之间的留白（V1 此处为纯空白 + V2 首轮曾放 chevron）
    V1.arrow(ax, (32.0, 49.4), (45.0, 49.4), color=hx(blend(FS.C["cyan"], 0.30)),
             lw=1.1, ms=9.0)
    head(ax, 58.4, 91.6, HEAD_Y["a"], "GWAS outcome", "40,024", ACC["a"])
    sub(ax, 74.1, SUB_Y["a"], "cases; 665,658 controls; female infertility; EUR")
    # ★ 子宫图标必须与图轴（含 -log10 P 轴标签，向左外扩 ~3.1 单位）拉开距离，
    #   否则会紧贴图轴（V2 首轮实测只余 0.45 单位）
    glyph_female(ax, 49.2, 48.9, r=1.15)
    V1.mini_regional(fig, 56.6, 45.6, 35.0, 7.0)
    sub(ax, 74.1, 43.4, "chr1p36.12 (Mb, GRCh38); lead rs56318008", "#7A7A7A")

    # ================================================== b1  工具筛选
    head(ax, COLS[0][0] + 2.0, COLS[0][1] - 2.0, HEAD_Y["b"],
         "Instrument selection", "26", ACC["b"])
    sub(ax, c1m, SUB_Y["b"], "instruments; one per gene-cell pair; all F > 10")
    V1.glyph_dna(ax, 6.0, 31.0, 8.4, 2.2)
    V1.arrow(ax, (14.8, 32.1), (18.2, 32.1), ms=5.2)
    V1.glyph_gene(ax, 18.6, 31.2, 8.4, 1.8, "Gene")
    chip_row(ax, COLS[0][0] + 2.6, COLS[0][1] - 2.6, 29.6,
             ["P < 5e-8", "r2 < 0.001", "10 Mb"], "#FFFFFF", EDGE, fs=4.8)
    # ★ x=7.0（不是 6.0）：小图足迹含 "tests" 轴标签会向左外扩 4.2 单位，
    #   x=6.0 时足迹左缘落到 1.8 < 列组左界 2.8（V2 首轮审计抓到）。
    _f = V1.mini_fstat(fig, 7.0, 23.6, 24.0, 4.8)
    # ★ V1 把 "F statistic" 画在轴内右下角，V2 更矮的托盘下会贴住 x 刻度标签 ->
    #   本图实例内隐藏（改由副题行 "all F > 10" 承担轴含义）；V1 源码零改动。
    for _t in _f.texts:
        if _t.get_text() == "F statistic":
            _t.set_visible(False)

    # ================================================== b2  两样本 MR
    head(ax, COLS[1][0] + 2.0, COLS[1][1] - 2.0, HEAD_Y["b"],
         "Two-sample MR", "1", ACC["b"])
    sub(ax, c2m, SUB_Y["b"], "SNP per pair; Wald ratio only; Steiger filtered")
    V1.glyph_dna(ax, 36.0, 30.7, 6.0, 2.2, n_rung=5)
    sub(ax, 39.0, 29.9, "eQTL", "#5A5A5A")
    V1.arrow(ax, (42.4, 31.8), (45.4, 31.8), ms=5.2)
    xx = np.linspace(0, 1, 60)
    ax.plot(45.8 + xx * 5.4, 31.8 + 0.95 * np.sin(xx * 5.4), color=FS.C["sig"], lw=0.75)
    for xd, yd in ((46.4, 32.75), (48.5, 30.85), (50.6, 32.65)):
        ax.add_patch(Circle((xd, yd), 0.24, fc=hx(blend(FS.C["moderate_sh"], 0.2)),
                            ec=FS.C["sig"], lw=0.4))
    sub(ax, 48.5, 29.9, "expression", "#5A5A5A")
    V1.arrow(ax, (51.4, 31.8), (54.2, 31.8), ms=5.2)
    glyph_female(ax, 56.8, 31.35, r=0.85)
    sub(ax, 56.8, 29.9, "infertility", "#5A5A5A")
    # 混杂：左侧虚线 + 中间标签（带底色压线）+ 右侧绿叉 + 向下虚线箭头
    ax.plot([43.4, 46.4], [27.9, 27.9], ls=(0, (2, 1.6)), color="#7A7A7A", lw=0.6)
    ax.plot([50.6, 53.8], [27.9, 27.9], ls=(0, (2, 1.6)), color="#7A7A7A", lw=0.6)
    ax.text(48.5, 27.9, "confounders", ha="center", va="center", fontsize=5.2,
            color=BODY, zorder=6,
            bbox=dict(boxstyle="square,pad=0.12", fc=WASH["b"], ec="none"))
    for xc in (43.4, 53.8):
        ax.plot([xc], [27.9], marker="x", ms=5.6, mew=1.5, color=FS.C["region"], zorder=6)
    V1.arrow(ax, (48.5, 27.6), (48.5, 26.2), color="#7A7A7A", lw=0.6, ms=4.6)
    chip_row(ax, COLS[1][0] + 2.6, COLS[1][1] - 2.6, 24.8,
             ["Wald ratio (primary)", "IVW / MR-Egger n/a"], "#FFFFFF", EDGE, fs=4.8)

    # ================================================== b3  共定位
    head(ax, COLS[2][0] + 2.0, COLS[2][1] - 1.8, HEAD_Y["b"],
         "Colocalization", "5 / 5 / 5", ACC["b"])
    sub(ax, c2h, SUB_Y["b"], "strong / moderate / no shared variants")
    V1.mini_coloc_curves(fig, 67.0, 26.6, 30.0, 7.2)
    chip_row(ax, COLS[2][0] + 2.6, COLS[2][1] - 2.6, 25.2,
             ["15 gene-cell pairs", "PP.H4 > 0.80"], "#FFFFFF", EDGE, fs=4.8)

    # ================================================== c1  区域标签
    head(ax, COLS[0][0] + 2.0, COLS[0][1] - 2.0, HEAD_Y["c"],
         "Region tagging", "12 / 14", ACC["c"])
    sub(ax, c1m, SUB_Y["c"], "cell types significant for both genes")
    _r = V1.mini_regtag(fig, 6.0, 5.8, 25.0, 6.2)
    # ★ V1 在轴内画了红色 "not in / panel" 两行字，位置与左上图例（in panel / significant）
    #   及 WNT4 的两个 "0" 标签互相叠压 —— V2 里隐藏它（信息已由下方注与 "0 0" 承担），
    #   只改本图实例，V1 源码零改动。
    for _t in _r.texts:
        if "not in" in _t.get_text():
            _t.set_visible(False)
    sub(ax, c1m, 3.6, "WNT4 not in the OneK1K panel (0 / 0)", "#7A7A7A")

    # ================================================== c2  SMR + HEIDI
    head(ax, COLS[1][0] + 2.0, COLS[1][1] - 2.0, HEAD_Y["c"],
         "SMR + HEIDI", "8 : 1", ACC["c"])
    sub(ax, c2m, SUB_Y["c"], "HEIDI-consistent cell types (CDC42 : LINC00339)")
    V1.mini_smr(fig, 37.6, 5.8, 24.6, 6.2)
    sub(ax, c2m, 3.6, "the two genes carry opposite effect directions", "#7A7A7A")

    # ================================================== c3  PheWAS 安全性
    head(ax, COLS[2][0] + 2.0, COLS[2][1] - 1.8, HEAD_Y["c"],
         "PheWAS safety", "52 -> 0", ACC["c"])
    sub(ax, c2h, SUB_Y["c"], "FDR < 0.05 signals before -> after conditioning")
    V1.mini_phewas(fig, 69.5, 6.6, 26.1, 5.0)
    sub(ax, c2h, 3.6, "2,469 endpoints x 2 instruments", "#7A7A7A")

    return fig


# ======================================================================
if __name__ == "__main__":
    # ★ 列组判据补正：V1.audit 的「越出列组」用的是模块级 COLS（= b/c 两带的 3 列组网格），
    #   而 a 带是 **2 组**（a1 占 2.8–30.4、a2 的图轴起于 56.6），不落在该网格上 ->
    #   把 a2 的真实归属范围补进判据列表（顺序：窄的在前，避免被宽范围抢先匹配而变松）。
    #   只改 V1 模块的这一个全局、仅在本进程内生效，跑完即还原；V1 源码零改动。
    _cols_orig = V1.COLS
    V1.COLS = list(_cols_orig) + [(44.0, 97.6)]
    s = 1.0
    try:
        for _ in range(8):
            f = fig1(s)
            w, _h = V1.measure_mm(f)
            if os.environ.get("FIG1_AUDIT") == "1":
                V1.audit(f)
            plt.close(f)
            print("   [fit] scale=%.4f -> %.2f mm" % (s, w))
            if _W_LO_MM <= w <= _W_HI_MM and abs(w - _W_TARGET_MM) <= 0.45:
                break
            s *= _W_TARGET_MM / w
        else:
            print("   [!! 标定未收敛]  最后一轮 %.2f mm" % w)
        paths, wmm, hmm = FS.save(fig1(s), "Fig1v2_study_design_graphical")
        print("   final scale=%.5f  %.1f x %.1f mm  %s" % (s, wmm, hmm, paths))
    finally:
        V1.COLS = _cols_orig
