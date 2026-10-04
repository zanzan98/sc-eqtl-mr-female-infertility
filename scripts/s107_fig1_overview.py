# -*- coding: utf-8 -*-
"""s107_fig1_overview.py —— Fig1「研究设计总览」重绘（分带式 study-design overview）

参照用户 2026-09-30 给出的 4 张参考图体例：
  ① 整幅由 3 个横向**分带**（band）构成，每带左上有 "a./b./c." 粗体小标题 + 细分隔线；
  ② 每带内 3 个**列组**，各有居中粗体组标题；
  ③ 组内 = 图标/示意 + **真实数据小图**（不使用伪造数据）；
  ④ 浅色底纹 + 细边框；带外不加箭头，流程由版面隐含。

★ 内容与数据源（全部现场从 tables/ 读取，图注可复核）：
  a1 单细胞 cis-eQTL（暴露）：OneK1K / 14 免疫细胞类型 / 980 donors；
     圆点 = 14 个细胞类型
  a2 结局 GWAS 汇总统计：`GCST90483463`（40,024 / 665,658, EUR）；敏感性 `GCST90483469`；
     小图 = chr1p36.12 区域关联（真实逐变异 -log10P，`_chr1_ld/chr1_plotdata.csv`）
  b1 工具变量筛选：P<5e-8、clumping r2<0.001 / 10 Mb；小图 = 真实 F 统计量分布
     （`10_discovery_MR_main.csv`，n=8,612，中位 61.0）
  b2 两样本 MR：cis-eQTL -> 基因表达 -> 女性不孕症；混杂（Confounders）用绿叉；
     单 SNP/对 ⇒ 只跑 Wald ratio，IVW/Egger 不可运行
  b3 共定位：示意密度曲线（概念图）+ 真实计数 `34_coloc_summary.csv`（15 对：5/5/5）
  c1 区域标签核对：`49b_task34_region_tag_summary.csv` —— 两工具在两基因上均 12/14 显著，
     而 `WNT4` 不在 OneK1K 面板内（0/0）⇒ 工具是**区域**标签
  c2 SMR + HEIDI：`46_smr_3p3_gene_summary.csv` —— SMR FDR<0.05 细胞数 11 / 8，
     HEIDI 一致细胞数 8 / 1（极不对称）
  c3 PheWAS 安全性扫描：`42b_phewas_safety_domains.csv` —— 52 个 FDR<0.05 信号按域
     （12 升高 / 40 降低）；条件化后全部消失

★★ 版面自检（硬纪律）：
  · 本脚本所有文字**不得越出画布** [0,100]x[0,63.2]，也**不得越出所在列组**；
  · 令 `FIG1_AUDIT=1` 运行可打印**每个文字对象**的逻辑坐标包围盒 + 越界/重叠报告。
  · 幅宽由 `_W_TARGET_MM` 反解 `FIGW`（见下），使 `bbox_inches='tight'` 的
    实测宽 = FIGW*25.4 + 2*pad = 182.4 mm，**不靠缩放标定硬凑**。

用法：
  python s107_fig1_overview.py                     # 渲染到 figures/（正式）
  set FIG_OUTDIR=<dir> 再运行                      # 渲染到预览目录（不动 figures/）
  set FIG1_AUDIT=1  再运行                         # 额外打印版面审计
"""
import csv
import io
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import _fig_style as FS                                   # noqa: E402
import matplotlib                                          # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt                            # noqa: E402
from matplotlib.patches import (Circle, Ellipse, FancyArrowPatch,   # noqa: E402
                                FancyBboxPatch, Polygon, Rectangle)
from matplotlib.path import Path                           # noqa: E402
from matplotlib.patches import PathPatch                    # noqa: E402
from matplotlib.text import Text                            # noqa: E402

T = os.path.join(FS.ROOT, "tables")
CW, CH = 100.0, 63.2                     # 画布逻辑坐标

# 幅宽反解：FS.save 用 bbox_inches='tight', pad_inches=0.02（见 _fig_style.py）。
# fig1() 里有一块铺满整幅的白底矩形 ⇒ tight bbox == 整幅画布，
#   实测宽(mm) = FIGW*25.4 + 2*pad_in*25.4
# 于是令 _W_TARGET_MM = 182.4 反解 FIGW，scale=1.0 时**直接**落标定值。
_SAVE_PAD_IN = 0.02
_W_TARGET_MM = 182.4
_W_LO_MM, _W_HI_MM = 181.5, 183.3
FIGW = _W_TARGET_MM / 25.4 - 2.0 * _SAVE_PAD_IN
FIGH = FIGW * CH / CW

BAND_X0, BAND_X1 = 1.2, 98.8             # 三个带的左右边界
COLS = [(2.8, 33.0), (34.4, 65.0), (66.4, 97.6)]   # 三个列组

ACC = {"a": FS.C["mhc"], "b": FS.C["sig"], "c": FS.C["region"]}
EDGE = "#B0B0B0"
BODY = "#3A3A3A"


def blend(h, t, base=(255, 255, 255)):
    c = tuple(int(h[i:i + 2], 16) for i in (1, 3, 5))
    return tuple(int(round(c[i] + (base[i] - c[i]) * t)) for i in range(3))


def hx(rgb):
    return "#%02X%02X%02X" % tuple(rgb)


WASH = {"a": hx(blend(FS.C["cyan"], 0.16)),
        "b": hx(blend(FS.C["risk_up"], 0.17)),
        "c": "#F3F3F3"}


def rd(fname):
    with io.open(os.path.join(T, fname), encoding="utf-8-sig", newline="") as fh:
        return list(csv.DictReader(fh))


def fnum(x, d=float("nan")):
    try:
        return float(x)
    except (TypeError, ValueError):
        return d


# ======================================================================
# 图标（全部用矢量图元手绘，无外部素材）
# ======================================================================
def glyph_dna(ax, x, y, w, h, c1=None, c2=None, n_rung=7):
    c1 = c1 or FS.C["cyan"]
    c2 = c2 or FS.C["sig"]
    t = np.linspace(0, 2.4 * math.pi, 160)
    xa = x + (t / t[-1]) * w
    ax.plot(xa, y + h / 2 + (h / 2) * np.sin(t), color=c1, lw=0.75, zorder=3)
    ax.plot(xa, y + h / 2 - (h / 2) * np.sin(t), color=c2, lw=0.75, zorder=3)
    for i in range(n_rung):
        tt = (i + 0.5) / n_rung * t[-1]
        xx = x + (tt / t[-1]) * w
        ax.plot([xx, xx], [y + h / 2 + (h / 2) * math.sin(tt),
                           y + h / 2 - (h / 2) * math.sin(tt)],
                color="#8A8A8A", lw=0.45, zorder=2)


def glyph_gene(ax, x, y, w, h, text="Gene", fc=None):
    fc = fc or hx(blend(FS.C["cyan"], 0.55))
    ax.add_patch(Rectangle((x, y), w, h, fc=fc, ec=FS.C["mhc"], lw=0.6, zorder=3))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center",
            fontsize=5.4, color="#1A1A1A", zorder=4)


def glyph_person(ax, cx, cy, s=1.0, fc=None, female=False, label=None):
    fc = fc or FS.C["mhc"]
    ax.add_patch(Circle((cx, cy + 0.92 * s), 0.34 * s, fc=fc, ec="none", zorder=3))
    ax.add_patch(Polygon([(cx - 0.42 * s, cy - 0.95 * s), (cx + 0.42 * s, cy - 0.95 * s),
                          (cx + 0.34 * s, cy + 0.42 * s), (cx - 0.34 * s, cy + 0.42 * s)],
                         closed=True, fc=fc, ec="none", zorder=3))
    ax.add_patch(Polygon([(cx - 0.34 * s, cy + 0.34 * s), (cx - 0.70 * s, cy - 0.55 * s),
                          (cx - 0.56 * s, cy - 0.60 * s), (cx - 0.24 * s, cy + 0.20 * s)],
                         closed=True, fc=fc, ec="none", zorder=3))
    ax.add_patch(Polygon([(cx + 0.34 * s, cy + 0.34 * s), (cx + 0.70 * s, cy - 0.55 * s),
                          (cx + 0.56 * s, cy - 0.60 * s), (cx + 0.24 * s, cy + 0.20 * s)],
                         closed=True, fc=fc, ec="none", zorder=3))
    if female:            # 裙摆
        ax.add_patch(Polygon([(cx - 0.30 * s, cy - 0.30 * s), (cx + 0.30 * s, cy - 0.30 * s),
                              (cx + 0.62 * s, cy - 0.98 * s), (cx - 0.62 * s, cy - 0.98 * s)],
                             closed=True, fc=fc, ec="none", zorder=4))
        ax.add_patch(Polygon([(cx - 0.34 * s, cy + 0.42 * s), (cx + 0.34 * s, cy + 0.42 * s),
                              (cx + 0.30 * s, cy - 0.30 * s), (cx - 0.30 * s, cy - 0.30 * s)],
                             closed=True, fc=fc, ec="none", zorder=4))
        ax.add_patch(Circle((cx, cy + 0.92 * s), 0.34 * s, fc=fc, ec="none", zorder=5))
    if label:
        ax.text(cx, cy - 1.35 * s, label, ha="center", va="top", fontsize=5.2)


def glyph_uterus(ax, cx, cy, s=1.0, fc=None, ec=None):
    """简化子宫：梨形宫体 + 宫颈 + 双侧输卵管与卵巢。"""
    fc = fc or hx(blend(FS.C["sig"], 0.30))
    ec = ec or FS.C["sig"]
    body = Path([(0.0, 1.00), (0.34, 0.94), (0.46, 0.62), (0.30, 0.05),
                 (-0.30, 0.05), (-0.46, 0.62), (-0.34, 0.94), (0.0, 1.00)],
                [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4,
                 Path.LINETO, Path.CURVE4, Path.CURVE4, Path.CURVE4])
    tr = matplotlib.transforms.Affine2D().scale(s).translate(cx, cy) + ax.transData
    ax.add_patch(PathPatch(body, transform=tr, fc=fc, ec=ec, lw=0.6, zorder=3))
    ax.add_patch(Polygon([(cx - 0.13 * s, cy + 0.06 * s), (cx + 0.13 * s, cy + 0.06 * s),
                          (cx + 0.13 * s, cy - 0.42 * s), (cx - 0.13 * s, cy - 0.42 * s)],
                         closed=True, fc=fc, ec=ec, lw=0.5, zorder=3))
    for sgn in (-1, 1):
        tube = Path([(sgn * 0.30 * s, 0.92 * s), (sgn * 0.62 * s, 1.10 * s),
                     (sgn * 0.86 * s, 1.06 * s), (sgn * 0.96 * s, 0.92 * s)],
                    [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4])
        ax.add_patch(PathPatch(tube, transform=tr, fc="none", ec=ec, lw=0.7, zorder=3))
        ax.add_patch(Ellipse((cx + sgn * 1.00 * s, cy + 0.86 * s), 0.34 * s, 0.24 * s,
                             fc=fc, ec=ec, lw=0.5, zorder=4))


def arrow(ax, p1, p2, style="-|>", color=None, lw=0.7, ls="-", ms=6.0):
    ax.add_patch(FancyArrowPatch(p1, p2, arrowstyle=style, mutation_scale=ms,
                                 lw=lw, color=color or "#4D4D4D", linestyle=ls,
                                 shrinkA=0, shrinkB=0, zorder=5))


# ======================================================================
# 小图（真实数据）
# ======================================================================
def iax(fig, x, y, w, h):
    return fig.add_axes([x / CW, y / CH, w / CW, h / CH])


def mini_regional(fig, x, y, w, h):
    """chr1p36.12 逐变异 -log10P（真实）。x 轴单位写在外层图注里，故不设 xlabel。"""
    d = rd(r"D:/endometriosis_project/_chr1_ld/chr1_plotdata.csv")
    pos = np.array([fnum(r["pos_38"]) for r in d]) / 1e6
    nl = np.array([fnum(r["neglogp"]) for r in d])
    ok = np.isfinite(pos) & np.isfinite(nl)
    pos, nl = pos[ok], nl[ok]
    ax = iax(fig, x, y, w, h)
    ax.scatter(pos, nl, s=0.7, c=FS.C["cyan"], lw=0, alpha=0.85, rasterized=True)
    top = int(np.argmax(nl))
    ax.scatter([pos[top]], [nl[top]], s=10, marker="D", c=FS.C["region"], lw=0.4,
               edgecolors="black", zorder=5)
    ax.axhline(-math.log10(5e-8), ls="--", lw=0.5, c=FS.C["grey"], zorder=2)
    ax.text(21.83, -math.log10(5e-8) + 0.4, "P = 5e-8", fontsize=5.2,
            color=FS.C["grey"], va="bottom", ha="left")
    ax.set_ylabel("-log10 P", fontsize=5.4)
    ax.set_xlim(21.80, 22.60)
    ymax = float(nl.max())
    ax.set_ylim(0.0, ymax * 1.16)
    # 基因位置：图内极少量方形标记 + 短点线（名称写在外层图注，避免与刻度冲突）
    for _g, p in (("WNT4", 22.448), ("LINC00339", 22.028), ("CDC42", 22.073)):
        ax.plot([p, p], [0.0, ymax * 0.055], ls=(0, (1.4, 1.4)), lw=0.55,
                color=FS.C["sig"], zorder=4)
        ax.plot([p], [ymax * 0.032], marker="s", ms=2.2, c=FS.C["sig"], zorder=5)
    ax.tick_params(labelsize=5.2, pad=1.2)
    ax.set_xticks([21.8, 22.0, 22.2, 22.4, 22.6])
    ax.set_yticks([0, 10, 20])          # ★ 显式刻度：否则自动定位器会造出视图外的
    FS.finalize(ax)                     #   幽灵刻度（如 30），污染版面审计
    return ax


def mini_fstat(fig, x, y, w, h):
    d = rd("10_discovery_MR_main.csv")
    f = np.array([fnum(r["F_stat"]) for r in d])
    f = f[np.isfinite(f)]
    ax = iax(fig, x, y, w, h)
    bins = np.logspace(np.log10(max(10.0, f.min())), np.log10(f.max()), 26)
    ax.hist(f, bins=bins, color=hx(blend(FS.C["mhc"], 0.28)), edgecolor=FS.C["mhc"], lw=0.35)
    ax.axvline(np.median(f), color=FS.C["region"], lw=0.8, ls="--")
    ax.text(np.median(f) * 1.35, ax.get_ylim()[1] * 0.78, "median F = %.0f" % np.median(f),
            fontsize=5.2, color=FS.C["region"])
    ax.set_xscale("log")
    ax.set_xticks([10, 100, 1000, 10000])   # ★ 显式刻度，避免出现视图外的 10^5 幽灵刻度
    ax.set_ylabel("tests", fontsize=5.4)
    # ★ 不用 set_xlabel：轴外 xlabel 会顶出 b 带下边界（实测足迹底 21.3 < 带底 22.1）；
    #   改为轴内右下角嵌入（该处柱高≈0，不会压柱）
    ax.text(0.985, 0.04, "F statistic", transform=ax.transAxes, ha="right", va="bottom",
            fontsize=5.0, color="#4D4D4D")
    ax.tick_params(labelsize=5.2, pad=1.2)
    FS.finalize(ax)
    return ax


def mini_coloc_curves(fig, x, y, w, h):
    """示意密度曲线：左 = 两性状 peak 不重合；右 = peak 重合。

    ★ 坐标约定（曾经写错并被抓到）：`tt` 是**以组中心 xc 为原点**的相对坐标，
      `d_a`/`d_b` 必须也是**相对偏移**，峰值落在 x = xc + d；
      绝不能把 xc 本身当成 d（那样峰会被推到 2*xc 处，左组的峰跑进右组、
      右组则整条退化成平线）。
    """
    ax = iax(fig, x, y, w, h)
    ax.set_axis_off()
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 1)
    tt = np.linspace(-2.2, 2.2, 300)

    def bump(t, d, sd=0.60, a=0.24):
        return a * np.exp(-0.5 * ((t - d) / sd) ** 2)

    for xc, kind in ((2.60, "left"), (7.65, "right")):
        d_a = -0.65 if kind == "left" else 0.0     # 实线（性状 1）
        d_b = +0.80 if kind == "left" else 0.0     # 虚线（eQTL）
        for yb in (0.68, 0.40):
            ax.plot([xc - 2.10, xc + 2.10], [yb, yb], color="#B0B0B0", lw=0.45, zorder=1)
            ax.plot(xc + tt, yb + bump(tt, d_a), color=FS.C["mhc"], lw=0.85, zorder=3)
            ax.plot(xc + tt, yb + bump(tt, d_b), color=FS.C["nominal"], lw=0.85,
                    ls=(0, (2.2, 1.4)), zorder=3)
        if kind == "right":
            for yb in (0.40, 0.68):
                ax.plot([xc, xc], [yb, yb + 0.28], color=FS.C["region"], lw=0.7,
                        ls=(0, (1.2, 1.2)), zorder=4)
    ax.text(2.60, 0.015, "no shared variants", ha="center", va="bottom", fontsize=5.0)
    ax.text(7.65, 0.015, "shared (PP.H4 > 0.80)", ha="center", va="bottom", fontsize=5.0)


def mini_regtag(fig, x, y, w, h):
    d = rd("49b_task34_region_tag_summary.csv")
    genes = ["WNT4", "LINC00339", "CDC42"]
    inst = sorted({r["SNP_rsID"] for r in d})
    get = {(r["gene"], r["SNP_rsID"]): (int(r["OneK1K_n_sig"]), int(r["OneK1K_n_available"]))
           for r in d}
    ax = iax(fig, x, y, w, h)
    ncol = len(inst)
    bw = 0.34
    for i, g in enumerate(genes):
        for j, s in enumerate(inst):
            sig, av = get.get((g, s), (0, 0))
            xb = i + (j - (ncol - 1) / 2.0) * (bw + 0.06)
            ax.bar(xb, av, width=bw, color="#DCDCDC", edgecolor="#9E9E9E", lw=0.35,
                   label="in panel" if (i == 0 and j == 0) else None)
            ax.bar(xb, sig, width=bw, color=FS.C["risk_up"] if g != "WNT4" else "#DCDCDC",
                   edgecolor="#4D4D4D", lw=0.35,
                   label="significant" if (i == 0 and j == 0) else None)
            ax.text(xb, max(av, 0.35) + 0.5, "%d" % sig if av else "0", ha="center",
                    va="bottom", fontsize=5.2,
                    color=BODY if av else FS.C["region"])
    ax.text(0, 2.4, "not in\npanel", ha="center", va="bottom", fontsize=5.0,
            color=FS.C["region"], linespacing=1.25)
    ax.set_xticks(range(len(genes)))
    ax.set_xticklabels(genes, fontsize=5.4)
    ax.set_ylabel("OneK1K cell types", fontsize=5.4)
    ax.set_ylim(0, 23.0)
    ax.set_yticks([0, 5, 10, 14])
    ax.tick_params(labelsize=5.2, pad=1.2)
    ax.legend(fontsize=5.2, loc="upper left", ncol=1, handlelength=1.0,
              handletextpad=0.35, labelspacing=0.25, borderpad=0.2, framealpha=0)
    FS.finalize(ax)
    return ax


def mini_smr(fig, x, y, w, h):
    d = [r for r in rd("46_smr_3p3_gene_summary.csv") if r["analysis"] == "main(5e-8)"]
    genes = [r["gene"] for r in d]
    smr = [int(r["n_smr_fdr05"]) for r in d]
    hei = [int(r["n_heidi_consistent"]) for r in d]
    ax = iax(fig, x, y, w, h)
    xp = np.arange(len(genes))
    ax.bar(xp - 0.19, smr, width=0.36, color=FS.C["cyan"], edgecolor="#4D4D4D", lw=0.35,
           label="SMR FDR < 0.05")
    ax.bar(xp + 0.19, hei, width=0.36, color=FS.C["mhc"], edgecolor="#4D4D4D", lw=0.35,
           label="HEIDI-consistent")
    for i in range(len(genes)):
        ax.text(xp[i] - 0.19, smr[i] + 0.3, "%d" % smr[i], ha="center", va="bottom", fontsize=5.2)
        ax.text(xp[i] + 0.19, hei[i] + 0.3, "%d" % hei[i], ha="center", va="bottom", fontsize=5.2)
    ax.set_xticks(xp)
    ax.set_xticklabels(genes, fontsize=5.4)
    ax.set_ylabel("number of cell types", fontsize=5.4)
    # ★ ylim 必须留到 17.5：否则右上图例（占 ~2.6 数据单位高）会压住 LINC00339 的 "8"
    ax.set_ylim(0, 17.5)
    ax.set_yticks([0, 5, 10, 15])
    ax.tick_params(labelsize=5.2, pad=1.2)
    # ★ 图例必须单行两列：图例尺寸由字号（磅）决定、与 ylim 无关，两行版会横跨
    #   LINC00339 的 "8" 标签所在位置（实测压字）；压成一行后纵向即让开。
    ax.legend(fontsize=5.2, loc="upper right", ncol=2, handlelength=1.0, handletextpad=0.35,
              columnspacing=0.7, labelspacing=0.25, borderpad=0.2, framealpha=0)
    FS.finalize(ax)
    return ax


def mini_phewas(fig, x, y, w, h):
    d = rd("42b_phewas_safety_domains.csv")
    lab = {"免疫/血液/肿瘤": "immune", "女性生殖/不孕": "reprod.",
           "妊娠/分娩/产褥": "preg.", "肌肉骨骼/结缔组织": "MSK", "其他": "other"}
    ntest = {"免疫/血液/肿瘤": "968", "女性生殖/不孕": "288", "妊娠/分娩/产褥": "222",
             "肌肉骨骼/结缔组织": "546", "其他": "2,914"}
    order = ["免疫/血液/肿瘤", "女性生殖/不孕", "妊娠/分娩/产褥", "肌肉骨骼/结缔组织", "其他"]
    idx = {r["domain"]: r for r in d}
    up = [int(idx[k]["n_risk_up"]) for k in order]
    dn = [int(idx[k]["n_risk_down"]) for k in order]
    ax = iax(fig, x, y, w, h)
    xp = np.arange(len(order))
    ax.bar(xp, up, width=0.62, color=FS.C["risk_up"], edgecolor="#4D4D4D", lw=0.35,
           label="risk-increasing")
    ax.bar(xp, dn, width=0.62, bottom=up, color=FS.C["risk_down"], edgecolor="#4D4D4D",
           lw=0.35, label="risk-decreasing")
    for i, k in enumerate(order):
        tot = up[i] + dn[i]
        ax.text(xp[i], tot + 0.7, "%d" % tot, ha="center", va="bottom", fontsize=5.2)
    # 每个域的检验数写在刻度标签第二行，省掉一条 xlabel（版面留给两条图注）
    ax.set_xticks(xp)
    ax.set_xticklabels([lab[k] + "\n(" + ntest[k] + ")" for k in order], fontsize=5.0,
                       linespacing=1.25)
    ax.set_ylabel("FDR < 0.05 signals", fontsize=5.4)
    ax.set_ylim(0, 40)
    ax.set_yticks([0, 10, 20, 30])
    ax.tick_params(labelsize=5.2, pad=1.2)
    ax.legend(fontsize=5.2, loc="upper right", handlelength=1.0, handletextpad=0.35,
              labelspacing=0.25, borderpad=0.2, framealpha=0)
    FS.finalize(ax)
    return ax


# ======================================================================
# 主图
# ======================================================================
BANDS = [("a", "a. Main data sources", 42.2, 62.0),
         ("b", "b. Primary analyses", 22.1, 40.7),
         ("c", "c. Follow-up analyses", 2.0, 20.6)]


def fig1(scale=1.0):
    FS.setup()
    fig = plt.figure(figsize=(FIGW * scale, FIGH * scale), facecolor="white")
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, CW)
    ax.set_ylim(0, CH)
    ax.axis("off")
    # 背景矩形：保证 tight bbox = 整幅（否则 inset axes 的 figure 分数坐标会错位）
    ax.add_patch(Rectangle((0, 0), CW, CH, fc="white", ec="none", zorder=-20))

    for key, title, y0, y1 in BANDS:
        ax.add_patch(FancyBboxPatch((BAND_X0, y0), BAND_X1 - BAND_X0, y1 - y0,
                                    boxstyle="round,pad=0,rounding_size=0.9",
                                    fc=WASH[key], ec=EDGE, lw=0.7, zorder=-10))
        ax.text(BAND_X0 + 1.8, y1 - 1.0, title, ha="left", va="top", fontsize=8.6,
                fontweight="bold", color=ACC[key], zorder=4)
        ax.plot([BAND_X0 + 1.2, BAND_X1 - 1.2], [y1 - 3.35, y1 - 3.35], color=ACC[key],
                lw=0.7, alpha=0.55, zorder=1)

    # ================================================== a1  列组 (2.8, 33.0) 邻域
    x0, x1 = COLS[0][0], 53.0
    xm = (x0 + x1) / 2
    ax.text(xm, 58.3, "Single-cell cis-eQTL datasets   (exposure)",
            ha="center", va="top", fontsize=7.6, fontweight="bold", color=BODY)
    ax.text(xm, 56.4, "OneK1K sc-eQTL  |  PBMC  |  980 donors  |  14 immune cell types",
            ha="center", va="top", fontsize=6.0)
    ax.text(xm, 55.3, "cis-eQTL re-run (TensorQTL), GRCh37",
            ha="center", va="top", fontsize=6.0)
    ax.text(xm, 54.2, "Replication: 1M-scBloodNL  |  9 merged cell types  |  N = 119  |  "
                      "nominal consistency only", ha="center", va="top", fontsize=5.6,
            color=FS.C["grey"])

    cells_r1 = [("B", "MEM"), ("B", "IN"), ("CD4", "NC"), ("CD4", "ET"),
                ("CD4", "SOX4"), ("CD8", "NC"), ("CD8", "ET")]
    cells_r2 = [("CD8", "S100B"), ("NK", ""), ("NK", "R"), ("Mono", "NC"),
                ("Mono", "C"), ("DC", ""), ("Plasma", "")]
    pal = [FS.C["region"], FS.C["nominal"], FS.C["sig"], FS.C["mhc"], FS.C["locus2"],
           FS.C["risk_down"], FS.C["green_dark"]]
    for row, cy in ((cells_r1, 51.7), (cells_r2, 46.3)):
        for i, (l1, l2) in enumerate(row):
            cx = 6.0 + i * 7.4
            ax.add_patch(Circle((cx, cy), 0.88, fc=hx(blend(pal[i], 0.55)),
                                ec=pal[i], lw=0.6, zorder=3))
            ax.text(cx, cy - 1.30, l1 + ("\n" + l2 if l2 else ""), ha="center", va="top",
                    fontsize=5.2, color=BODY, linespacing=1.25, zorder=4)

    # ================================================== a2  列组 (54.6, 97.4)
    ax0, ax1 = 54.6, 97.4
    ax.text((ax0 + ax1) / 2, 58.3, "GWAS summary statistics   (outcome)",
            ha="center", va="top", fontsize=7.6, fontweight="bold", color=BODY)
    cols = [ax0 + 0.4, ax0 + 16.6, ax0 + 28.8]
    ax.text(cols[0], 56.6, "Phenotype", ha="left", va="top", fontsize=5.6,
            fontweight="bold", color=BODY)
    ax.text(cols[1], 56.6, "Source", ha="left", va="top", fontsize=5.6,
            fontweight="bold", color=BODY)
    ax.text(cols[2], 56.6, "cases / controls", ha="left", va="top", fontsize=5.6,
            fontweight="bold", color=BODY)
    ax.plot([ax0 + 0.2, ax1 - 0.2], [55.9, 55.9], color=EDGE, lw=0.5)
    ax.text(cols[0], 55.3, "Female infertility", ha="left", va="top", fontsize=5.6)
    ax.text(cols[1], 55.3, "GCST90483463", ha="left", va="top", fontsize=5.6)
    ax.text(cols[2], 55.3, "40,024 / 665,658", ha="left", va="top", fontsize=5.6)
    ax.text(cols[0], 54.1, "sensitivity outcome", ha="left", va="top", fontsize=5.6,
            color=FS.C["grey"])
    ax.text(cols[1], 54.1, "GCST90483469", ha="left", va="top", fontsize=5.6,
            color=FS.C["grey"])
    ax.text(cols[2], 54.1, "(EUR, GRCh37)", ha="left", va="top", fontsize=5.6,
            color=FS.C["grey"])
    glyph_uterus(ax, 59.6, 47.4, s=1.80)
    mini_regional(fig, 65.0, 45.4, 32.0, 7.2)
    ax.text(81.0, 43.9, "chr1p36.12 regional association (Mb, GRCh38); lead rs56318008",
            ha="center", va="top", fontsize=5.4, color=FS.C["grey"])

    # ================================================== b1
    c1m = (COLS[0][0] + COLS[0][1]) / 2
    ax.text(c1m, 36.6, "Instrument selection",
            ha="center", va="top", fontsize=7.4, fontweight="bold", color=BODY)
    ax.text(c1m, 34.6,
            "cis-eQTL P < 5e-8  |  LD clumping r2 < 0.001, 10 Mb",
            ha="center", va="top", fontsize=5.6)
    ax.text(c1m, 33.5, "one instrument per gene-cell pair (26 / 26)",
            ha="center", va="top", fontsize=5.6)
    glyph_dna(ax, 5.6, 28.95, 11.0, 3.1)
    glyph_gene(ax, 20.6, 29.3, 9.6, 2.1, "Gene")
    arrow(ax, (17.1, 30.35), (20.2, 30.35), ms=5.2)
    mini_fstat(fig, 7.0, 24.8, 24.4, 4.1)

    # ================================================== b2
    c2m = (COLS[1][0] + COLS[1][1]) / 2
    ax.text(c2m, 36.6, "Two-sample Mendelian randomization",
            ha="center", va="top", fontsize=7.4, fontweight="bold", color=BODY)
    glyph_dna(ax, 36.2, 29.0, 6.2, 2.6, n_rung=5)
    ax.text(39.3, 28.5, "cis-eQTL", ha="center", va="top", fontsize=5.2)
    arrow(ax, (42.9, 30.3), (46.6, 30.3))
    xx = np.linspace(0, 1, 60)
    ax.plot(47.0 + xx * 5.6, 30.3 + 1.0 * np.sin(xx * 5.4), color=FS.C["sig"], lw=0.75)
    for xd, yd in ((47.6, 31.3), (49.8, 29.2), (52.1, 31.0)):
        ax.add_patch(Circle((xd, yd), 0.24, fc=hx(blend(FS.C["moderate_sh"], 0.2)),
                            ec=FS.C["sig"], lw=0.4))
    ax.text(49.8, 28.5, "gene expression", ha="center", va="top", fontsize=5.2)
    arrow(ax, (53.0, 30.3), (56.4, 30.3))
    glyph_person(ax, 58.4, 30.2, s=1.5, fc=FS.C["mhc"], female=True)
    glyph_person(ax, 61.6, 30.2, s=1.5, fc=FS.C["risk_down"])
    ax.text(60.0, 28.5, "female infertility", ha="center", va="top", fontsize=5.2)
    # 混杂：左右绿叉 + 中间断开（标签用带底色的小白块压住虚线，避免线穿字）
    ax.plot([43.0, 46.9], [34.3, 34.3], ls=(0, (2, 1.6)), color="#7A7A7A", lw=0.6)
    ax.plot([52.7, 56.0], [34.3, 34.3], ls=(0, (2, 1.6)), color="#7A7A7A", lw=0.6)
    ax.text(49.8, 34.3, "Confounders", ha="center", va="center", fontsize=5.4,
            color=BODY, zorder=6,
            bbox=dict(boxstyle="square,pad=0.12", fc=WASH["b"], ec="none"))
    ax.plot([49.8, 49.8], [34.1, 32.0], ls=(0, (2, 1.6)), color="#7A7A7A", lw=0.6)
    arrow(ax, (49.8, 34.1), (49.8, 32.0), style="-|>", color="#7A7A7A", lw=0.6, ms=5.0)
    for xc in (43.0, 56.0):
        ax.plot([xc], [34.3], marker="x", ms=6.0, mew=1.6, color=FS.C["region"], zorder=6)
    ax.text(c2m, 26.9, "Wald ratio (primary)", ha="center", va="top", fontsize=5.6)
    ax.text(c2m, 25.9, "IVW / MR-Egger not runnable (1 SNP per pair)", ha="center",
            va="top", fontsize=5.6)
    ax.text(c2m, 24.8, "Steiger filtering applied to every pair", ha="center", va="top",
            fontsize=5.6, color=FS.C["grey"])
    ax.text(c2m, 23.7, "instrument strength: F = 61 (median)  |  all F > 10", ha="center",
            va="top", fontsize=5.6, color=FS.C["grey"])

    # ================================================== b3
    c3m = (COLS[2][0] + COLS[2][1]) / 2
    ax.text(c3m, 36.6, "Colocalization validation",
            ha="center", va="top", fontsize=7.4, fontweight="bold", color=BODY)
    mini_coloc_curves(fig, 66.8, 27.2, 30.4, 7.6)
    ax.text(c3m, 26.5, "coloc.abf (primary), coloc.susie (diagnostic)",
            ha="center", va="top", fontsize=5.6)
    ax.text(c3m, 25.5, "15 gene-cell pairs", ha="center", va="top", fontsize=5.6)
    ax.text(c3m, 24.5, "5 strong / 5 moderate / 5 no shared", ha="center", va="top",
            fontsize=5.6)
    ax.text(c3m, 23.5, "threshold PP.H4 > 0.80", ha="center", va="top", fontsize=5.6,
            color=FS.C["grey"])

    # ================================================== c1
    ax.text(c1m, 16.6, "Region tagging (fine-mapping)",
            ha="center", va="top", fontsize=7.4, fontweight="bold", color=BODY)
    ax.text(c1m, 14.7, "susie_rss credible set (external LD reference)",
            ha="center", va="top", fontsize=5.6)
    ax.text(c1m, 13.7, "the instrument tags the REGION, not a single gene",
            ha="center", va="top", fontsize=5.6)
    mini_regtag(fig, 6.0, 5.8, 25.4, 6.5)
    ax.text(c1m, 4.4, "12 / 14 cell types significant for both genes",
            ha="center", va="top", fontsize=5.4, color=FS.C["grey"])

    # ================================================== c2
    ax.text(c2m, 16.6, "SMR + HEIDI", ha="center", va="top",
            fontsize=7.4, fontweight="bold", color=BODY)
    ax.text(c2m, 14.7, "cis-SMR + HEIDI per gene and cell type",
            ha="center", va="top", fontsize=5.6)
    ax.text(c2m, 13.7, "HEIDI-consistent: CDC42 8 / 11, LINC00339 1 / 8",
            ha="center", va="top", fontsize=5.6, color=FS.C["grey"])
    mini_smr(fig, 37.6, 5.8, 24.6, 6.5)
    ax.text(c2m, 4.4, "the two genes carry OPPOSITE effect directions", ha="center",
            va="top", fontsize=5.4, color=FS.C["grey"])

    # ================================================== c3
    ax.text(c3m, 16.6, "PheWAS safety scan of the chr1p36.12 region",
            ha="center", va="top", fontsize=7.4, fontweight="bold", color=BODY)
    ax.text(c3m, 14.7, "2,469 FinnGen R12 endpoints x 2 instruments",
            ha="center", va="top", fontsize=5.6)
    ax.text(c3m, 13.7, "4,938 tests; 52 signals at FDR < 0.05",
            ha="center", va="top", fontsize=5.6)
    mini_phewas(fig, 69.5, 7.2, 26.1, 5.2)
    ax.text(c3m, 4.6, "conditioning on the region lead removes all 52 signals",
            ha="center", va="top", fontsize=5.4, color=FS.C["region"])
    ax.text(c3m, 3.6, "no risk-increasing immune-domain signal within power",
            ha="center", va="top", fontsize=5.4, color=FS.C["grey"])

    return fig


# ======================================================================
# 版面审计 / 尺寸实测
# ======================================================================
def _lg(bb, w_px, h_px):
    """display 包围盒 -> 逻辑坐标 (lx0, lx1, ly0, ly1)。"""
    return (bb.x0 / w_px * CW, bb.x1 / w_px * CW,
            bb.y0 / h_px * CH, bb.y1 / h_px * CH)


def _ax_footprint(a, rend, w_px, h_px):
    """小图托盘足迹（含刻度 / 刻度标签 / 轴标签）—— 直接取 matplotlib 自有口径。

    ★ 比对过的三种写法：
      · `a.get_children()` 递归不到 tick 标签（它们挂在 XAxis/YAxis 之下）-> 严重低估；
      · `a.findobj()` 递归会把**视图外的幽灵刻度**（如 log 轴的 10^5、ylim 外的 30）
        以及 `set_axis_off()` 后仍 visible 的刻度一起算进来 -> 高估、误报；
      · `a.get_tightbbox(render)` 正确排除以上两者 —— 故采用它。
      残余偏差只有一处：matplotlib 恒建的 `title/_left_title/_right_title` 三个**空**
      Text 会把足迹顶端抬高 ~1.2 逻辑单位；本审计只在**横向**比对列组边界、
      纵向只比画布边界，故不受影响（且横向为空文本不贡献宽度）。
    """
    return _lg(a.get_tightbbox(rend), w_px, h_px)


def audit(fig, verbose=True):
    """双口径版面审计：
      ① 逐 Text（只取我手写的注释；tick 标签的 window_extent 会退化为空，故跳过）
      ② 逐小图 `get_tightbbox`（含刻度/轴标签在内的真实足迹）—— 这条才是判「托盘越界」的权威口径
    返回 (越出画布的文字, 越出所在列组的小图, 越出画布的小图)。
    """
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    w_px, h_px = fig.get_size_inches() * fig.dpi
    skip = set()
    axlist = [a for a in fig.axes]
    for a in axlist:
        skip.update(a.get_xticklabels())
        skip.update(a.get_yticklabels())
        skip.update(a.get_xticklabels(minor=True))
        skip.update(a.get_yticklabels(minor=True))
    rows = []
    for t in fig.findobj(match=Text):
        s = t.get_text()
        if not s.strip() or t in skip:
            continue
        bb = t.get_window_extent(renderer=rend)
        # ★ tick 标签（含 axis("off") 后不可见的主轴刻度）在本环境下 window_extent
        #   会退化成 ~1 px 的空框 -> 必须按像素阈值剔除，否则整份审计全是假阳性。
        if bb.width < 3.0 or bb.height < 3.0:
            continue
        rows.append((s,) + _lg(bb, w_px, h_px))
    bad_txt = [r for r in rows
               if r[1] < BAND_X0 - 0.4 or r[2] > BAND_X1 + 0.4
               or r[3] < 0.6 or r[4] > CH - 0.6]
    bad_ax, off_ax, off_band = [], [], []
    for a in axlist:
        pos = a.get_position()
        if pos.width >= 0.999 and pos.height >= 0.999:   # 主画布轴
            continue
        lx0, lx1, ly0, ly1 = _ax_footprint(a, rend, w_px, h_px)
        tag = "小图@[%.1f,%.1f]x[%.1f,%.1f]" % (pos.x0 * CW, (pos.x0 + pos.width) * CW,
                                                pos.y0 * CH, (pos.y0 + pos.height) * CH)
        rec = (tag, lx0, lx1, ly0, ly1)
        if lx0 < BAND_X0 - 0.2 or lx1 > BAND_X1 + 0.2 or ly0 < 0.6 or ly1 > CH - 0.6:
            bad_ax.append(rec)
        # ★ 托盘必须留在**所属横带**内（刻度/轴标签常把足迹顶出带边界，是易漏的一类缺陷）
        ycen = (pos.y0 + pos.height / 2) * CH
        for _k, _t, by0, by1 in BANDS:
            if by0 <= ycen <= by1:
                if ly0 < by0 - 0.25 or ly1 > by1 + 0.25:
                    off_band.append(rec + ("band %s [%.1f,%.1f]" % (_k, by0, by1),))
                break
        for cx0, cx1 in COLS:
            if (pos.x0 * CW) >= cx0 - 0.1 and (pos.x0 + pos.width) * CW <= cx1 + 0.1:
                if lx0 < cx0 - 0.35 or lx1 > cx1 + 0.35:
                    off_ax.append(rec + ("col %.1f-%.1f" % (cx0, cx1),))
                break
    if verbose:
        rows2 = sorted(rows, key=lambda r: -(r[2] - r[1]))
        print("   [audit] 手写注释 n=%d ；最宽 14 条：" % len(rows))
        for s, a, b, c, d in rows2[:14]:
            print("      w=%5.1f  x[%5.1f,%5.1f]  y[%5.1f,%5.1f]  %s"
                  % (b - a, a, b, c, d, s[:54].replace("\n", "/")))
        print("   [audit] 小图托盘足迹（含刻度/轴标签）：")
        for a in axlist:
            pos = a.get_position()
            if pos.width >= 0.999 and pos.height >= 0.999:
                continue
            lx0, lx1, ly0, ly1 = _ax_footprint(a, rend, w_px, h_px)
            lg = ""
            if a.get_legend() is not None:
                g = _lg(a.get_legend().get_window_extent(rend), w_px, h_px)
                lg = "  legend x[%.1f,%.1f] y[%.1f,%.1f]" % g
            print("      轴@[%5.1f,%5.1f]x[%5.1f,%5.1f] -> 足迹 x[%5.1f,%5.1f] y[%5.1f,%5.1f]%s"
                  % (pos.x0 * CW, (pos.x0 + pos.width) * CW, pos.y0 * CH,
                     (pos.y0 + pos.height) * CH, lx0, lx1, ly0, ly1, lg))
        print("   [audit] 越出画布: 文字 %d / 小图 %d ；越出列组 %d ；越出所属横带 %d"
              % (len(bad_txt), len(bad_ax), len(off_ax), len(off_band)))
        for r in bad_txt:
            print("      !! 文字 x[%5.1f,%5.1f] y[%5.1f,%5.1f]  %s"
                  % (r[1], r[2], r[3], r[4], r[0][:54].replace("\n", "/")))
        for r in bad_ax:
            print("      !! 小图 x[%5.1f,%5.1f] y[%5.1f,%5.1f]  %s"
                  % (r[1], r[2], r[3], r[4], r[0]))
        for r in off_ax:
            print("      !! 列组外 x[%5.1f,%5.1f] y[%5.1f,%5.1f]  %s  (%s)"
                  % (r[1], r[2], r[3], r[4], r[0], r[5]))
        for r in off_band:
            print("      !! 横带外 x[%5.1f,%5.1f] y[%5.1f,%5.1f]  %s  (%s)"
                  % (r[1], r[2], r[3], r[4], r[0], r[5]))
    return bad_txt, bad_ax, off_ax, off_band


def measure_mm(fig, dpi=600):
    import io as _io
    from PIL import Image
    buf = _io.BytesIO()
    fig.savefig(buf, format="png", dpi=dpi, bbox_inches="tight", pad_inches=_SAVE_PAD_IN,
                facecolor="white")
    buf.seek(0)
    im = Image.open(buf)
    return im.size[0] / dpi * 25.4, im.size[1] / dpi * 25.4


if __name__ == "__main__":
    s = 1.0
    for _ in range(8):
        f = fig1(s)
        w, _h = measure_mm(f)
        if os.environ.get("FIG1_AUDIT") == "1":
            audit(f)
        plt.close(f)
        print("   [fit] scale=%.4f -> %.2f mm" % (s, w))
        if _W_LO_MM <= w <= _W_HI_MM and abs(w - _W_TARGET_MM) <= 0.45:
            break
        s *= _W_TARGET_MM / w
    else:
        print("   [!! 标定未收敛]  最后一轮 %.2f mm" % w)
    paths, wmm, hmm = FS.save(fig1(s), "Fig1_study_design")
    print("   final scale=%.5f  %.1f x %.1f mm  %s" % (s, wmm, hmm, paths))
