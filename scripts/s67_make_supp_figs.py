#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""s67_make_supp_figs.py -- 补充图 Fig S1 / Fig S2 重做

图注依据（已定稿，逐字）：
  39b_技术文档_SM与SN_v8.md §7 补充信息 -> ### 补充图
    "图 S1 各细胞类型的工具变量效力与分析门控（各细胞类型的工具数、名义与 FDR 显著检验数，
      以及效力标记）。图 S2 与敏感性结局 GCST90483469 的对照（15/15 一致）。"
  38_英文正文_论文体重写_v5.md -> ### Supplementary figures
    "Fig. S1 Instrument power and analytical gating by cell type (instrument counts,
      nominally and FDR-significant test counts, and power flags).
     Fig. S2 Comparison with the sensitivity outcome GCST90483469 (15/15 concordant)."

数据源（不改动任何已有表）：
  Fig S1  tables/20_celltype_power.csv        (14 细胞类型: n_instruments/n_nominal/n_FDR05/power_flag)
  Fig S2  tables/39_sensitivity_vs_main_GCST90483469.csv  (b_main vs b_sens, both_FDR05)
          + tables/15_discovery_significant_with_locus.csv (主分析 se / locus_id)

规格：Nature 双栏，标定到 182.4 mm（硬窗口 181.5-183.3）；Arial；图内文字全 ASCII；
      600 dpi；FS.save 双重字形核验（告警 + cmap 逐字符）。
"""
import csv
import io
import json
import math
import os
import sys
from collections import OrderedDict

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import _fig_style as FS                       # noqa: E402
import matplotlib                             # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt               # noqa: E402
import matplotlib.gridspec as gridspec        # noqa: E402

FS.setup()
BASE = os.path.join(FS.ROOT, "tables")
LOGS = os.path.join(FS.ROOT, "logs")
os.makedirs(LOGS, exist_ok=True)
LOG = []
def W(s=""):
    LOG.append(str(s))
    print(s)

TARGET_MM, LO_MM, HI_MM = 182.4, 181.5, 183.3
# 单栏（Nature single column = 88.9 mm）。Fig. S4 是单面板图，按历史 FigD1 的紧凑比例出图，
# 不用整栏宽 —— 整栏宽会把单面板火山图拉成 2.66:1 的横条（用户 2026-09-29 指出）。
SC_TARGET_MM, SC_LO_MM, SC_HI_MM = 88.9, 88.0, 90.0
# 每张图各自的栏宽规格（供标定留下可审计的痕迹）
FIG_SPEC = {
    "FigS1_instrument_power_gating": (TARGET_MM, LO_MM, HI_MM),
    "FigS2_sensitivity_outcome_concordance": (TARGET_MM, LO_MM, HI_MM),
    "FigS3_manhattan_by_celltype": (TARGET_MM, LO_MM, HI_MM),
    "FigS4_discovery_volcano": (SC_TARGET_MM, SC_LO_MM, SC_HI_MM),
}
CAL = OrderedDict()

FLAG_C = {"signal": FS.C["mhc"], "LOW_POWER_DO_NOT_READ_AS_NULL": FS.C["nominal"]}
GENE_C = {"CDC42": FS.C["locus1"], "LINC00339": FS.C["locus2"],
          "YME1L1": FS.C["gene_ym"], "ANXA4": FS.C["gene_an"]}
LOCUS_ORDER = ["L1", "L2", "L3"]
S2GEO = {}                      # build_s2 回填的几何量（供日志核对）
CENTER_RESID = {}               # 2026-10-01：各行面板组居中后的残差(mm)（本文件暂未用）
ROW_GAP_RESID = {}              # 2026-10-01：各相邻行带的实测可见空白(mm)，供核验脚本读取
S3ROWS = []                     # 2026-10-01：build_s3 回填的逐行 axes 组（供调试/核验脚本用）


# ---------------------------------------------------------------- helpers
def load(f):
    with io.open(os.path.join(BASE, f), encoding="utf-8-sig", newline="") as fh:
        return list(csv.DictReader(fh))


def measure_mm(fig, dpi=600):
    from PIL import Image
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=dpi, bbox_inches="tight", pad_inches=0.02,
                facecolor="white")
    buf.seek(0)
    im = Image.open(buf)
    return im.size[0] / dpi * 25.4, im.size[1] / dpi * 25.4


def fit_save(build, name, target=TARGET_MM, tol=0.45, iters=10, lo=LO_MM, hi=HI_MM):
    s = 1.0
    w_mm = float("nan")

    def _good(x):
        return lo <= x <= hi and abs(x - target) <= tol

    for _ in range(iters):
        fig = build(s)
        w_mm, _h = measure_mm(fig)
        plt.close(fig)
        if _good(w_mm):
            break
        s *= target / w_mm
    paths, wmm, hmm = FS.save(build(s), name)
    tries = 0
    while not _good(wmm) and tries < 3:
        tries += 1
        s *= target / wmm
        paths, wmm, hmm = FS.save(build(s), name)
        W("   [re-fit %d] %s -> %.2f mm" % (tries, name, wmm))
    CAL[name] = round(s, 5)
    if not _good(wmm):
        W("   [!! WIDTH OUT OF RANGE] %s = %.2f mm (want %.1f-%.1f)"
          % (name, wmm, lo, hi))
    return paths, wmm, hmm


def fmt_int(v):
    return "%d" % v


def text_width_in(s, fontsize, family="Arial"):
    """用 TextPath 精确测量字符串宽度（英寸）。

    不用渲染器、与 dpi / figure 尺寸无关 —— 用于**程序化**决定多面板间距，
    避免「y 轴长标签伸进左面板造成文字压文字」（本图 Fig S2 首版实际踩中）。
    """
    from matplotlib.textpath import TextPath
    from matplotlib.font_manager import FontProperties
    tp = TextPath((0, 0), s, size=fontsize, prop=FontProperties(family=family))
    return tp.get_extents().width / 72.0


# ---------------------------------------------------------------- data
def get_s1():
    rows = load("20_celltype_power.csv")
    assert len(rows) == 14, "S1: expect 14 cell types, got %d" % len(rows)
    rows.sort(key=lambda r: -int(r["n_instruments"]))
    flags = sorted(set(r["power_flag"] for r in rows))
    assert flags == ["LOW_POWER_DO_NOT_READ_AS_NULL", "signal"], flags
    low = [r["cell_type"] for r in rows if r["power_flag"].startswith("LOW")]
    assert set(low) == {"DC", "NK_R", "Plasma", "CD4_SOX4"}, low
    W("[S1] 14 cell types OK ; low-power flagged: %s" % ", ".join(low))
    W("     instruments range = %s..%s ; nominal range = %s..%s ; FDR range = %s..%s"
      % (rows[-1]["n_instruments"], rows[0]["n_instruments"],
         min(r["n_nominal"] for r in rows), max(r["n_nominal"] for r in rows),
         min(r["n_FDR05"] for r in rows), max(r["n_FDR05"] for r in rows)))
    return rows


def get_s2():
    r39 = load("39_sensitivity_vs_main_GCST90483469.csv")
    both = [r for r in r39 if r["both_FDR05"].strip().lower() == "true"]
    sign_ok = [r for r in both if r["sign_consistent"].strip().lower() == "true"]
    assert len(r39) == 8598, "S2: expect 8598 rows, got %d" % len(r39)
    assert len(both) == 15 and len(sign_ok) == 15, (len(both), len(sign_ok))
    # main-analysis se + locus from the closed table
    r15 = load("15_discovery_significant_with_locus.csv")
    assert len(r15) == 15
    idx = {(r["gene"], r["cell_type"]): r for r in r15}
    out = []
    for r in both:
        k = (r["gene"], r["cell_type"])
        assert k in idx, "S2: no main-analysis row for %s" % (k,)
        m = idx[k]
        db = abs(float(m["b"]) - float(r["b_main"]))
        assert db < 1e-12, "S2: b_main mismatch %.3g for %s" % (db, k)
        out.append({
            "locus": m["locus_id"], "gene": r["gene"], "cell_type": r["cell_type"],
            "b_main": float(r["b_main"]), "se_main": float(m["se"]),
            "b_sens": float(r["b_sens"]), "se_sens": float(r["se_MR"]),
            "p_main": float(r["p_main"]), "qval_main": float(r["qval"]),
            "p_sens": float(r["p_sens"]), "p_FDR_sens": float(r["p_FDR"]),
            "F_stat": float(r["F_stat"]),
        })
    out.sort(key=lambda d: -d["b_main"])
    xs = [d["b_main"] for d in out]
    ys = [d["b_sens"] for d in out]
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    sxy = sum((a - mx) * (b - my) for a, b in zip(xs, ys))
    sxx = sum((a - mx) ** 2 for a in xs)
    syy = sum((b - my) ** 2 for b in ys)
    r = sxy / (sxx * syy) ** 0.5
    slope = sxy / sxx
    assert all(a * b > 0 for a, b in zip(xs, ys)), "S2: sign not consistent"
    W("[S2] rows=%d ; both_FDR05=15 ; sign-consistent=15/15 ; b_main==b(15_) OK" % len(r39))
    W("     Pearson r = %.6f ; slope = %.6f ; b_main range %.6f..%.6f"
      % (r, slope, min(xs), max(xs)))
    return out, r, slope


# ---------------------------------------------------------------- Fig S1
def build_s1(rows, qmin, scale):
    n = len(rows)
    fig = plt.figure(figsize=(7.55 * scale, 3.25 * scale))
    # ★ 2026-10-01 本轮：**保持 wspace=0.17**。
    #   面板字母改用全局标准后走 `anchor="corner"`（见文件末注释）：字母右缘距本面板绘图区
    #   左缘 6.41 mm、字母总占宽 ≈ 8.5 mm < 列间隙 6.60 mm 之外的部分落在**轴框上方**的
    #   空白带里（本图面板无 set_title，那里没有任何墨迹）⇒ 无需重排列间距。
    gs = gridspec.GridSpec(1, 4, width_ratios=[1.60, 1.0, 0.92, 1.04], wspace=0.17)
    y = list(range(n))[::-1]          # index 0 (largest) on top

    def style_ax(ax, labels=None):
        # ★ 必须先 set_ticks 再给 labels —— 反序会触发
        #   "set_ticklabels() should only be used with a fixed number of ticks"
        ax.set_yticks(y, labels=labels)
        ax.set_ylim(-0.7, n - 0.3)
        FS.finalize(ax)

    # (a) instruments, coloured by power flag
    axa = fig.add_subplot(gs[0, 0])
    vals = [int(r["n_instruments"]) for r in rows]
    cols = [FLAG_C[r["power_flag"]] for r in rows]
    axa.barh(y, vals, height=0.66, color=cols, edgecolor="white", linewidth=0.35, zorder=3)
    for yy, v in zip(y, vals):
        axa.text(v + max(vals) * 0.02, yy, fmt_int(v), va="center", ha="left", fontsize=5.6)
    axa.set_xlim(0, max(vals) * 1.16)
    axa.set_xlabel("Instruments (cis-eQTL genes)")
    style_ax(axa, [r["cell_type"] for r in rows])
    # 面板字母统一移到 `subplots_adjust` 之后（2026-10-01）：锚定版必须在**最终**布局下量。
    from matplotlib.patches import Patch
    axa.legend(handles=[Patch(facecolor=FLAG_C["signal"], edgecolor="none",
                              label="power: signal"),
                        Patch(facecolor=FLAG_C["LOW_POWER_DO_NOT_READ_AS_NULL"],
                              edgecolor="none", label="low power (not a null result)")],
               loc="lower right", bbox_to_anchor=(1.02, -0.02), fontsize=5.4,
               handlelength=1.0, handleheight=0.7, borderaxespad=0.1)

    # (b) nominal
    axb = fig.add_subplot(gs[0, 1])
    vn = [int(r["n_nominal"]) for r in rows]
    axb.barh(y, vn, height=0.66, color=FS.C["nominal"], edgecolor="white", linewidth=0.35, zorder=3)
    for yy, v in zip(y, vn):
        axb.text(v + max(vn) * 0.02, yy, fmt_int(v), va="center", ha="left", fontsize=5.6)
    axb.set_xlim(0, max(vn) * 1.20)
    axb.set_xlabel("Nominally significant tests")
    style_ax(axb, [])

    # (c) FDR < 0.05
    axc = fig.add_subplot(gs[0, 2])
    vf = [int(r["n_FDR05"]) for r in rows]
    axc.barh(y, vf, height=0.66, color=FS.C["sig"], edgecolor="white", linewidth=0.35, zorder=3)
    for yy, v in zip(y, vf):
        axc.text(v + 0.09, yy, fmt_int(v), va="center", ha="left", fontsize=5.6)
    axc.set_xlim(-0.12, max(vf) + 0.95)
    axc.set_xticks([0, 1, 2, 3])
    axc.set_xlabel("FDR < 0.05 tests")
    style_ax(axc, [])

    # (d) -log10(minimum q)  —— 2026-09-29 扩板（原 FigD4 独有量，正文其它图未给）
    axd = fig.add_subplot(gs[0, 3])
    vq = [-math.log10(max(qmin[r["cell_type"]], 1e-320)) for r in rows]
    axd.barh(y, vq, height=0.66, color=FS.C["mhc"], edgecolor="white", linewidth=0.35, zorder=3)
    axd.axvline(-math.log10(0.05), ls="--", lw=0.6, color=FS.C["sig"], zorder=2)
    axd.text(3.0, -math.log10(0.05) - 0.20, "FDR = 0.05", fontsize=5.4,
             color=FS.C["sig"], va="top", ha="left")
    axd.set_xlim(0, max(vq) * 1.10)
    axd.set_xlabel("-log10(minimum q)")
    style_ax(axd, [])

    fig.subplots_adjust(left=0.082, right=0.995, top=0.90, bottom=0.155)
    # ★★ 全局标准（2026-10-01）：「面板字母 ↔ 本面板 Y 轴最上面那个词」= Fig2 标准
    #    （Δx 6.41 / Δy 2.42 mm）。必须在 `subplots_adjust` **之后**锚定（要量最终位置）。
    #    ★ 本图是**单行** 4 面板 ⇒ 无"行带间距"可调；整幅墨迹中心 = 本行中心 ⇒ 无需居中。
    for _ax, _lb in ((axa, "a"), (axb, "b"), (axc, "c"), (axd, "d")):
        FS.panel_label_anchored(_ax, _lb, anchor="corner")
    return fig


# ---------------------------------------------------------------- Fig S2
def build_s2(rows, r, slope, scale):
    n = len(rows)
    fig_w_in, fig_h_in = 7.2 * scale, 3.35 * scale
    fig = plt.figure(figsize=(fig_w_in, fig_h_in))

    # ---- 程序化定间距：panel b 的 y 标签向左伸出，必须整段落在两面板之间 ----
    TICK_FS = 6.8                                   # panel b y 刻度字号
    row_labels = ["%s | %s" % (d["gene"], d["cell_type"]) for d in rows]
    widths = [text_width_in(t, TICK_FS) for t in row_labels]
    lab_w_in = max(widths)
    longest = row_labels[widths.index(lab_w_in)]
    # matplotlib rcParams: ytick.major.size=2.2pt + ytick.major.pad(默认3.5pt) + 呼吸
    gap_in = lab_w_in + (2.2 + 3.5) / 72.0 + 0.075
    LM, RM, BOT, HGT = 0.118, 0.006, 0.150, 0.745
    gap_f = gap_in / fig_w_in
    span = 1.0 - LM - RM - gap_f
    w1 = span / (1.0 + 1.42)
    w2 = span - w1
    assert gap_in > lab_w_in, "S2: panel gap smaller than the longest tick label"
    S2GEO.clear()
    S2GEO.update({"longest_label": longest, "label_w_mm": round(lab_w_in * 25.4, 2),
                  "gap_mm": round(gap_in * 25.4, 2),
                  "panelA_mm": round(w1 * fig_w_in * 25.4, 1),
                  "panelB_mm": round(w2 * fig_w_in * 25.4, 1)})

    # (a) concordance scatter
    axa = fig.add_axes([LM, BOT, w1, HGT])
    xs = [d["b_main"] for d in rows]
    ys = [d["b_sens"] for d in rows]
    lo, hi = min(xs + ys), max(xs + ys)
    pad = (hi - lo) * 0.10
    lo, hi = lo - pad, hi + pad
    axa.plot([lo, hi], [lo, hi], ls=(0, (3.2, 2.2)), lw=0.7, color="#9A9A9A", zorder=1)
    seen = []
    # ★ 全局「点的大小与风格」标准（_fig_style.py）：面积 ×2 + 同色系深 2 度描边。
    #   ★ 原 `mec="white"` 是"白边做分隔"；按新标准统一换成**同色系深 2 度**轮廓
    #     （白底会被浅色填充吃掉分隔作用，深轮廓在浅底上同样能把相邻点分开）。
    _msa = FS.pt_up(3.6)
    for d in rows:
        g = d["gene"]
        mk = {"L1": "o", "L2": "s", "L3": "^"}[d["locus"]]
        axa.plot([d["b_main"]], [d["b_sens"]], marker=mk, ms=_msa,
                 color=GENE_C[g], mec=FS.edge_color(GENE_C[g]), mew=FS.mew_for(_msa),
                 ls="none", zorder=3, label=(g if g not in seen else None))
        if g not in seen:
            seen.append(g)
    axa.set_xlim(lo, hi)
    axa.set_ylim(lo, hi)
    axa.set_aspect("equal", adjustable="box")
    axa.set_xlabel("MR effect, main outcome (beta)")
    axa.set_ylabel("MR effect, sensitivity outcome (beta)")
    axa.text(0.045, 0.955,
             "15/15 direction-consistent\nr = %.4f ; slope = %.3f\nall 15 pairs FDR < 0.05\nin both outcomes"
             % (r, slope),
             transform=axa.transAxes, va="top", ha="left", fontsize=5.6, linespacing=1.45)
    # ★ 图例标记**绝对尺寸不变**（图上点放大了 √2 ⇒ markerscale 同除 √2），见 _fig_style.py
    axa.legend(loc="lower right", fontsize=5.4, markerscale=FS.markerscale_keep(1.15),
               handletextpad=0.35, labelspacing=0.22, borderaxespad=0.25)
    FS.finalize(axa)

    # (b) paired forest: main vs sensitivity
    axb = fig.add_axes([LM + w1 + gap_f, BOT, w2, HGT])
    yy = list(range(n))[::-1]
    off = 0.17
    # ★ 全局「点的大小与风格」标准（_fig_style.py）：面积 ×2。
    #   实心（Main outcome）→ 加同色系深 2 度描边；空心（Sensitivity outcome）
    #   是**纯描边型**标记 ⇒ 保留原 mew=0.6，只把边线换成同色系深 2 度。
    _msb2 = FS.pt_up(2.9)
    for y0, d in zip(yy, rows):
        axb.plot([d["b_main"], d["b_sens"]], [y0, y0], lw=0.6, color="#C8C8C8", zorder=1)
        axb.errorbar([d["b_main"]], [y0 + off],
                     xerr=[[1.96 * d["se_main"]], [1.96 * d["se_main"]]],
                     fmt="o", ms=_msb2, color="#2B2B2B", mfc="#2B2B2B",
                     mec=FS.edge_color("#2B2B2B"), mew=FS.mew_for(_msb2),
                     elinewidth=0.6, capsize=1.5, capthick=0.6, zorder=3,
                     label=("Main outcome" if y0 == yy[0] else None))
        axb.errorbar([d["b_sens"]], [y0 - off],
                     xerr=[[1.96 * d["se_sens"]], [1.96 * d["se_sens"]]],
                     fmt="o", ms=_msb2, mew=0.6, color=FS.edge_color(FS.C["nominal"]),
                     mfc="white",
                     elinewidth=0.6, capsize=1.5, capthick=0.6, zorder=3,
                     label=("Sensitivity outcome" if y0 == yy[0] else None))
    axb.axvline(0, lw=0.6, color="#7F7F7F", zorder=0)
    axb.set_yticks(yy)
    axb.set_yticklabels(row_labels)
    axb.tick_params(axis="y", labelsize=TICK_FS)
    axb.set_ylim(-0.75, n - 0.25)
    axb.set_xlabel("MR effect (beta, 95% CI)")
    # ★ 图例标记**绝对尺寸不变**（图上点放大了 √2 ⇒ markerscale 同除 √2），见 _fig_style.py
    axb.legend(loc="lower right", fontsize=5.4, markerscale=FS.markerscale_keep(1.1),
               handletextpad=0.35, labelspacing=0.25, borderaxespad=0.25)
    FS.finalize(axb)

    # ★★ 全局标准（2026-10-01）：面板字母距离 = Fig2 标准（Δx 6.41 / Δy 2.42 mm）。
    #   ★ 本图两个面板都用 **anchor="corner"**：面板 (b) 的 y 刻度是**长分类标签**
    #     （"CDC42 | Mono_NC" ≈ 19 mm），按"词左缘 −6.41 mm"锚会把字母推到面板 (a)
    #     的数据区里 ⇒ 改锚到**本面板绘图区左上角**（口径不变，仍为 6.41 / 2.42 mm）。
    #   ★ 本图是**单行** 2 面板 ⇒ 无行带间距可调、整幅墨迹中心 = 本行中心 ⇒ 无需居中。
    for _ax, _lb in ((axa, "a"), (axb, "b")):
        FS.panel_label_anchored(_ax, _lb, anchor="corner")
    return fig


# ---------------------------------------------------------------- Fig S3 / S4
def get_discovery():
    """Discovery 全表（8,612 检验）+ 与 s13 完全一致的染色体内归一化位置 xpos。

    xpos 复刻 scripts/s13_discovery_figures.py 的算法：
        frac = tss_pos_grch38 / max(tss_pos_grch38 | 同一染色体)
        xpos = chr + frac * 0.86
    以保证图 S3 与历史诊断图逐点一致（可追溯）。
    """
    rows = load("10_discovery_MR_main.csv")
    out = []
    for r in rows:
        try:
            c = int(float(r["chr"]))
            pos = float(r["tss_pos_grch38"])
        except (ValueError, TypeError):
            continue
        out.append({"chr": c, "pos": pos, "p": float(r["p"]), "qval": float(r["qval"]),
                    "cell_type": r["cell_type"], "gene": r["gene"], "b": float(r["b"]),
                    "is_MHC": r["is_MHC"].strip().lower() == "true"})
    assert len(out) == 8612, "S3/S4: expect 8612 discovery tests, got %d" % len(out)
    maxpos = {}
    for x in out:
        maxpos[x["chr"]] = max(maxpos.get(x["chr"], 0.0), x["pos"])
    for x in out:
        x["neglogp"] = -math.log10(max(x["p"], 1e-320))
        x["xpos"] = x["chr"] + (x["pos"] / maxpos[x["chr"]]) * 0.86
    n_mhc = sum(1 for x in out if x["is_MHC"])
    n_sig = sum(1 for x in out if x["qval"] < 0.05)
    n_sig_mhc = sum(1 for x in out if x["is_MHC"] and x["qval"] < 0.05)
    bonf = 0.05 / len(out)                                  # 与正文 Fig. 2a 同一 Bonferroni 线
    bh = max(x["p"] for x in out if x["qval"] < 0.05)       # BH 阈值 = FDR<0.05 组内最大 p
    W("[S3/S4] tests=%d ; FDR<0.05=%d ; MHC=%d (%.2f%%) ; MHC_and_FDR=%d"
      % (len(out), n_sig, n_mhc, 100.0 * n_mhc / len(out), n_sig_mhc))
    W("          Bonferroni p = %.6g (-log10=%.4f) ; BH p = %.6g (-log10=%.4f)"
      % (bonf, -math.log10(bonf), bh, -math.log10(bh)))
    # 口径守卫：若 MHC 检验中出现 FDR 显著对，则"不区分 MHC"的画法就站不住了
    assert n_sig_mhc == 0, "S3/S4: MHC pairs are FDR-significant -> MHC colouring must be revisited"
    return out, bonf, bh


def build_s3(d, bonf, scale):
    """图 S3 · 各细胞类型分层的 Discovery 曼哈顿图（14 小倍数）。

    口径与正文 Fig. 2a **严格一致**：同一 8,612 个检验、同一 Bonferroni 虚线（p = 0.05/8612）、
    同一 FDR<0.05 判定（红菱形）。**不引入 MHC 着色**——正文未讨论 MHC，且 293/8612 = 3.40% 的
    MHC 检验中无一是 FDR 显著对（见 get_discovery 的守卫断言）。
    """
    cells = sorted(set(x["cell_type"] for x in d))
    assert len(cells) == 14, cells
    ymax = int(math.ceil(max(x["neglogp"] for x in d))) + 1
    fig = plt.figure(figsize=(7.40 * scale, 5.60 * scale))
    # ★ 2026-10-01 本轮：wspace 0.17 -> 0.22。
    #   原因：面板字母改用全局标准（anchor="corner"，字母右缘距本面板绘图区左缘 6.41 mm，
    #   字母总占宽 ≈ 8.5 mm）。原 wspace=0.17 时列间隙只有 6.64 mm ⇒ (b)(c)(d) 列的字母会
    #   侵入左邻面板绘图区 1.9 mm，而本图每个面板**都有 set_title**（居中对齐、右端距
    #   绘图区右缘仅约 2 mm）⇒ 余量只剩 0.1 mm，实测不可接受。
    #   0.22 使列间隙 ≈ 8.4 mm（面板宽由 39.1 降到 37.8 mm，-3.2%；数据/口径不变）。
    gs = gridspec.GridSpec(4, 4, left=0.058, right=0.995, top=0.955, bottom=0.062,
                           wspace=0.22, hspace=0.46)
    letters = "abcdefghijklmn"
    _rows = [[], [], [], []]                    # 逐行的 axes 组（供行带间距调整用）
    thr_y = -math.log10(bonf)
    for i in range(16):
        ax = fig.add_subplot(gs[i // 4, i % 4])
        _rows[i // 4].append(ax)
        if i >= len(cells):
            ax.axis("off")
            continue
        cell = cells[i]
        g = [x for x in d if x["cell_type"] == cell]
        ng = [x for x in g if x["qval"] >= 0.05]
        sg = [x for x in g if x["qval"] < 0.05]
        # ★ 全局「点的大小与风格」标准（_fig_style.py）：面积 ×2 + 同色系深 2 度描边。
        #   `scatter` 的 `s` 即"直径(pt)²" ⇒ `s_up()` 直接把它 ×2（等价直径 ×√2）；
        #   `lw` 用 `mew_for(等价直径)` 换算（`ms_of_s`），统一走同一条描边宽度规则。
        _s_ng = FS.s_up(2.0)
        ax.scatter([x["xpos"] for x in ng], [x["neglogp"] for x in ng],
                   s=_s_ng, c=FS.C["ns"], lw=FS.mew_for(FS.ms_of_s(_s_ng)),
                   edgecolors=FS.edge_color(FS.C["ns"]),
                   alpha=0.65, rasterized=True, zorder=2)
        if sg:
            _s_sg = FS.s_up(13)
            ax.scatter([x["xpos"] for x in sg], [x["neglogp"] for x in sg],
                       s=_s_sg, c=FS.C["sig"], lw=FS.mew_for(FS.ms_of_s(_s_sg)),
                       edgecolors=FS.edge_color(FS.C["sig"]), marker="D", zorder=5)
        ax.axhline(thr_y, ls="--", lw=0.5, color=FS.C["region"], zorder=1)
        ax.set_xlim(0.4, 23.2)
        ax.set_ylim(0, ymax)
        ax.set_yticks([0, 4, 8, 12])
        ax.set_xticks([1, 6, 12, 18, 22])
        ax.set_title("%s  (n = %d%s)" % (cell, len(g),
                                         ", FDR < 0.05 = %d" % len(sg) if sg else ""),
                     fontsize=5.9, pad=1.6)
        ax.set_xticklabels(["1", "6", "12", "18", "22"] if i >= 12 else [])
        ax.set_yticklabels(["0", "4", "8", "12"] if i % 4 == 0 else [])
        ax.tick_params(labelsize=5.3, length=1.8, pad=1.4)
        FS.finalize(ax)
        # 面板字母统一移到循环后（2026-10-01）：锚定版必须在**最终**布局下量。
    # ★★ 2026-10-01：x 0.010 -> 0.0059（左移 0.75 mm）。
    #   原因：本条是**跨 4 行的共用旋转 y 轴总标题**，它的墨迹（旋转后"水平"宽度 = 字高）
    #   占 x 2.19–3.69 mm，而第 1 列面板字母（`anchor="corner"`：字母右缘 = 绘图区左缘 − 6.41 mm）
    #   占 x 1.90–4.27 mm ⇒ **两者抢同一条左侧留白带**。逐字实测（`_s173`）：字母 I（窄，
    #   盒宽 0.85 mm ⇒ 墨迹从 3.6 mm 起）的行 y 恰好落在总标题的 y 跨度 55.8–89.2 mm 内
    #   ⇒ 与总标题的 ")" 字形**相碰 0.26 mm**（`_s145` 以 bbox 口径报 0.61 × 2.65 mm）。
    #   左移 0.75 mm 后总标题墨迹变为 1.44–2.94 mm ⇒ 对字母 I 留出 0.66 mm 净空。
    #   ★ 幅面影响：tight 左缘由 1.84 -> 1.09 mm（+0.75 mm），`fit_save` 会把宽度重新标定回
    #     182.4 mm（scale 略降），行带由 `set_row_ink_gaps` 重新收敛 ⇒ 无需手工补偿。
    fig.text(0.0059, 0.52, "-log10(P)   (MR, single-SNP Wald)", rotation=90,
             ha="left", va="center", fontsize=6.2)
    # 图例放在第 15 个空槽内 —— 放在图外会与 (d) 面板标题发生文字压文字（首版实测）
    # ★ `gs[3, 3]` 已在上面循环的 i=15 轮里建过同一个 axes（`add_subplot` 对同一
    #   GridSpec 槽位返回同一对象）⇒ 此处**不再**重复 append 到 `_rows[3]`。
    axl = fig.add_subplot(gs[3, 3])
    axl.axis("off")
    axl.text(0.02, 0.92, "Thresholds (identical to Fig. 2a):", transform=axl.transAxes,
             fontsize=5.8, va="top", ha="left", color="#7F7F7F")
    axl.plot([0.03, 0.21], [0.66, 0.66], ls="--", lw=0.6, color=FS.C["region"],
             transform=axl.transAxes, clip_on=False)
    axl.text(0.25, 0.66, "Bonferroni  p = 0.05/8612", transform=axl.transAxes,
             fontsize=5.8, va="center", ha="left", color=FS.C["dark"])
    axl.plot([0.12], [0.38], marker="D", ms=3.4, ls="", color=FS.C["sig"],
             transform=axl.transAxes, clip_on=False)
    axl.text(0.25, 0.38, "FDR < 0.05", transform=axl.transAxes,
             fontsize=5.8, va="center", ha="left", color=FS.C["dark"])

    # ★★ 全局标准（2026-10-01）：「面板字母 ↔ 本面板图形」= Fig2 标准（Δx 6.41 / Δy 2.42 mm）。
    #   ★ 本图用 **anchor="corner"**：左列虽有数值刻度（"0/4/8/12"），但 4 列若混用两种锚
    #     会让同一行的字母横向差 4.5 mm（**同行不齐**），且词锚会把第 1 列字母推出画布
    #     2.1 mm。统一锚到绘图区左上角 ⇒ 同行同列完全对齐，且都在画布内。
    #   ★★ 调用时机：必须在**所有** axes（含末行 legend）都建好之后。
    for _i, _lb in enumerate(letters):
        FS.panel_label_anchored(_rows[_i // 4][_i % 4], _lb, anchor="corner")

    # ★★ 全局标准（2026-10-01）：「图与图之间的距离」= Fig2 的可见空白带 4.74 mm。
    #   ★ 4x4 网格 ⇒ 3 个行边界，**逐边界测量并竖直平移下方各行**（单一 hspace 不可能同时
    #     满足，见 `_fig_style.set_row_ink_gaps` 的说明）。
    ROW_GAP_RESID["figS3"] = FS.set_row_ink_gaps(fig, _rows, [4.74, 4.74, 4.74])
    S3ROWS.clear()
    S3ROWS.extend([[a for a in g] for g in _rows])

    # 全局 x 标签：必须排在行带调整**之后** —— 末行被竖直平移后，其刻度标签位置已变。
    _y4 = min(a.get_position().y0 for a in _rows[3])
    fig.text(0.5, max(0.012, _y4 - 0.032), "chromosome", ha="center", va="bottom",
             fontsize=6.2)
    return fig


def build_s4(d, bonf, bh, scale):
    """图 S4 · Discovery MR 效应量分布（火山图）。

    x = MR 效应量（beta），y = -log10(P)；两条阈值线对应正文所用的两种口径
    （Bonferroni p = 0.05/8612，与 Fig. 2a 一致；BH p 阈值对应 FDR < 0.05）。
    与 Fig. 2a（x = 基因组位置）互补：本图展示效应量的**方向分布**。

    ★ 版式（2026-09-29 按用户意见改）：**单栏、近方形**，沿用历史 FigD1 的画布比例
      （4.60 x 3.50 in，tight bbox 后约 1.22:1）。原先按整栏（182.4 mm）出图会把
      单面板拉成 2.66:1 的横条，与 D1 的观感不符。

    ★ **不引入 MHC 着色**（历史 FigD1 有）：正文全文对 "MHC" 零次提及（大小写不敏感实查 0 命中），
      且 293/8,612 = 3.40% 的 MHC 检验中无一是 FDR 显著对（get_discovery 的守卫断言），
      故沿用正文 Fig. 2a 的口径 —— 只用 Bonferroni 线，不做区域着色。
    """
    fig = plt.figure(figsize=(4.60 * scale, 3.50 * scale))
    ax_rect = [0.135, 0.145, 0.850, 0.840]
    ax = fig.add_axes(ax_rect)
    ng = [x for x in d if x["p"] >= 0.05 and x["qval"] >= 0.05]
    nm = [x for x in d if x["p"] < 0.05 and x["qval"] >= 0.05]
    sg = sorted([x for x in d if x["qval"] < 0.05], key=lambda x: x["p"])
    # ★ 全局「点的大小与风格」标准（_fig_style.py）：面积 ×2 + 同色系深 2 度描边。
    #   `scatter` 的 `s` 即"直径(pt)²" ⇒ `s_up()` 把它 ×2（等价直径 ×√2）；
    #   `lw` 用 `mew_for(等价直径)`（`ms_of_s`）统一换算。
    _s_ng = FS.s_up(2.0)
    ax.scatter([x["b"] for x in ng], [x["neglogp"] for x in ng], s=_s_ng, c=FS.C["ns"],
               lw=FS.mew_for(FS.ms_of_s(_s_ng)), edgecolors=FS.edge_color(FS.C["ns"]),
               alpha=0.45, rasterized=True, zorder=2)
    _s_nm = FS.s_up(2.4)
    ax.scatter([x["b"] for x in nm], [x["neglogp"] for x in nm], s=_s_nm,
               c=FS.C["nominal"], lw=FS.mew_for(FS.ms_of_s(_s_nm)),
               edgecolors=FS.edge_color(FS.C["nominal"]),
               alpha=0.72, rasterized=True, zorder=3)
    _s_sg = FS.s_up(15)
    ax.scatter([x["b"] for x in sg], [x["neglogp"] for x in sg], s=_s_sg, c=FS.C["sig"],
               lw=FS.mew_for(FS.ms_of_s(_s_sg)), edgecolors=FS.edge_color(FS.C["sig"]),
               marker="D", zorder=6)
    yb, yf = -math.log10(bonf), -math.log10(bh)
    ax.axhline(yb, ls="--", lw=0.6, color=FS.C["region"], zorder=1)
    ax.axhline(yf, ls=":", lw=0.6, color=FS.C["grey"], zorder=1)
    ax.axvline(0, ls="-", lw=0.5, color=FS.C["light"], zorder=1)
    bmin, bmax = min(x["b"] for x in d), max(x["b"] for x in d)
    pad = (bmax - bmin) * 0.06
    ax.set_xlim(bmin - pad, bmax + pad)
    ax.set_ylim(0, 13.0)
    ax.text(ax.get_xlim()[0] + pad * 0.3, yb + 0.18, "Bonferroni  p = 0.05/8612",
            fontsize=5.5, color=FS.C["region"], ha="left", va="bottom")
    ax.text(ax.get_xlim()[1] - pad * 0.3, yf + 0.18, "FDR = 0.05",
            fontsize=5.5, color="#7F7F7F", ha="right", va="bottom")
    # 基因标注：**程序化避让** —— 首版把 "ANXA4 (Mono_NC)" 向左排，文字横跨到
    # LINC00339 的三个红菱形上（原分辨率裁切才看得出）。改为候选位逐个探测：
    # 既不压任何红菱形，也不与已放置的标注重叠（宽度用 TextPath 实测）。
    LAB_FS = 5.4
    x0, x1 = ax.get_xlim()
    y0, y1 = ax.get_ylim()
    ax_w_in = ax_rect[2] * 4.60 * scale
    ax_h_in = ax_rect[3] * 3.50 * scale
    to_dx = lambda pt: pt / 72.0 * (x1 - x0) / ax_w_in     # noqa: E731
    to_dy = lambda pt: pt / 72.0 * (y1 - y0) / ax_h_in     # noqa: E731
    th_data = to_dy(LAB_FS * 1.18)
    # ★ 2026-10-01：红菱形按全局「点大小」标准放大后（ms 3.873 -> 5.477，半径 1.94 -> 2.74 pt），
    #   避让半径必须同步放大，否则标注会重新压到菱形上。2.6 = 旧半径 1.94 + 呼吸 0.66
    #   ⇒ 新值 = 2.74 + 0.66 = **3.40 pt**。
    mk_x, mk_y = to_dx(3.40), to_dy(3.40)                  # 菱形半径 + 呼吸
    pts = [(q["b"], q["neglogp"]) for q in sg]
    placed = []
    CANDS = [(4, 2, "left"), (-4, 2, "right"), (4, 9, "left"), (-4, 9, "right"),
             (0, -13, "center"), (0, 12, "center")]
    seen = set()
    for q in sg:
        if q["gene"] in seen:
            continue
        seen.add(q["gene"])
        txt = "%s (%s)" % (q["gene"], q["cell_type"])
        w = text_width_in(txt, LAB_FS) * (x1 - x0) / ax_w_in
        chosen = None
        for dx_pt, dy_pt, ha in CANDS:
            axx = q["b"] + to_dx(dx_pt)
            ayy = q["neglogp"] + to_dy(dy_pt)
            lo, hi = (axx, axx + w) if ha == "left" else \
                     ((axx - w, axx) if ha == "right" else (axx - w / 2, axx + w / 2))
            blo, bhi = ayy, ayy + th_data
            if lo < x0 + 0.004 or hi > x1 - 0.004:
                continue
            if any(not (hi <= a or lo >= b or bhi <= c or blo >= d) for a, b, c, d in placed):
                continue
            if any(not (hi <= mx - mk_x or lo >= mx + mk_x or bhi <= my - mk_y or blo >= my + mk_y)
                   for mx, my in pts):
                continue
            chosen = (lo, hi, blo, bhi, dx_pt, dy_pt, ha)
            break
        if chosen is None:                       # 兜底：仍放在锚点右上方
            axx = q["b"] + to_dx(4)
            ayy = q["neglogp"] + to_dy(2)
            chosen = (axx, axx + w, ayy, ayy + th_data, 4, 2, "left")
            W("   [!! S4 label no collision-free slot] %s" % txt)
        lo, hi, blo, bhi, dx_pt, dy_pt, ha = chosen
        placed.append((lo, hi, blo, bhi))
        ax.annotate(txt, (q["b"], q["neglogp"]), textcoords="offset points",
                    xytext=(dx_pt, dy_pt), ha=ha, va="bottom",
                    fontsize=LAB_FS, color=FS.C["dark"])
    ax.set_xlabel("MR effect (beta per SD of expression)", fontsize=6.8)
    ax.set_ylabel("-log10(P)", fontsize=6.8)
    ax.tick_params(labelsize=6.0, length=2.0)
    # 图例用显式 handle —— markerscale 会把红菱形放大到与数据点不成比例（首版实测）
    from matplotlib.lines import Line2D
    handles = [
        Line2D([], [], marker="o", ls="", ms=2.8, mfc=FS.C["ns"], mec="none",
               label="not significant (n = %d)" % len(ng)),
        Line2D([], [], marker="o", ls="", ms=3.2, mfc=FS.C["nominal"], mec="none",
               label="nominal P < 0.05 (n = %d)" % len(nm)),
        Line2D([], [], marker="D", ls="", ms=4.2, mfc=FS.C["sig"], mec="none",
               label="FDR < 0.05 (n = %d)" % len(sg)),
    ]
    ax.legend(handles=handles, loc="upper left", fontsize=5.5, handletextpad=0.45,
              labelspacing=0.30, borderaxespad=0.25)
    FS.finalize(ax)
    return fig


# ---------------------------------------------------------------- source data
def get_s1_qmin(s1):
    """图 S1 面板 (d) 的最小 q 值。源 = tables/17_discovery_celltype_summary.csv。

    ★ 按 cell_type **键**连接（不按位置拼接），并用 n_FDR_main 与 20_ 表的 n_FDR05 交叉核对。
    """
    r17 = load("17_discovery_celltype_summary.csv")
    assert len(r17) == 14, len(r17)
    byc = {r["cell_type"]: r for r in r17}
    assert set(byc) == set(r["cell_type"] for r in s1), \
        "S1(d): cell-type sets differ between 17_ and 20_"
    for r in s1:
        assert int(r["n_FDR05"]) == int(byc[r["cell_type"]]["n_FDR_main"]), \
            "S1(d): n_FDR05 mismatch for %s" % r["cell_type"]
    m = {c: float(r["qval_min"]) for c, r in byc.items()}
    order = sorted(m.items(), key=lambda kv: kv[1])
    W("[S1d] min q: best = %s (%.3g) ; worst = %s (%.3g)"
      % (order[0][0], order[0][1], order[-1][0], order[-1][1]))
    return m


def write_sources(s1, qmin, s2, disc, bonf, bh):
    p1 = os.path.join(BASE, "64a_FigS1_celltype_power_gating.csv")
    with io.open(p1, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("cell_type,n_instruments,n_genes,n_tested,n_nominal,n_FDR05,"
                 "pct_nominal,qval_min,neglog10_qval_min,power_flag\n")
        for r in s1:
            q = qmin[r["cell_type"]]
            fh.write("%s,%s,%s,%s,%s,%s,%s,%.10g,%.6f,%s\n"
                     % (r["cell_type"], r["n_instruments"], r["n_genes"], r["n_tested"],
                        r["n_nominal"], r["n_FDR05"], r["pct_nominal"], q,
                        -math.log10(max(q, 1e-320)), r["power_flag"]))
    p2 = os.path.join(BASE, "64b_FigS2_sensitivity_vs_main_15pairs.csv")
    with io.open(p2, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("locus_id,gene,cell_type,b_main,se_main,b_sens,se_sens,p_main,"
                 "qval_main,p_sens,p_FDR_sens,F_stat,sign_consistent,both_FDR05\n")
        for d in s2:
            fh.write("%s,%s,%s,%.10g,%.10g,%.10g,%.10g,%.6g,%.6g,%.6g,%.6g,%.6g,True,True\n"
                     % (d["locus"], d["gene"], d["cell_type"], d["b_main"], d["se_main"],
                        d["b_sens"], d["se_sens"], d["p_main"], d["qval_main"],
                        d["p_sens"], d["p_FDR_sens"], d["F_stat"]))
    # 图 S3 / S4 的阈值与每细胞类型计数（点级源数据见 10_discovery_MR_main.csv，
    # 已作为 Fig2a 的 panel source 导出到 supplementary_data，此处不重复 1 MB 大表）
    p4 = os.path.join(BASE, "64d_FigS3S4_discovery_thresholds.csv")
    with io.open(p4, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("scope,cell_type,n_tests,n_nominal_p05,n_FDR05,bonferroni_p,bh05_p,"
                 "neglog10_bonferroni,neglog10_bh05,point_source\n")
        fh.write("ALL,,%d,%d,%d,%.10g,%.10g,%.6f,%.6f,10_discovery_MR_main.csv\n"
                 % (len(disc), sum(1 for x in disc if x["p"] < 0.05),
                    sum(1 for x in disc if x["qval"] < 0.05), bonf, bh,
                    -math.log10(bonf), -math.log10(bh)))
        for c in sorted(set(x["cell_type"] for x in disc)):
            g = [x for x in disc if x["cell_type"] == c]
            fh.write("cell_type,%s,%d,%d,%d,,,,,10_discovery_MR_main.csv\n"
                     % (c, len(g), sum(1 for x in g if x["p"] < 0.05),
                        sum(1 for x in g if x["qval"] < 0.05)))
    return p1, p2, p4


def main():
    W("s67_make_supp_figs @ %s" % __import__("datetime").datetime.now().isoformat(timespec="seconds"))
    W("")
    s1 = get_s1()
    qmin = get_s1_qmin(s1)
    s2, r, slope = get_s2()
    disc, bonf, bh = get_discovery()
    W("")
    W("--- rendering ---")
    p1, w1, h1 = fit_save(lambda s: build_s1(s1, qmin, s), "FigS1_instrument_power_gating")
    p2, w2, h2 = fit_save(lambda s: build_s2(s2, r, slope, s), "FigS2_sensitivity_outcome_concordance")
    p3, w3, h3 = fit_save(lambda s: build_s3(disc, bonf, s), "FigS3_manhattan_by_celltype")
    p4, w4, h4 = fit_save(lambda s: build_s4(disc, bonf, bh, s), "FigS4_discovery_volcano",
                          target=SC_TARGET_MM, lo=SC_LO_MM, hi=SC_HI_MM)
    W("")
    W("[S2 layout] longest y-label %r = %.1f mm ; inter-panel gap = %.1f mm"
      % (S2GEO.get("longest_label"), S2GEO.get("label_w_mm", float("nan")),
         S2GEO.get("gap_mm", float("nan"))))
    W("            panel A = %.1f mm ; panel B = %.1f mm"
      % (S2GEO.get("panelA_mm", float("nan")), S2GEO.get("panelB_mm", float("nan"))))
    sp1, sp2, sp4 = write_sources(s1, qmin, s2, disc, bonf, bh)
    for sp in (sp1, sp2, sp4):
        W("[sources] %s" % os.path.basename(sp))
    with io.open(os.path.join(FS.FIGDIR, "_fig_scale_supp.json"), "w",
                 encoding="utf-8", newline="\n") as fh:
        fh.write(json.dumps({"target_mm": TARGET_MM, "window_mm": [LO_MM, HI_MM],
                             "per_figure": {k: {"target_mm": v[0],
                                                "window_mm": [v[1], v[2]]}
                                            for k, v in FIG_SPEC.items()},
                             "calibration": CAL}, ensure_ascii=False, indent=1) + "\n")
    W("")
    W("VERDICT: FigS1 %.1fx%.1f ; FigS2 %.1fx%.1f ; FigS3 %.1fx%.1f (double col %.1f mm) ; "
      "FigS4 %.1fx%.1f (single col %.1f mm)"
      % (w1, h1, w2, h2, w3, h3, TARGET_MM, w4, h4, SC_TARGET_MM))


# ★★ 2026-10-01：原来 `main()` 直接写在模块层 ⇒ **任何 `import s67_make_supp_figs` 都会
#   顺带重出 S1–S4 并覆写 `tables/64{a,b,d}_*.csv`、`figures/_fig_scale_supp.json`（探针踩过）。
#   现收进 `if __name__ == "__main__":` —— import 变为**零副作用**，
#   且允许按图分进程重出（连续多图 fit_save 在同一进程里会累积性硬崩，见 `_s163_regen_one.py`）。
#   ★ 回退：把本 guard 去掉即恢复原行为（不推荐）。
def _run_main():
    with io.open(os.path.join(LOGS, "_s67_make_supp_figs.log"), "w",
                 encoding="utf-8", newline="\n") as fh:
        try:
            main()
        except BaseException:
            import traceback
            W("FATAL:\n" + traceback.format_exc())
            fh.write("\n".join(LOG) + "\n")
            raise
        fh.write("\n".join(LOG) + "\n")


if __name__ == "__main__":
    _run_main()
