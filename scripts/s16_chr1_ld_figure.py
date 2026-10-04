# -*- coding: utf-8 -*-
"""FigD5 —— chr1p36.12 基因座 LD 结构图
Panel A: 结局区域关联图（按与结局 top 变异的 r2 着色）+ 6 个工具变量位置 + 基因注释
Panel B: 索引变异 + 6 个工具变量的 r2 矩阵（OneK1K 980 人真实基因型）
"""
import os, sys
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _fig_style as S
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrow
from matplotlib.lines import Line2D

ROOT = r"D:\endometriosis_project\11_sc_eqtl_mr_project"
T = os.path.join(ROOT, "tables")
WORK = r"D:\endometriosis_project\_chr1_ld"
INDEX_38 = 22139327

S.setup()

region = pd.read_csv(os.path.join(WORK, "chr1_plotdata.csv"))
mat = pd.read_csv(os.path.join(WORK, "chr1_matrix.csv"), index_col=0)
sig = pd.read_csv(os.path.join(T, "15_discovery_significant_with_locus.csv"))
plan = pd.read_csv(os.path.join(T, "22_replication_test_plan.csv"),
                   dtype={"variant_id_grch37": str, "variant_id_grch38": str})

# 工具变量 -> (gene, 细胞群列表, 方向)
meta = {}
for _, r in plan.iterrows():
    key = r["variant_id_grch37"]
    d = meta.setdefault(key, dict(gene=set(), cells=[], signs=set(), pos38=r["variant_id_grch38"]))
    d["gene"].add(r["gene"])
    d["cells"].append(r["cell_type_scBloodNL"])
    d["signs"].add("+" if r["discovery_b"] > 0 else "-")

# GRCh37 id -> GRCh38 pos（仅 chr1 窗口内；Tier1 的 chr10/chr2 工具不属于本图）
p37_to_38 = {r["variant_id_grch37"]: int(r["variant_id_grch38"].split(":")[1]) for _, r in plan.iterrows()}
p37_to_38["1:22465820"] = INDEX_38   # 索引变异
XLIM = (21.90, 22.50)
CHR1_VARS = [v for v in meta
             if v.startswith("1:") and XLIM[0] * 1e6 <= p37_to_38.get(v, -1) <= XLIM[1] * 1e6]
CHR1_VARS.append("1:22465820")
meta = {k: v for k, v in meta.items() if k in CHR1_VARS}
assert len(CHR1_VARS) == 7, "chr1 工具+索引应为 7 个，实得 %d" % len(CHR1_VARS)

# rsID 查询
rs_by_pos = dict(zip(region["pos_38"], region["rsid"]))

TSS = {"LINC00339": 22028206, "CDC42": 22072786}

fig = plt.figure(figsize=(7.05, 5.25))
gs = fig.add_gridspec(2, 1, height_ratios=[1.45, 1.0], hspace=0.42)
axA = fig.add_subplot(gs[0])
axB = fig.add_subplot(gs[1])


def _plabel(ax, s):
    """面板标签固定于画布左缘，避免落在 axes 外撑爆 tight bbox。

    ★ P1（2026-10-01）：字号 9 -> **8 pt**（Nature 面板字母统一 8 pt 粗体）。
    ★★ 2026-10-01（图件细改 · 问题3）：**强制大写**。
       本图不走 `_fig_style.panel_label`（面板标签要钉在画布左缘），故需单独改这一处。
       回退 = 删掉 `.upper()`。
    """
    fig.text(0.002, ax.get_position().y1 + 0.010, str(s).upper(), fontsize=8,
             fontweight="bold", va="bottom", ha="left", color="black")

# ---------------- Panel A ----------------
BINS = [(-0.001, 0.2, "#D9D9D9", "0.0-0.2"),
        (0.2, 0.4, "#4DBBD5", "0.2-0.4"),
        (0.4, 0.6, "#91D1C2", "0.4-0.6"),
        (0.6, 0.8, "#F39B7F", "0.6-0.8"),
        (0.8, 1.001, "#DC0000", "0.8-1.0")]

axA.scatter(region.loc[region["r2_index"].isna(), "pos_38"] / 1e6,
            region.loc[region["r2_index"].isna(), "neglogp"],
            s=3.2, c="#EDEDED", lw=0, zorder=1, label="not in OneK1K")
for lo, hi, col, lab in BINS:
    m = region["r2_index"].notna() & (region["r2_index"] > lo) & (region["r2_index"] <= hi)
    axA.scatter(region.loc[m, "pos_38"] / 1e6, region.loc[m, "neglogp"],
                s=5.5, c=col, lw=0, zorder=3, label="r2 %s" % lab)

idxrow = region[region["pos_38"] == INDEX_38]
if len(idxrow):
    axA.scatter(idxrow["pos_38"] / 1e6, idxrow["neglogp"], s=42, marker="D",
                c="#7E2F8E", lw=0.6, edgecolors="black", zorder=6)

# 工具变量位置（细竖线，按基因着色）
GENECOL = {"CDC42": "#E64B35", "LINC00339": "#3C5488"}
ymax = region["neglogp"].max()
instr_mb = {}
for v37, d in sorted(meta.items(), key=lambda kv: p37_to_38[kv[0]]):
    x = p37_to_38[v37] / 1e6
    gene = "/".join(sorted(d["gene"]))
    instr_mb.setdefault(gene, []).append(x)
    if v37 != "1:22465820":
        axA.axvline(x, ymin=0.0, ymax=0.62, color=GENECOL.get(gene, "#000000"),
                    lw=0.85, ls="-", zorder=4, alpha=0.85)

# 工具变量清单（右上空白区，替代易拥挤的旋转标签）
_tx, _ty = 22.195, ymax * 0.70
for gene in ["CDC42", "LINC00339"]:
    if gene not in instr_mb:
        continue
    lst = " / ".join("%.3f" % v for v in sorted(instr_mb[gene]))
    axA.text(_tx, _ty, "%s instruments (%s):  %s Mb"
             % (gene, "+" if gene == "CDC42" else "-", lst),
             color=GENECOL[gene], fontsize=5.7, va="top", ha="left", zorder=7)
    _ty -= ymax * 0.072
axA.text(_tx, _ty, "index variant = %s (outcome lead)" % rs_by_pos.get(INDEX_38, "--"),
         color="#7E2F8E", fontsize=5.7, va="top", ha="left", zorder=7)

# 基因 TSS 标记
yb = -ymax * 0.045
for g, p in TSS.items():
    axA.plot([p / 1e6], [yb], marker="s", ms=3.0, color=GENECOL[g], clip_on=False, zorder=5)
    axA.text(p / 1e6, yb - ymax * 0.022, g, fontsize=6, color=GENECOL[g],
             ha="center", va="top")

axA.axhline(-np.log10(5e-8), ls="--", lw=0.55, c=S.C["grey"], zorder=2)
axA.text(21.905, -np.log10(5e-8) + 0.25, "P = 5e-8", fontsize=6, color=S.C["grey"],
         va="bottom", ha="left", bbox=dict(facecolor="white", edgecolor="none", pad=0.6))

axA.set_xlim(21.90, 22.50)
axA.set_ylim(-ymax * 0.125, ymax * 1.13)
axA.set_ylabel("-log10 P  (female infertility, GCST90483463)")
axA.set_xlabel("Position on chr1 (Mb, GRCh38)")
axA.set_xticks(np.arange(21.9, 22.51, 0.1))
for t in axA.get_xticklabels():
    t.set_rotation(0)
axA.legend(loc="upper left", ncol=3, fontsize=5.5, columnspacing=0.7, handletextpad=0.35,
           labelspacing=0.25, bbox_to_anchor=(0.0, 1.005))
S.finalize(axA)
_plabel(axA, "a")

# ---------------- Panel B ----------------
ids = list(mat.index)
labels = []
for v in ids:
    p38 = p37_to_38.get(v, np.nan)
    rs = rs_by_pos.get(p38, "--")
    gene = "/".join(sorted(meta.get(v, {}).get("gene", {"index"})))
    if v == "1:22465820":
        gene = "INDEX"
    labels.append("%s\n%s" % (rs, gene))
labels_full = []
for v in ids:
    p38 = p37_to_38.get(v, np.nan)
    rs = rs_by_pos.get(p38, "--")
    gene = "/".join(sorted(meta.get(v, {}).get("gene", {"index"})))
    if v == "1:22465820":
        gene = "INDEX (outcome lead)"
    labels_full.append("%s / %.6f Mb / %s" % (rs, p38 / 1e6 if p38 == p38 else float("nan"), gene))

M = mat.values.astype(float)
im = axB.imshow(M, cmap="RdYlBu_r", vmin=0, vmax=1, aspect="auto")
for i in range(len(ids)):
    for j in range(len(ids)):
        v = M[i, j]
        if v == v:
            axB.text(j, i, "%.2f" % v, ha="center", va="center", fontsize=5.2,
                     color=("white" if (v > 0.62 or v < 0.18) else "black"))
axB.set_xticks(range(len(ids)))
axB.set_yticks(range(len(ids)))
axB.set_xticklabels(labels, fontsize=5.3)
axB.set_yticklabels(labels, fontsize=5.3)
axB.set_xticks(np.arange(-0.5, len(ids), 1), minor=True)
axB.set_yticks(np.arange(-0.5, len(ids), 1), minor=True)
axB.grid(which="minor", color="white", lw=0.5)
axB.tick_params(which="minor", length=0)
# 索引行列加框
ki = ids.index("1:22465820")
axB.add_patch(Rectangle((-0.5, ki - 0.5), len(ids), 1, fill=False, edgecolor="black", lw=0.8))
axB.add_patch(Rectangle((ki - 0.5, -0.5), 1, len(ids), fill=False, edgecolor="black", lw=0.8))
cb = fig.colorbar(im, ax=axB, fraction=0.030, pad=0.015, aspect=13)
cb.set_label("r2 (OneK1K 980 donors, GRCh37 genotypes)", fontsize=6)
cb.ax.tick_params(labelsize=6, width=0.5, length=2)
cb.outline.set_linewidth(0.4)
_plabel(axB, "b")

# ---- 诊断：定位撑爆 tight bbox 的 artist ----
fig.canvas.draw()
_r = fig.canvas.get_renderer()
_tb = fig.get_tightbbox(_r)
print("   [dbg] fig inches = %s ; tight bbox inches = %.2f x %.2f"
      % (fig.get_size_inches(), _tb.width, _tb.height))
from matplotlib.text import Text as _T
for _t in fig.findobj(match=_T):
    try:
        _bx = _t.get_window_extent(_r)
    except Exception:
        continue
    if _bx.width > 1500 or abs(_bx.x0) > 20000 or abs(_bx.x1) > 30000:
        print("   [dbg WIDE] %r x0=%.0f x1=%.0f w=%.0f" % (_t.get_text()[:34], _bx.x0, _bx.x1, _bx.width))
for _ax in fig.axes:
    _b = _ax.get_window_extent(_r)
    print("   [dbg AX] pos=%s bbox: x0=%.0f x1=%.0f w=%.0f" % (_ax.get_label(), _b.x0, _b.x1, _b.width))
    for _t in _ax.get_xticklabels() + _ax.get_yticklabels():
        _bx = _t.get_window_extent(_r)
        if _bx.width > 800 or abs(_bx.x1) > 30000:
            print("      [dbg TICK] %r x0=%.0f x1=%.0f w=%.0f" % (_t.get_text()[:34], _bx.x0, _bx.x1, _bx.width))

paths, wmm, hmm = S.save(fig, "FigD5_chr1_LD_structure")

# 同步输出带标签的矩阵表
out = mat.copy()
out.index = labels_full
out.columns = [l.split(" / ")[0] for l in labels_full]
out.to_csv(os.path.join(T, "23_chr1_LD_matrix.csv"), encoding="utf-8-sig")
print("   width=%.1f mm height=%.1f mm" % (wmm, hmm))
