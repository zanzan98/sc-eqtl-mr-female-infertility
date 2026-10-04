#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""_fig_style.py —— 本子项目（11_sc_eqtl_mr_project）绘图样式模块

规格（沿用主项目 10_pain_infertility_project/scripts/_fig_style.py 的既定纪律）：
  Arial 8-10 pt ; 刻度朝内 ; 轴线 0.5 pt ; 600 dpi ; 全部文字 ASCII

★ 字形纪律（通用，血泪教训）：Arial 不含 U+2205(∅) / U+2227(∧) / U+2605(★) / U+221D(∝)
  等字符（含 □■→×≤≥▲●ρμ±≈）。matplotlib 缺字形时**不报错**，只发 warning 并在
  LastResortHE 下渲染成**空豆腐块**。故 save() 做双重核验：
    ① 捕获 "Glyph NNNN missing from font(s)" 告警（权威判据）
    ② 用 rcParams['font.sans-serif'][0] 的真实 cmap 逐字符比对

============================================================================
调色板机制（2026-09-30 指令 S：把 GB12 写进项目 PALETTES）
----------------------------------------------------------------------------
  · `FIGPAL` 环境变量选调色板；**默认 "gb12"**（用户 2026-09-30 给定的 12 色，
    已固化进 skill `nature-figure` 作默认盘）。`FIGPAL=npg` 可逐字节回到历史配色。
  · `FIG_OUTDIR` 可把输出重定向到别处（试色 / 回归比对时不污染 figures/）。
  · `C` 的构造：`C = dict(_C_BASE) ; C.update(PALETTES[PALETTE]["C"])`
      - `_C_BASE` = **历史语义色**（勿改！决定 npg 默认出图的字节一致性）
      - `PALETTES["npg"]["C"] = {}` ⇒ npg 下 `C is _C_BASE` 的等值副本，逐字节不变
      - `PALETTES["gb12"]["C"]` 来自同目录 `_gb12_palette.py`（GB12 映射的唯一真源）
  · `cmap_for(kind)`：热图色图。**npg 必须原样返回内置色图名字符串**，否则默认出图
    不再逐字节可复现；gb12 时由 `C["mhc"]`（seqA）/ `C["sig"]`（seqB）与黑白插值派生。
  · 任何作图脚本只用 `NPG` / `C` 两个名字取色 ⇒ 换调色板无需改动作图逻辑。

★ 回归纪律（本轮实测）：改完本文件 + 把 s37/s67 里残留的硬编码语义 hex 换成 `FS.C[...]`
  之后，`FIGPAL=npg` 必须仍能**逐字节**复现 figures/ 现行 9 件（Fig1-5 / FigS1-4，9/9 IDENTICAL）。
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager

ROOT = "D:/endometriosis_project/11_sc_eqtl_mr_project"
# ★ FIG_OUTDIR 可把输出重定向到别处（试色 / 回归时不污染 figures/）
FIGDIR = os.environ.get("FIG_OUTDIR") or os.path.join(ROOT, "figures")
os.makedirs(FIGDIR, exist_ok=True)

_NPG_BASE = ["#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F",
             "#8491B4", "#91D1C2", "#DC0000", "#7E6148", "#B09C85"]

# npg 默认语义色 = 历史值（勿改！决定默认出图的字节一致性）
_C_BASE = {
    "sig":      "#DC0000",   # FDR 显著
    "ns":       "#B8B8B8",   # 不显著
    "mhc":      "#3C5488",   # MHC 区
    "nominal":  "#F39B7F",   # 仅名义显著
    "locus1":   "#E64B35",   # chr1 WNT4/CDC42
    "locus2":   "#00A087",   # chr10 YME1L1
    "locus3":   "#3C5488",   # chr2 ANXA4
    "grey":     "#7F7F7F",
    "dark":     "#2B2B2B",
    "light":    "#D9D9D9",
    # P0 审稿意见修正：coloc strong_shared 二级命名（与正文 §3.4 / §D2 完全一致）
    "robust_strong":  "#8B0000",   # 稳健 strong（45 组敏感性 9/9 均 PP.H4 > 0.8）
    "suggest_strong": "#E87070",   # 提示性 strong（仅 p12 = 1e-6 时跌破 0.8）
    "moderate_sh":    "#F7C4A0",   # moderate_shared
    "no_shared":      "#BFBFBF",
    # P0：区域归因 / 不可报告
    "region":  "#B22222",
    "nr":      "#B22222",
    # ↓↓ 2026-09-30 收编键：把 s37/s67 里原本硬编码的语义色并入真源
    #    （默认值 = 原硬编码值，故 npg 默认出图不受影响）
    "risk_up":     "#C0392B",   # PheWAS 风险升高
    "risk_down":   "#2E86C1",   # PheWAS 风险降低
    "gene_ym":     "#7E6148",   # YME1L1（+ ANXA4_Mono_NC 的 unstable 组）
    "gene_an":     "#B09C85",   # ANXA4
    # ★★ 2026-10-02（指令 S34）：Fig5(b) 谱系组 E / F 原为**硬编码**在 s37 里的
    #    "#8C8C8C" / "#C7C7C7"，本轮收编为语义键。这里刻意保留**原硬编码灰值**
    #    ⇒ FIGPAL=npg 出图与历史逐字节一致；GB12 侧的真值见 `_gb12_palette.GB12_C`
    #    （"lineage_e": "#9ED17B" / "lineage_f": "#BAD2E1"）。
    "lineage_e":   "#8C8C8C",
    "lineage_f":   "#C7C7C7",
    "cyan":        "#4DBBD5",   # 替代 lead / 二线强调
    "green_dark":  "#1A7F37",   # Fig5 五级色序第 1
    "amber":       "#F0AD4E",   # Fig5 五级色序第 3
    "lineage_b":   "#E8834A",   # 谱系组 B（胎盘附着/产科出血）
    "lineage_c":   "#7B52A1",   # 谱系组 C（早产）
    "tint_cool":   "#EAF1F8",   # 冷色浅底（说明框）
    "tint_warm":   "#FDEEEA",   # 暖色浅底（说明框）
}


def _load_gb12():
    """按**路径**载入同目录的 `_gb12_palette.py`（避免 sys.path 依赖）。"""
    import importlib.util as _ilu
    p = os.path.join(os.path.dirname(os.path.abspath(__file__)), "_gb12_palette.py")
    spec = _ilu.spec_from_file_location("_gb12_palette", p)
    m = _ilu.module_from_spec(spec)
    spec.loader.exec_module(m)
    return list(m.GB12), dict(m.GB12_C)


try:
    _GB12, _GB12_C = _load_gb12()
except Exception as _e:                                       # pragma: no cover
    import sys as _sys
    _sys.stderr.write("WARN gb12 palette not loaded: %r\n" % (_e,))
    _GB12, _GB12_C = [], {}

PALETTES = {"npg": {"NPG": _NPG_BASE, "C": {}},
            "gb12": {"NPG": _GB12, "C": _GB12_C}}

PALETTE = os.environ.get("FIGPAL", "gb12")
if PALETTE not in PALETTES:
    raise SystemExit("unknown FIGPAL=%r ; available=%s" % (PALETTE, sorted(PALETTES)))

NPG = list(PALETTES[PALETTE]["NPG"])
C = dict(_C_BASE)
C.update(PALETTES[PALETTE]["C"])


def _mix(hex_a, hex_b, t):
    a = [int(hex_a[i:i + 2], 16) for i in (1, 3, 5)]
    b = [int(hex_b[i:i + 2], 16) for i in (1, 3, 5)]
    return "#%02X%02X%02X" % tuple(int(round(a[i] + (b[i] - a[i]) * t)) for i in range(3))


def cmap_for(kind):
    """热图色图：kind='seqA'（原 YlGnBu，PP.H4 热图）/ 'seqB'（原 OrRd，chr1p36 热图）。

    ★ 默认调色板 npg 必须**原样返回内置色图名字符串**，否则默认出图不再逐字节可复现。
    """
    if PALETTE == "npg":
        return "YlGnBu" if kind == "seqA" else "OrRd"
    from matplotlib.colors import LinearSegmentedColormap
    base = C["mhc"] if kind == "seqA" else C["sig"]
    lo = _mix(base, "#FFFFFF", 0.90)
    mid = _mix(base, "#FFFFFF", 0.45)
    return LinearSegmentedColormap.from_list(
        "pal_%s_%s" % (PALETTE, kind), [lo, mid, base, _mix(base, "#000000", 0.35)])


def setup():
    have = {f.name for f in font_manager.fontManager.ttflist}
    fam = "Arial" if "Arial" in have else ("Liberation Sans" if "Liberation Sans" in have else "DejaVu Sans")
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": [fam, "Arial", "Liberation Sans", "DejaVu Sans"],
        "font.size": 8, "axes.labelsize": 8, "axes.titlesize": 8,
        "xtick.labelsize": 7, "ytick.labelsize": 7,
        "legend.fontsize": 6.5, "figure.titlesize": 9.5,
        "mathtext.fontset": "dejavusans",
        "axes.linewidth": 0.5, "lines.linewidth": 0.7, "lines.markersize": 2.4,
        "patch.linewidth": 0.5,
        "xtick.major.width": 0.5, "ytick.major.width": 0.5,
        "xtick.major.size": 2.2, "ytick.major.size": 2.2,
        "xtick.direction": "in", "ytick.direction": "in",
        "xtick.top": False, "ytick.right": False,
        "legend.frameon": False, "legend.handlelength": 1.2,
        "legend.handletextpad": 0.5, "legend.labelspacing": 0.3,
        "figure.dpi": 150, "savefig.dpi": 600, "savefig.pad_inches": 0.02,
        "ps.fonttype": 42, "pdf.fonttype": 42, "svg.fonttype": "none",
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.grid": False, "errorbar.capsize": 2,
        # ★ 图文字一律 ASCII：matplotlib 默认用 U+2212(MINUS SIGN) 排负号刻度，
        #   虽然 Arial 有该字形不会触发豆腐块，但违反「全 ASCII」纪律（且换字体会翻车）
        #   -> 关掉后刻度负号走 ASCII '-'。
        "axes.unicode_minus": False,
    })
    return fam


def check_glyphs(fig):
    from matplotlib.font_manager import FontProperties, findfont
    from matplotlib.text import Text
    fam = plt.rcParams["font.sans-serif"][0]
    missing = set()
    try:
        from fontTools.ttLib import TTFont
    except Exception:
        return missing
    try:
        path = findfont(FontProperties(family=fam))
        cm = set(TTFont(path, fontNumber=0, lazy=True).getBestCmap().keys())
    except Exception:
        return missing
    for t in fig.findobj(match=Text):
        s = t.get_text()
        if not s:
            continue
        for ch in s:
            if ord(ch) < 128:
                continue
            if ord(ch) not in cm:
                missing.add("%s(U+%04X) in %r" % (ch, ord(ch), s[:28]))
    return missing


def save(fig, name, dpi=600, exts=("pdf", "png")):
    """保存并做双重字形核验；返回 (paths, 实测宽度mm, 实测高度mm)。"""
    import warnings
    cm_bad = check_glyphs(fig)
    with warnings.catch_warnings(record=True) as wl:
        warnings.simplefilter("always")
        paths = []
        for ext in exts:
            p = os.path.join(FIGDIR, "%s.%s" % (name, ext))
            fig.savefig(p, dpi=dpi, bbox_inches="tight", pad_inches=0.02, facecolor="white")
            paths.append(p)
    warn_bad = set()
    for w in wl:
        m = str(w.message)
        if "missing from font" in m:
            warn_bad.add(m.split("UserWarning: ")[-1][:110])
    # 实测尺寸（bbox_inches='tight' 会改变画布，必须实测而非读 rcParams）
    try:
        from PIL import Image
        im = Image.open(paths[-1])
        wmm, hmm = im.size[0] / dpi * 25.4, im.size[1] / dpi * 25.4
    except Exception:
        wmm = hmm = float("nan")
    plt.close(fig)
    print("   [saved] %s  (%s @ %d dpi)  实测 %.1f x %.1f mm"
          % (name, "/".join(exts), dpi, wmm, hmm))
    if warn_bad:
        print("   [!! GLYPH MISSING -> LastResort tofu] %s" % "; ".join(sorted(warn_bad)))
    if cm_bad:
        print("   [!! CMAP MISSING in %s] %s"
              % (plt.rcParams["font.sans-serif"][0], "; ".join(sorted(cm_bad))))
    return paths, wmm, hmm


# ============================================================================
# ★★ 「同色系加深 k 度」（2026-10-01 用户定案：k = 2，1 度 = CIE L* 降 10）
#   口径：**保色相 H、保 HSL 饱和度 S，只把 CIE L* 降 k×10**。
#   ★ 与 skill `cvd-safe-palette-design/scripts/cvd_sim.py` 的 `hsl_hex_L(h, L*, S)`
#     同口径；本函数**自包含实现**（不依赖外部 skill 路径），
#     `scripts/_s140_darken_crosscheck.py` 做了逐色对照（<1/255）。
#   ★ 为什么不用「每个 RGB 通道乘系数」：那样会同时改色相与饱和度，
#     浅色（如 #FDDAEC）乘完会明显偏色；锁 L* 是唯一能保证"只是更深"的做法。
# ============================================================================
def _hex2rgb01(h):
    h = h.lstrip("#")
    return tuple(int(h[i:i + 2], 16) / 255.0 for i in (0, 2, 4))


def _rgb01_to_hex(rgb):
    return "#" + "".join("%02X" % int(round(max(0.0, min(1.0, c)) * 255)) for c in rgb)


def _rgb01_to_lab(rgb):
    """sRGB(0-1) -> CIE L*a*b*（D65，2°）。"""
    def f(u):
        return u / 12.92 if u <= 0.04045 else ((u + 0.055) / 1.055) ** 2.4
    r, g, b = (f(c) for c in rgb)
    X = 0.4124564 * r + 0.3575761 * g + 0.1804375 * b
    Y = 0.2126729 * r + 0.7151522 * g + 0.0721750 * b
    Z = 0.0193339 * r + 0.1191920 * g + 0.9503041 * b
    Xn, Yn, Zn = 0.95047, 1.0, 1.08883

    def g_(t):
        return t ** (1.0 / 3.0) if t > 0.008856 else (7.787 * t + 16.0 / 116.0)
    fx, fy, fz = g_(X / Xn), g_(Y / Yn), g_(Z / Zn)
    return (116.0 * fy - 16.0, 500.0 * (fx - fy), 200.0 * (fy - fz))


def lstar(hexv):
    """CIE L*（0-100）。"""
    return _rgb01_to_lab(_hex2rgb01(hexv))[0]


def darken(hexv, dl=20.0):
    """同色系加深：**保色相 H、保 HSL 饱和度 S**，把 CIE L* 降 `dl`。

    ★ 返回 hex。`dl = 20` 即用户的「加深 2 个度」（1 度 = L* 降 10）。
    ★ 已到 L*=0（纯黑）时返回 #000000。
    """
    import colorsys
    r, g, b = _hex2rgb01(hexv)
    h, l, s = colorsys.rgb_to_hls(r, g, b)
    tgt = _rgb01_to_lab((r, g, b))[0] - float(dl)
    if tgt <= 0.0:
        return "#000000"
    lo, hi = 0.0, l
    for _ in range(40):                     # 二分：HSL 明度单调映射到 L*
        mid = 0.5 * (lo + hi)
        if _rgb01_to_lab(colorsys.hls_to_rgb(h, mid, s))[0] > tgt:
            hi = mid
        else:
            lo = mid
    return _rgb01_to_hex(colorsys.hls_to_rgb(h, 0.5 * (lo + hi), s))


def panel_label(ax, label, dx=-0.11, dy=1.03, size=8):
    """面板字母。Nature 规定：**8 pt、粗体、直立（非斜体）、小写 a,b,c**。

    ★ P1（2026-10-01）：默认 size 由 9 → **8**（Nature 硬指标）。
      凡未显式传 size 的调用点（s37 的 Fig2-5）自动随此改动。
    ★★ 2026-10-01（图件细改 · 问题3）：**强制大写**。
      用户明确要求 Fig1-Fig5 与 FigS1-FigS4 的面板字母一律改为大写 A,B,C,...
      —— 注意这与 Nature 官方规范（小写 a,b,c）**相反**，是用户对目标期刊体例的明确选择。
      ★ 回退方式：删掉下面这一行 `label = str(label).upper()` 即整体恢复小写，
        **无需改动任何调用点**（20+ 处调用传的都是小写单字母）。
    """
    label = str(label).upper()
    ax.text(dx, dy, label, transform=ax.transAxes, fontsize=size,
            fontweight="bold", va="bottom", ha="left", color="black")


# ============================================================================
# ★★★ 全局标准（2026-10-01 用户定案）：「面板字母 ↔ 该面板图形」的距离
# ----------------------------------------------------------------------------
#   用户原话：「记住现在字母与图的位置，**后面所有的图都改成这个距离**」。
#   口径（Fig2 第二~六轮逐轮收敛出来的，别再另立）：
#       Δx = 字母**右缘** → 该面板 **Y 轴最上面那个词**的**左缘**
#       Δy = 字母**下缘** → 该词的**上缘**
#   数值（Fig2 第六轮定案并已复核到 ±0.01 mm）：
LETTER_GAP_DX_MM = 6.41
LETTER_GAP_DY_MM = 2.42
#   ★★ 两个必须遵守的细节（都踩过坑）：
#     ① 「Y 轴最上面那个词」**必须过滤到绘图区内**。`ax.get_yticklabels()` 会返回
#        落在 ylim 之外、**根本不会被绘制**的标签 —— Fig2 面板 A 的 '14' 与 '-2' 即是；
#        不剔除就会选中看不见的 '14'，把 A 的竖直间隙算成 −4.69 mm（字母跑到词"下方"）。
#     ② 该词默认**右对齐**到绘图区左缘，各面板"最宽那行"不同 ⇒ **dx 必然逐面板不同**，
#        不要试图用一个 dx 套所有面板。
#   ★ dx / dy 是 **axes 分数**：只动「行距 / 幅高 / 画布比例」而**不动面板尺寸**时，
#     已定值继续精确有效；一旦面板尺寸（axes 的 mm 宽高）变了，必须用下面的函数重算。
#   ★★ **scale 依赖（务必注意）**：文字尺寸以 pt 固定、不随图幅缩放，而 axes 的 mm 尺寸随
#      `fit_save` 的标定 scale 等比缩放 ⇒ **同一个 dx 在不同 scale 下对应不同的 mm**。
#      实测（Fig2）：scale 1.00000 → dx = −0.46478；scale 1.02351 → dx = −0.45428（差 0.53 mm）。
#      ⇒ ① 函数必须在**最终 scale** 下调用；
#         ② 新图第一轮建议**跑两遍**（第一遍 fit_save 收敛并把 scale 落到 `_fig_scale.json`，
#            第二遍 import 期读回 CAL 后才是最终 scale），或直接显式传入标定后的 scale。
#         ③ ★ 另一个陷阱：`FIG_OUTDIR` 一旦被设置，`FIGDIR` 就指向别处，CAL 会**从那个目录**
#            读 `_fig_scale.json` —— 读不到就静默退化成 1.0。做探针/自检时要直接读
#            `figures/_fig_scale.json`，不要用 `CAL.get(name, 1.0)`。
# ============================================================================
def letter_offset_for(ax, gap_dx_mm=None, gap_dy_mm=None, size=8, label="X",
                      anchor="word", shift_mm=None):
    """求 `panel_label` 应传的 (dx, dy)，使字母与「本面板 Y 轴最上面那个词」的位移 = 全局标准。

    anchor : ``"word"``（默认，Fig2 口径）或 ``"corner"``。
      ★★★ 2026-10-01 补充图统一时新增 `"corner"`，用途与判据如下 ——
      对**补充图 S1/S2/S3** 这类面板，"Y 轴最上面那个词"往往是**长分类标签**
      （如 S2 面板 b 的 ``"CDC42 | Mono_NC"`` ≈ 19 mm）或**长行名**（S1 面板 a），
      此时按 ``"word"`` 锚会把字母推到绘图区外 12–20 mm、**压到相邻面板的数据/标题上**，
      已不再是"字母与本面板的距离"。
      ⇒ 这类面板改用 **``"corner"``**：锚 = **本面板绘图区左上角**，
         口径不变，仍为 Δx = 6.41 mm（字母右缘 → 本面板绘图区左缘）、
         Δy = 2.42 mm（字母下缘 → 本面板绘图区上缘）。
      ★ 只有**数值刻度**（Fig2–5 与 S3 左列）才用 ``"word"``：那才是 Fig2 定义口径时的场景。

    返回 (dx, dy)，单位 = **axes 分数**（与 `panel_label` 同口径，dx/dy 即字母左下角）。

    ★ 必须在**该面板全部文字都布置完成之后**调用（即原本调用 `panel_label` 的那一行），
      因为它要真实渲染才能量到文字的 bbox。

    ★★ 必须在**最终 scale** 下调用（见上方「scale 依赖」条）：文字尺寸不随图幅缩放，
      而 axes 的 mm 尺寸随 scale 等比缩放 ⇒ 同一 dx 在不同 scale 下对应的 mm 不同。

    用法（替代 `panel_label`）::

        dx, dy = FS.letter_offset_for(ax, label="c")
        FS.panel_label(ax, "c", dx=dx, dy=dy)
    """
    import numpy as _np

    if gap_dx_mm is None:
        gap_dx_mm = LETTER_GAP_DX_MM
    if gap_dy_mm is None:
        gap_dy_mm = LETTER_GAP_DY_MM

    fig = ax.figure
    fig.canvas.draw()                                   # 需要真实 renderer 才能量 bbox
    r = fig.canvas.get_renderer()
    dpi = float(fig.dpi)
    px_per_mm = dpi / 25.4

    pw, ph = (float(v) for v in _np.asarray(fig.get_size_inches()).ravel()[:2])
    p = ax.get_position()
    axx0 = p.x0 * pw * dpi
    axy0 = p.y0 * ph * dpi
    axy1 = p.y1 * ph * dpi
    axw = (p.x1 - p.x0) * pw * dpi
    axh = (p.y1 - p.y0) * ph * dpi

    # ① 本面板"最上面那个词" —— 仅取**完整落在绘图区内**的 y 刻度标签（见上文坑①）
    best = None
    _best_txt = "<none>"
    _n_used = _n_all = 0
    _force_corner = str(anchor).lower() in ("corner", "axes", "ax_corner")
    if not _force_corner:
        for t in ax.get_yticklabels():
            if not t.get_text().strip():
                continue
            _n_all += 1
            b = t.get_window_extent(r)
            if b.y0 < axy0 - 1e-6 or b.y1 > axy1 + 1e-6:
                continue
            _n_used += 1
            if best is None or b.y1 > best.y1:
                best = b
                _best_txt = t.get_text()
    if best is None:
        # ★★ 回退（2026-10-01 全局统一时新增）：该面板**没有**落在绘图区内的 y 刻度标签
        #    （典型：Fig4 面板 (b) 的轨道图 `set_yticks([])`）。此时"Y 轴最上面那个词"
        #    根本不存在，口径无定义 ⇒ 改以**绘图区左上角**为锚：
        #        字母右缘 → 绘图区左缘 = gap_dx_mm
        #        字母下缘 → 绘图区上缘 = gap_dy_mm
        #    ★ 视觉上与"距词"口径基本一致：有词时词右对齐到绘图区左缘（差一个 tick pad
        #      ~1.2 mm），故两种口径给出的字母落点相差不超过约 1 mm。
        class _Anchor(object):
            __slots__ = ("x0", "y0", "x1", "y1")

        best = _Anchor()
        best.x0 = axx0
        best.x1 = axx0
        best.y0 = axy1
        best.y1 = axy1

    # ② 字母自身 bbox —— 用临时 text 量，量完即删（字形/字号与 panel_label 完全一致）
    tmp = ax.text(0.0, 0.0, str(label).upper(), transform=ax.transAxes, fontsize=size,
                  fontweight="bold", va="bottom", ha="left")
    fig.canvas.draw()
    lb = tmp.get_window_extent(fig.canvas.get_renderer())
    tmp.remove()

    tgt_lx0 = best.x0 - gap_dx_mm * px_per_mm - lb.width     # 目标字母左缘
    tgt_ly0 = best.y1 + gap_dy_mm * px_per_mm                # 目标字母下缘
    # ★★ 2026-10-01（逐图细改 · Fig3 用户微调）：在**标准锚点之上**再叠加一个毫米级位移。
    #    用户原话：「D 再往上移一点，不要跟 Y 轴标题连在一起」「B 往上移，跟 C 在同一水平线」
    #    「C 往右移一点点」—— 这类"就一点点"的手工微调不能改全局标准，只能逐字母叠加。
    #    shift_mm = (Δx_mm, Δy_mm)；**+Δy = 向上**（与 matplotlib 同向）。
    if shift_mm:
        tgt_lx0 += float(shift_mm[0]) * px_per_mm
        tgt_ly0 += float(shift_mm[1]) * px_per_mm
    if os.environ.get("LETTER_ANCHOR_DIAG"):
        ANCHOR_DIAG.append(dict(label=str(label), anchor=str(anchor),
                                axx0_mm=round(axx0 / px_per_mm, 3),
                                word=_best_txt, n_used=_n_used, n_all=_n_all,
                                word_x0_mm=round(best.x0 / px_per_mm, 3),
                                probe_x_mm=round(tgt_lx0 / px_per_mm, 3),
                                lw_mm=round(lb.width / px_per_mm, 3)))
    return ((tgt_lx0 - axx0) / axw, (tgt_ly0 - axy0) / axh)


LETTER_NUDGE_LOG = {}
# ★ 列对齐记录：{(字母1,字母2,...): 实际采用的公共左缘 mm}；None = 无可用候选
LETTER_COLUMN_LOG = {}
# ★ 列对齐诊断：[{候选左缘 mm, 被否原因 or None}, ...]（每个候选一条；每次调用清空追加）
_COL_DIAG = []
# ★ 列对齐诊断：各面板的标准探针左缘（mm）与其 axes 的 [x0, w]（mm）
_COL_PROBE = []
# ★ 锚点诊断（仅当环境变量 LETTER_ANCHOR_DIAG 置位时写入）
ANCHOR_DIAG = []


def _nudge_letter_clear(ax, label, dx, dy, size=8, thresh_mm2=0.30):
    """把刚放好的面板字母**微移**到"不与任何其他文字相碰"的第一个位置。

    ★★ 为什么必须做（2026-10-01 Fig5 实测，本轮回归）：
      标准锚定式是「字母右缘 = 本面板 Y 轴最上面那个词左缘 − 6.41 mm」。对**有旋转 y 轴标题**
      的面板，那块 6.41 mm 的留白区**正好被 y 轴标题占着**（Fig5 面板 (e)：旋转标题
      `number of tests (of 968)` 占 x 111.8–118.0 mm，锚词左缘在 119.7 mm ⇒ 字母落到
      111.4–113.3 mm，与标题**压在一起 1.51 × 2.60 mm**）。用户对字母的硬要求是
      「**不与图的字重叠**，都在图的左上角但距离适中」⇒ 锚定之后必须**做一次碰撞消解**。

    做法：在 (dx, dy) 附近按「上 → 右 → 下 → 左 → 斜」的顺序穷举小位移（0.5–6 mm），
      取**第一个**与任何其他可见 `Text`（含本面板刻度标签/轴标题、其他面板的一切文字）
      重叠面积 ≤ `thresh_mm2`、且**仍在画布内**的位置；找不到就原样返回并记入
      `LETTER_NUDGE_LOG`（构建日志会打印）。
    ★ 阈值 0.30 mm²：文字的逻辑包围盒之间常有无害的 0.05–0.08 mm 细缝重叠
      （多行标签的兄弟行、`Text` 盒的 leading），按面积卡掉它们、只处理真碰。
    """
    from matplotlib.text import Text          # ★ 局部导入：本模块顶层没有 Text 名（曾漏而报 NameError）

    fig = ax.figure
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    fd = float(fig.dpi)
    px_per_mm = fd / 25.4
    pw, ph = (float(v) for v in fig.get_size_inches())

    mine = None
    for t in ax.texts:
        if t.get_text() == str(label).upper() and abs(float(t.get_fontsize()) - size) < 1e-6:
            mine = t
    if mine is None:
        return dx, dy, 0.0
    others = []
    for t in fig.findobj(match=Text):
        if t is mine or not t.get_visible():
            continue
        if not t.get_text().strip():
            continue
        try:
            b = t.get_window_extent(r)
        except Exception:
            continue
        if b.width > 0 and b.height > 0:
            others.append(b)

    p = ax.get_position()
    axw_mm = (p.x1 - p.x0) * pw * 25.4
    axh_mm = (p.y1 - p.y0) * ph * 25.4

    def _ok(ox_mm, oy_mm):
        mine.set_position((dx + ox_mm / axw_mm, dy + oy_mm / axh_mm))
        b = mine.get_window_extent(r)
        if b.x0 < 0 or b.y0 < 0 or b.x1 > pw * fd or b.y1 > ph * fd:
            return False                       # 出画布 ⇒ 会撑大 tight bbox
        for o in others:
            ox = min(b.x1, o.x1) - max(b.x0, o.x0)
            oy = min(b.y1, o.y1) - max(b.y0, o.y0)
            if ox > 0 and oy > 0 and (ox * oy) / (px_per_mm ** 2) > thresh_mm2:
                return False
        return True

    if _ok(0.0, 0.0):
        mine.set_position((dx, dy))
        return dx, dy, 0.0
    cands = [(0.0, 1.0), (0.0, 2.0), (0.0, 2.5), (0.0, 3.5), (0.0, 4.5),
             (1.0, 0.0), (2.0, 0.0), (3.0, 0.0), (4.0, 0.0), (5.0, 0.0), (6.0, 0.0),
             (0.0, -1.0), (0.0, -2.0), (0.0, -3.0), (0.0, -4.5),
             (-1.5, 0.0), (-3.0, 0.0),
             (1.5, 1.5), (3.0, 1.5), (1.5, 3.0)]
    for ox, oy in cands:
        if _ok(ox, oy):
            LETTER_NUDGE_LOG[str(label).upper()] = (ox, oy)
            return dx + ox / axw_mm, dy + oy / axh_mm, (ox ** 2 + oy ** 2) ** 0.5
    mine.set_position((dx, dy))
    LETTER_NUDGE_LOG[str(label).upper()] = None
    return dx, dy, 0.0


def align_letter_column(plants, size=8, anchor="word", gap_dx_mm=None, gap_dy_mm=None,
                        thresh_mm2=0.30, edge_mm=0.90, shifts_mm=None, verbose=False):
    """把**同一列**的若干面板字母对齐到同一条竖线（用户 2026-10-01：「D 要跟 AB 对齐」）。

    为什么需要它：
      `letter_offset_for(anchor="word")` 的口径是「字母右缘 = **本面板** Y 轴最上面那个词的左缘
      − 6.41 mm」。各面板"最上面那个词"**宽度不同** ⇒ **各面板字母的 x 天生不同**。
      Fig2 就是这样（A x0=5.68 / B x0=13.60），用户接受；但 Fig3 的 A、B、D 是**同一列**
      （A、D 通栏，B 左列），横向不齐就非常显眼（v10 里 D 在 0.72 mm、A/B 在 18.0 mm）。

    规则（三条，按序）：
      ① 先对每个面板按**标准口径**各算一个锚点 x（保留各自的标准 dy —— 竖直方向仍逐面板）；
      ② 候选公共 x = 各锚点去重后**从大到小**（= 尽量贴近标签列、远离画布左缘）。
      ③ 取**第一个**满足「三处字母都落在画布内(距边 >= edge_mm) 且与任何其他可见 Text
         重叠面积 <= thresh_mm2」的候选。全不满足则退化为各自标准锚点并记录。
      ★ 用"最大"而非"最小"：任何一个面板的最上面那个词只要够宽，就能把整列字母拖到画布边缘
        （Fig3 面板 B 的 `CDC42_B_MEM  (robust)` 会把字母拖到 x0 = 0.51 mm）。

    shifts_mm : ``{"B": (0.0, 2.9), ...}`` 在共线结果之上再叠加的**逐字母手工微调**（mm，
       +Δy = 向上）；只动被点名的那个字母（Fig3 用户要求 B 与 C 齐平、D 抬高离开 y 轴标题）。
       收集进碰撞核验（改后的位置同样不许压字）。

    返回：{label: (dx, dy, x_px)}；并把实际用的公共 x 记入 `LETTER_COLUMN_LOG`。

    ★★ 与 `_nudge_letter_clear` 互斥：对齐之后**不能再逐字母微移**（会把对齐破坏掉）。
    ★★ 必须在**该行所有文字布置完成之后**、且**最终 scale** 下调用。

    用法::

        FS.align_letter_column([(axa, "a"), (axb, "b"), (axd, "d")])
    """
    from matplotlib.text import Text

    if not plants:
        return {}
    del _COL_DIAG[:]                 # 只保留「最后一次构建」的诊断（fit_save 会构建十余次）
    fig = plants[0][0].figure
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    fd = float(fig.dpi)
    px_per_mm = fd / 25.4
    pw, ph = (float(v) for v in fig.get_size_inches())
    P = [tuple(ax.get_position().bounds) for ax, _ in plants]

    def _dx_for(ax_pos, x_px):
        return (x_px - ax_pos[0] * pw * fd) / (ax_pos[2] * pw * fd)

    # ① 各面板的标准锚点（dx, dy）
    offs = [letter_offset_for(ax, gap_dx_mm=gap_dx_mm, gap_dy_mm=gap_dy_mm,
                              size=size, label=lb, anchor=anchor)
            for ax, lb in plants]
    probe_x = [pp[0] * pw * fd + offs[i][0] * pp[2] * pw * fd for i, pp in enumerate(P)]

    # ② 其它文字（碰撞核验用；位置不随字母移动而变，量一次即可）
    ours = {str(lb).upper() for _, lb in plants}
    others = []
    for t in fig.findobj(match=Text):
        if not t.get_visible() or not t.get_text().strip():
            continue
        tb = t.get_window_extent(r)
        if tb.width > 0 and tb.height > 0:
            others.append(tb)

    # ③ 字母 bbox 的真实尺寸：临时各建一个同字号/同字重的字母量出来（量完即删）
    tmp = []
    for ax, lb in plants:
        tmp.append(ax.text(0.0, 0.0, str(lb).upper(), transform=ax.transAxes,
                           fontsize=size, fontweight="bold", va="bottom", ha="left"))
    fig.canvas.draw()
    _tbb = [t.get_window_extent(r) for t in tmp]
    lw = max(b.width for b in _tbb)
    lh = max(b.height for b in _tbb)
    for t in tmp:
        t.remove()
    fig.canvas.draw()

    # ★ 逐字母手工微调（mm，+Δy = 向上）：见 `letter_offset_for` 的 `shift_mm` 说明。
    shifts_mm = shifts_mm or {}

    def _sh(lb):
        s = shifts_mm.get(str(lb).upper())
        return (float(s[0]), float(s[1])) if s else (0.0, 0.0)

    dy_px = [P[i][1] * ph * fd + offs[i][1] * P[i][3] * ph * fd + _sh(lb)[1] * px_per_mm
             for i, (ax, lb) in enumerate(plants)]

    def _boxes_at(x_px):
        out = []
        for i, pp in enumerate(P):
            sx = _sh(plants[i][1])[0] * px_per_mm
            out.append((x_px + sx, dy_px[i], x_px + sx + lw, dy_px[i] + lh))
        return out

    def _clear(boxes):
        """返回 None = 通过；否则返回一句话说明被谁挡住（供诊断）。"""
        # ★ 边界判据（曾写成 0.5 英寸的上下边距 ⇒ 上排字母被判出局、三个候选全否）：
        #   左右：距画布边 >= edge_mm（避免字母贴边；也避免它成为最左墨迹把内容拖宽）
        #   上下：只需**落在画布内**（竖直位置由各自的标准 dy 决定，本函数不负责）
        for j, b in enumerate(boxes):
            if b[0] < edge_mm * px_per_mm:
                return "盒%d 左缘 %.3f mm < 边距 %.2f mm" % (j, b[0] / px_per_mm, edge_mm)
            if b[2] > pw * fd or b[1] < 0.0 or b[3] > ph * fd:
                return "盒%d 出画布" % j
            for o in others:
                ox = min(b[2], o.x1) - max(b[0], o.x0)
                oy = min(b[3], o.y1) - max(b[1], o.y0)
                if ox > 0 and oy > 0 and (ox * oy) / (px_per_mm ** 2) > thresh_mm2:
                    return ("盒%d 压住 x=[%.2f,%.2f] y=[%.2f,%.2f] 面积 %.2f mm^2"
                            % (j, o.x0 / px_per_mm, o.x1 / px_per_mm,
                               (ph * fd - o.y1) / px_per_mm, (ph * fd - o.y0) / px_per_mm,
                               (ox * oy) / (px_per_mm ** 2)))
        return None

    del _COL_PROBE[:]
    for i, (ax, lb) in enumerate(plants):
        _COL_PROBE.append((str(lb).upper(), round(probe_x[i] / px_per_mm, 3),
                           [round(P[i][0] * pw * 25.4, 3), round(P[i][2] * pw * 25.4, 3)]))
    cands = sorted({round(x, 3) for x in probe_x}, reverse=True)
    chosen = None
    for c in cands:
        why = _clear(_boxes_at(c))
        _COL_DIAG.append((round(c / px_per_mm, 3), why))
        if why is None:
            chosen = c
            break
    if chosen is None:
        chosen = min(probe_x)
        LETTER_COLUMN_LOG[tuple(sorted(ours))] = None
    else:
        LETTER_COLUMN_LOG[tuple(sorted(ours))] = chosen / px_per_mm

    out = {}
    for i, (ax, lb) in enumerate(plants):
        sx, sy = _sh(lb)
        dx = _dx_for(P[i], chosen + sx * px_per_mm)
        dy = (dy_px[i] - P[i][1] * ph * fd) / (P[i][3] * ph * fd)
        panel_label(ax, lb, dx=dx, dy=dy, size=size)
        out[str(lb).upper()] = (dx, dy, chosen + sx * px_per_mm)
    if verbose:
        print("[align_letter_column] %s -> 公共左缘 %.3f mm (候选 %s)"
              % (sorted(ours), chosen / px_per_mm,
                 ", ".join("%.3f" % (c / px_per_mm) for c in cands)))
    return out


def align_letter_left_edges(plants, size=8, anchor="word", gap_dx_mm=None, gap_dy_mm=None,
                            x_px=None, verbose=False):
    """把若干面板字母的**左缘**对齐到同一条竖线（用户 2026-10-01 Fig4 第九轮：
       「并且三个字母也要最左端对齐」）。

    与 `align_letter_column` 的关系与差别（两者都对齐**左缘**，差别在"取哪条竖线"）：
      · `align_letter_column`：候选 = 各面板标准探针 x 去重后**从大到小**，取**最靠右**且
        不撞字、不出画布的那条 —— 用于 Fig3「D 要跟 AB 对齐」（字母都在画布内部，且有
        y 轴标题可让位）。
      · 本函数：取**最靠左**的那条标准探针 x（= 最宽的那个字母的标准位置），其余字母
        统一左移补差。理由：全局标准「字母**右缘** → 该面板 Y 轴最上面那个词的左缘
        = 6.41 mm」是**逐面板**的定义，而各面板"最上面那个词"的宽度不同 ⇒ 标准位置
        天生不同。取最靠左者 ⇒ 至少对**最宽的那个字母**仍精确满足 6.41 mm，其余字母
        只会**离标签更远**（不会更近）⇒ 不会因对齐而压到标签文字。
      ★ 副作用：字母可能越过画布左缘（本项目统一走 `bbox_inches='tight'` 裁剪，无碍）。
        **这一点与 `align_letter_column` 的 `edge_mm` 判据相反**，故不能复用那个函数。
      ★ 竖直方向仍逐面板用各自的标准 dy（本函数只统一 x）。
      ★ 与 `_nudge_letter_clear` 互斥：对齐后**不得**再逐字母微移（会把对齐破坏掉）。
      ★ 必须在**该行所有文字布置完成之后**、且**最终 scale** 下调用。

    x_px : 指定的公共左缘（px，fig.dpi）。None = 取各面板标准探针里最靠左的那个。

    用法::

        FS.align_letter_left_edges([(axa, "a"), (axc, "c"), (axe, "e")])
    """
    if not plants:
        return {}
    fig = plants[0][0].figure
    fig.canvas.draw()
    fd = float(fig.dpi)
    pw = float(fig.get_size_inches()[0])
    P = [tuple(ax.get_position().bounds) for ax, _ in plants]
    offs = [letter_offset_for(ax, gap_dx_mm=gap_dx_mm, gap_dy_mm=gap_dy_mm,
                              size=size, label=lb, anchor=anchor)
            for ax, lb in plants]
    probe_x = [P[i][0] * pw * fd + offs[i][0] * P[i][2] * pw * fd
               for i in range(len(plants))]
    chosen = float(x_px) if x_px is not None else min(probe_x)

    out = {}
    for i, (ax, lb) in enumerate(plants):
        dx = (chosen - P[i][0] * pw * fd) / (P[i][2] * pw * fd)
        dy = offs[i][1]
        panel_label(ax, lb, dx=dx, dy=dy, size=size)
        out[str(lb).upper()] = (dx, dy, chosen)
    LETTER_COLUMN_LOG[tuple(sorted(str(lb).upper() for _, lb in plants))] = \
        chosen / (fd / 25.4)
    if verbose:
        print("[align_letter_left_edges] 公共左缘 %.3f mm  (各标准探针 %s)"
              % (chosen / (fd / 25.4),
                 ", ".join("%.3f" % (p / (fd / 25.4)) for p in probe_x)))
    return out


def _avoid_default():
    """碰撞消解的**总开关**：默认开；环境变量 `LETTER_AVOID=0` 可全局关掉（做对照实验用）。

    ★★ 为什么要留这个开关（2026-10-01 实测）：消解会**移动字母** ⇒ 改变 tight bbox ⇒
      `fit_save` 重新标定 scale ⇒ **整幅几何都变**。Fig3 面板 B 被右移 5 mm 后内容变窄
      1.9%，scale 从 1.03449 顶到 1.05452，幅高 168.19 → **170.56 mm**，越过 P0 的
      ≤170 mm 红线。所以要能一键回到"纯锚定、不消解"的状态做对照。
    """
    return os.environ.get("LETTER_AVOID", "1").strip().lower() not in ("0", "false", "no", "off")


def panel_label_anchored(ax, label, size=8, gap_dx_mm=None, gap_dy_mm=None, anchor="word",
                         avoid_texts=None, shift_mm=None):
    """`panel_label` 的**锚定版**：按全局标准摆放字母（一步到位，推荐给所有图用）。

    anchor : ``"word"``（默认；锚到本面板 Y 轴最上面那个词 —— Fig2 口径）
             或 ``"corner"``（锚到本面板绘图区左上角 —— 供 S1/S2/S3 这类长分类标签面板用，
             见 `letter_offset_for` 的说明）。
    avoid_texts : 摆好后做一次**碰撞消解**（见 `_nudge_letter_clear`）。默认由
             `LETTER_AVOID` 环境变量决定（不设 = 开）。
    shift_mm : ``(Δx_mm, Δy_mm)`` 在标准锚点之上再叠加的**逐字母手工微调**（+Δy = 向上）。
    """
    if avoid_texts is None:
        avoid_texts = _avoid_default()
    dx, dy = letter_offset_for(ax, gap_dx_mm=gap_dx_mm, gap_dy_mm=gap_dy_mm,
                               size=size, label=label, anchor=anchor, shift_mm=shift_mm)
    panel_label(ax, label, dx=dx, dy=dy, size=size)
    if avoid_texts:
        dx, dy, _ = _nudge_letter_clear(ax, label, dx, dy, size=size)
    return dx, dy


# ============================================================================
# ★★★ 全局标准（2026-10-01 用户定案）：「图与图之间的距离」
# ----------------------------------------------------------------------------
#   用户原话：「记住现在图与图之间的距离，**后面所有图的距离与图2一致**」。
#   以 Fig2 的行带间距为基准（2026-10-01 第八版实测）：
#     · 相邻面板块之间的**可见空白带**（墨迹到墨迹）  = **4.74 mm**
#         —— 这是"眼睛看到的距离"，口径 = 600 dpi PNG 逐行扫「最长全白横带」
#            （`scripts/_s136_fig2_row_gap.py`），不用先决定"哪些文字算上一行"。
#     · 上面这条空白带对应的**轴框到轴框**距离            = **14.23 mm**
#         —— 两者相差 9.49 mm = 下行块的"向上外溢"（如面板字母）+ 上行块的"向下外溢"
#            （x 轴刻度数字 + x 轴标题）。**外溢量只与文字尺寸有关，与行距无关**。
#   ★ 换到别的图时**用哪一条**：
#     · 想跟用户看到的一致 ⇒ 对准 **可见空白带 4.74 mm**（推荐）；
#     · 若该图上下两块都几乎没有文字外溢 ⇒ 两条几乎相等，随便用哪条。
#   ★★ 换算公式（n 行网格，height_ratios = ratios，网格占图高比例 f = top − bottom）：
#           u      = f·H / (S + h·(S/n)·(n−1))          # u = 行高单位(mm)，S = sum(ratios)
#           行高_i  = ratios[i] · u
#           行间距  = h · (S/n) · u
#     ⇒ 已知目标行间距 gap 与"希望保持的行高"⇒ u = 行高/ratio ⇒ **h = gap·n / (S·u)**
#     ★ 关键纪律（Fig2 踩过）：**H 与 hspace 不能只动一个**。H 固定时改 hspace 只在
#       「行高 ↔ 行间距」之间重新分配 —— 会把面板整体拉高/压扁。要"只改距离、不改面板"，
#       必须**同时**改 hspace 与 H（Fig2 即：hspace 0.42→0.2521 且 图高比 7.3→6.7938）。
#   ★ Fig2 复现参数：`FIG2_H_RATIO = 6.7938`、`FIG2_HSPACE = 0.2521`（在 s37_final_figures.py）。
# ============================================================================
FIG_ROW_INK_GAP_MM = 4.74        # 可见空白带（墨迹↔墨迹），Fig2 实测
FIG_ROW_AXES_GAP_MM = 14.23      # 对应的轴框↔轴框距离，Fig2 实测


def hspace_for_row_gap(ratios, gap_mm, u_mm):
    """求令"轴框行间距 = gap_mm(mm)"的 GridSpec `hspace`（在保持既定行高的前提下）。

    参数
    ----
    ratios : 序列，GridSpec 的 height_ratios（如 [0.92, 1.25]）
    gap_mm : 目标**轴框到轴框**的行间距（mm）。若要对准"可见空白带"，先加上文字外溢量
             （= 下块上溢出 + 上块下溢出，见上方 FIG_ROW_INK_GAP_MM 注释）
    u_mm   : 行高单位（mm）= 任一现成行高 / 它的 ratio。**用它把"行高不变"这个前提固定住**。

    返回 `hspace`（相对平均行高的比例，matplotlib 口径）。
    """
    n = len(ratios)
    S = float(sum(ratios))
    return gap_mm * n / (S * float(u_mm))


# ============================================================================
# ★★★ 全局标准（2026-10-01 用户定案）：「各行面板组居中于整幅」
# ----------------------------------------------------------------------------
#   用户原话「图要整体居中」；口径经确认 = **各行面板组居中于整幅**。
#   定义：把图内**每一行**的面板组（含其刻度标签 / 轴标题 / colorbar）的**墨迹并集中心**
#         对齐到**整幅图的墨迹中心**（= 全图 tight bbox 的水平中心）。
#   ★ 为什么基准取"整幅墨迹中心"而不是"画布中心"：出图走 `bbox_inches='tight'`，成品 PNG
#     的左右边界就是墨迹边界 ⇒ 肉眼可辨的只有「行与行之间是否对齐」。以整幅墨迹中心为共同
#     基准，等价于"各行落在同一条竖直中心线上"；改用画布中心反而可能把内容推出画布边界。
#   ★★ 调用时机（硬性）：必须在 `fig.subplots_adjust(...)` **之后**、面板字母**之前**调用，
#      否则会被 gridspec 重新布局覆盖（Fig2 的 `FIG2_A_SHIFT` 踩过同一个坑）。
#      面板字母写在它之后，`letter_offset_for()` 才量得到最终的 axes 位置。
#   ★ colorbar 是**独立 axes**（`cb.ax`）⇒ 必须与宿主面板放进同一组一起平移。
# ============================================================================
def ink_center_x(axes, renderer=None):
    """一组 axes 的墨迹并集中心（inches，画布坐标）。"""
    fig = axes[0].figure
    if renderer is None:
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
    fd = float(fig.dpi)
    x0 = min(ax.get_tightbbox(renderer).x0 for ax in axes) / fd
    x1 = max(ax.get_tightbbox(renderer).x1 for ax in axes) / fd
    return 0.5 * (x0 + x1)


def center_row_groups(fig, groups, iters=3, verbose=False):
    """把每组 axes 的墨迹并集中心平移到**整幅墨迹中心**。

    groups : 可迭代；每项是一个 axes 序列（同一"行"的面板，含其 colorbar 的 `.ax`）。
    返回  : 每组的最终残差（mm，正 = 该组中心在整幅中心右侧）。

    ★ 必须在 `subplots_adjust` 之后调用（见上方注释）。
    """
    groups = [[a for a in g] for g in groups]

    def _allx(r):
        fd = float(fig.dpi)
        out = []
        for ax in fig.axes:
            b = ax.get_tightbbox(r)
            out.append((b.x0 / fd, b.x1 / fd))
        return out

    for _ in range(int(iters)):
        fig.canvas.draw()
        r = fig.canvas.get_renderer()
        fw = float(fig.get_size_inches()[0])
        allx = _allx(r)
        uc = 0.5 * (min(t[0] for t in allx) + max(t[1] for t in allx))
        shift = 0.0
        for g in groups:
            d = (uc - ink_center_x(g, r)) / fw          # figure fraction
            mm = abs(d) * fw * 25.4
            if mm > 1e-3:
                shift = max(shift, mm)
                for ax in g:
                    p = ax.get_position()
                    ax.set_position([p.x0 + d, p.y0, p.width, p.height])
        if verbose:
            print("      [center_row_groups] max shift = %.3f mm" % shift)
        if shift < 1e-3:
            break

    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    allx = _allx(r)
    uc = 0.5 * (min(t[0] for t in allx) + max(t[1] for t in allx))
    return [(ink_center_x(g, r) - uc) * 25.4 for g in groups]


def align_axes_left(ref_ax, mov_axes, fig=None, verbose=False):
    """把 `mov_axes` 的**水平位置**整体搬到与 `ref_ax` **左缘对齐**（各自宽度、纵向不动）。

    为什么需要它（2026-10-01 实测）：
      `GridSpec(3, 2)` 下 `gs[0, :]` 与 `gs[2, :]` 本是**同一列**、左右缘完全相同；
      但 `center_row_groups()` 为让它那一"行"的墨迹居中，会**逐行独立平移** ——
      Fig3 实测：A 绘图区被推到 [28.700, 178.350]、D 被推到 [14.140, 163.790]，
      左缘相差 **14.560 mm**（用户看到的就是"D 的 Y 轴没跟 A 对齐"）。
      ⇒ 居中之后用本函数把目标面板强制拉回参考面板的列上。

    ★ 调用时机（硬性）：`center_row_groups()` **之后**、面板字母锚定 **之前** ——
      字母锚定读的是 axes 的最终位置。
    ★ 副作用提醒：被搬动的行**不再满足**"该行墨迹居中"（本函数是那个标准的显式例外）。
    """
    fig = fig if fig is not None else ref_ax.figure
    fw = float(fig.get_size_inches()[0])
    pa = ref_ax.get_position()
    mm = []
    for ax in mov_axes:
        p = ax.get_position()
        mm.append((pa.x0 - p.x0) * fw * 25.4)
        ax.set_position([pa.x0, p.y0, p.width, p.height])
    if verbose:
        print("      [align_axes_left] ref x0=%.4f（%d 个 axes 位移 mm：%s）"
              % (pa.x0, len(mov_axes), ", ".join("%+.3f" % v for v in mm)))
    return mm


# ============================================================================
# ★★★ 全局标准（2026-10-01 用户定案 · Fig4 第九轮）：「Y 轴刻度文字左缘对齐」
# ----------------------------------------------------------------------------
#   用户原话「使图A、C、E的Y轴图标最左端对齐」；口径经 AskUserQuestion 裁定 =
#   **Y 轴刻度文字（tick labels）的左缘**，不是绘图区轴线、也不是两者都要。
#   ★ 与 `align_axes_left` 是**互斥的两件事**，必须先想清要哪一条：
#     · `align_axes_left` 对齐的是**绘图区左缘** ⇒ 各面板轴线共线，但文字列仍错位
#       （各面板刻度文字宽度不同：Fig4 实测 A 11.13 / C 4.40 / E 21.29 mm）；
#     · 本函数对齐的是**文字列左缘** ⇒ 文字列共线，但轴线必然按各自文字宽度**错开**
#       （x0_i = L + w_i）。用户本轮明确要后者的视觉整齐感。
#   ★ 各面板刻度文字的宽度只由**字号 + 文字内容**决定，与 `fit_save` 的 scale 无关
#     （文字以 pt 固定）⇒ 本对齐在任意 scale 下都成立，不必等收敛。
#   ★ 调用时机：与 `align_axes_left` 相同 = `center_row_groups()` **之后**、字母锚定
#     **之前**；同样是「各行组居中」的显式例外（须同步更新 `_s168` 的 MARGIN_BASE）。
#   ★ 只做**水平平移**（宽度、纵向不动）；面板的附属 axes（如 colorbar）必须一起搬。
# ============================================================================
def ytick_label_left_px(ax, renderer=None, center_test=True):
    """本面板 **y 刻度文字列的左缘**（px，fig.dpi）。无可用标签时返回 None。

    ★ 标签筛选（沿用 `letter_offset_for` 坑①的精神，但判据更宽容一档）：
      `ax.get_yticklabels()` 会返回**落在 ylim 之外、根本不会被绘制**的标签
      （Fig2 面板 A 的 '14' / '-2' 即是），必须剔除。
      · `center_test=True`（默认）：要求标签 bbox 的**竖直中心**落在绘图区内 ——
        比 `letter_offset_for` 的"整框落入"宽松 0.5 个行高，避免最靠边的那一行
        （如 Fig4 (c) 的 '0.00'，只高出轴底 0.14 mm）被误剔除。
      · `center_test=False`：整框必须落入（与 `letter_offset_for` 完全一致的口径）。
    """
    if renderer is None:
        ax.figure.canvas.draw()
        renderer = ax.figure.canvas.get_renderer()
    b = ax.get_window_extent(renderer)
    lx = None
    for t in ax.get_yticklabels():
        if not t.get_text().strip():
            continue
        tb = t.get_window_extent(renderer)
        if center_test:
            cy = 0.5 * (tb.y0 + tb.y1)
            if cy < b.y0 - 1e-6 or cy > b.y1 + 1e-6:
                continue
        else:
            if tb.y0 < b.y0 - 1e-6 or tb.y1 > b.y1 + 1e-6:
                continue
        if lx is None or tb.x0 < lx:
            lx = tb.x0
    return lx


def align_ytick_label_left(groups, iters=2, verbose=False):
    """把若干面板的 **y 刻度文字列左缘** 对齐到同一条竖线（取最靠左者为准）。

    groups : ``[(axes_seq, ref_ax), ...]``；`ref_ax` 提供本组的刻度文字列（通常就是
             该组的主面板），`axes_seq` = 需要**一起平移**的全部 axes（主面板 +
             它的 colorbar 等附属 axes）。
    返回   : 每组的最终左缘残差（mm，正 = 该组文字列仍在全组基准右侧）。

    ★ 时机：`center_row_groups()` 之后、面板字母锚定之前（见上方全局标准块）。
    ★ 与 `align_axes_left` 一样是「各行组居中」的显式例外。
    """
    groups = [([a for a in g[0]], g[1]) for g in groups]
    fig = groups[0][1].figure
    for _ in range(int(iters)):
        fig.canvas.draw()
        r = fig.canvas.get_renderer()
        fw = float(fig.get_size_inches()[0])
        lefts = [ytick_label_left_px(ref, r) for _, ref in groups]
        if any(v is None for v in lefts):
            if verbose:
                print("      [align_ytick_label_left] 有面板取不到刻度标签，跳过")
            break
        tgt = min(lefts)
        move = False
        for i, (axs, _ref) in enumerate(groups):
            d = (tgt - lefts[i]) / (fw * float(fig.dpi))       # figure fraction
            if abs(d) * fw * 25.4 > 1e-3:
                move = True
                for ax in axs:
                    p = ax.get_position()
                    ax.set_position([p.x0 + d, p.y0, p.width, p.height])
        if not move:
            break
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    lefts = [ytick_label_left_px(ref, r) for _, ref in groups]
    if verbose:
        print("      [align_ytick_label_left] 左缘 mm = %s"
              % ", ".join("n/a" if v is None else "%.3f" % (v / (float(fig.dpi) / 25.4))
                          for v in lefts))
    base = lefts[0] if lefts[0] is not None else float("nan")
    return [(v - base) / (float(fig.dpi) / 25.4) if v is not None else float("nan")
            for v in lefts]


# ============================================================================
# ★★★ 全局标准（2026-10-01 用户定案）：「图与图之间的距离」= Fig2 的可见空白带
# ----------------------------------------------------------------------------
#   用户原话「图与图之间的距离也要与模板 fig2 一致」。Fig2 实测可见空白带 = **4.74 mm**
#   （口径 = 600 dpi PNG 逐行扫「最长全白横带」，见 `_s136_fig2_row_gap.py` / `_s158`）。
#   ★★ 为什么必须"测量 + 平移"而不是"调 hspace"：
#     · 相邻两行之间**眼睛看到的距离** = 轴框间隙 − （上行向下外溢的 x 刻度/轴标题 +
#       下行向上外溢的面板字母）。**外溢量由文字尺寸决定，与 hspace 无关**；
#     · 更关键：**同一张图里两个边界的外溢量并不相等**（Fig3 实测 14.80 vs 17.55 mm，
#       因为 B/C 行下面是两行刻度 + 轴标题，而 D 行上面只有图例）⇒ **单一 hspace 数学上
#       不可能让两个边界同时等于 4.74 mm**。
#     ⇒ 本函数对**每个边界单独**测量并竖直平移其下方的所有行，直到可见空白 = 目标。
#   ★ 测量用 `fig.canvas` 的渲染缓冲（复用上一次 draw，**零额外渲染成本**）。
#   ★★ 调用时机（硬性）：与 `center_row_groups` 同为 `subplots_adjust` 之后；但本函数应
#      排在 **面板字母之后** —— 因为字母本身是"下行最上方的墨迹"，先放字母再量才准。
#      字母挂在 axes 上（axes 分数坐标）⇒ 之后竖直平移行不会破坏字母与面板的相对位移。
# ============================================================================
def ink_mask(fig):
    """渲染后的墨迹掩码（bool，fig.dpi 像素，y 向下）。复用 canvas buffer。"""
    import numpy as _np
    fig.canvas.draw()
    a = _np.asarray(fig.canvas.buffer_rgba())
    return (a[..., :3] < 250).any(axis=2)


def _longest_white_band(mask, y0, y1):
    """mask[y0:y1] 内最长的连续全白行段，返回 (start, end)（end 开区间）或 None。"""
    y0 = max(0, int(y0))
    y1 = min(int(mask.shape[0]), int(y1))
    if y1 - y0 < 2:
        return None
    rh = mask[y0:y1].any(axis=1)
    best = (0, 0, 0)
    i, n = 0, len(rh)
    while i < n:
        if not rh[i]:
            j = i
            while j + 1 < n and not rh[j + 1]:
                j += 1
            if (j - i + 1) > best[0]:
                best = (j - i + 1, i, j + 1)
            i = j + 1
        else:
            i += 1
    return (y0 + best[1], y0 + best[2]) if best[0] else None


def group_ink_y(g, renderer=None):
    """一组 axes 的**墨迹竖直范围**（inches，figure 坐标，y 向上）⇒ (y_lo, y_hi)。

    = 该组所有 axes 的 `get_tightbbox()` 并集（含标题 / 刻度标签 / 轴标题 / 面板字母 /
      colorbar）。返回 None 表示该组没有任何有效墨迹。
    """
    if renderer is None:
        g[0].figure.canvas.draw()
        renderer = g[0].figure.canvas.get_renderer()
    fd = float(g[0].figure.dpi)
    lo, hi = None, None
    for a in g:
        if not a.get_visible():
            continue
        try:
            b = a.get_tightbbox(renderer)
        except Exception:
            continue
        if b is None or b.width <= 0 or b.height <= 0:
            continue
        lo = b.y0 / fd if lo is None else min(lo, b.y0 / fd)
        hi = b.y1 / fd if hi is None else max(hi, b.y1 / fd)
    return None if lo is None else (lo, hi)


def figlevel_crossrow_artists(fig, groups, renderer=None, margin_mm=1.0):
    """找出**跨越行边界**的 figure 级 artist（典型：多行网格共用的旋转 y 轴总标题）。

    ★ 为什么必须把它们排除出「行带」的测量：全宽白带的定义要求该像素行**整幅无墨**。
      这类总标题在 figure 坐标里是**固定**的（不随行平移），纵向一跨就是好几行 ⇒ 它会把
      行边界处的白带**切成碎片**，实测值退化成「该文字自身的字间空白」，而且**对行平移完全
      不响应**（FigS3 实测：边界 1 恒读 0.762 / 1.947 mm，而真值 4.74）。
    ★ 它不属于任何一行，也不是"图与图之间的距离"的一部分 ⇒ 测量时**暂时隐藏**，交付时照常显示。

    口径：artist 的竖直范围与任一「行边界区间」（下行墨迹顶 → 上行墨迹底，各外扩 margin）
      相交 ⇒ 判为跨行。返回要隐藏的 artist 列表（调用方负责恢复 `set_visible(True)`）。
    """
    if renderer is None:
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
    fd = float(fig.dpi)
    yy = [group_ink_y(g, renderer) for g in groups]
    ivs = []
    for i in range(len(groups) - 1):
        if yy[i] is None or yy[i + 1] is None:
            continue
        ivs.append((yy[i + 1][1] - margin_mm / 25.4, yy[i][0] + margin_mm / 25.4))
    if not ivs:
        return []
    cand = (list(getattr(fig, "texts", [])) + list(getattr(fig, "legends", []))
            + list(getattr(fig, "artists", [])))
    out = []
    for a in cand:
        if not a.get_visible():
            continue
        try:
            b = a.get_window_extent(renderer)
        except Exception:
            continue
        if b is None or (b.height <= 0 and b.width <= 0):
            continue
        y0, y1 = b.y0 / fd, b.y1 / fd
        if any(not (y1 < lo or y0 > hi) for (lo, hi) in ivs):
            out.append(a)
    return out


def set_row_ink_gaps(fig, groups, targets_mm, iters=8, dpi=600, verbose=False):
    """把**相邻行组的可见墨迹间隙**逐个调到 `targets_mm`（长度 = len(groups) − 1）。

    groups     : 从上到下的行组（每组为 axes 序列）。
    targets_mm : 每个边界的可见空白带目标（mm）。Fig2 标准 = 4.74。
    返回       : 各边界的最终实测间隙（mm）。

    ★ 必须在 `subplots_adjust` 之后、且**面板字母放置之后**调用（见上方注释）。

    ★★★ 2026-10-01 四轮重写，前三次都踩到坑，最终口径 = **交付同款渲染 + 全幅白带清单 + 就近分配**：

      · 初版：`fig.canvas.buffer_rgba()` 光栅扫白带。**缺陷**：光栅只覆盖**画布矩形**，
        而本函数的竖直平移会把下方的行推出画布（FigS3 实测 row3 的 y0 被推到 −0.072），
        溢出的墨迹不在 buffer 里 ⇒ 测量静默失真（返回 2.032，真实 3.725）。
      · 二版：`get_tightbbox()` 并集（不受画布限制）。**缺陷**：`get_tightbbox` 是 artists 的
        **逻辑包围盒**，比实际墨迹大 0.1–0.8 mm ⇒ 与用户看到的 PNG **口径不一致**
        （明明调到 4.740，PNG 上量出来是 4.57 / 5.29 / 5.04）。
      · 三版：光栅 + 「窗口 = 两行 bbox 之间各外扩 ±win_mm」。**两个残留缺陷**：
        ① 外扩的 3 mm 若在行内恰好全白 ⇒ 白带被**虚增**（Fig4 实测虚增 0.42 mm）；
        ② **跨行 figure 级文字**（S3 的 4 行共用旋转 y 轴总标题）会把白带切成碎片，
           实测退化成该文字的字间空白且**对行平移零响应**（恒读 0.762 / 1.947 mm，
           `bias` 被限幅累加到 8.0 也毫无效果）。
      ⇒ 最终版：**按交付口径渲染**（`bbox_inches='tight'` + `pad_inches=0.02` + `dpi=600`，
        与 `_fig_style.save()` 完全同款 ⇒ 实测值 = 用户在成品 PNG 上量的值），扫**全幅白带清单**，
        再按「白带中心落在本边界区间 [下行墨迹顶, 上行墨迹底] 内」**就近分配**给各边界
        （区间由 bbox 给出、必为真白带的超集 ⇒ 既不会截断也不会伸进行内）。
        测量时**暂时隐藏跨行的 figure 级文字**（`figlevel_crossrow_artists`）。
      ★ Fig2 的 4.74 mm 基准本身就是**光栅口径**（`_s136` / `_s158` 在 600 dpi PNG 上扫
        最长全白横带）⇒ 用光栅口径才能"与模板 fig2 一致"。
    """
    groups = [[a for a in g] for g in groups]
    n = len(groups)
    if len(targets_mm) != n - 1:
        raise ValueError("set_row_ink_gaps: targets 数(%d) != 边界数(%d)"
                         % (len(targets_mm), n - 1))
    fh_mm = float(fig.get_size_inches()[1]) * 25.4
    targets = [float(t) for t in targets_mm]

    def _bbox_gaps():
        fig.canvas.draw()
        r = fig.canvas.get_renderer()
        yy = [group_ink_y(g, r) for g in groups]
        return [float("nan") if (yy[i] is None or yy[i + 1] is None)
                else (yy[i][0] - yy[i + 1][1]) * 25.4 for i in range(n - 1)]

    def _converge_bbox(goal):
        """用 bbox 口径（快，无需光栅）把各边界调到 goal。

        ★★ 每个边界调整前**重新测量**：调整边界 i 时平移的是 groups[i+1:]，会同时改变
          **下方所有**边界的现状；若沿用同一轮开头的测量值（stale），多边界图会过调震荡。
          FigS3（3 个边界）实测：旧写法下 `bias` 一路涨到 8.0 而第 2 个边界的 raster 读数
          卡在 0.76–1.95 mm 不动。边界数 ≥3 时这是**必现**问题。
        """
        for _ in range(int(iters)):
            moved = 0.0
            for i in range(n - 1):
                gaps = _bbox_gaps()              # ★ 每边界重测，勿用 stale
                if gaps[i] != gaps[i]:
                    continue
                d = gaps[i] - goal[i]
                if abs(d) > 0.02:
                    dy = d / fh_mm               # 间隙偏大 => 下方行上移 (+dy, y 向上)
                    for g in groups[i + 1:]:
                        for a in g:
                            p = a.get_position()
                            a.set_position([p.x0, p.y0 + dy, p.width, p.height])
                    moved = max(moved, abs(d))
            if moved < 0.02:
                break
        return _bbox_gaps()

    # ① bbox 口径快速收敛到目标；② 用**光栅**（= 用户所见）标定 bbox 与光栅的口径偏置。
    #   ★ dpi=600 = 交付口径 ⇒ 1 px = 0.042 mm，实测值与成品 PNG 逐像素可比。
    bias = [0.0] * (n - 1)
    raster = [float("nan")] * (n - 1)
    for _round in range(6):
        _converge_bbox([targets[i] + bias[i] for i in range(n - 1)])
        raster = row_ink_gaps_raster(fig, groups, dpi=dpi)
        err = [raster[i] - targets[i] for i in range(n - 1)]
        if verbose:
            print("      [set_row_ink_gaps] round %d  raster=%s  bias=%s"
                  % (_round, ["%.3f" % v for v in raster], ["%.3f" % v for v in bias]))
        if all(e != e or abs(e) <= 0.05 for e in err):      # 全部收敛（NaN 视为已够宽容）
            break
        # ★★ 符号：光栅口径与 bbox 口径只差一个**近似常数** δ = raster − bbox_gap
        #   （bbox 是 artists 的逻辑包围盒，比实际墨迹大 ⇒ 光栅 ≥ bbox）。
        #   要 raster = target ⇒ 需 bbox_gap = target − δ
        #   ⇒ 而本函数的 bbox 目标是 `target + bias` ⇒ **bias 必须取 −δ = −err**。
        #   ★ 若写成 `bias += err`，递推化为 raster_{k+1} = 2·raster_k − target（斜率 2）
        #     ⇒ **发散**（2026-10-01 实测：4.826 → 4.953 → 5.207 → 5.588）。
        #   限幅 ±2.5 mm：防止读数异常时把行一次性压到重叠。
        bias = [bias[i] - (0.0 if err[i] != err[i] else max(-2.5, min(2.5, err[i])))
                for i in range(n - 1)]
    return raster


def row_ink_gaps_raster(fig, groups, dpi=600, pad_inches=0.02, hide_crossrow=True,
                        min_band_mm=0.3):
    """各相邻行组的**可见墨迹间隙**（mm），口径 = **交付同款光栅**（= 用户在成品 PNG 上看到的）。

    实现三步：
      ① 暂时隐藏**跨行的 figure 级文字**（如 4x4 网格共用的旋转 y 轴总标题，见
         `figlevel_crossrow_artists`）—— 否则它会把白带切碎；
      ② 用**与 `_fig_style.save()` 完全同款**的参数渲染：`dpi=600, bbox_inches='tight',
         pad_inches=0.02` ⇒ 得到的像素与交付 PNG 逐像素可比；
      ③ 扫**全幅白带清单**（连续全白像素行，长度 ≥ `min_band_mm`），再按
         「白带中心 ∈ [下行墨迹顶, 上行墨迹底]（bbox 口径）」就近分配给各边界。
         ★ bbox 是墨迹的**超集** ⇒ 该区间必包含真白带、又不会伸进行内 ⇒ 既不截断也不虚增。
    """
    import io as _io
    import numpy as _np
    from PIL import Image as _Image

    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    yy = [group_ink_y(g, r) for g in groups]
    n = len(groups)

    hidden = figlevel_crossrow_artists(fig, groups, r) if hide_crossrow else []
    for a in hidden:
        a.set_visible(False)
    try:
        fb = fig.get_tightbbox(r)                    # ★ inches（figure 坐标）
        buf = _io.BytesIO()
        fig.savefig(buf, format="png", dpi=dpi, bbox_inches="tight",
                    pad_inches=pad_inches, facecolor="white")
    finally:
        for a in hidden:
            a.set_visible(True)
    buf.seek(0)
    arr = _np.asarray(_Image.open(buf).convert("RGB"))
    buf.close()
    mask = (arr < 250).any(axis=2)
    ppm = dpi / 25.4
    top_in = float(fb.y1) + pad_inches                # PNG 第 0 行对应的 figure y（inch）

    def _px(y_in):
        return (top_in - y_in) * dpi

    # 全幅白带清单
    rh = mask.any(axis=1)
    bands = []
    i, nh = 0, len(rh)
    while i < nh:
        if not rh[i]:
            j = i
            while j + 1 < nh and not rh[j + 1]:
                j += 1
            if (j - i + 1) / ppm >= min_band_mm:
                bands.append((i, j + 1))
            i = j + 1
        else:
            i += 1

    out = []
    for k in range(n - 1):
        if yy[k] is None or yy[k + 1] is None:
            out.append(float("nan"))
            continue
        lo, hi = _px(yy[k + 1][1]), _px(yy[k][0])
        if lo > hi:
            lo, hi = hi, lo
        cand = [b for b in bands if lo <= 0.5 * (b[0] + b[1]) <= hi]
        if not cand:                                  # 兜底：取与区间相交的最长白带
            cand = [b for b in bands if b[1] > lo and b[0] < hi]
        if not cand:
            out.append(float("nan"))
            continue
        best = max(cand, key=lambda b: b[1] - b[0])
        out.append((best[1] - best[0]) / ppm)
    return out


# ============================================================================
# ★★★ 全局标准（2026-10-01 用户定案）：「点的大小与风格」—— 以 Fig2 为模板
# ----------------------------------------------------------------------------
#   用户原话：「fig2 现在不改动了，**记住现在 fig2 点的大小与风格**，后面的图风格与其一致」。
#   三条口径（Fig2 定案时用户逐条选定，此处原样推广）：
#     ① **放大**：**面积 ×2 ⇒ 直径 ×√2 = 1.41421**（数值直径 `ms` 直接乘 √2）。
#        ★ 对 `scatter` 的 `s`（matplotlib 的 `s` 就是"直径(pt)的平方" ⇒ 面积口径）：
#           **`s` 直接 ×2** 即等价于直径 ×√2。`s_up()` 负责这层换算。
#     ② **描边**：同色系**深 2 个度** = CIE L* 降 20
#        （`darken()`：保色相 H、保 HSL 饱和度 S，只降 L*）。
#     ③ **描边宽**：`mew_for(ms)` = 数值直径的 10%，clamp 到 **[0.30, 0.70] pt**
#        （落在 Nature 描边 0.25–1 pt 窗口内）。
#   ★★ 两个例外（Fig2 已采用，推广时同样适用）：
#     · **"纯描边型"标记保留原 `mew`** —— 指 `mfc` 为白/空的标记（Fig2 的 (c) 白底复核
#       方块 `mew=0.7`、`[NR]` 叉 `mew=1.0`）。那已经是它们唯一的线重，按比例缩放会改变
#       可辨性；只把它们**边线色**换成同色系深 2 度、**尺寸**照常放大。
#     · **图例标记的绝对尺寸保持不变** ⇒ 图上点放大 k 倍后 `markerscale` 必须**同除 k**
#       （用 `markerscale_keep()`）。用户要放大的是**图上的点**，不是图例。
#   ★ 回退：把 `P_UP` 与 `P_DL` 当作唯一开关 —— `P_UP=1.0, P_DL=0` 即回到"未放大、无描边"。
# ============================================================================
P_UP = 1.41421      # 面积 ×2 ⇒ 直径 ×√2
P_DL = 20.0         # 同色系加深 2 个度（1 度 = CIE L* −10）


def pt_up(ms, k=P_UP):
    """把点的**数值直径**按「面积 ×k²」放大后返回新直径（默认 k=√2 ⇒ 面积 ×2）。"""
    return float(ms) * float(k)


def s_up(s, k=P_UP):
    """把 `scatter` 的 **`s`**（面积口径，单位 pt²）按「面积 ×k²」放大后返回。

    ★ `matplotlib` 的 `s` 就等于"标记直径(pt)²" ⇒ **面积 ×2 ⇔ `s` ×2 ⇔ 直径 ×√2**。
      所以本函数体内的乘法是 `k*k`（k=√2 时正好 ×2），与 `pt_up` 口径自洽。
    """
    return float(s) * float(k) * float(k)


def ms_of_s(s):
    """`scatter` 的 `s`（pt²）反推**等价直径 `ms`**（pt），供 `mew_for` 使用。"""
    import math as _m
    return _m.sqrt(float(s))


def mew_for(ms):
    """描边宽度 = 数值直径的 10%，clamp 到 [0.30, 0.70] pt。"""
    return min(0.70, max(0.30, 0.10 * float(ms)))


def edge_color(hexv, dl=P_DL):
    """同色系描边色 = 该点色按 CIE L* 降 `dl`（保色相 H、保 HSL 饱和度 S）。"""
    return darken(hexv, dl)


def markerscale_keep(orig, k=P_UP):
    """图上点放大 k 倍后，要让**图例标记绝对尺寸不变** ⇒ `markerscale` 同除 k。"""
    return float(orig) / float(k)


def finalize(ax):
    for sp in ax.spines.values():
        sp.set_linewidth(0.5)
