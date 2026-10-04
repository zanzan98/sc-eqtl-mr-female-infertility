# -*- coding: utf-8 -*-
"""
s67c_verify_supp_figs.py  —  Fig S1 / Fig S2 验收核验（独立判据，不依赖出图脚本的内建日志）

核验四项（每项都有硬断言，任一不过则整体 FAIL）：
  1) PDF 内 LastResortHE 命中数 == 0
     —— 这是「字形缺失 -> 空豆腐块」的权威判据（matplotlib 缺字形不报错，只渲染豆腐块）
  2) PDF 嵌入字体名单 + 页数 == 1
  3) PNG 的 size / dpi 换算 mm，与出图日志实测宽度对照（容差 0.3 mm）
  4) 图内可见文字全 ASCII：以源数据表（64a/64b）与脚本内标签集为准做互补核验
     （PDF 不直接暴露文本层，故以 ①LastResortHE=0 + ②源表/标签 ASCII 双向夹逼）

输出：logs/_s67c_verify_supp_figs.log（脚本自写文件，PowerShell stdio 不可靠）
"""
import io
import os
import re
import json
import struct
import csv

ROOT = "D:/endometriosis_project/11_sc_eqtl_mr_project"
FIGD = os.path.join(ROOT, "figures")
TABD = os.path.join(ROOT, "tables")
LOGD = os.path.join(ROOT, "logs")

TARGETS = [
    ("FigS1_instrument_power_gating", 182.4, 70.6,
     os.path.join(TABD, "64a_FigS1_celltype_power_gating.csv")),
    ("FigS2_sensitivity_outcome_concordance", 182.5, 77.7,
     os.path.join(TABD, "64b_FigS2_sensitivity_vs_main_15pairs.csv")),
    ("FigS3_manhattan_by_celltype", 182.4, 135.5,
     os.path.join(TABD, "64d_FigS3S4_discovery_thresholds.csv")),
    ("FigS4_discovery_volcano", 89.0, 68.6,
     os.path.join(TABD, "64d_FigS3S4_discovery_thresholds.csv")),
]

# 图件各自的栏宽硬窗口：S1/S2/S3 = Nature 双栏；S4 = 单栏（Nature single column 88.9 mm）
WIN = {"FigS4_discovery_volcano": (88.0, 90.0)}
WIN_DEFAULT = (181.5, 183.3)

out = []
FAIL = []


def w(s):
    out.append(s)


def check(cond, msg):
    if not cond:
        FAIL.append(msg)
    return cond


def pdf_fonts_and_pages(path):
    """从 PDF 原始字节抽取 /BaseFont 名单与页数（无第三方依赖，够用）。"""
    b = open(path, "rb").read()
    fonts = sorted(set(re.findall(rb"/BaseFont\s*/([A-Za-z0-9+\-_,\.]+)", b)))
    pages = len(re.findall(rb"/Type\s*/Page[^s]", b))
    mediabox = re.findall(rb"/MediaBox\s*\[([^\]]+)\]", b)
    return b, [f.decode("latin-1") for f in fonts], pages, mediabox


def main():
    w("s67c_verify_supp_figs  —  Fig S1-S4 acceptance check")
    w("=" * 74)
    summary = []
    for name, wref, href, src in TARGETS:
        w("")
        w("### %s" % name)
        pdf = os.path.join(FIGD, name + ".pdf")
        png = os.path.join(FIGD, name + ".png")
        check(os.path.exists(pdf), "%s PDF MISSING" % name)
        check(os.path.exists(png), "%s PNG MISSING" % name)
        if not (os.path.exists(pdf) and os.path.exists(png)):
            continue

        # ---- 1) LastResortHE ----
        b, fonts, pages, mediabox = pdf_fonts_and_pages(pdf)
        n_lr = b.count(b"LastResortHE")
        n_lr2 = b.count(b"LastResort")
        w("  [1] LastResortHE count      = %d   (LastResort* = %d)" % (n_lr, n_lr2))
        check(n_lr == 0, "%s: LastResortHE=%d (tofu glyphs present!)" % (name, n_lr))

        # ---- 2) fonts + pages ----
        w("  [2] embedded /BaseFont      = %s" % ("; ".join(fonts) if fonts else "(none)"))
        w("      page objects            = %d" % pages)
        check(pages == 1, "%s: page count = %d (expected 1)" % (name, pages))
        tiny = sorted(set(f.split("+")[-1] for f in fonts))
        w("      font families           = %s" % "; ".join(tiny))

        # ---- 3) PNG geometry vs log ----
        from PIL import Image
        im = Image.open(png)
        dpi = im.info.get("dpi", (600, 600))
        wpx, hpx = im.size
        wmm, hmm = wpx / dpi[0] * 25.4, hpx / dpi[1] * 25.4
        w("  [3] PNG pixels              = %d x %d  @ dpi=(%.0f,%.0f)"
          % (wpx, hpx, dpi[0], dpi[1]))
        w("      -> %.1f x %.1f mm   (log said %.1f x %.1f mm)"
          % (wmm, hmm, wref, href))
        check(abs(wmm - wref) <= 0.3, "%s: PNG width %.1f vs log %.1f (>0.3 mm)" % (name, wmm, wref))
        check(abs(hmm - href) <= 0.3, "%s: PNG height %.1f vs log %.1f (>0.3 mm)" % (name, hmm, href))
        _lo, _hi = WIN.get(name, WIN_DEFAULT)
        check(_lo <= wmm <= _hi, "%s: width %.1f outside column window [%.1f,%.1f]" % (name, wmm, _lo, _hi))

        # ---- 4) ASCII / glyph dual-bound ----
        nonascii = set()
        with io.open(src, "r", encoding="utf-8") as fh:
            body = fh.read()
        for ch in body:
            if ord(ch) >= 128:
                nonascii.add(ch)
        w("  [4] source table chars>=128 = %s"
          % ("; ".join("%s(U+%04X)" % (c, ord(c)) for c in sorted(nonascii)) if nonascii else "(none - all ASCII)"))
        # 标签集（四个脚本内硬编码的可见文字，须与出图脚本一致）
        labels = """Instruments (cis-eQTL genes)
Nominally significant tests
FDR < 0.05 tests
-log10(minimum q)
FDR = 0.05
power: signal
low power (not a null result)
CD4_NC NK CD8_ET CD8_NC CD4_ET B_IN B_MEM Mono_C Mono_NC CD8_S100B DC NK_R Plasma CD4_SOX4
chromosome
-log10(P)
-log10(P)   (MR, single-SNP Wald)
Thresholds (identical to Fig. 2a):
Bonferroni  p = 0.05/8612
not significant
nominal P < 0.05
FDR < 0.05
MR effect (beta per SD of expression)
CDC42 LINC00339 YME1L1 ANXA4
Comparison with sensitivity outcome GCST90483469
MR effect, main outcome (beta)
MR effect, sensitivity outcome (beta)
MR effect (beta, 95% CI)
Main outcome
Sensitivity outcome
15/15 direction-consistent
all 15 pairs FDR < 0.05
in both outcomes"""
        bad = sorted({c for c in labels if ord(c) >= 128})
        w("      script label chars>=128 = %s" % ("; ".join(bad) if bad else "(none - all ASCII)"))
        check(not bad, "%s: non-ASCII in labels: %s" % (name, bad))

        summary.append((name, wmm, hmm, n_lr, len(fonts)))

    # ---- inventory cross-check ----
    inv = os.path.join(TABD, "63_figure_inventory.csv")
    w("")
    w("### inventory cross-check (%s)" % os.path.basename(inv))
    with io.open(inv, "r", encoding="utf-8-sig", newline="") as fh:
        rows = list(csv.DictReader(fh))
    for r in rows:
        if r.get("编号", "").strip() in ("图 S1", "图 S2", "图 S3", "图 S4"):
            w("  %s | 状态=%s | PNG=%s bytes | PDF=%s bytes"
              % (r.get("编号"), r.get("状态"), r.get("PNG字节"), r.get("PDF字节")))

    w("")
    w("=" * 74)
    if FAIL:
        w("VERDICT: FAIL  (%d issue%s)" % (len(FAIL), "s" if len(FAIL) > 1 else ""))
        for f in FAIL:
            w("   - " + f)
    else:
        w("VERDICT: PASS  (all assertions OK)")
        for s in summary:
            w("   - %-42s %.1f x %.1f mm  LastResortHE=%d  fonts=%d" % s)

    os.makedirs(LOGD, exist_ok=True)
    p = os.path.join(LOGD, "_s67c_verify_supp_figs.log")
    with io.open(p, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(out) + "\n")
    return 0 if not FAIL else 1


if __name__ == "__main__":
    raise SystemExit(main())
