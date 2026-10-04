# -*- coding: utf-8 -*-
"""s111_fig_audit_nature.py —— 对 figures/ 全部图件做「Nature 规范符合性」实测盘点。

测什么（全部从**已落盘的 PDF** 读，不读 rcParams —— 只有产物是 ground truth）：
  ① 幅宽 / 幅高（mm）：对照 89 / 183 / ≤170 mm
  ② 正文文字字号直方图（pt，按字符数加权）+ 低于 5 pt / 高于 7 pt 的清单
  ③ 字体清单（须全为 Arial 族子集；出现 LastResort* = 豆腐块）
  ④ 描边宽度清单（Nature 要求最终尺寸下 0.25–1 pt）
  ⑤ 内嵌位图数量（超大量点云须 rasterized，否则 PDF 膨胀）

★ 只读不改：不写 figures/，不重出任何图。
"""
import collections
import glob
import io
import os

import pymupdf

ROOT = r"D:/endometriosis_project/11_sc_eqtl_mr_project"
FIGD = os.path.join(ROOT, "figures")
LOG = os.path.join(ROOT, "logs", "_s111_fig_audit_nature.txt")

IN_PAPER = {"Fig1_study_design", "Fig2_discovery_replication", "Fig3_coloc_sensitivity",
            "Fig4_chr1p36_12_finemap", "Fig5_opentargets_phewas_safety",
            "FigS1_instrument_power_gating", "FigS2_sensitivity_outcome_concordance",
            "FigS3_manhattan_by_celltype", "FigS4_discovery_volcano"}

L = []


def W(s=""):
    L.append(str(s))


def audit(pdf):
    d = pymupdf.open(pdf)
    pg = d[0]
    r = pg.rect
    wmm, hmm = r.width / 72.0 * 25.4, r.height / 72.0 * 25.4

    sizes = collections.Counter()          # pt -> 字符数
    fonts = collections.Counter()
    small, big = [], []
    for blk in pg.get_text("dict")["blocks"]:
        for ln in blk.get("lines", []):
            for sp in ln["spans"]:
                t = sp["text"]
                if not t.strip():
                    continue
                sz = round(sp["size"], 2)
                sizes[sz] += len(t)
                fonts[sp["font"]] += len(t)
                if sz < 5.0:
                    small.append((sz, t[:34]))
                elif sz > 7.05:
                    big.append((sz, t[:34]))

    widths = collections.Counter()
    nfill = 0
    try:
        for dr in pg.get_drawings():
            w = dr.get("width")
            if w is None or w == 0:
                nfill += 1
            else:
                widths[round(w, 3)] += 1
    except Exception as e:                                   # noqa: BLE001
        widths[("ERR", repr(e)[:40])] += 1

    nimg = len(pg.get_images(full=True))
    d.close()
    return dict(wmm=wmm, hmm=hmm, sizes=sizes, fonts=fonts, widths=widths,
                nfill=nfill, nimg=nimg, small=small, big=big)


targets = sorted(glob.glob(os.path.join(FIGD, "Fig*.pdf")))
W("=" * 104)
W("figures/ 全库图件 · Nature 规范符合性实测（数据来源 = 已落盘 PDF，1 pt = 最终尺寸 1 pt）")
W("=" * 104)
W("%-34s %-13s %-6s %-9s %-26s %-8s %s"
  % ("file", "mm (W x H)", "入论文", "字号区间pt", "字体", "描边pt", "内嵌位图"))
W("-" * 104)

summary = []
for p in targets:
    stem = os.path.basename(p)[:-4]
    r = audit(p)
    wl = [k for k in r["widths"] if isinstance(k, float)]
    wl_s = "%.2f-%.2f" % (min(wl), max(wl)) if wl else "-"
    szs = sorted(k for k in r["sizes"] if isinstance(k, float))
    sz_s = "%.1f-%.1f" % (szs[0], szs[-1]) if szs else "-"
    fonts = ",".join(sorted(r["fonts"]))
    W("%-34s %-13s %-6s %-9s %-26s %-8s %d"
      % (stem, "%.1f x %.1f" % (r["wmm"], r["hmm"]),
         "Y" if stem in IN_PAPER else "-", sz_s, fonts[:26], wl_s, r["nimg"]))
    summary.append((stem, r))

# ---- 逐图细目 ----
W("")
W("=" * 104)
W("细目")
W("=" * 104)
for stem, r in summary:
    W("")
    W("### %s   %.1f x %.1f mm   %s"
      % (stem, r["wmm"], r["hmm"], "(in paper)" if stem in IN_PAPER else ""))
    W("   字号直方图（pt: 字符数）：%s"
      % ", ".join("%.2f:%d" % (k, v) for k, v in sorted(r["sizes"].items())))
    lo = [k for k in r["sizes"] if isinstance(k, float) and k < 5.0]
    hi = [k for k in r["sizes"] if isinstance(k, float) and k > 7.05]
    W("   低于 5 pt：%s" % (("量=%d（字号 %s）" % (sum(r["sizes"][k] for k in lo),
                                              ",".join("%.2f" % k for k in sorted(lo))))
                          if lo else "无"))
    if hi:
        W("   高于 7 pt：量=%d（字号 %s）"
          % (sum(r["sizes"][k] for k in hi), ",".join("%.2f" % k for k in sorted(hi))))
    if r["small"]:
        ex = collections.Counter(t for _s, t in r["small"])
        W("   小字样例：%s" % "; ".join("%.2fpt %r" % (s, t)
                                    for s, t in r["small"][:6]))
    W("   字体：%s" % ", ".join("%s(%d)" % (k, v) for k, v in r["fonts"].most_common()))
    W("   描边宽度(pt: 次数)：%s ；无描边填充路径 %d 条"
      % (", ".join("%s:%d" % (k, v) for k, v in sorted(r["widths"].items(), key=str)),
         r["nfill"]))

# ---- 判定 ----
W("")
W("=" * 104)
W("判定（Nature 官方硬指标）")
W("=" * 104)
W("· 幅宽须 ∈ {89, 183} mm（或 120–136 mm 的 1.5 栏）；最大高 170 mm")
W("· 正文文字 5–7 pt（面板字母 8 pt 粗体另计）；描边 0.25–1 pt")
W("· 字体须 Arial/Helvetica 族；不得出现 LastResort*（= 缺字形豆腐块）")
W("· 主图须矢量（.pdf/.eps 首选），禁 .png/.tiff/.jpeg")

bad_w = [s for s, r in summary if not (88.0 <= r["wmm"] <= 90.0 or 180.5 <= r["wmm"] <= 183.5
                                       or 119.0 <= r["wmm"] <= 137.0)]
bad_h = [s for s, r in summary if r["hmm"] > 170.0]
bad_fs = [s for s, r in summary if any(isinstance(k, float) and k < 5.0 for k in r["sizes"])]
bad_fw = [s for s, r in summary
          if any(isinstance(k, float) and (k < 0.25 or k > 1.0) for k in r["widths"])]
bad_ft = [s for s, r in summary if any("LastResort" in k or "Arial" not in k
                                       and "Helvetica" not in k for k in r["fonts"])]
W("")
W("幅宽不合规：%s" % (bad_w or "无"))
W("幅高超 170 mm：%s" % (bad_h or "无"))
W("存在 <5 pt 文字：%s" % (bad_fs or "无"))
W("存在 <0.25 或 >1 pt 描边：%s" % (bad_fw or "无"))
W("字体非 Arial/Helvetica 族：%s" % (bad_ft or "无"))

os.makedirs(os.path.dirname(LOG), exist_ok=True)
with io.open(LOG, "w", encoding="utf-8", newline="\n") as fh:
    fh.write("\n".join(L) + "\n")
print("audit done ->", LOG)
