#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""s114_p0_shrink_height.py —— P0：把 Fig3 / Fig4 / Fig5 压到幅高 <= 170 mm。

★ 只重出 **Fig3/4/5**（不动 Fig1/Fig2 的任何字节）。
★ 硬闸门：出图后逐图断言 h <= 170 mm，否则脚本退出码 1 并如实报告。
★ 输出目录随 FIG_OUTDIR（干跑 = 临时目录；正式 = figures/）。
★ 结果写文件再读；不 print 非 ASCII。

用法（必须显式设环境变量，防误覆盖 figures/ 正本）：
  FIGPAL=gb12  FIG_NOTES=0  FIG_OUTDIR=<dir>  python s114_p0_shrink_height.py
"""
import io
import json
import os
import sys
import traceback

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

assert os.environ.get("FIGPAL") == "gb12", "必须显式 FIGPAL=gb12"
assert os.environ.get("FIG_NOTES") == "0", "必须显式 FIG_NOTES=0"

import s37_final_figures as S      # noqa: E402
import _fig_style as FS            # noqa: E402

HARD_H = 170.0
JOBS = [(S.fig3, "Fig3_coloc_sensitivity"),
        (S.fig4, "Fig4_chr1p36_12_finemap"),
        (S.fig5, "Fig5_opentargets_phewas_safety")]

LOGDIR = os.environ.get("FIG_OUTDIR") or os.path.join(FS.ROOT, "logs")
os.makedirs(LOGDIR, exist_ok=True)
LOG = os.path.join(LOGDIR, "_s114_p0_shrink.txt")
SQ = os.path.join(FS.FIGDIR, "_fig_scale.json")
L = []


def w(s=""):
    L.append(str(s))


def main():
    w("=" * 96)
    w("P0 幅高压缩 · Fig3/4/5   FIGDIR = %s" % FS.FIGDIR)
    w("=" * 96)
    w("LAY = %s" % json.dumps(S.LAY, ensure_ascii=False, sort_keys=True))
    w("")
    w("%-34s %-12s %-14s %s" % ("name", "scale", "实测 mm", "h<=170?"))
    w("-" * 96)
    ok = 0
    bad = []
    for fn, name in JOBS:
        try:
            _p, wmm, hmm = S.fit_save(fn, name)
            good = hmm <= HARD_H
            ok += 1 if good else 0
            if not good:
                bad.append((name, hmm))
            w("%-34s %-12.5f %-14s %s" %
              (name, S.CAL[name], "%.2f x %.2f" % (wmm, hmm), "OK" if good else "FAIL"))
        except Exception as e:                                   # noqa: BLE001
            w("%-34s  !! FAILED: %s" % (name, e))
            traceback.print_exc()
    # 合并 scale（保留未重出图的既有值）
    old = {}
    if os.path.exists(SQ):
        try:
            with io.open(SQ, encoding="utf-8") as fh:
                old = json.load(fh)
        except Exception:
            old = {}
    old.setdefault("scales", {})
    old["target_mm"] = S.TARGET_MM
    old["scales"].update(S.CAL)
    old["p0_height_cap_mm"] = HARD_H
    with io.open(SQ, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(json.dumps(old, indent=1, ensure_ascii=False) + "\n")
    w("")
    w("scale 已合并 -> %s" % SQ)
    w("通过 %d/%d ；超限 %s" % (ok, len(JOBS), bad if bad else "无"))
    with io.open(LOG, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(L) + "\n")
    print("done ok=%d/%d -> %s" % (ok, len(JOBS), LOG))
    return 0 if not bad else 1


if __name__ == "__main__":
    sys.exit(main())
