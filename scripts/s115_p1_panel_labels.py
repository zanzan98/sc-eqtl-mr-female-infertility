#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""s115_p1_panel_labels.py —— P1：重出 s37 的 Fig2/Fig3/Fig4/Fig5。

改动来源（均在本次 P1）：
  · `_fig_style.panel_label` 默认 size 9 -> 8  ⇒ Fig2/3/4/5 的面板字母随之 8 pt
  · `s37.fig2` 红叉标记 mew 1.1 -> 1.0 pt
★ 不动 Fig1（s37.fig1 是流程示意图，无面板字母、无描边改动）。
★ 硬闸门：幅宽须落在 181.5-183.3 mm；幅高须 <= 170 mm。
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

JOBS = [(S.fig2, "Fig2_discovery_replication"),
        (S.fig3, "Fig3_coloc_sensitivity"),
        (S.fig4, "Fig4_chr1p36_12_finemap"),
        (S.fig5, "Fig5_opentargets_phewas_safety")]

LOGDIR = os.environ.get("FIG_OUTDIR") or os.path.join(FS.ROOT, "logs")
os.makedirs(LOGDIR, exist_ok=True)
LOG = os.path.join(LOGDIR, "_s115_p1_s37.txt")
SQ = os.path.join(FS.FIGDIR, "_fig_scale.json")
L = []


def w(s=""):
    L.append(str(s))


def main():
    w("=" * 96)
    w("P1 面板字母 8 pt + Fig2 描边<=1.0 pt · 重出 Fig2/3/4/5   FIGDIR = %s" % FS.FIGDIR)
    w("=" * 96)
    w("%-34s %-12s %-16s %s" % ("name", "scale", "实测 mm", "验收"))
    w("-" * 96)
    bad = []
    for fn, name in JOBS:
        try:
            _p, wmm, hmm = S.fit_save(fn, name)
            good = (181.5 <= wmm <= 183.3) and (hmm <= 170.0)
            if not good:
                bad.append((name, wmm, hmm))
            w("%-34s %-12.5f %-16s %s" %
              (name, S.CAL[name], "%.2f x %.2f" % (wmm, hmm), "OK" if good else "FAIL"))
        except Exception as e:                                   # noqa: BLE001
            w("%-34s  !! FAILED: %s" % (name, e))
            traceback.print_exc()
    old = {}
    if os.path.exists(SQ):
        try:
            with io.open(SQ, encoding="utf-8") as fh:
                old = json.load(fh)
        except Exception:
            old = {}
    old.setdefault("scales", {})
    old["scales"].update(S.CAL)
    old["p1_panel_label_pt"] = 8
    with io.open(SQ, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(json.dumps(old, indent=1, ensure_ascii=False) + "\n")
    w("")
    w("通过 %d/%d ；不合格 %s" % (len(JOBS) - len(bad), len(JOBS), bad if bad else "无"))
    with io.open(LOG, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(L) + "\n")
    print("done bad=%d -> %s" % (len(bad), LOG))
    return 0 if not bad else 1


if __name__ == "__main__":
    sys.exit(main())
