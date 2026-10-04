# -*- coding: utf-8 -*-
"""s119_regen_panels_upper.py —— 面板字母大写后重出 s37 的 Fig3/4/5。

背景：`_fig_style.panel_label` 已加 `.upper()`（问题3），故凡走该函数的图都需重出。
      Fig2 已由 `s118_fig2_revise.py` 单独重出（含布局与配色改动），此处不重复。

★ FIGPAL=gb12、FIG_NOTES=0，不设 FIG_OUTDIR ⇒ 落盘到项目 figures/。
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
import s37_final_figures as S           # noqa: E402

LOG = r"D:/endometriosis_project/11_sc_eqtl_mr_project/logs/_s119_upper_regen.txt"
L = []


def w(s=""):
    L.append(str(s))


JOBS = [
    (S.fig3, "Fig3_coloc_sensitivity"),
    (S.fig4, "Fig4_chr1p36_12_finemap"),
    (S.fig5, "Fig5_opentargets_phewas_safety"),
]

w("=" * 96)
w("面板字母大写后重出 Fig3/4/5    FIGDIR = %s" % S.FS.FIGDIR)
w("=" * 96)
w("%-34s %-10s %-20s %s" % ("name", "scale", "实测 mm", "验收"))
w("-" * 96)
bad = 0
for fn, name in JOBS:
    paths, wmm, hmm = S.fit_save(fn, name)
    ok = (181.5 <= wmm <= 183.3) and (hmm <= 170.0)
    if not ok:
        bad += 1
    w("%-34s %-10.5f %-20s %s" % (name, S.CAL.get(name, float("nan")),
                                  "%.2f x %.2f" % (wmm, hmm), "OK" if ok else "FAIL"))
w("")
w("通过 %d/%d ；不合格 %s" % (len(JOBS) - bad, len(JOBS), "无" if not bad else str(bad)))

with io.open(LOG, "w", encoding="utf-8", newline="\n") as fh:
    fh.write("\n".join(L) + "\n")
print("done -> %s" % LOG)
