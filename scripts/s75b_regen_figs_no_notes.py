# -*- coding: utf-8 -*-
"""s75b · FIG_NOTES=0 重出 Fig2–Fig5（**不重出 Fig1**）

用户裁定 2026-09-29：仅删除 Fig2–Fig5 的图内脚注（panel notes）。
做法：设 FIG_NOTES=0 -> s37 内 4 处 fig.text 脚注跳过；subplots_adjust 不动，
      故面板绝对位置/尺寸不变，仅 tight 裁剪后整图变矮。
"""
import io, os, sys, json, hashlib

ROOT = r"D:/endometriosis_project/11_sc_eqtl_mr_project"
os.environ["FIG_NOTES"] = "0"
sys.path.insert(0, os.path.join(ROOT, "scripts"))
LOG = os.path.join(ROOT, "logs", "_s75b_regen_no_notes.txt")

out = []
def w(s=""):
    out.append(str(s)); print(s)

OLD_SCALE = {}
try:
    OLD_SCALE = json.load(io.open(os.path.join(ROOT, "figures", "_fig_scale.json"), encoding="utf-8")).get("scales", {})
except Exception as e:
    w("读旧 scale 失败: %s" % e)

import s37_final_figures as S37
w("SHOW_NOTES = %s" % S37.SHOW_NOTES)
assert S37.SHOW_NOTES is False, "FIG_NOTES=0 未生效"

JOBS = [(S37.fig2, "Fig2_discovery_replication", "Fig2_discovery_replication"),
        (S37.fig3, "Fig3_coloc_sensitivity", "Fig3_coloc_sensitivity"),
        (S37.fig4, "Fig4_chr1p36_12_finemap", "Fig4_chr1p36_12_finemap"),
        (S37.fig5, "Fig5_opentargets_phewas_safety", "Fig5_opentargets_phewas_safety")]

from PIL import Image
res = {}
for fn, name, png in JOBS:
    _p, wmm, hmm = S37.fit_save(fn, name)
    png_p = os.path.join(ROOT, "figures", png + ".png")
    pdf_p = os.path.join(ROOT, "figures", png + ".pdf")
    px = Image.open(png_p).size
    md = hashlib.md5(open(png_p, "rb").read()).hexdigest()
    old_s = OLD_SCALE.get(name, float("nan"))
    new_s = S37.CAL.get(name, float("nan"))
    res[name] = dict(w_mm=wmm, h_mm=hmm, px=px, md5=md,
                     scale_new=new_s, scale_old=old_s,
                     bytes_png=os.path.getsize(png_p), bytes_pdf=os.path.getsize(pdf_p))
    w("  [OK] %-32s  %.2f x %.2f mm   px=%s   scale %.5f -> %.5f" % (name, wmm, hmm, px, old_s, new_s))

with io.open(os.path.join(ROOT, "logs", "_s75b_new_sizes.json"), "w", encoding="utf-8", newline="\n") as fh:
    json.dump(res, fh, indent=1, ensure_ascii=False)
w("")
w("尺寸/哈希表 -> logs/_s75b_new_sizes.json")
w("DONE %d/4" % len(res))

with io.open(LOG, "w", encoding="utf-8", newline="\n") as fh:
    fh.write("\n".join(out) + "\n")
