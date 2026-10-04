# -*- coding: utf-8 -*-
"""s121_update_deliverables.py —— 更新交付清单 + 核验源数据完整性。

① 更新 `tables/63_figure_inventory.csv`：对本次**实际重绘过**的图件，刷新
   PNG/PDF 字节、出图日期、状态列里的尺寸备注，并把面板数一栏的字母改为大写。
② 只读核验 `supplementary_data/_manifest.csv` 里 49 条产物的 sha256
   —— 证明本次图件改动**没有碰任何数据文件**（红线 1）。
"""
import csv
import hashlib
import io
import os

ROOT = r"D:/endometriosis_project/11_sc_eqtl_mr_project"
INV = os.path.join(ROOT, "tables", "63_figure_inventory.csv")
MAN = os.path.join(ROOT, "supplementary_data", "_manifest.csv")
SUP = os.path.join(ROOT, "supplementary_data")
LOG = os.path.join(ROOT, "logs", "_s121_deliverables.txt")
L = []


def w(s=""):
    L.append(str(s))


TODAY = "2026-10-01"

# 本次实际重绘过的图件（s118 / s119 / s67 / s16）
REGEN = {
    "figures/Fig2_discovery_replication.png": "Fig2_discovery_replication",
    "figures/Fig3_coloc_sensitivity.png": "Fig3_coloc_sensitivity",
    "figures/Fig4_chr1p36_12_finemap.png": "Fig4_chr1p36_12_finemap",
    "figures/Fig5_opentargets_phewas_safety.png": "Fig5_opentargets_phewas_safety",
    "figures/FigS1_instrument_power_gating.png": "FigS1_instrument_power_gating",
    "figures/FigS2_sensitivity_outcome_concordance.png": "FigS2_sensitivity_outcome_concordance",
    "figures/FigS3_manhattan_by_celltype.png": "FigS3_manhattan_by_celltype",
    "figures/FigS4_discovery_volcano.png": "FigS4_discovery_volcano",
    "figures/FigD5_chr1_LD_structure.png": "FigD5_chr1_LD_structure",
}

PANEL_FIX = {"3 (a\u2013c)": "3 (A\u2013C)", "4 (a\u2013d)": "4 (A\u2013D)",
             "6 (a\u2013f)": "6 (A\u2013F)", "4 (a-d)": "4 (A-D)",
             "2 (a-b)": "2 (A-B)", "14 (a-n)": "14 (A-N)"}

# ---------------------------------------------------------------- ① 更新 inventory
with io.open(INV, encoding="utf-8-sig", newline="") as fh:
    rows = list(csv.DictReader(fh))
    fields = list(rows[0].keys())

from PIL import Image       # noqa: E402

w("=" * 100)
w("① tables/63_figure_inventory.csv 更新")
w("=" * 100)
changed = []
for r in rows:
    png = r.get("PNG文件", "")
    pdf = (r.get("PDF文件", "") or "").strip()
    if png in REGEN:
        name = REGEN[png]
        ppng = os.path.join(ROOT, png)
        ppdf = os.path.join(ROOT, pdf)
        if os.path.exists(ppng):
            r["PNG字节"] = str(os.path.getsize(ppng))
            im = Image.open(ppng)
            wmm = im.size[0] / 600 * 25.4
            hmm = im.size[1] / 600 * 25.4
            r["出图日期"] = TODAY
            r["状态"] = "已交付 %.1fx%.1fmm" % (wmm, hmm)
            changed.append((r["编号"], name, r["PNG字节"], "%.1fx%.1f" % (wmm, hmm)))
        if os.path.exists(ppdf):
            r["PDF字节"] = str(os.path.getsize(ppdf))
    old = (r.get("面板数", "") or "").strip()
    if old in PANEL_FIX:
        r["面板数"] = PANEL_FIX[old]
        for c in changed:
            pass

for num, name, b, sz in changed:
    w("   %-10s %-34s PNG %9s B    %s mm" % (num, name, b, sz))

with io.open(INV, "w", encoding="utf-8-sig", newline="") as fh:
    wr = csv.DictWriter(fh, fieldnames=fields, lineterminator="\n")
    wr.writeheader()
    wr.writerows(rows)
w("   已写回 %d 行" % len(rows))
w("")

# ---------------------------------------------------------------- ② 核验 manifest sha256
w("=" * 100)
w("② supplementary_data/_manifest.csv 完整性核验（只读；证明数据未被本次改动触碰）")
w("=" * 100)
with io.open(MAN, encoding="utf-8-sig", newline="") as fh:
    mrows = list(csv.DictReader(fh))


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


okn = missn = badn = 0
bad = []
for m in mrows:
    f = m["file"]
    p = os.path.join(SUP, f)
    if not os.path.exists(p):
        missn += 1
        bad.append((f, "文件不存在"))
        continue
    if os.path.getsize(p) != int(m["bytes"]):
        badn += 1
        bad.append((f, "字节不符 %d vs %s" % (os.path.getsize(p), m["bytes"])))
        continue
    if sha256(p) != m["sha256"]:
        badn += 1
        bad.append((f, "sha256 不符"))
        continue
    okn += 1
w("   核验通过 %d / %d ；缺失 %d ；不符 %d" % (okn, len(mrows), missn, badn))
for f, why in bad[:20]:
    w("     !! %s  %s" % (f, why))
w("   >>> 结论：%s" % ("全部源数据与补充表未变（本次仅改图件与图注）"
                       if (missn == 0 and badn == 0) else "存在变动，需人工复核"))

with io.open(LOG, "w", encoding="utf-8", newline="\n") as fh:
    fh.write("\n".join(L) + "\n")
print("done -> %s" % LOG)
