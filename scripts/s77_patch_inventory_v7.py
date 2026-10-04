# -*- coding: utf-8 -*-
"""s77 · 更新 tables/63_figure_inventory.csv 到 v7

内容：
  ① Fig2–Fig5 的 PNG字节 / PDF字节 / 出图日期 更新（重出后）
  ② Fig2–Fig5 的「状态」补图件实测尺寸（与补充图体例一致）
  ③ 新增列「在v7预览PDF页」，插在「状态」之前（与 v5/v6 两列相邻）
纪律：BOM / LF 保持；写后回读复核；备份
"""
import io, os, csv, hashlib

ROOT = r"D:/endometriosis_project/11_sc_eqtl_mr_project"
CSV = os.path.join(ROOT, "tables", "63_figure_inventory.csv")
BAK = r"D:/_transfer_logs/_figs_backup_20260929_nonotes/63_figure_inventory.csv.bak3_v7"
LOG = os.path.join(ROOT, "logs", "_s77_patch_inventory_v7.txt")

out = []
def w(s=""):
    out.append(str(s)); print(s)

raw = open(CSV, "rb").read()
bom = raw.startswith(b"\xef\xbb\xbf")
text = raw[3:].decode("utf-8") if bom else raw.decode("utf-8")
crlf = text.count("\r\n")
lf = text.count("\n") - crlf
nl = "\r\n" if crlf > 0 and lf == 0 else "\n"
w("before bytes=%d  BOM=%s  CRLF=%d LF=%d  nl=%r" % (len(raw), bom, crlf, lf, nl))
w("before md5=%s" % hashlib.md5(raw).hexdigest())

rows = [r for r in csv.reader(io.StringIO(text)) if r]
hdr = rows[0]
w("header before = %s" % (hdr,))
assert "在v7预览PDF页" not in hdr, "已含 v7 列，幂等中止"

# ① 数据更新
UPD = {
    "Fig2": dict(png=623880,  pdf=167933, size="182.8x145.9mm"),
    "Fig3": dict(png=775725,  pdf=57734,  size="182.5x196.4mm"),
    "Fig4": dict(png=882327,  pdf=85276,  size="182.3x242.9mm"),
    "Fig5": dict(png=905835,  pdf=62186,  size="181.9x202.7mm"),
}
V7PAGE = {"Fig1": "4", "Fig2": "8", "Fig3": "10", "Fig4": "13", "Fig5": "15",
          "图 S1": "22", "图 S2": "23", "图 S3": "23", "图 S4": "24",
          "FigD1": "-", "FigD2": "-", "FigD3": "-", "FigD4": "-", "FigD5": "-"}

i_png = hdr.index("PNG字节")
i_pdf = hdr.index("PDF字节")
i_date = hdr.index("出图日期")
i_stat = hdr.index("状态")
i_num = hdr.index("编号")

NEW_COL = "在v7预览PDF页"
new_hdr = list(hdr)
new_hdr.insert(new_hdr.index("状态"), NEW_COL)
w("header after  = %s" % (new_hdr,))

outrows = [new_hdr]
changed = []
for r in rows[1:]:
    if len(r) != len(hdr):
        w("!! 列数异常行，原样保留: %s" % (r[:3],)); outrows.append(r); continue
    num = r[i_num]
    rr = list(r)
    if num in UPD:
        u = UPD[num]
        changed.append("%s: png %s->%d  pdf %s->%d  date %s->2026-09-29  status '%s'->'已交付 %s'"
                       % (num, r[i_png], u["png"], r[i_pdf], u["pdf"], r[i_date], r[i_stat], u["size"]))
        rr[i_png] = str(u["png"]); rr[i_pdf] = str(u["pdf"])
        rr[i_date] = "2026-09-29"; rr[i_stat] = "已交付 %s" % u["size"]
    rr.insert(i_stat, V7PAGE.get(num, "-"))
    outrows.append(rr)

w("")
w("-- 变更明细 --")
for c in changed:
    w("  " + c)
w("  新增列 %s：%d 行全部赋值" % (NEW_COL, len(outrows) - 1))

buf = io.StringIO()
csv.writer(buf, lineterminator=nl).writerows(outrows)
data = (b"\xef\xbb\xbf" if bom else b"") + buf.getvalue().encode("utf-8")
with io.open(BAK, "wb") as fh:
    fh.write(raw)
with io.open(CSV, "wb") as fh:
    fh.write(data)

sz = os.path.getsize(CSV)
w("")
w("after bytes=%d (Δ%+d)  md5=%s" % (sz, sz - len(raw), hashlib.md5(open(CSV, "rb").read()).hexdigest()))

# 回读
raw2 = open(CSV, "rb").read()
b2 = raw2.startswith(b"\xef\xbb\xbf")
rows2 = [r for r in csv.reader(io.StringIO(raw2[3:].decode("utf-8") if b2 else raw2.decode("utf-8"))) if r]
h2 = rows2[0]
w("recheck BOM=%s rows=%d" % (b2, len(rows2)))
w("recheck header = %s" % (h2,))
tab = {r[h2.index("编号")]: r for r in rows2[1:]}
ok = True
ok &= (h2[-1] == "状态") and (h2[-2] == NEW_COL)
for k, v in V7PAGE.items():
    got = tab[k][h2.index(NEW_COL)]
    ok &= (got == v)
    w("  %-8s v7 页 = %-3s (期望 %s)" % (k, got, v))
for k, u in UPD.items():
    ok &= (tab[k][h2.index("PNG字节")] == str(u["png"]))
    ok &= (tab[k][h2.index("状态")] == "已交付 %s" % u["size"])
w("")
w("VERDICT = %s" % ("PASS" if ok else "FAIL"))

with io.open(LOG, "w", encoding="utf-8", newline="\n") as fh:
    fh.write("\n".join(out) + "\n")
