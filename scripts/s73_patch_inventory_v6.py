# -*- coding: utf-8 -*-
"""s73_patch_inventory_v6.py —— 给 tables/63_figure_inventory.csv 增列「在v6预览PDF页」

依据：`44_论文初稿_带图_v6_预览.pdf`（24 页）的页内定位实测（verify_pdf_draft_v6 输出）。
  正文 5 图 = 4 / 8 / 10 / 13 / 15（与 v5 相同，未移动）
  图 S1–S4 = 22 / 23 / 23 / 24（本轮新入 docx 附录 A）
  D 系列 = `-`（不上传、不引用）

硬纪律：BOM 保留、LF 保留、写后 getsize + 计数复核。
"""
import io
import os
import shutil
import hashlib

ROOT = r"D:/endometriosis_project/11_sc_eqtl_mr_project"
CSV = os.path.join(ROOT, "tables", "63_figure_inventory.csv")
BACKUP_DIR = r"D:/_transfer_logs/_ms_backup_20260929_s1s4"
LOG = os.path.join(ROOT, "logs", "_s73_patch_inventory_v6.log")

V6PAGE = {
    "Fig1": "4", "Fig2": "8", "Fig3": "10", "Fig4": "13", "Fig5": "15",
    "图 S1": "22", "图 S2": "23", "图 S3": "23", "图 S4": "24",
    "FigD1": "-", "FigD2": "-", "FigD3": "-", "FigD4": "-", "FigD5": "-",
}

L = []
def p(*a):
    s = " ".join(str(x) for x in a); L.append(s)


raw = open(CSV, "rb").read()
p("bytes before = %d  BOM=%s  LF=%d  CRLF=%d" % (
    len(raw), raw.startswith(b"\xef\xbb\xbf"), raw.count(b"\n"), raw.count(b"\r\n")))
assert raw.count(b"\r\n") == 0, "该 CSV 现为 CRLF，脚本假设 LF"

text = raw.decode("utf-8-sig")
lines = [ln for ln in text.split("\n")]
p("lines = %d" % len(lines))
hdr = lines[0].split(",")
p("header before = %r" % hdr)
assert "在v6预览PDF页" not in lines[0], "已存在该列（幂等保护）"
assert hdr[-1] == "状态", "末列不是 状态，实际 %r" % hdr[-1]

out = [lines[0] + ",在v6预览PDF页"]
n = 0
for ln in lines[1:]:
    if not ln.strip():
        out.append(ln)
        continue
    cells = ln.split(",")
    tag = cells[1]
    assert tag in V6PAGE, "未在映射表内的编号：%r" % tag
    out.append(ln + "," + V6PAGE[tag])
    n += 1
    p("  %-7s -> v6 页 %s" % (tag, V6PAGE[tag]))
assert n == len(V6PAGE), "改写行数 %d != %d" % (n, len(V6PAGE))

new_text = "\n".join(out)
new_raw = new_text.encode("utf-8")
if raw.startswith(b"\xef\xbb\xbf"):
    new_raw = b"\xef\xbb\xbf" + new_raw

os.makedirs(BACKUP_DIR, exist_ok=True)
shutil.copy2(CSV, os.path.join(BACKUP_DIR, "63_figure_inventory.csv.bak_v6"))
with io.open(CSV, "wb") as fh:
    fh.write(new_raw)

chk = open(CSV, "rb").read()
p("bytes after = %d (delta %+d)" % (len(chk), len(chk) - len(raw)))
p("BOM preserved = %s ; CRLF = %d" % (chk.startswith(b"\xef\xbb\xbf"), chk.count(b"\r\n")))
t2 = chk.decode("utf-8-sig")
p("header after = %r" % t2.split("\n")[0])
for tag, pg in V6PAGE.items():
    hit = [ln for ln in t2.split("\n") if ln.startswith("") and (","+tag+",") in ln + ","]
    assert any(ln.endswith("," + pg) for ln in hit), "行 %s 未追加上 %s" % (tag, pg)
p("14 行全部追加上 v6 页码 ✅")
p("VERDICT = PASS")
with io.open(LOG, "w", encoding="utf-8", newline="\n") as fh:
    fh.write("\n".join(L) + "\n")
print("\n".join(L))
