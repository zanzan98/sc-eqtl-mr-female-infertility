# -*- coding: utf-8 -*-
"""
s74 · 修正 tables/63_figure_inventory.csv 的列序
问题：s73 用「行末追加」得到 ...,在v5预览PDF页,状态,在v6预览PDF页
目标：把「在v6预览PDF页」移到「状态」之前 => ...,在v5预览PDF页,在v6预览PDF页,状态
纪律：with io.open 写盘 + 写完 getsize + 回读复核；BOM / 换行风格保持
"""
import io, os, csv, shutil, hashlib

ROOT = r"D:/endometriosis_project/11_sc_eqtl_mr_project"
CSV = os.path.join(ROOT, "tables", "63_figure_inventory.csv")
BAK = r"D:/_transfer_logs/_ms_backup_20260929_s1s4/63_figure_inventory.csv.bak2_reorder"
LOG = os.path.join(ROOT, "logs", "_s74_reorder_inventory.log")

_lines = []
def p(s=""):
    _lines.append(str(s))
    with io.open(LOG, "a", encoding="utf-8", newline="\n") as fh:
        fh.write(str(s) + "\n")

# 清空日志
with io.open(LOG, "w", encoding="utf-8", newline="\n") as fh:
    fh.write("")

p("==== s74 修正 inventory 列序 ====")
raw = open(CSV, "rb").read()
bom = raw.startswith(b"\xef\xbb\xbf")
body = raw[3:] if bom else raw
crlf = body.count(b"\r\n")
lf = body.count(b"\n") - crlf
nl = "\r\n" if crlf > 0 and crlf == body.count(b"\r\n") else "\n"
p("bytes before = %d  BOM=%s  LF=%d  CRLF=%d  nl=%r" % (len(raw), bom, lf, crlf, nl))
p("md5 before = %s" % hashlib.md5(raw).hexdigest())

text = body.decode("utf-8")
rows = list(csv.reader(io.StringIO(text)))
rows = [r for r in rows if r]          # 去空行
hdr = rows[0]
p("header before = %s" % (hdr,))

COL_NEW = "在v6预览PDF页"
COL_STATE = "状态"
assert COL_NEW in hdr, "缺新列"
assert COL_STATE in hdr, "缺状态列"
assert hdr[-1] == COL_NEW, "新列不在行末，无需重排"

i_new = hdr.index(COL_NEW)
i_state = hdr.index(COL_STATE)
p("index: 状态=%d  在v6预览PDF页=%d" % (i_state, i_new))

# 重排 header：先移除末列，再插到 状态 之前
new_hdr = list(hdr)
new_hdr.pop(i_new)
i_state2 = new_hdr.index(COL_STATE)
new_hdr.insert(i_state2, COL_NEW)
p("header after  = %s" % (new_hdr,))

out = [new_hdr]
for r in rows[1:]:
    if len(r) != len(hdr):
        p("!! 列数异常行（跳过重排、原样保留）: %s" % (r[:3],))
        out.append(r)
        continue
    v = r[i_new]
    rr = list(r)
    rr.pop(i_new)
    rr.insert(i_state2, v)      # i_state2 = pop 后「状态」列在新表中的位置
    out.append(rr)
    p("  %-8s -> v6 页 %s" % (r[1], v))

# 写回
buf = io.StringIO()
w = csv.writer(buf, lineterminator=nl)
w.writerows(out)
s = buf.getvalue()
data = (b"\xef\xbb\xbf" if bom else b"") + s.encode("utf-8")
with io.open(BAK, "wb") as fh:
    fh.write(raw)
with io.open(CSV, "wb") as fh:
    fh.write(data)

sz = os.path.getsize(CSV)
p("bytes after  = %d (delta %+d)" % (sz, sz - len(raw)))
p("md5 after    = %s" % hashlib.md5(open(CSV, "rb").read()).hexdigest())

# ---- 回读复核 ----
raw2 = open(CSV, "rb").read()
b2 = raw2.startswith(b"\xef\xbb\xbf")
rows2 = [r for r in csv.reader(io.StringIO(raw2[3:].decode("utf-8") if b2 else raw2.decode("utf-8"))) if r]
h2 = rows2[0]
p("recheck: BOM=%s  rows=%d" % (b2, len(rows2)))
p("recheck header = %s" % (h2,))
ok = True
ok &= (h2[-2:] == [COL_NEW, COL_STATE])
ok &= (h2[9] == "在v5预览PDF 页" or h2[9].startswith("在v5预览PDF"))
# 数据搬运正确性：每一行的 v6 页值应与重排前同值
old_rows = {r[1]: r[i_new] for r in rows[1:] if len(r) == len(hdr)}
new_rows = {r[1]: r[-2] for r in rows2[1:] if len(r) == len(h2)}
p("values match = %s" % (old_rows == new_rows))
ok &= (old_rows == new_rows)
p("行数一致 = %s" % (len(rows2) == len(rows)))
ok &= (len(rows2) == len(rows))

p("")
p("VERDICT = %s" % ("PASS" if ok else "FAIL"))
print("\n".join(_lines))
