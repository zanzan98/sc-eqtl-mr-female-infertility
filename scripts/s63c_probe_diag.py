# -*- coding: utf-8 -*-
"""s63c_probe_diag.py —— 定位 verify 中唯一 miss 的探针在 PDF 文本层的实际形态"""
import fitz, io, os
PDF = r"D:\_wb_docx_build\_exp_draft_v4.pdf"
OUT = r"D:\endometriosis_project\11_sc_eqtl_mr_project\logs\s63c_probe_diag.log"
buf = []
def p(*a):
    s = " ".join(str(x) for x in a); buf.append(s)
    try: print(s)
    except UnicodeEncodeError: print(s.encode("ascii", "backslashreplace").decode("ascii"))

doc = fitz.open(PDF)
raw = "".join(pg.get_text() for pg in doc)
flat = raw.replace(" ", "").replace("\n", "")
key = "工具并集内最大两两"
i = flat.find(key)
p("find '%s' -> %d" % (key, i))
if i >= 0:
    seg = flat[i - 20:i + 80]
    p("SEGMENT = %r" % seg)
    p("CODEPOINTS = %s" % " ".join("U+%04X(%s)" % (ord(c), c) for c in seg))
# 变体探测
for v in ["r²达1.0000", "r²达", "1.0000", "r2达", "r² 达 1.0000"]:
    p("contains %-14s = %s" % (v, v.replace(" ", "") in flat))
# 也看看 I4 段落所在页
for pi, pg in enumerate(doc):
    t = pg.get_text()
    if "MVMR" in t or "多变量孟德尔随机化" in t:
        p("PAGE %d 命中 MVMR/多变量" % (pi + 1))
        f = t.replace(" ", "").replace("\n", "")
        j = f.find("工具并集")
        if j >= 0:
            p("   PAGE_SEG = %r" % f[max(0, j - 10):j + 70])
io.open(OUT, "w", encoding="utf-8", newline="\n").write("\n".join(buf) + "\n")
print("DIAG_DONE")
