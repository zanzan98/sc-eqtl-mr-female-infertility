# -*- coding: utf-8 -*-
"""s63d_lnnum_diag.py —— 核查 Word 底本是否启用了行号（w:lnNumType），并 dump 相关页原始文本"""
import io, os, zipfile, re
import fitz
ROOT = r"D:/endometriosis_project/11_sc_eqtl_mr_project"
OUT = os.path.join(ROOT, "logs", "s63d_lnnum_diag.log")
buf = []
def p(*a):
    s = " ".join(str(x) for x in a); buf.append(s)
    try: print(s)
    except UnicodeEncodeError: print(s.encode("ascii", "backslashreplace").decode("ascii"))

for name in ("41_论文初稿_带图_v3.docx", "42_论文初稿_带图_v4.docx"):
    f = os.path.join(ROOT, name)
    p("== %s (%d B) ==" % (name, os.path.getsize(f)))
    z = zipfile.ZipFile(f)
    doc = z.read("word/document.xml").decode("utf-8")
    p("  document.xml len = %d" % len(doc))
    for tag in ("lnNumType", "pgNumType", "w:sectPr", "w:footerReference"):
        p("  count(%s) = %d" % (tag, doc.count(tag)))
    m = re.search(r"<w:lnNumType[^>]*/>", doc)
    p("  lnNumType element = %r" % (m.group(0) if m else None))
    # 该 I4 段的邻近片段
    i = doc.find("工具并集内最大两两")
    if i >= 0:
        seg = re.sub(r"<[^>]+>", "", doc[i - 200:i + 260])
        p("  XML_TEXT_AROUND = %r" % seg)
    else:
        p("  '工具并集内最大两两' 在 document.xml 中未找到（可能被 run 切分）")
        # 退一步：按 1.0000 找
        j = doc.find("1.0000")
        p("  find '1.0000' = %d" % j)

pdf = r"D:\_wb_docx_build\_exp_draft_v4.pdf"
d = fitz.open(pdf)
p("== PDF 原始页文本（18、19 页，逐行带行首行尾）==")
for pi in (17, 18):
    p("---- page %d ----" % (pi + 1))
    for ln in d[pi].get_text().splitlines():
        if ln.strip():
            p("   |%s|" % ln)
io.open(OUT, "w", encoding="utf-8", newline="\n").write("\n".join(buf) + "\n")
print("DIAG2_DONE")
