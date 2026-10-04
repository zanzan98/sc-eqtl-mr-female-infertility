# -*- coding: utf-8 -*-
"""s72a_recon_body.py —— 侦察 v4 docx body 的 XML 子元素序列（含空段落）与关键段落的 pPr/style

为 v6 构建提供：①插图段与图注段的确切 pPr / style；②块间空段落约定；③插入锚点的 XML 位置。
输出 logs/_s72a_recon_body.txt
"""
import io, os, re
from docx import Document
from docx.oxml.ns import qn

ROOT = r"D:/endometriosis_project/11_sc_eqtl_mr_project"
DOCX = os.path.join(ROOT, "42_论文初稿_带图_v4.docx")
LOG = os.path.join(ROOT, "logs", "_s72a_recon_body.txt")
NS_WP = "http://schemas.openxmlformats.org/drawingml/2006/wordprocessingDrawing"

L = []
def p(*a):
    L.append(" ".join(str(x) for x in a))

doc = Document(DOCX)
body = doc.element.body
kids = list(body)

p("body 子元素数 = %d" % len(kids))
p("")
p("========== body 序列（idx | tag | style | 有无图 | 文本前 60 字） ==========")
para_pos = []       # (body_idx, Paragraph)
for i, k in enumerate(kids):
    tag = k.tag.split("}")[-1]
    if tag != "p":
        p("[%03d] <%s>" % (i, tag)); continue
    from docx.text.paragraph import Paragraph
    par = Paragraph(k, doc)
    t = par.text.strip()
    has_img = bool(k.findall(".//{%s}inline" % NS_WP))
    st = par.style.name if par.style is not None else None
    p("[%03d] p style=%-12s img=%s | %s" % (i, st, "Y" if has_img else "-", t[:60]))
    para_pos.append((i, par))

p("")
p("========== 关键段落 XML（pPr 摘录） ==========")


def dump(par, label):
    p("--- %s | text=%r" % (label, par.text.strip()[:50]))
    pPr = par._p.find(qn("w:pPr"))
    p("    pPr = %s" % (pPr.xml.replace("\n", " ") if pPr is not None else None))
    runs = par._p.findall(qn("w:r"))
    p("    n_runs = %d" % len(runs))
    if runs:
        rPr = runs[0].find(qn("w:rPr"))
        p("    run0.rPr = %s" % (rPr.xml.replace("\n", " ") if rPr is not None else None))
    p("")


allp = [Paragraph(k, doc) for k in kids if k.tag.split("}")[-1] == "p"]

def first_where(fn, label):
    for par in allp:
        if fn(par):
            dump(par, label)
            return par
    p("--- %s NOT FOUND" % label); p("")
    return None

first_where(lambda x: x._p.findall(".//{%s}inline" % NS_WP), "插图段（第 1 张图）")
first_where(lambda x: x.text.strip().startswith("图 1 |"), "图注段 图 1")
first_where(lambda x: x.text.strip() == "摘要", "标题段 摘要")
first_where(lambda x: x.text.strip() == "附录 A 补充信息", "标题段 附录 A 补充信息")
first_where(lambda x: x.text.strip() == "参考文献", "标题段 参考文献")
first_where(lambda x: "Supplementary_Data.zip" in x.text and "补充数据见" in x.text, "附录 A 指针句")
first_where(lambda x: x.text.strip() == "4 讨论", "标题段 4 讨论")
first_where(lambda x: x.text.strip().startswith("图 5 |"), "图注段 图 5")

p("========== 5 张图的 inline 结构（第 1 张详列） ==========")
n = 0
for par in allp:
    inls = par._p.findall(".//{%s}inline" % NS_WP)
    if not inls:
        continue
    n += 1
    if n == 1:
        p(inls[0].xml[:3000])
        p("")
p("图片段总数 = %d" % n)

os.makedirs(os.path.dirname(LOG), exist_ok=True)
with io.open(LOG, "w", encoding="utf-8", newline="\n") as fh:
    fh.write("\n".join(L) + "\n")
print("recon -> %s (%d lines)" % (LOG, len(L)))
