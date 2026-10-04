# -*- coding: utf-8 -*-
"""s71_recon_v4.py —— v4 docx 结构侦察（为 v6 构建定位锚点）

输出 logs/_s71_recon_v4.txt：
  · 全部非空段落 idx + 前 100 字
  · 图片 inline 所在段 idx + extent（cm）+ 图片关系 ID
  · 关键锚点定位：1000 Genomes / 4 讨论 / 附录 A / 性别分层
  · 段落内 run 数分布（判断能否做 run 级外科）
"""
import io, os, re
from docx import Document
from docx.oxml.ns import qn

ROOT = r"D:/endometriosis_project/11_sc_eqtl_mr_project"
DOCX = os.path.join(ROOT, "42_论文初稿_带图_v4.docx")
LOG = os.path.join(ROOT, "logs", "_s71_recon_v4.txt")
NS_WP = "http://schemas.openxmlformats.org/drawingml/2006/wordprocessingDrawing"
NS_A = "http://schemas.openxmlformats.org/drawingml/2006/main"
NS_R = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"

L = []
def p(*a):
    s = " ".join(str(x) for x in a); L.append(s)

doc = Document(DOCX)
paras = doc.paragraphs
np = [x for x in paras if x.text.strip()]
p("底本 %s" % DOCX)
p("总段落 = %d ；非空段落 = %d" % (len(paras), len(np)))
p("")
p("========== 全部非空段落（idx = 非空序号） ==========")
for i, x in enumerate(np):
    if i < 10 or True:
        p("[%02d] runs=%2d | %s" % (i, len(x._p.findall(qn("w:r"))), x.text.strip()[:110]))
p("")
p("========== 图片 inline ==========")
for i, x in enumerate(np):
    for inl in x._p.findall(".//{%s}inline" % NS_WP):
        e = inl.find("{%s}extent" % NS_WP)
        docPr = inl.find("{%s}docPr" % NS_WP)
        blip = inl.find(".//{%s}blip" % NS_A)
        rid = blip.get("{%s}embed" % NS_R) if blip is not None else None
        p("[para %02d] extend=%.2f x %.2f cm  docPr.name=%r  rid=%s  | %s" % (
            i, int(e.get("cx"))/360000, int(e.get("cy"))/360000,
            docPr.get("name") if docPr is not None else None, rid, x.text.strip()[:50]))
p("")
p("========== 关键锚点 ==========")
def find(kw, label):
    hits = [(i, x.text.strip()) for i, x in enumerate(np) if kw in x.text]
    p("--- %-22s kw=%r  hits=%d" % (label, kw, len(hits)))
    for i, t in hits[:6]:
        # 上下文
        s = t.index(kw) if kw in t else 0
        p("    [%02d] …%s…" % (i, t[max(0, s-90):s+120]))
for kw, lab in [("1000 Genomes", "Fig4a 图注"),
                ("OneK1K", "OneK1K"),
                ("附录 A", "附录 A"),
                ("补充信息", "补充信息"),
                ("性别分层", "待删句"),
                ("Supplementary", "Supplementary"),
                ("补充图", "补充图"),
                ("讨论", "讨论标题"),
                ("质量控制与分层视图", "引用句"),
                ("图 1", "图注1"),
                ("Fig. 4", "Fig.4 图注"),
                ("Fig. 5", "Fig.5 图注"),
                ("参考文献", "参考文献"),
                ("声明", "声明")]:
    find(kw, lab)
p("")
p("========== 末尾 12 段 ==========")
for i, x in enumerate(np[-12:], start=len(np)-12):
    p("[%02d] %s" % (i, x.text.strip()[:150]))
p("")
p("========== 图注段落定位（以「图 4」等开头的段） ==========")
for i, x in enumerate(np):
    t = x.text.strip()
    if re.match(r"^图\s*[1-5]\s*\|", t) or re.match(r"^图\s*S", t) or t.startswith("附录"):
        p("[%02d] %s" % (i, t[:150]))

os.makedirs(os.path.dirname(LOG), exist_ok=True)
with io.open(LOG, "w", encoding="utf-8", newline="\n") as fh:
    fh.write("\n".join(L) + "\n")
print("recon -> %s (%d lines)" % (LOG, len(L)))
