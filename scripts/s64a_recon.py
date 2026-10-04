# -*- coding: utf-8 -*-
"""s64a_recon.py —— 只读侦察
(A) 档务：枚举项目根「前缀号重复」的文件，并对每个文件统计全项目（含 memory）的**引用次数**
(B) 文本：定位局限段中「性别分层」句在 v4.docx 的段落与 run 结构（判断可否整段重写）
日志：logs/s64a_recon.log
"""
import io, os, re, sys
from collections import defaultdict

ROOT = r"D:/endometriosis_project/11_sc_eqtl_mr_project"
MEM = r"D:/endometriosis_project/.workbuddy/memory"
DOCX = os.path.join(ROOT, "42_论文初稿_带图_v4.docx")
LOG = os.path.join(ROOT, "logs", "s64a_recon.log")
buf = []
def p(*a):
    s = " ".join(str(x) for x in a); buf.append(s)
    try: print(s)
    except UnicodeEncodeError: print(s.encode("ascii", "backslashreplace").decode("ascii"))

# ---------------- (A) 档务 ----------------
p("== (A) 项目根前缀号重复枚举 ==")
files = sorted(os.listdir(ROOT))
pref = defaultdict(list)
for f in files:
    if not os.path.isfile(os.path.join(ROOT, f)):
        continue
    m = re.match(r"^(\d+[a-z]?)_", f)
    if m:
        pref[m.group(1)].append(f)
dup = {k: v for k, v in pref.items() if len(v) > 1}
p("重复前缀组数 = %d" % len(dup))
for k in sorted(dup, key=lambda x: (len(x), x)):
    p("  [%s] %d 份" % (k, len(dup[k])))
    for f in dup[k]:
        p("      %s" % f)

# 引用统计（扫项目内 *.md + memory/*.md，排除自身）
scan_dirs = [ROOT, MEM]
corpus = []
for d in scan_dirs:
    for dp, dn, fn in os.walk(d):
        if any(seg in dp for seg in ("00_data_raw", "02_data_processed", "_")):
            continue
        for f in fn:
            if f.endswith((".md", ".txt", ".py")):
                fp = os.path.join(dp, f)
                try:
                    corpus.append((fp, io.open(fp, encoding="utf-8", errors="replace").read()))
                except Exception:
                    pass
p("")
p("== 重复前缀文件的被引用次数（扫 %d 个文本文件）==" % len(corpus))
for k in sorted(dup, key=lambda x: (len(x), x)):
    for f in dup[k]:
        n = 0
        who = []
        for fp, txt in corpus:
            if os.path.basename(fp) == f:
                continue
            c = txt.count(f)
            if c:
                n += c; who.append(os.path.relpath(fp, r"D:/endometriosis_project"))
        p("  %-52s ref=%d %s" % (f, n, ("| " + ", ".join(who[:4])) if who else ""))

# ---------------- (B) 局限段 run 结构 ----------------
p("")
p("== (B) v4.docx 中「性别分层」句所在段落的 run 结构 ==")
from docx import Document
from docx.oxml.ns import qn
doc = Document(DOCX)
hits = [x for x in doc.paragraphs if "性别分层" in x.text]
p("含「性别分层」的段落数 = %d" % len(hits))
for hi, para in enumerate(hits):
    p("---- 段落 %d ----" % hi)
    p("  文本 = %r" % para.text)
    runs = para._p.findall(qn("w:r"))
    p("  run 数 = %d" % len(runs))
    for ri, r in enumerate(runs):
        rpr = r.find(qn("w:rPr"))
        flags = []
        if rpr is not None:
            for tag, name in (("w:b", "B"), ("w:i", "I"), ("w:u", "U"),
                              ("w:color", "COLOR"), ("w:vertAlign", "VA"),
                              ("w:rFonts", "FONT")):
                e = rpr.find(qn(tag))
                if e is not None:
                    flags.append(name + (":" + (e.get(qn("w:val")) or e.get(qn("w:ascii")) or "1")))
        t = "".join(x.text or "" for x in r.findall(qn("w:t")))
        p("    [%02d] %-28s %r" % (ri, ",".join(flags) or "-", t))
    # 段落级 pPr 摘要
    ppr = para._p.find(qn("w:pPr"))
    if ppr is not None:
        p("  pPr = %s" % re.sub(r"\s+", " ", ppr.xml)[:300])
io.open(LOG, "w", encoding="utf-8", newline="\n").write("\n".join(buf) + "\n")
print("RECON_DONE")
