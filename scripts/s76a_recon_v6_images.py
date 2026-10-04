# -*- coding: utf-8 -*-
"""s76a · 侦察 v6 docx 的图片结构：media 文件 / rId / 每个 inline 的 extent"""
import io, os, re, zipfile, hashlib

ROOT = r"D:/endometriosis_project/11_sc_eqtl_mr_project"
V6 = os.path.join(ROOT, "44_论文初稿_带图_v6.docx")
FIGD = os.path.join(ROOT, "figures")
LOG = os.path.join(ROOT, "logs", "_s76a_recon_v6_images.txt")

out = []
def w(s=""):
    out.append(str(s)); print(s)

def md5b(b):
    return hashlib.md5(b).hexdigest()

w("v6 docx = %s  (%d B)" % (V6, os.path.getsize(V6)))
z = zipfile.ZipFile(V6)
names = z.namelist()

# figures/ 目录 PNG 的 md5 索引（用于反查 media 对应哪张图）
fig_md5 = {}
for f in os.listdir(FIGD):
    if f.lower().endswith(".png"):
        p = os.path.join(FIGD, f)
        fig_md5[md5b(open(p, "rb").read())] = f
w("figures/ 下 PNG 数 = %d" % len(fig_md5))

w("")
w("== word/media/ ==")
media = sorted([n for n in names if n.startswith("word/media/")])
media_info = {}
for m in media:
    b = z.read(m)
    h = md5b(b)
    media_info[m] = h
    w("  %-32s %9d B  md5 %s  -> %s" % (m, len(b), h, fig_md5.get(h, "?? 未匹配 figures/")))

w("")
w("== word/_rels/document.xml.rels ==")
rels = z.read("word/_rels/document.xml.rels").decode("utf-8")
rid2tgt = dict(re.findall(r'Id="([^"]+)"[^>]*Target="([^"]+)"', rels))
img_rid = {k: v for k, v in rid2tgt.items() if "image" in v}
for k, v in sorted(img_rid.items()):
    w("  %-8s -> %s" % (k, v))

doc = z.read("word/document.xml").decode("utf-8")
w("")
w("document.xml 长度 = %d 字符" % len(doc))
w("wp:inline 数 = %d ; wp:anchor 数 = %d"
  % (len(re.findall(r"<wp:inline\b", doc)), len(re.findall(r"<wp:anchor\b", doc))))

w("")
w("== 每个 inline 块 ==")
blocks = list(re.finditer(r"<wp:inline\b.*?</wp:inline>", doc, re.S))
for i, m in enumerate(blocks):
    seg = m.group(0)
    emb = re.search(r'r:embed="([^"]+)"', seg)
    rid = emb.group(1) if emb else "?"
    tgt = img_rid.get(rid, "?")
    h = media_info.get("word/" + tgt.lstrip("/"), "?")
    tag = fig_md5.get(h, "?")
    exts = re.findall(r"<(wp:extent|a:ext|a:chExt|a:off)\b([^>]*)>", seg)
    w("  [#%d] pos=%d..%d  r:embed=%s  target=%s  => %s" % (i, m.start(), m.end(), rid, tgt, tag))
    for name, attrs in exts:
        cx = re.search(r'cx="(\d+)"', attrs)
        cy = re.search(r'cy="(\d+)"', attrs)
        x = re.search(r'\bx="(-?\d+)"', attrs)
        y = re.search(r'\by="(-?\d+)"', attrs)
        w("        <%s> cx=%s cy=%s x=%s y=%s" % (name,
          cx.group(1) if cx else "-", cy.group(1) if cy else "-",
          x.group(1) if x else "-", y.group(1) if y else "-"))
    # cx/cy 全量
    allcxcy = re.findall(r'cx="(\d+)"\s+cy="(\d+)"', seg)
    w("        cx/cy 对 = %s" % (allcxcy,))

w("")
w("== media 对应结论（顺序） ==")
for i, m in enumerate(blocks):
    seg = m.group(0)
    emb = re.search(r'r:embed="([^"]+)"', seg)
    rid = emb.group(1) if emb else "?"
    tgt = img_rid.get(rid, "?")
    h = media_info.get("word/" + tgt.lstrip("/"), "?")
    w("  inline[%d] -> %s -> %s" % (i, tgt, fig_md5.get(h, "??")))

with io.open(LOG, "w", encoding="utf-8", newline="\n") as fh:
    fh.write("\n".join(out) + "\n")
