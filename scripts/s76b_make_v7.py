# -*- coding: utf-8 -*-
"""s76b · 构建 45_论文初稿_带图_v7.docx

做法（ZIP 级外科手术，最可控）：
  ① media/image2..5.png 换为无脚注版新 PNG（image1=Fig1 与 image6..9=S1..S4 不动）
  ② document.xml 中 inline[1..4] 的 cy 由旧值改为按新宽高比计算的值
     （cx 恒 5760000 = 16.0 cm 不变）
  ③ 每个旧 cy 在 document.xml 中应恰好出现 2 次（wp:extent + a:ext）=> 全局替换安全
  ④ 除上述 8 处数字外，document.xml 逐字不变

自证：先用【备份的旧 PNG】复现 v6 中的旧 cy，确认算法无误后再算新 cy。
"""
import io, os, re, zipfile, hashlib

ROOT = r"D:/endometriosis_project/11_sc_eqtl_mr_project"
FIGD = os.path.join(ROOT, "figures")
BAK = r"D:/_transfer_logs/_figs_backup_20260929_nonotes"
V6 = os.path.join(ROOT, "44_论文初稿_带图_v6.docx")
V7 = os.path.join(ROOT, "45_论文初稿_带图_v7.docx")
LOG = os.path.join(ROOT, "logs", "_s76b_make_v7.txt")

out = []
def w(s=""):
    out.append(str(s)); print(s)

from docx.image.image import Image

def emu_size(p):
    im = Image.from_file(p)
    return im.width, im.height      # python-docx: 已是 EMU

CX = 5760000                        # 16.0 cm
TARGETS = [
    (1, "word/media/image2.png", "Fig2_discovery_replication.png", 5268537),
    (2, "word/media/image3.png", "Fig3_coloc_sensitivity.png", 7021880),
    (3, "word/media/image4.png", "Fig4_chr1p36_12_finemap.png", 8314354),
    (4, "word/media/image5.png", "Fig5_opentargets_phewas_safety.png", 7588040),
]

w("==== v7 构建 ====")
w("v6 = %d B  md5 %s" % (os.path.getsize(V6), hashlib.md5(open(V6, "rb").read()).hexdigest()))

w("")
w("== ① 自证 cy 算法：用备份的旧 PNG 复现 v6 中的 cy ==")
w("   容差 <= 2 EMU（= 5.6e-7 cm，纯浮点末位舍入；Fig4 实测差 1 EMU）")
TOL = 2
ALGO_OK = True
for idx, media, png, v6cy in TARGETS:
    ow, oh = emu_size(os.path.join(BAK, png))
    calc = int(round(CX * oh / ow))
    d = calc - v6cy
    m = abs(d) <= TOL
    ALGO_OK &= m
    w("  inline[%d] %-34s EMU=(%d,%d)  calc=%d  v6=%d  Δ=%+d  ok=%s"
      % (idx, png, ow, oh, calc, v6cy, d, m))
w("  算法自证 = %s" % ("PASS" if ALGO_OK else "FAIL"))
assert ALGO_OK, "cy 算法自证失败，停止（禁止带着不确定的算法改稿）"

w("")
w("== ② 计算新 cy（新 PNG 的 EMU 宽高比） ==")
NEW = {}
for idx, media, png, v6cy in TARGETS:
    nw, nh = emu_size(os.path.join(FIGD, png))
    ncy = int(round(CX * nh / nw))
    NEW[idx] = (media, png, v6cy, ncy, nw, nh)
    w("  inline[%d] %-34s 新EMU=(%d,%d)  cy %d -> %d  高度 %.3f -> %.3f cm"
      % (idx, png, nw, nh, v6cy, ncy, v6cy / 360000.0, ncy / 360000.0))

zin = zipfile.ZipFile(V6)
members = zin.infolist()
doc = zin.read("word/document.xml").decode("utf-8")
w("")
w("== ③ document.xml 中旧 cy 的唯一性检查 ==")
for idx, media, png, v6cy in TARGETS:
    c = doc.count('cy="%d"' % v6cy)
    w("  cy=\"%d\" 出现 %d 次 (应为 2)" % (v6cy, c))
    assert c == 2, "旧 cy 不唯一，停止"

doc2 = doc
for idx, media, png, v6cy in TARGETS:
    ncy = NEW[idx][3]
    doc2 = doc2.replace('cy="%d"' % v6cy, 'cy="%d"' % ncy)

# 除 cy 数字外逐字不变（把所有 cy="数字" 归一化后比较）
norm_old = re.sub(r'cy="\d+"', 'cy="#"', doc)
norm_new = re.sub(r'cy="\d+"', 'cy="#"', doc2)
w("  归一化后 doc 相同 = %s" % (norm_old == norm_new))
assert norm_old == norm_new, "除 cy 外还有改动，停止"

w("")
w("== ④ 重打包 ==")
new_media = {}
for idx, media, png, v6cy in TARGETS:
    b = open(os.path.join(FIGD, png), "rb").read()
    new_media[media] = b

zout = zipfile.ZipFile(V7, "w", zipfile.ZIP_DEFLATED)
for zi in members:
    if zi.filename == "word/document.xml":
        data = doc2.encode("utf-8")
    elif zi.filename in new_media:
        data = new_media[zi.filename]
    else:
        data = zin.read(zi.filename)
    zout.writestr(zi, data)
zout.close()
zin.close()

w("")
w("== ⑤ 复核 ==")
w("v7 = %d B  md5 %s" % (os.path.getsize(V7), hashlib.md5(open(V7, "rb").read()).hexdigest()))

z2 = zipfile.ZipFile(V7)
fig_md5 = {}
for f in os.listdir(FIGD):
    if f.lower().endswith(".png"):
        p = os.path.join(FIGD, f)
        fig_md5[hashlib.md5(open(p, "rb").read()).hexdigest()] = f

media_names = sorted([n for n in z2.namelist() if n.startswith("word/media/")])
w("  media 数 = %d" % len(media_names))
for m in media_names:
    b = z2.read(m)
    h = hashlib.md5(b).hexdigest()
    w("    %-24s %9d B  md5 %s  -> %s" % (m, len(b), h, fig_md5.get(h, "(旧图/未匹配)")))

# 逐字节确认：未在替换清单中的 media 与 v6 一致
zin2 = zipfile.ZipFile(V6)
w("")
w("  -- 与 v6 的 media 逐字节比对 --")
for m in media_names:
    a = zin2.read(m)
    b = z2.read(m)
    w("    %-24s 相同=%s" % (m, a == b))

w("")
w("  -- document.xml 差异行数（应只涉及 cy 数字）--")
d1 = zin2.read("word/document.xml").decode("utf-8")
d2 = z2.read("word/document.xml").decode("utf-8")
w("    长度 %d -> %d" % (len(d1), len(d2)))
i = 0
for idx, media, png, v6cy in TARGETS:
    w("    %s: cy %d -> %d" % (media, v6cy, NEW[idx][3]))

z2.close(); zin2.close()

with io.open(LOG, "w", encoding="utf-8", newline="\n") as fh:
    fh.write("\n".join(out) + "\n")
w("")
w("VERDICT = %s" % ("PASS" if os.path.exists(V7) else "FAIL"))
