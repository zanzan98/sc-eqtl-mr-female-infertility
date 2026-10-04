# -*- coding: utf-8 -*-
"""s109_make_v8_fig1v2.py —— 中文版论文 v8：把 Fig1 换成【V2 图形化版】

底本：`45_论文初稿_带图_v7.docx`（项目根；md5 记录在下）
产物：`46_中文版论文_带图_v8.docx`

做法（ZIP 级外科手术，与 s76b_make_v7.py 同一套路数）：
  ① `word/media/image1.png`  <- `D:/_transfer_logs/_fig1_draft_v2/Fig1v2_study_design_graphical.png`
  ② `document.xml` 中 **第 1 个内联**（Fig1）的 `cy` 按新 PNG 宽高比重算；
     `cx` 恒 5760000（= 16.0 cm）不动。其余 8 个内联尺寸一字不动。
  ③ 除该 cy 数字外，`document.xml` **逐字不变**；其余 8 个 media **逐字节不变**。

★ 自证顺序（不通过就停，禁止带着不确定的算法改稿）：
  A. 确认 v7 的 `image1.png` md5 == `figures/Fig1_study_design.png`（即「旧 Fig1」）；
  B. 用**旧** Fig1 PNG 的像素宽高复现 v7 中的 `cy=3612702`（容差 ≤2 EMU）；
  C. 再用**新** V2 PNG 算出新 cy；断言新旧不相等（确实需要改）。
"""
import hashlib
import io
import os
import re
import zipfile

ROOT = r"D:/endometriosis_project/11_sc_eqtl_mr_project"
FIGD = os.path.join(ROOT, "figures")
OLD_FIG1 = os.path.join(FIGD, "Fig1_study_design.png")
NEW_FIG1 = r"D:/_transfer_logs/_fig1_draft_v2/Fig1v2_study_design_graphical.png"
V7 = os.path.join(ROOT, "45_论文初稿_带图_v7.docx")
V8 = os.path.join(ROOT, "46_中文版论文_带图_v8.docx")
LOG = os.path.join(ROOT, "logs", "_s109_make_v8.txt")

CX = 5760000                      # 16.0 cm，全部 9 个内联共用
MEDIA1 = "word/media/image1.png"
TOL = 2                           # EMU

out = []


def w(s=""):
    out.append(str(s))
    print(s)


def md5(b):
    return hashlib.md5(b).hexdigest()


def png_px(p):
    """读 PNG IHDR 得 (w, h) 像素 —— 不依赖任何第三方库。"""
    with open(p, "rb") as fh:
        head = fh.read(33)
    assert head[:8] == b"\x89PNG\r\n\x1a\n", "不是 PNG: %s" % p
    assert head[12:16] == b"IHDR", "IHDR 位置异常: %s" % p
    return (int.from_bytes(head[16:20], "big"), int.from_bytes(head[20:24], "big"))


w("==== 中文版论文 v8 构建（Fig1 -> V2 图形化版）====")
w("底本 v7 = %d B  md5 %s" % (os.path.getsize(V7), md5(open(V7, "rb").read())))

zin = zipfile.ZipFile(V7)
members = zin.infolist()
doc = zin.read("word/document.xml").decode("utf-8")
img1_v7 = zin.read(MEDIA1)

# ---------------------------------------------------------------- A
w("")
w("== A. image1 身份核对 ==")
h_old = md5(open(OLD_FIG1, "rb").read())
h_img1 = md5(img1_v7)
w("   v7 image1  md5 %s  %d B" % (h_img1, len(img1_v7)))
w("   figures/Fig1_study_design.png md5 %s  %d B" % (h_old, os.path.getsize(OLD_FIG1)))
assert h_img1 == h_old, "v7 的 image1 不是 figures/Fig1_study_design.png，停止"
w("   -> SAME（image1 确为旧版 Fig1，可安全替换）")

# ---------------------------------------------------------------- B
w("")
w("== B. cy 算法自证（用旧 PNG 复现 v7 的 cy）==")
ow, oh = png_px(OLD_FIG1)
calc_old = int(round(CX * oh / ow))
m_first = re.search(r'<wp:extent cx="(\d+)" cy="(\d+)"/>', doc)
v7_cy1 = int(m_first.group(2))
w("   旧 Fig1 像素 = (%d, %d)   calc_cy = %d   v7 第1内联 cy = %d   Δ = %+d"
  % (ow, oh, calc_old, v7_cy1, calc_old - v7_cy1))
assert abs(calc_old - v7_cy1) <= TOL, "cy 算法自证失败（Δ>%d EMU），停止" % TOL
w("   -> PASS（容差 ≤%d EMU）" % TOL)

# ---------------------------------------------------------------- C
w("")
w("== C. 新 cy（V2 PNG 宽高比）==")
nw, nh = png_px(NEW_FIG1)
cy_new = int(round(CX * nh / nw))
w("   新 Fig1 像素 = (%d, %d)   cy %d -> %d   （显示高度 %.3f -> %.3f cm）"
  % (nw, nh, v7_cy1, cy_new, v7_cy1 / 360000.0, cy_new / 360000.0))
assert cy_new != v7_cy1, "新旧 cy 相同，无需改 document.xml（预期应不同）"

cnt = doc.count('cy="%d"' % v7_cy1)
w('   cy="%d" 在 document.xml 中出现 %d 次（应为 2：wp:extent + a:ext）' % (v7_cy1, cnt))
assert cnt == 2, "旧 cy 不唯一，停止"

doc2 = doc.replace('cy="%d"' % v7_cy1, 'cy="%d"' % cy_new)
norm_old = re.sub(r'cy="\d+"', 'cy="#"', doc)
norm_new = re.sub(r'cy="\d+"', 'cy="#"', doc2)
w("   归一化（把所有 cy 数字抹成 #）后逐字相同 = %s" % (norm_old == norm_new))
assert norm_old == norm_new, "除 cy 数字外还有改动，停止"

# ---------------------------------------------------------------- 重打包
w("")
w("== D. 重打包 ==")
new_img1 = open(NEW_FIG1, "rb").read()
# ★ 坑（本轮实测）：不能把**源包的 ZipInfo 对象**直接交给 `zout.writestr(zi, data)` ——
#   writestr 会把 zi.header_offset 改写为**新包**里的偏移，而 zi 就是 zin.NameToInfo 里的
#   同一个对象，于是之后任何 `zin.read(name)` 都会按被污染的偏移去源包里定位，
#   抛 `zipfile.BadZipFile: Bad magic number for file header`（复核段实测命中）。
#   正解：为每个成员**新建一份 ZipInfo 副本**（保留日期/压缩类型/属性），源对象零污染。
zout = zipfile.ZipFile(V8, "w", zipfile.ZIP_DEFLATED)
for zi in members:
    if zi.filename == "word/document.xml":
        data = doc2.encode("utf-8")
    elif zi.filename == MEDIA1:
        data = new_img1
    else:
        data = zin.read(zi.filename)
    ni = zipfile.ZipInfo(zi.filename, date_time=zi.date_time)
    ni.compress_type = zi.compress_type
    ni.external_attr = zi.external_attr
    ni.internal_attr = zi.internal_attr
    ni.create_system = zi.create_system
    zout.writestr(ni, data)
zout.close()

# ---------------------------------------------------------------- 复核
w("")
w("== E. 复核 ==")
w("v8 = %d B  md5 %s" % (os.path.getsize(V8), md5(open(V8, "rb").read())))
z2 = zipfile.ZipFile(V8)
w("   部件数 v7=%d  v8=%d  单侧独有=%d"
  % (len(members), len(z2.namelist()),
     len(set(zin.namelist()) ^ set(z2.namelist()))))
w("   image1 md5 %s -> %s   == V2 PNG ? %s"
  % (h_img1, md5(z2.read(MEDIA1)), md5(z2.read(MEDIA1)) == md5(new_img1)))
same = diff = 0
# 预期唯二改动：image1.png（换图）与 document.xml（仅 cy 数字，下面单独判定）
SPECIAL = ("word/media/image1.png", "word/document.xml")
for n in sorted(set(zin.namelist()) & set(z2.namelist())):
    if n in SPECIAL:
        continue
    a, b = zin.read(n), z2.read(n)
    if a == b:
        same += 1
    else:
        diff += 1
        w("   !! 与 v7 不同：%s" % n)
w("   除 image1 与 document.xml 外，与 v7 逐字节相同 = %d 件；不同 = %d 件" % (same, diff))
assert diff == 0, "存在意外改动，停止"

d1 = zin.read("word/document.xml").decode("utf-8")
d2 = z2.read("word/document.xml").decode("utf-8")
w("   document.xml 长度 %d -> %d；差异 = %s"
  % (len(d1), len(d2), "仅 cy 数字" if norm_old == norm_new else "异常"))
w("")
w("VERDICT = %s" % ("PASS" if (diff == 0 and os.path.exists(V8)) else "FAIL"))

with io.open(LOG, "w", encoding="utf-8", newline="\n") as fh:
    fh.write("\n".join(out) + "\n")
print("WROTE", LOG)
