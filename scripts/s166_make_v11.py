# -*- coding: utf-8 -*-
"""s166_make_v11.py —— v10 → v11：把论文里的 Fig2–Fig5 / FigS1–S4 换成**三项全局标准统一后**的版本。

用户 2026-10-01 指令：「每张图图序字母与图的距离要跟模板fig2的一致，图要整体居中，
图与图之间的距离也要与模板fig2一致，请修改，修改后同样插入论文给我预览」。
★ 本轮改的是 Fig3/4/5 与 FigS1/S2/S3（Fig2 是模板、冻结；Fig1 无面板字母、单面板）。

★ 换图清单（**6 张**；被排除的 image1 / image2 / image9 见下方 MEDIA_MAP 的逐条说明）：
  image2 ← figures/Fig2_discovery_replication.png     （v12：背景点再放大 + 图例分开）
  image3 ← figures/Fig3_coloc_sensitivity.png         （(c) 点放大 + 同色系描边）
  image4 ← figures/Fig4_chr1p36_12_finemap.png        （(b)(c)(e) 点放大 + 描边）
  image5 ← figures/Fig5_opentargets_phewas_safety.png （(a)(b) 点放大 + 描边）
  image6 ← figures/FigS1_instrument_power_gating.png
  image7 ← figures/FigS2_sensitivity_outcome_concordance.png
  image8 ← figures/FigS3_manhattan_by_celltype.png
  image9 ← figures/FigS4_discovery_volcano.png

★★ 关键差别（**与 v9 那次不同**）：v9 换图时新旧图**像素尺寸逐张相等**，所以 `document.xml`
   逐字节不变。**本轮不是** —— v9 内嵌的 Fig2/3/4/5 是**幅高压缩之前**的旧版，纵横比差很多
   （如 Fig4：旧 4307×5737 AR=0.751，新 4307×3961 AR=1.087）。
   ⇒ 内联尺寸 `cy` **必须按新纵横比重算**，否则 Word 会把图**拉伸变形**。

   Word 的显示尺寸 = `wp:extent` 的 (cx, cy)（EMU，1 mm = 36000 EMU）。v9 里 9 张图
   **一律 cx = 5760000 EMU = 160.00 mm**（正文栏宽）⇒ 本脚本**只改 cy**，令
   `cy = round(cx × h_px / w_px)`，保持 160 mm 宽、按真实纵横比缩放。
   ★ `wp:extent` 与 `a:ext` 两处 cy 必须同步改（只改一处，Word 里会出现"框对图不对"）。

★ 实现纪律（沿用 v9）：
  ④ 自证：`image1`/`image2`/`image9` **不在**换图清单内（各自的 md5 判定理由见 MEDIA_MAP）。
  ① 重打包**为每个成员新建 ZipInfo 副本**（复用源 `zi` 会污染 `header_offset` ⇒ BadZipFile）；
  ② 自证：除「8 件 media + document.xml」外其余部件与 v9 **逐字节相同**；
  ③ 自证：新 media 的 md5 必须 ≠ 旧 md5；回读后 v11 内 media 的 md5 == 源文件 md5。
"""
import hashlib
import io
import os
import re
import zipfile

ROOT = r"D:/endometriosis_project/11_sc_eqtl_mr_project"
FIGD = os.path.join(ROOT, "figures")
V9 = os.path.join(ROOT, "48_中文版论文_带图_v10_GB12.docx")
V10 = os.path.join(ROOT, "49_中文版论文_带图_v11_GB12.docx")
MD9 = os.path.join(ROOT, "48_中文版论文_带图_v10_GB12.md")
MD10 = os.path.join(ROOT, "49_中文版论文_带图_v11_GB12.md")
LOG = os.path.join(ROOT, "logs", "_s166_make_v11.txt")

EMU_CX = 5760000           # = 160.00 mm，v10 里 9 张图统一的显示宽度

# media 名 -> 新图源文件
# ★★ 2026-10-01 v11 换图清单 = **6 张**（不是 v10 的 8 张）。逐条说明被排除的 3 张：
#   · image1 (Fig1) —— **绝不能换**：v10 内嵌的是用户选定的 **V2 / 背景 α=0.35** 版
#     （md5 fffbf2a5…），而 `figures/Fig1_study_design.png` 仍是**未提升的正本**
#     （md5 2ce5ac81…）。按 media 名机械换图会把 V2 降级回旧版。
#   · image2 (Fig2) —— 本轮 Fig2 是**模板、冻结不改**；v10 内嵌 PNG 与 figures/ 当前
#     PNG **逐字节相同**（md5 3509e6a4…），列入只会触发"md5 未变"的假警报。
#   · image9 (FigS4) —— 本轮未改（单面板、无面板字母、无行带；
#     上一轮点风格统一时已改并已进 v10），md5 bbd53af7… 前后一致。
MEDIA_MAP = [
    ("word/media/image3.png", os.path.join(FIGD, "Fig3_coloc_sensitivity.png")),
    ("word/media/image4.png", os.path.join(FIGD, "Fig4_chr1p36_12_finemap.png")),
    ("word/media/image5.png", os.path.join(FIGD, "Fig5_opentargets_phewas_safety.png")),
    ("word/media/image6.png", os.path.join(FIGD, "FigS1_instrument_power_gating.png")),
    ("word/media/image7.png", os.path.join(FIGD, "FigS2_sensitivity_outcome_concordance.png")),
    ("word/media/image8.png", os.path.join(FIGD, "FigS3_manhattan_by_celltype.png")),
]

L = []


def W(s=""):
    L.append(str(s))


def md5b(b):
    return hashlib.md5(b).hexdigest()


def png_wh(b):
    """读 PNG IHDR 得 (w, h) 像素；不依赖第三方库。"""
    assert b[:8] == b"\x89PNG\r\n\x1a\n", "不是 PNG"
    assert b[12:16] == b"IHDR", "IHDR 位置异常"
    return (int.from_bytes(b[16:20], "big"), int.from_bytes(b[20:24], "big"))


zin = zipfile.ZipFile(V9)
members = zin.infolist()

W("=" * 100)
W("s166  v10 -> v11（换 Fig3/4/5 + FigS1/2/3 共 6 张；并按新纵横比重算内联 cy）")
W("=" * 100)
W("底本 : %s  (%d 件部件)" % (V9, len(members)))
W("产物 : %s" % V10)
W("")

# ---------------------------------------------------------------- ① rId -> media 映射
rels = zin.read("word/_rels/document.xml.rels").decode("utf-8")
rmap = dict(re.findall(r'Id="([^"]+)"[^>]*?Target="([^"]+)"', rels))
r2m = {k: ("word/" + v if not v.startswith("word/") else v) for k, v in rmap.items()}

doc = zin.read("word/document.xml").decode("utf-8")

# ---------------------------------------------------------------- ② 读新图 + 算新 cy
newdata = {}
new_cy = {}          # media 全名 -> (old_cy, new_cy)
W("---- 换图清单（%d 张）/ 新旧像素尺寸 / 重算后的 cy ----" % len(MEDIA_MAP))
W("%-22s %-14s %-14s %-11s %-11s %s" % ("media", "旧 px", "新 px", "旧 cy", "新 cy", "判定"))
W("-" * 100)
for name, src in MEDIA_MAP:
    ob = zin.read(name)
    nb = open(src, "rb").read()
    ow, nw = png_wh(ob), png_wh(nb)
    newdata[name] = nb
    # 该 media 在 document.xml 里对应哪一个 drawing（按 rId 反查）
    rid = [k for k, v in r2m.items() if v == name]
    assert len(rid) == 1, "%s: rId 反查得到 %d 个" % (name, len(rid))
    pat = r'<w:drawing>(?:(?!</w:drawing>).)*?r:embed="%s"(?:(?!</w:drawing>).)*?</w:drawing>' % re.escape(rid[0])
    blk = re.search(pat, doc, re.S)
    assert blk, "%s: document.xml 里找不到对应 drawing" % name
    cys = re.findall(r'<(?:wp:extent|a:ext) cx="%d" cy="(\d+)"' % EMU_CX, blk.group(0))
    assert len(cys) == 2 and cys[0] == cys[1], "%s: wp:extent/a:ext 的 cy 不是 2 处一致 -> %s" % (name, cys)
    o_cy = int(cys[0])
    n_cy = int(round(EMU_CX * nw[1] / float(nw[0])))
    new_cy[name] = (o_cy, n_cy)
    W("%-22s %-14s %-14s %-11d %-11d %s"
      % (os.path.basename(name), "%dx%d" % ow, "%dx%d" % nw, o_cy, n_cy,
         ("尺寸同" if ow == nw else "尺寸变") + ("/md5变" if md5b(ob) != md5b(nb) else "/!!md5未变")))
W("")

# ---------------------------------------------------------------- ③ 改 document.xml
W("---- document.xml：逐 drawing 替换 cy（wp:extent 与 a:ext 同步） ----")
out_parts, pos, nblk = [], 0, 0
for m in re.finditer(r"<w:drawing>.*?</w:drawing>", doc, re.S):
    out_parts.append(doc[pos:m.start()])
    blk = m.group(0)
    rid = re.search(r'r:embed="([^"]+)"', blk)
    tgt = r2m.get(rid.group(1)) if rid else None
    if tgt in new_cy:
        o_cy, n_cy = new_cy[tgt]
        nb2, nrep = re.subn(r'(<(?:wp:extent|a:ext) cx="%d" cy=")%d(")' % (EMU_CX, o_cy),
                            lambda mm: mm.group(1) + str(n_cy) + mm.group(2), blk)
        assert nrep == 2, "%s: 期望替换 2 处 cy，实际 %d" % (tgt, nrep)
        W("   %-22s cy %d -> %d   (%.2f -> %.2f mm)"
          % (os.path.basename(tgt), o_cy, n_cy, o_cy / 36000.0, n_cy / 36000.0))
        blk = nb2
        nblk += 1
    out_parts.append(blk)
    pos = m.end()
out_parts.append(doc[pos:])
doc2 = "".join(out_parts)
assert doc2.count("<w:drawing>") == doc.count("<w:drawing>"), "drawing 个数变了"
assert nblk == len(MEDIA_MAP), "实际改动的 drawing 数 %d != %d" % (nblk, len(MEDIA_MAP))
W("   合计改动 %d 个 drawing（期望 %d）" % (nblk, len(MEDIA_MAP)))
W("")

# ---------------------------------------------------------------- ④ 重打包
newparts = dict(newdata)
newparts["word/document.xml"] = doc2.encode("utf-8")
zout = zipfile.ZipFile(V10, "w", zipfile.ZIP_DEFLATED)
for zi in members:
    data = newparts.get(zi.filename)
    if data is None:
        data = zin.read(zi.filename)
    ni = zipfile.ZipInfo(zi.filename, date_time=zi.date_time)   # ★ 新建副本，勿复用 zi
    ni.compress_type = zi.compress_type
    ni.external_attr = zi.external_attr
    ni.internal_attr = zi.internal_attr
    ni.create_system = zi.create_system
    zout.writestr(ni, data)
zout.close()

# ---------------------------------------------------------------- ⑤ 自证 B
z2 = zipfile.ZipFile(V10)
n_a, n_b = set(zin.namelist()), set(z2.namelist())
W("---- 自证 B（部件级） ----")
W("部件数 v10=%d  v11=%d ；仅 v10 有=%s ；仅 v11 有=%s"
  % (len(n_a), len(n_b), sorted(n_a - n_b) or "无", sorted(n_b - n_a) or "无"))
assert n_a == n_b, "部件集合不一致，停止"

same = diff = 0
changed = []
for n in sorted(n_a):
    if zin.read(n) == z2.read(n):
        same += 1
    else:
        diff += 1
        changed.append(n)
W("与 v10 逐字节相同的部件 = %d ；不同的 = %d" % (same, diff))
W("不同的部件清单：%s" % changed)
expected = sorted(list(newdata.keys()) + ["word/document.xml"])
assert sorted(changed) == expected, "出现了预期之外的改动，停止 \n  实际 %s\n  期望 %s" % (sorted(changed), expected)
W("   >>> 仅 media×8 + document.xml 变化 —— 与预期一致")
W("")

# ---------------------------------------------------------------- ⑥ 回读校验
W("---- 回读校验（v11 内 media == 源文件；cy 已生效） ----")
ok = True
d2 = z2.read("word/document.xml").decode("utf-8")
for name, src in MEDIA_MAP:
    a = md5b(z2.read(name))
    b = md5b(open(src, "rb").read())
    o_cy, n_cy = new_cy[name]
    hit = ('cx="%d" cy="%d"' % (EMU_CX, n_cy)) in d2
    ok &= (a == b) and hit
    W("  %-22s md5 %s   cy=%d 在文档中 %s" % (os.path.basename(name), "OK " if a == b else "!!", n_cy, "OK" if hit else "!! 缺失"))
assert ok, "回读不一致，停止"
W("")

# ---------------------------------------------------------------- ⑦ 页宽核算
W("---- 显示尺寸（160 mm 宽，按真实纵横比） ----")
tot = 0.0
for name, src in MEDIA_MAP:
    w, h = png_wh(newdata[name])
    mm_h = EMU_CX * h / w / 36000.0
    W("  %-22s 160.00 x %7.2f mm" % (os.path.basename(name), mm_h))
W("")

# ---------------------------------------------------------------- ⑧ 伴生 md
with open(MD9, "rb") as fh:
    mdb = fh.read()
with open(MD10, "wb") as fh:
    fh.write(mdb)
W("伴生 markdown：%s（与 v10 的 md **逐字节相同**，本轮只换图、不改正文）" % os.path.basename(MD10))
W("")

W("文件大小 v10=%d B  v11=%d B  (Δ=%+d)"
  % (os.path.getsize(V9), os.path.getsize(V10), os.path.getsize(V10) - os.path.getsize(V9)))
W("v11 docx md5 = %s" % md5b(open(V10, "rb").read()))
W("")
W("VERDICT = PASS")

os.makedirs(os.path.dirname(LOG), exist_ok=True)
with io.open(LOG, "w", encoding="utf-8", newline="\n") as fh:
    fh.write("\n".join(L) + "\n")
print("v11 done -> %s" % V10)
