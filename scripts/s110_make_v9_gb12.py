# -*- coding: utf-8 -*-
"""s110_make_v9_gb12.py —— v8 → v9：把论文里**全部 9 张图**换成 GB12／α=0.35 版。

用户 2026-10-01 指令：
  「我最终选择35%的，帮我插入论文中，顺便把之前用GB12色盘画的图也一并插入论文中」

★ 换图清单（**9 张全换**）：
  image1  ← `_fig1_draft_v2/band_alpha/Fig1v2_bandalpha_a035.png`（V2 图形化版，a/b 带 α=0.35）
  image2  ← figures/Fig2_discovery_replication.png         （GB12）
  image3  ← figures/Fig3_coloc_sensitivity.png             （GB12）
  image4  ← figures/Fig4_chr1p36_12_finemap.png            （GB12）
  image5  ← figures/Fig5_opentargets_phewas_safety.png     （GB12）
  image6  ← figures/FigS1_instrument_power_gating.png      （GB12）
  image7  ← figures/FigS2_sensitivity_outcome_concordance.png（GB12）
  image8  ← figures/FigS3_manhattan_by_celltype.png        （GB12）
  image9  ← figures/FigS4_discovery_volcano.png            （GB12）

★★ 关键前提（本脚本自证 A）：**新旧图逐张像素尺寸完全相同**（w 与 h 都相等）
   ⇒ docx 内联的 `cy` 一律不用重算 ⇒ **`document.xml` 逐字节不变**。
   这是本项目里最干净的一种换图：只替换 `word/media/*` 的字节。

★ 实现纪律：
  ① 重打包必须**为每个成员新建 ZipInfo 副本**（上一轮踩过：把源包 zi 直接交给 writestr
     会污染 zi.header_offset，导致之后 zin.read 抛 BadZipFile）；
  ② 自证 B：除这 9 件 media 外，其余部件与 v8 **逐字节相同**；`document.xml` **逐字节相同**；
  ③ 自证 C：新 media 的 md5 必须 ≠ 旧 md5（否则等于没换）。
"""
import hashlib
import io
import os
import zipfile

ROOT = r"D:/endometriosis_project/11_sc_eqtl_mr_project"
FIGD = os.path.join(ROOT, "figures")
V8 = os.path.join(ROOT, "46_中文版论文_带图_v8.docx")
V9 = os.path.join(ROOT, "47_中文版论文_带图_v9_GB12.docx")
LOG = os.path.join(ROOT, "logs", "_s110_make_v9.txt")

MEDIA_MAP = [
    ("word/media/image1.png",
     r"D:/_transfer_logs/_fig1_draft_v2/band_alpha/Fig1v2_bandalpha_a035.png"),
    ("word/media/image2.png", os.path.join(FIGD, "Fig2_discovery_replication.png")),
    ("word/media/image3.png", os.path.join(FIGD, "Fig3_coloc_sensitivity.png")),
    ("word/media/image4.png", os.path.join(FIGD, "Fig4_chr1p36_12_finemap.png")),
    ("word/media/image5.png", os.path.join(FIGD, "Fig5_opentargets_phewas_safety.png")),
    ("word/media/image6.png", os.path.join(FIGD, "FigS1_instrument_power_gating.png")),
    ("word/media/image7.png", os.path.join(FIGD, "FigS2_sensitivity_outcome_concordance.png")),
    ("word/media/image8.png", os.path.join(FIGD, "FigS3_manhattan_by_celltype.png")),
    ("word/media/image9.png", os.path.join(FIGD, "FigS4_discovery_volcano.png")),
]

L = []


def W(s=""):
    L.append(str(s))


def md5b(b):
    return hashlib.md5(b).hexdigest()


def png_wh(b):
    """读 PNG IHDR 得 (w, h) 像素；不依赖任何第三方库。"""
    assert b[:8] == b"\x89PNG\r\n\x1a\n", "不是 PNG"
    assert b[12:16] == b"IHDR", "IHDR 位置异常"
    return (int.from_bytes(b[16:20], "big"), int.from_bytes(b[20:24], "big"))


zin = zipfile.ZipFile(V8)
members = zin.infolist()

W("=" * 96)
W("s110  v8 -> v9（换全部 9 张图）")
W("=" * 96)
W("底本 : %s  (%d 件部件)" % (V8, len(members)))
W("产物 : %s" % V9)
W("")

# ---------------- 读新图 ----------------
newdata = {}
W("---- 换图清单与自证 A（新旧像素尺寸必须逐张相等） ----")
W("%-24s %-22s %-22s %s" % ("media", "旧 px", "新 px", "判定"))
W("-" * 96)
okA = True
old_wh = {}
for name, src in MEDIA_MAP:
    ob = zin.read(name)
    nb = open(src, "rb").read()
    ow, nw = png_wh(ob), png_wh(nb)
    same = (ow == nw)
    okA &= same
    old_wh[name] = ow
    newdata[name] = nb
    W("%-24s %-22s %-22s %s   %s"
      % (os.path.basename(name), "%dx%d" % ow, "%dx%d" % nw,
         "SAME" if same else "!! DIFFER",
         "" if md5b(ob) != md5b(nb) else "!! 新旧 md5 相同（等于没换）"))
assert okA, "自证 A 失败：有图的像素尺寸发生变化 ⇒ document.xml 的 cy 必须重算，脚本停止"
W("")
W("自证 A: PASS —— 9/9 像素尺寸一致 ⇒ document.xml 无需改动")
W("")

# 逐张 md5 变化表
W("---- 自证 C（md5 必变） ----")
W("%-24s %-34s -> %-34s" % ("media", "旧 md5", "新 md5"))
W("-" * 96)
okC = True
for name, _src in MEDIA_MAP:
    o = md5b(zin.read(name))
    n = md5b(newdata[name])
    okC &= (o != n)
    W("%-24s %-34s -> %-34s" % (os.path.basename(name), o, n))
assert okC, "自证 C 失败：存在未变化的图"
W("")

# ---------------- 重打包 ----------------
zout = zipfile.ZipFile(V9, "w", zipfile.ZIP_DEFLATED)
for zi in members:
    if zi.filename in newdata:
        data = newdata[zi.filename]
    else:
        data = zin.read(zi.filename)
    ni = zipfile.ZipInfo(zi.filename, date_time=zi.date_time)   # ★ 新建副本，勿复用 zi
    ni.compress_type = zi.compress_type
    ni.external_attr = zi.external_attr
    ni.internal_attr = zi.internal_attr
    ni.create_system = zi.create_system
    zout.writestr(ni, data)
zout.close()

# ---------------- 自证 B ----------------
z2 = zipfile.ZipFile(V9)
n_v8, n_v9 = set(zin.namelist()), set(z2.namelist())
W("---- 自证 B（部件级） ----")
W("部件数 v8=%d  v9=%d ；仅 v8 有=%s ；仅 v9 有=%s"
  % (len(n_v8), len(n_v9), sorted(n_v8 - n_v9) or "无", sorted(n_v9 - n_v8) or "无"))
assert n_v8 == n_v9, "部件集合不一致，停止"

same = diff = 0
changed = []
for n in sorted(n_v8):
    a, b = zin.read(n), z2.read(n)
    if a == b:
        same += 1
    else:
        diff += 1
        changed.append(n)
W("与 v8 逐字节相同的部件 = %d ；不同的 = %d" % (same, diff))
W("不同的部件清单：%s" % changed)
expected = sorted(n for n, _ in MEDIA_MAP)
assert sorted(changed) == expected, "出现了预期之外的改动，停止"

dm = zin.read("word/document.xml")
dm2 = z2.read("word/document.xml")
W("document.xml  v8 md5=%s  v9 md5=%s  ->  %s"
  % (md5b(dm), md5b(dm2), "逐字节相同" if dm == dm2 else "!! 已变更"))
assert dm == dm2, "document.xml 发生了变化，停止"

# 回读校验：v9 内 media 的 md5 必须等于源文件 md5
W("")
W("---- 回读校验（v9 内 media == 源文件） ----")
okD = True
for name, src in MEDIA_MAP:
    a = md5b(z2.read(name))
    b = md5b(open(src, "rb").read())
    okD &= (a == b)
    W("  %-24s %s  %s" % (os.path.basename(name), "OK " if a == b else "!! ", a))
assert okD, "回读不一致，停止"

W("")
W("文件大小 v8=%d B  v9=%d B  (Δ=%+d)"
  % (os.path.getsize(V8), os.path.getsize(V9), os.path.getsize(V9) - os.path.getsize(V8)))
W("产物 md5 = %s" % md5b(open(V9, "rb").read()))
W("")
W("VERDICT = PASS")

# ---- 顺带盘点：正文里是否还有会与新配色冲突的颜色词 ----
try:
    import re
    txt = re.sub(r"<[^>]+>", "", dm.decode("utf-8"))
    hits = []
    for kw in ("红色", "深红", "浅红", "红叉", "红框", "橙色", "深蓝", "浅蓝", "绿色",
               "灰色", "紫色", "黄色"):
        c = txt.count(kw)
        if c:
            hits.append("%s x%d" % (kw, c))
    W("")
    W("★ 正文颜色词盘点（供图注一致性核查，本脚本**不改文本**）：%s" % ("；".join(hits) or "无"))
except Exception as e:                                          # noqa: BLE001
    W("颜色词盘点失败：%r" % (e,))

os.makedirs(os.path.dirname(LOG), exist_ok=True)
with io.open(LOG, "w", encoding="utf-8", newline="\n") as fh:
    fh.write("\n".join(L) + "\n")
print("v9 done -> %s" % V9)
