# -*- coding: utf-8 -*-
"""s72c_crosscheck.py —— v6 与 v5（docx + md）的交叉核验

目的（独立于构建脚本的自证）：
  1. v6 docx 的「局限段」必须与 **v5 docx** 的同一段**逐字相同**
     （v5 docx 是 s64b 依裁定 5 做纯删除后的既成产物 → 是删除动作的ground truth；
      若不一致说明我的删除与 v5 口径有出入）
  2. v6 docx 的「图 4 图注」必须与 **v5 md** 的同段（去 markdown 标记后）逐字相同
  3. v6 docx 的 4 条 S 图注必须与 **v6 md**（去标记后）逐字相同
  4. 顺带核实 v4→v5 登记册里「432 → 383 字」的字符数记录是否与实测一致

日志：logs/_s72c_crosscheck.txt
"""
import io
import os
import re
from docx import Document

ROOT = r"D:/endometriosis_project/11_sc_eqtl_mr_project"
LOG = os.path.join(ROOT, "logs", "_s72c_crosscheck.txt")
L = []
def p(*a):
    L.append(" ".join(str(x) for x in a))


def norm(s):
    s = re.sub(r"[*`_>#\[\]]", "", s)
    return re.sub(r"\s+", "", s)


def paras_text(path):
    d = Document(path)
    return [x.text.strip() for x in d.paragraphs if x.text.strip()]


docx_v4 = paras_text(os.path.join(ROOT, "42_论文初稿_带图_v4.docx"))
docx_v5 = paras_text(os.path.join(ROOT, "43_论文初稿_带图_v5.docx"))
docx_v6 = paras_text(os.path.join(ROOT, "44_论文初稿_带图_v6.docx"))
md_v5 = io.open(os.path.join(ROOT, "43_论文初稿_带图_v5.md"), encoding="utf-8").read()
md_v6 = io.open(os.path.join(ROOT, "44_论文初稿_带图_v6.md"), encoding="utf-8").read()

p("非空段落数： v4=%d  v5=%d  v6=%d" % (len(docx_v4), len(docx_v5), len(docx_v6)))


def one(lst, kw, label):
    h = [t for t in lst if kw in t]
    assert len(h) == 1, "%s 命中 %d" % (label, len(h))
    return h[0]


# ---- 1. 局限段三版对照 ----
lim_v4 = one(docx_v4, "本研究将因果结论限定在区域层级", "v4 局限段")
lim_v5 = one(docx_v5, "本研究将因果结论限定在区域层级", "v5 局限段")
lim_v6 = one(docx_v6, "本研究将因果结论限定在区域层级", "v6 局限段")
p("")
p("---- 1. 局限段（§4）----")
p("v4 长度 = %d 字；以「性别分层」结尾段：…%s" % (len(lim_v4), lim_v4[-70:]))
p("v5 长度 = %d 字" % len(lim_v5))
p("v6 长度 = %d 字" % len(lim_v6))
p("v4 → v5 净减 = %d 字（登记册记「432 → 383，净减 49」）" % (len(lim_v4) - len(lim_v5)))
p("v5 == v6 ? %s" % (lim_v5 == lim_v6))
assert lim_v5 == lim_v6, "★v6 局限段与 v5 docx 不一致"
p("✅ v6 的删除结果与 v5 docx 逐字相同（删除口径一致）")
p("v6 局限段尾部：…%s" % lim_v6[-70:])
p("v6 含「性别分层」= %s ；含「。在安全性评估一侧」= %s" % (
    "性别分层" in lim_v6, "。在安全性评估一侧" in lim_v6))
assert "性别分层" not in lim_v6 and "。在安全性评估一侧" in lim_v6
# md 侧
lim_md_v5 = one([ln.strip() for ln in md_v5.split("\n") if ln.strip()],
                "本研究将因果结论限定在区域层级", "v5 md 局限段")
lim_md_v6 = one([ln.strip() for ln in md_v6.split("\n") if ln.strip()],
                "本研究将因果结论限定在区域层级", "v6 md 局限段")
p("md 侧 v5/v6 一致 = %s ；md == docx（归一化）= %s" % (
    lim_md_v5 == lim_md_v6, norm(lim_md_v6) == norm(lim_v6)))
assert norm(lim_md_v6) == norm(lim_v6)

# ---- 2. 图 4 图注 ----
p("")
p("---- 2. 图 4 图注 ----")
cap4_v4 = one(docx_v4, "图 4 | chr1p36.12", "v4 图 4 图注")
cap4_v6 = one(docx_v6, "图 4 | chr1p36.12", "v6 图 4 图注")
cap4_md5 = one([ln.strip() for ln in md_v5.split("\n") if ln.strip()], "**图 4 |", "v5 md 图 4 图注")
cap4_md6 = one([ln.strip() for ln in md_v6.split("\n") if ln.strip()], "**图 4 |", "v6 md 图 4 图注")
p("v4 图注含 '1000 Genomes' = %s" % ("1000 Genomes" in cap4_v4))
p("v6 图注含 'OneK1K 参考面板，980 名供者' = %s" % ("OneK1K 参考面板，980 名供者" in cap4_v6))
p("md v5 == md v6（图 4 图注） = %s" % (cap4_md5 == cap4_md6))
p("docx v6 == md v6（归一化） = %s" % (norm(cap4_v6) == norm(cap4_md6)))
assert norm(cap4_v6) == norm(cap4_md6), "★docx 与 md 的图 4 图注不一致"
p("✅ docx 图 4 图注与 v6 md 逐字一致（含 §3.3 图注更正）")

# ---- 3. S1–S4 图注 ----
p("")
p("---- 3. 补充图 S1–S4 图注（docx vs md）----")
for k in ["S1", "S2", "S3", "S4"]:
    cd = one(docx_v6, "图 %s |" % k, "v6 docx 图 %s" % k)
    cm = one([ln.strip() for ln in md_v6.split("\n") if ln.strip()], "**图 %s |" % k, "v6 md 图 %s" % k)
    ok = norm(cd) == norm(cm)
    p("  图 %s：docx %d 字 / md %d 字 / 归一化一致=%s" % (k, len(cd), len(cm), ok))
    assert ok, "图 %s 的 docx 与 md 图注不一致" % k
    assert ("（a" in cd) or (k == "S4"), "图 %s 图注缺面板标签" % k
p("✅ 4 条 S 图注在 docx 与 md 中逐字一致")

# ---- 4. 稳定性复核 ----
p("")
p("---- 4. 参考串在位 ----")
t6 = "\n".join(docx_v6)
for s, exp in [("质量控制与分层视图另见补充图 S1–S4。", 1), ("1000 Genomes", 0),
               ("性别分层", 0), ("MHC", 0), ("补充图 S1", 1)]:
    c = t6.count(s)
    p("  count(%s) = %d（期望 %d）" % (s, c, exp))
    assert c == exp, "参考串计数不符：%s" % s

p("")
p("VERDICT = PASS")
with io.open(LOG, "w", encoding="utf-8", newline="\n") as fh:
    fh.write("\n".join(L) + "\n")
print("WROTE", LOG)
