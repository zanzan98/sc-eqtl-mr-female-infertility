# -*- coding: utf-8 -*-
"""s69b_verify_manuscript.py — 复核 s69 补丁结果（裁定 2/3/5 的判据）"""
import io, os, re
ROOT = r"D:/endometriosis_project/11_sc_eqtl_mr_project"
FILES = ["38_英文正文_论文体重写_v5.md", "39_英文正文_论文体v6.md",
         "38_中文版_论文体重写_v5_正文与图注.md", "39_中文版_论文体v7_正文与图注.md",
         "43_论文初稿_带图_v5.md", "39b_技术文档_SM与SN_v8.md"]
CAP_EN_NEW = b"(r\xc2\xb2, OneK1K reference panel, 980 donors)"
CAP_ZH_NEW = u"（r²，OneK1K 参考面板，980 名供者）".encode("utf-8")
CITE_EN = "Additional quality-control and stratified views are provided in Supplementary Figs. S1\u2013S4."
CITE_ZH = u"质量控制与分层视图另见补充图 S1–S4。"
old_leg_en = b"Fig. S1 Instrument power and analytical gating by cell type (instrument counts"
old_leg_zh = u"图 S1 各细胞类型的工具变量效力与分析门控（各细胞类型的工具数".encode("utf-8")

out = []
for f in FILES:
    p = os.path.join(ROOT, f)
    d = open(p, "rb").read()
    t = d.decode("utf-8")
    out.append("=" * 80)
    out.append("FILE %s  bytes=%d" % (f, len(d)))
    out.append("  1000 Genomes count      = %d" % d.count(b"1000 Genomes"))
    out.append("  cite EN count           = %d" % t.count(CITE_EN))
    out.append("  cite ZH count           = %d" % t.count(CITE_ZH))
    out.append("  cap EN new count        = %d" % d.count(CAP_EN_NEW))
    out.append("  cap ZH new count        = %d" % d.count(CAP_ZH_NEW))
    out.append("  old legend EN count     = %d" % d.count(old_leg_en))
    out.append("  old legend ZH count     = %d" % d.count(old_leg_zh))
    out.append("  new legend EN (S4) count= %d" % d.count(b"Fig. S4 Distribution of discovery"))
    out.append("  new legend ZH (S4) count= %d" % d.count(u"图 S4 Discovery MR 效应量分布".encode("utf-8")))
    # 引用上下文
    for mark, name in [(CITE_EN, "CITE_EN"), (CITE_ZH, "CITE_ZH")]:
        j = t.find(mark)
        if j >= 0:
            out.append("  [%s] @%d  ...%s..." % (name, j, t[max(0, j - 90):j + len(mark) + 90].replace("\r\n", "\\r\\n").replace("\n", "\\n")))
    j = t.find("OneK1K")
    if j >= 0:
        out.append("  [OneK1K] @%d  ...%s..." % (j, t[max(0, j - 90):j + 90].replace("\r\n", "\\r\\n").replace("\n", "\\n")))
    # 全局重复相邻行守卫
    lines = t.split("\n")
    dup = [i for i in range(1, len(lines)) if lines[i].strip() and lines[i] == lines[i - 1] and len(lines[i].strip()) > 20]
    out.append("  adjacent-dup line idx    = %s" % (dup if dup else "none"))

p = os.path.join(ROOT, "logs", "_s69b_verify.txt")
with io.open(p, "w", encoding="utf-8", newline="\n") as fh:
    fh.write("\n".join(out) + "\n")
print("wrote", p)
