# -*- coding: utf-8 -*-
"""s63e_locate_inserts.py —— 定位 4 处插入与 5 张图在新 v4 预览 PDF 中的页码"""
import fitz, io, os
PDF = r"D:\_wb_docx_build\_exp_draft_v4.pdf"
OUT = r"D:\endometriosis_project\11_sc_eqtl_mr_project\logs\s63e_locate_inserts.log"
buf = []
def p(*a):
    s = " ".join(str(x) for x in a); buf.append(s)
    try: print(s)
    except UnicodeEncodeError: print(s.encode("ascii", "backslashreplace").decode("ascii"))

doc = fitz.open(PDF)
pages = []
for pg in doc:
    t = pg.get_text()
    t = "\n".join(ln for ln in t.splitlines() if not __import__("re").fullmatch(r"\s*\d{1,3}\s*", ln))
    pages.append(t.replace(" ", "").replace("\n", ""))

KEYS = [
    ("I1 插入(三位点比较)", "这构成本研究把该区域单独展开的理由"),
    ("I2 插入(GTEx+OneK1K)", "在三个信号座中的特殊性"),
    ("I4 插入(MVMR)", "上述分解应视为探索性证据"),
    ("I3 插入(方法学校准)", "而不是两套标准"),
    ("I1 锚点", "记忆 B 细胞的复核状态因而有待更大规模的单细胞队列确认"),
    ("I2 锚点", "而 CDC42 更可能是这一区域效应在特定免疫细胞中的"),
    ("I4 锚点", "而非被检验基因自身的 eQTL"),
    ("I3 锚点", "两层共同解释了 WNT4 在全血可检出却在单细胞血液面板中缺席这一格局"),
    ("图1", "研究设计与分析流程"),
    ("图2", "细胞类型特异性 MR 信号的发现与复核"),
    ("图3", "共定位后验概率与各位点的稳健性"),
    ("图4", "多基因连锁不平衡结构与基因归属边界"),
    ("图5", "Open Targets 关联谱与 FinnGen R12"),
]
for name, k in KEYS:
    hits = [i + 1 for i, t in enumerate(pages) if k.replace(" ", "") in t]
    p("%-12s -> pages %s" % (name, hits))
io.open(OUT, "w", encoding="utf-8", newline="\n").write("\n".join(buf) + "\n")
print("LOCATE_DONE")
