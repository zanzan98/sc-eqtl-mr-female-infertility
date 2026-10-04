# -*- coding: utf-8 -*-
"""s72d_render_pages.py —— 渲染 v6 预览 PDF 的关键页为 PNG（目视复核用）"""
import fitz, os

PDF = r"D:/endometriosis_project/11_sc_eqtl_mr_project/44_论文初稿_带图_v6_预览.pdf"
OUT = r"D:/endometriosis_project/11_sc_eqtl_mr_project/logs/_s72d_pages"
PAGES = [20, 21, 22, 23, 24]
os.makedirs(OUT, exist_ok=True)
doc = fitz.open(PDF)
for n in PAGES:
    if n > doc.page_count:
        continue
    doc[n - 1].get_pixmap(dpi=110).save(os.path.join(OUT, "v6_p%02d.png" % n))
print("pages rendered", PAGES, "->", OUT)
