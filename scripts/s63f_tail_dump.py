# -*- coding: utf-8 -*-
import fitz, io, os, re
PDF = r"D:\_wb_docx_build\_exp_draft_v4.pdf"
OUT = r"D:\endometriosis_project\11_sc_eqtl_mr_project\logs\s63f_tail_dump.log"
buf = []
def p(*a):
    s = " ".join(str(x) for x in a); buf.append(s)
    try: print(s)
    except UnicodeEncodeError: print(s.encode("ascii","backslashreplace").decode("ascii"))
d = fitz.open(PDF)
for pi in range(15, d.page_count):
    p("========== page %d ==========" % (pi + 1))
    for ln in d[pi].get_text().splitlines():
        p("   |%s|" % ln)
io.open(OUT, "w", encoding="utf-8", newline="\n").write("\n".join(buf) + "\n")
print("TAIL_DUMP_DONE")
