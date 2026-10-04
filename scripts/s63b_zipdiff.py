# -*- coding: utf-8 -*-
"""s63b_zipdiff.py —— 对底本 docx 与输出 docx 做 zip 部件级比对（只读）
确认 python-docx 往返未丢件（尤其 word/media/*）。
日志：logs/s63b_zipdiff.log
"""
import io, os, zipfile, hashlib
ROOT = r"D:/endometriosis_project/11_sc_eqtl_mr_project"
A = os.path.join(ROOT, "41_论文初稿_带图_v3.docx")
B = os.path.join(ROOT, "42_论文初稿_带图_v4.docx")
LOG = os.path.join(ROOT, "logs", "s63b_zipdiff.log")
L = []
def p(*a):
    s = " ".join(str(x) for x in a); L.append(s)
    try: print(s)
    except UnicodeEncodeError: print(s.encode("ascii","backslashreplace").decode("ascii"))

za, zb = zipfile.ZipFile(A), zipfile.ZipFile(B)
na = {i.filename: i for i in za.infolist()}
nb = {i.filename: i for i in zb.infolist()}
p("底本部件数 = %d / 输出部件数 = %d" % (len(na), len(nb)))
only_a = sorted(set(na) - set(nb)); only_b = sorted(set(nb) - set(na))
p("-- 仅存于底本的部件（%d）--" % len(only_a))
for k in only_a: p("   A_ONLY %-52s %9d B" % (k, na[k].file_size))
p("-- 仅存于输出的部件（%d）--" % len(only_b))
for k in only_b: p("   B_ONLY %-52s %9d B" % (k, nb[k].file_size))

tot_a = sum(na[k].file_size for k in na)
tot_b = sum(nb[k].file_size for k in nb)
p("解压后总字节：底本 %d / 输出 %d（差 %+d）" % (tot_a, tot_b, tot_b - tot_a))

p("-- 共有部件中内容有变的（按解压后 sha1）--")
chg = []
for k in sorted(set(na) & set(nb)):
    ha = hashlib.sha1(za.read(k)).hexdigest()
    hb = hashlib.sha1(zb.read(k)).hexdigest()
    if ha != hb:
        chg.append((k, na[k].file_size, nb[k].file_size))
for k, sa, sb in chg:
    p("   CHANGED %-48s %9d -> %9d" % (k, sa, sb))

p("-- word/media 逐项 --")
for k in sorted(x for x in set(na) & set(nb) if x.startswith("word/media/")):
    p("   %-32s %9d -> %9d  %s" % (k, na[k].file_size, nb[k].file_size,
      "SAME" if hashlib.sha1(za.read(k)).hexdigest() == hashlib.sha1(zb.read(k)).hexdigest() else "DIFF"))
p("MEDIA_A=%d MEDIA_B=%d" % (sum(1 for k in na if k.startswith("word/media/")),
                             sum(1 for k in nb if k.startswith("word/media/"))))
io.open(LOG, "w", encoding="utf-8", newline="\n").write("\n".join(L) + "\n")
print("ZIPDIFF_DONE")
