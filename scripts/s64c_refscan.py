# -*- coding: utf-8 -*-
"""s64c_refscan.py —— 引用扫描（修正版：按**目录名分量**排除，不再用子串 '_'）
目的：确认对候选待重命名文件的引用是否存在，避免重命名后产生悬空引用。
日志：logs/s64c_refscan.log
"""
import io, os, re

ROOT = r"D:/endometriosis_project/11_sc_eqtl_mr_project"
MEM = r"D:/endometriosis_project/.workbuddy/memory"
LOG = os.path.join(ROOT, "logs", "s64c_refscan.log")

SKIP_DIRS = {"00_data_raw", "02_data_processed", "figures", "_onek1k_plink",
             "_chr1_ld", "_panel", "_tools", "__pycache__", ".git"}
INCLUDE_EXT = (".md", ".txt", ".py", ".csv", ".json")

buf = []
def p(*a):
    s = " ".join(str(x) for x in a); buf.append(s)
    try: print(s)
    except UnicodeEncodeError: print(s.encode("ascii", "backslashreplace").decode("ascii"))

corpus = []
for base in (ROOT, MEM):
    for dp, dn, fn in os.walk(base):
        dn[:] = [d for d in dn if d not in SKIP_DIRS and not d.startswith("_fig")]
        for f in fn:
            if f.endswith(INCLUDE_EXT):
                fp = os.path.join(dp, f)
                try:
                    corpus.append((fp, io.open(fp, encoding="utf-8", errors="replace").read()))
                except Exception:
                    pass
p("扫到文本文件 = %d" % len(corpus))

CAND = [
    "28b_任务A_三位点比较_暂停点12.md",
    "28c_正文框架_审稿意见修正版.md",
    "28_论文正文初稿与图注定稿_暂停点12.md",
    "31_任务三3.2条件PheWAS_暂停点14.md",
    "31b_裁定1GTEx补充与任务B_MVMR_暂停点14.md",
]
p("")
p("== 候选文件被引用情况 ==")
for name in CAND:
    stem = name[:-3]                     # 去 .md
    who = []
    for fp, txt in corpus:
        if os.path.basename(fp) == name:
            continue
        c = txt.count(name) + txt.count(stem)
        if c:
            who.append((os.path.relpath(fp, r"D:/endometriosis_project"), c))
    p("  %-48s ref=%d" % (name, sum(c for _, c in who)))
    for w, c in who:
        p("        <- %s (x%d)" % (w, c))

p("")
p("== 全根目录「同前缀且同暂停点标签」的硬冲突 ==")
files = [f for f in os.listdir(ROOT) if os.path.isfile(os.path.join(ROOT, f))]
from collections import defaultdict
g = defaultdict(list)
for f in files:
    m = re.match(r"^(\d+[a-z]?)_", f)
    if not m:
        continue
    mp = re.search(r"暂停点(\d+)", f)
    g[(m.group(1), mp.group(1) if mp else "-")].append(f)
hard = {k: v for k, v in g.items() if len(v) > 1}
for k in sorted(hard, key=lambda x: (len(x[0]), x[0])):
    p("  prefix=%-5s 暂停点=%-3s → %d 份：%s" % (k[0], k[1], len(hard[k]), " | ".join(hard[k])))
p("硬冲突组数 = %d" % len(hard))
io.open(LOG, "w", encoding="utf-8", newline="\n").write("\n".join(buf) + "\n")
print("REFSCAN_DONE")
