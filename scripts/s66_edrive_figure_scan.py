# -*- coding: utf-8 -*-
"""
s66_edrive_figure_scan.py  ——  E 盘图件资产检索 + 本项目快照差异比对（只读）

用途：核查「补充图 Fig S1/S2 是否存在于 E 盘」，并给出一份可复核的图件清单。
产出：
  tables/63b_Edrive_images.csv       E 盘全部图件/矢量文件（path, bytes, mtime）
  tables/63b_Edrive_supp_named.csv   E 盘文件名命中「补充图式命名」的文件（含非图件）
  logs/_s66_edrive_scan.txt          运行摘要（判定用）
纪律：
  * 只读；不移动、不改名、不删除任何文件（外置盘，用户数据）。
  * 一律 Python io.open(..., encoding='utf-8', newline='\\n') 写盘；写完核 size。
  * 匹配「补充图」必须用精确正则 + 只匹配文件名；短 token（som/supp）会命中
    Eigen 头文件里的 Support/Symmetry（实测 443 条噪声）。
"""
import io
import os
import re
import sys
import traceback
from collections import Counter, defaultdict
from datetime import datetime

PRJ = "D:/endometriosis_project/11_sc_eqtl_mr_project"
OUT_CSV_IMG = os.path.join(PRJ, "tables", "63b_Edrive_images.csv")
OUT_CSV_SUPP = os.path.join(PRJ, "tables", "63b_Edrive_supp_named.csv")
OUT_LOG = os.path.join(PRJ, "logs", "_s66_edrive_scan.txt")

TARGET = "E:\\"
OURS_ON_E = os.path.join(TARGET, "生信2", "endometriosis_project", "11_sc_eqtl_mr_project")
OURS_ON_D = PRJ

IMG_EXT = {".png", ".pdf", ".svg", ".tif", ".tiff", ".jpg", ".jpeg", ".eps", ".ai", ".psd", ".bmp", ".gif"}
SKIP_DIRS = {"$RECYCLE.BIN", "System Volume Information", "$WinREAgent"}

SUPP_FILE_RE = re.compile(
    r"(?<![a-z0-9])("
    r"fig(ure)?[\s_\-]*s\d*"                       # Fig S1 / Figure S1 / fig_s1 / FigS
    r"|figs\d*"                                     # figs1 / figs
    r"|s(upp)?[\s_\-]*fig(ure)?"                    # supp fig
    r"|supp(lementary|lement)?[\s_\-]*(fig|figure|data|material)"
    r"|补充图|附图|附图件|补充材料"
    r"|fig[\s_\-]*s[12]\b"
    r")", re.I)

buf = []
def P(s=""):
    buf.append(str(s))

def human(n):
    if n is None or n < 0:
        return "?"
    for u in ("B", "KB", "MB", "GB"):
        if n < 1024:
            return "%.1f%s" % (n, u)
        n /= 1024.0
    return "%.1fTB" % n

def snap(d):
    m = {}
    if os.path.isdir(d):
        for f in os.listdir(d):
            fp = os.path.join(d, f)
            if os.path.isfile(fp):
                m[f] = os.path.getsize(fp)
    return m

def main():
    P("s66 @ %s" % datetime.now().isoformat(timespec="seconds"))
    P("target = %s" % TARGET)
    if not os.path.exists(TARGET):
        P("VERDICT: E: NOT PRESENT -> abort")
        return
    P("E: present = True")
    P("")

    images, supp = [], []
    n_files = n_dirs = 0
    for dirpath, dirnames, filenames in os.walk(TARGET, onerror=lambda e: None):
        dirnames[:] = [d for d in dirnames if d not in SKIP_DIRS]
        n_dirs += len(dirnames)
        for fn in filenames:
            n_files += 1
            fp = os.path.join(dirpath, fn)
            ext = os.path.splitext(fn.lower())[1]
            if ext in IMG_EXT:
                try:
                    st = os.stat(fp)
                    images.append((fp, st.st_size,
                                   datetime.fromtimestamp(st.st_mtime).strftime("%Y-%m-%d %H:%M")))
                except Exception:
                    images.append((fp, -1, "?"))
            if SUPP_FILE_RE.search(fn):
                try:
                    sz = os.path.getsize(fp)
                except Exception:
                    sz = -1
                supp.append((fp, sz))

    P("=== SUMMARY ===")
    P("dirs=%d files=%d image_or_vector=%d supp_named=%d" % (n_dirs, n_files, len(images), len(supp)))
    P("")

    P("=== SUPP-NAMED FILES GROUPED BY TOP-2 PATH ===")
    grp = defaultdict(list)
    for fp, sz in supp:
        parts = fp.split(os.sep)
        key = os.sep.join(parts[:4]) if len(parts) >= 4 else fp
        grp[key].append((fp, sz))
    for k in sorted(grp):
        P("  [%d] %s" % (len(grp[k]), k))
    P("")

    P("=== OUR PROJECT SNAPSHOT ON E: (figures/) ===")
    se, sd = snap(os.path.join(OURS_ON_E, "figures")), snap(os.path.join(OURS_ON_D, "figures"))
    P("E: figures files = %d ; D: figures files = %d" % (len(se), len(sd)))
    P("only on E: %s" % (sorted(set(se) - set(sd)) or "none"))
    P("only on D: %s" % (sorted(set(sd) - set(se)) or "none"))
    diff = [f for f in sorted(set(se) & set(sd)) if se[f] != sd[f]]
    P("same name / different size = %d : %s" % (len(diff), diff))
    figs = [f for f in se if re.search(r"figs\d", f, re.I)]
    P("figures named FigS* on E: copy = %s" % (figs or "NONE"))
    P("")

    P("=== OUR PROJECT SNAPSHOT ON E: (tables/) ===")
    te, td = snap(os.path.join(OURS_ON_E, "tables")), snap(os.path.join(OURS_ON_D, "tables"))
    P("E: tables = %d ; D: tables = %d" % (len(te), len(td)))
    P("only on E: %s" % (sorted(set(te) - set(td)) or "none"))
    P("only on D count = %d" % len(set(td) - set(te)))
    P("")

    with io.open(OUT_CSV_IMG, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("path,bytes,mtime\n")
        for fp, sz, mt in sorted(images):
            fh.write('"%s",%d,%s\n' % (fp.replace('"', '""'), sz, mt))
    with io.open(OUT_CSV_SUPP, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("path,bytes\n")
        for fp, sz in sorted(supp):
            fh.write('"%s",%d\n' % (fp.replace('"', '""'), sz))
    with io.open(OUT_LOG, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(buf) + "\n")
    P("wrote: %s (%d B)" % (OUT_CSV_IMG, os.path.getsize(OUT_CSV_IMG)))
    P("wrote: %s (%d B)" % (OUT_CSV_SUPP, os.path.getsize(OUT_CSV_SUPP)))
    P("wrote: %s (%d B)" % (OUT_LOG, os.path.getsize(OUT_LOG)))
    with io.open(OUT_LOG, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(buf) + "\n")

try:
    main()
except BaseException:
    buf.append("FATAL:\n" + traceback.format_exc())
    with io.open(OUT_LOG, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(buf) + "\n")
    raise

sys.stdout.write("DONE\n")
