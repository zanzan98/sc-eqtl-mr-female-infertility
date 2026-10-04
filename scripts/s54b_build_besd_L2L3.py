# -*- coding: utf-8 -*-
"""s54b_build_besd_L2L3.py — 任务A：L2/L3 SMR BESD 构建（复刻 s43b）

--qfile <query> --make-besd --out <prefix> --add-n <N>
N 按细胞类型真实有效样本量给出（CD4_NC=980, Mono_NC=690），不再一律取 980。
（注：单 SNP SMR 统计量不使用 eQTL N；N 仅进 .esi 元数据。）
"""
import os, io, sys, glob, subprocess, ctypes

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

PROJ = r"D:\endometriosis_project\11_sc_eqtl_mr_project"
SMRD = os.path.join(PROJ, "00_data_raw", "smr_L2L3")
EQTL = os.path.join(SMRD, "eqtl")
BESD = os.path.join(SMRD, "besd")
SMR_EXE = r"D:\endometriosis_project\_tools\smr-1.3.1-win-x86_64\smr-1.3.1-win-x86_64\smr-1.3.1-win.exe"
LOG = r"D:\endometriosis_project\_s54b_build_besd_L2L3.log"

NEFF = {"L2_YME1L1_CD4_NC": 980, "L3_ANXA4_Mono_NC": 690}

buf = []
def w(s=""):
    buf.append(str(s)); print(s)

def flush():
    io.open(LOG, "w", encoding="utf-8", newline="\n").write("\n".join(buf) + "\n")

class MS(ctypes.Structure):
    _fields_ = [("dwLength", ctypes.c_ulong), ("dwMemoryLoad", ctypes.c_ulong),
                ("ullTotalPhys", ctypes.c_ulonglong), ("ullAvailPhys", ctypes.c_ulonglong),
                ("ullTotalPageFile", ctypes.c_ulonglong), ("ullAvailPageFile", ctypes.c_ulonglong),
                ("ullTotalVirtual", ctypes.c_ulonglong), ("ullAvailVirtual", ctypes.c_ulonglong),
                ("ullAvailExtendedVirtual", ctypes.c_ulonglong)]

def mem():
    st = MS(); st.dwLength = ctypes.sizeof(MS)
    ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(st))
    return st.ullAvailPhys / 1048576.0, st.dwMemoryLoad

os.makedirs(BESD, exist_ok=True)
av, load = mem()
w("=" * 78)
w("s54b — L2/L3 SMR BESD 构建（修订版 N：按细胞类型真实 N）")
w("起始内存: 可用 %.0f MB (占 %d%%)" % (av, load))
w("=" * 78)

tags = sorted(os.path.basename(p).replace(".query.txt", "")
              for p in glob.glob(os.path.join(EQTL, "*.query.txt")))
w("数据集 (%d): %s" % (len(tags), ", ".join(tags)))
w()

summary = []
for i, tag in enumerate(tags, 1):
    qf = os.path.join(EQTL, "%s.query.txt" % tag)
    out = os.path.join(BESD, tag)
    n = NEFF.get(tag, 980)
    if all(os.path.exists(out + e) for e in (".besd", ".esi", ".epi")):
        w("[%d/%d] %-20s 已存在，跳过" % (i, len(tags), tag))
        summary.append((tag, "exists", 0, n)); continue
    cmd = [SMR_EXE, "--qfile", qf, "--make-besd", "--out", out, "--add-n", str(n)]
    p = subprocess.run(cmd, capture_output=True, text=True, encoding="utf-8",
                       errors="replace", cwd=BESD)
    so = (p.stdout or "").strip().replace("\r\n", " | ")
    se = (p.stderr or "").strip().replace("\r\n", " | ")
    ok = all(os.path.exists(out + e) for e in (".besd", ".esi", ".epi"))
    sz = os.path.getsize(out + ".besd") if ok else 0
    w("[%d/%d] %-20s N=%-4d rc=%d besd=%s (%.1f KB)" % (i, len(tags), tag, n, p.returncode, ok, sz / 1024.0))
    if so: w("      stdout: %s" % so[:300])
    if se: w("      stderr: %s" % se[:300])
    summary.append((tag, "ok" if ok else "FAIL", sz, n))

w()
w("-" * 78)
w("汇总: %d/%d 成功" % (sum(1 for s in summary if s[1] in ("ok", "exists")), len(tags)))
for tag in tags:
    epi = os.path.join(BESD, tag + ".epi"); esi = os.path.join(BESD, tag + ".esi")
    if os.path.exists(epi):
        ls = io.open(epi, encoding="utf-8", errors="replace").read().splitlines()
        w("  %s.epi (%d 行): %s" % (tag, len(ls), " || ".join(ls[:3])))
    if os.path.exists(esi):
        ls = io.open(esi, encoding="utf-8", errors="replace").read().splitlines()
        w("  %s.esi: %s" % (tag, " || ".join(ls)))
av2, load2 = mem()
w(); w("VERDICT = %s" % ("PASS" if all(s[1] in ("ok", "exists") for s in summary) else "FAIL"))
w("结束内存: 可用 %.0f MB (占 %d%%)" % (av2, load2))
flush()
