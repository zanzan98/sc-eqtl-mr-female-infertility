# -*- coding: utf-8 -*-
"""s43b_build_besd.py — 任务三 3.3 步骤 1：把 14 个 per-cell eQTL query 文件转成 SMR BESD 格式。

依据 SMR 1.3.1 官方手册（D:\\endometriosis_project\\_smr_manual.txt）：
  --qfile <query>  --make-besd --out <prefix> [--add-n N]
  query 格式 = SNP Chr BP A1 A2 Freq Probe Probe_Chr Probe_bp Gene Orientation b se p
输出 <prefix>.besd / .esi / .epi

样本量 N：LD 参考与 OneK1K 各细胞类型同源（980 供者）。eQTL N 在单 SNP SMR 中不被使用，
此处仅为 .esi 元数据完备性，透明取 980（已在本日志中标注）。
"""
import os, io, sys, glob, subprocess, json, ctypes

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

PROJ = r"D:\endometriosis_project\11_sc_eqtl_mr_project"
SMR3 = os.path.join(PROJ, "00_data_raw", "smr_3p3")
EQTL = os.path.join(SMR3, "eqtl")
BESD = os.path.join(SMR3, "besd")
SMR_EXE = r"D:\endometriosis_project\_tools\smr-1.3.1-win-x86_64\smr-1.3.1-win-x86_64\smr-1.3.1-win.exe"
LOG = os.path.join(PROJ, "scripts", "_s43b_build_besd.log")

N_SAMPLE = 980

buf = []
def w(s=""):
    buf.append(str(s))
    print(s)

def flush():
    os.makedirs(os.path.dirname(LOG), exist_ok=True)
    with io.open(LOG, "w", encoding="utf-8") as f:
        f.write("\n".join(buf) + "\n")

class MEMSTATUSEX(ctypes.Structure):
    _fields_ = [("dwLength", ctypes.c_ulong), ("dwMemoryLoad", ctypes.c_ulong),
                ("ullTotalPhys", ctypes.c_ulonglong), ("ullAvailPhys", ctypes.c_ulonglong),
                ("ullTotalPageFile", ctypes.c_ulonglong), ("ullAvailPageFile", ctypes.c_ulonglong),
                ("ullTotalVirtual", ctypes.c_ulonglong), ("ullAvailVirtual", ctypes.c_ulonglong),
                ("ullAvailExtendedVirtual", ctypes.c_ulonglong)]

def mem():
    st = MEMSTATUSEX(); st.dwLength = ctypes.sizeof(MEMSTATUSEX)
    ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(st))
    return (st.ullAvailPhys / 1048576.0, st.dwMemoryLoad, st.ullTotalPhys / 1048576.0)

os.makedirs(BESD, exist_ok=True)

w("=" * 78)
w("s43b_build_besd.py — SMR BESD 构建（任务三 3.3 步骤 1）")
w("SMR_EXE = %s" % SMR_EXE)
w("N_SAMPLE (输入 .esi 的样本量，取 LD 参考供者数) = %d" % N_SAMPLE)
av, load, tot = mem()
w("起始内存: 可用 %.0f MB / 总 %.0f MB (占用 %d%%)" % (av, tot, load))
w("=" * 78)

cells = sorted(os.path.basename(p).replace(".query.txt", "")
               for p in glob.glob(os.path.join(EQTL, "*.query.txt")))
w("细胞类型 (%d): %s" % (len(cells), ", ".join(cells)))
w()

summary = []
for i, cell in enumerate(cells, 1):
    qf = os.path.join(EQTL, "%s.query.txt" % cell)
    out = os.path.join(BESD, cell)
    if os.path.exists(out + ".besd") and os.path.exists(out + ".esi") and os.path.exists(out + ".epi"):
        w("[%2d/%d] %-12s 已存在，跳过" % (i, len(cells), cell))
        summary.append((cell, "exists", 0, ""))
        continue
    cmd = [SMR_EXE, "--qfile", qf, "--make-besd", "--out", out, "--add-n", str(N_SAMPLE)]
    p = subprocess.run(cmd, capture_output=True, text=True,
                       encoding="utf-8", errors="replace", cwd=BESD)
    so = (p.stdout or "").strip().replace("\r\n", " | ")
    se = (p.stderr or "").strip().replace("\r\n", " | ")
    ok = all(os.path.exists(out + e) for e in (".besd", ".esi", ".epi"))
    sz = os.path.getsize(out + ".besd") if ok else 0
    w("[%2d/%d] %-12s rc=%d  besd=%s (%.1f KB)" % (i, len(cells), cell, p.returncode, ok, sz / 1024.0))
    if so:
        w("      stdout: %s" % so[:400])
    if se:
        w("      stderr: %s" % se[:400])
    summary.append((cell, "ok" if ok else "FAIL", sz, (so + " " + se)[:200]))

w()
w("-" * 78)
w("汇总: %d/%d 成功" % (sum(1 for s in summary if s[1] in ("ok", "exists")), len(cells)))
for s in summary:
    w("  %-12s %-8s %s" % (s[0], s[1], ("%.1f KB" % (s[2] / 1024.0)) if s[2] else ""))

# 核验 .epi 行数与 .esi 内容
w()
w("核验 .epi / .esi（取 B_MEM 与前 2 个细胞）")
for cell in cells[:3]:
    epi = os.path.join(BESD, cell + ".epi")
    esi = os.path.join(BESD, cell + ".esi")
    if os.path.exists(epi):
        with io.open(epi, "r", encoding="utf-8", errors="replace") as f:
            lines = f.read().splitlines()
        w("  %s.epi (%d 行): %s" % (cell, len(lines), " || ".join(lines[:4])))
    if os.path.exists(esi):
        with io.open(esi, "r", encoding="utf-8", errors="replace") as f:
            lines = f.read().splitlines()
        w("  %s.esi (%d 行): %s" % (cell, len(lines), " || ".join(lines)))

av2, load2, _ = mem()
w()
w("VERDICT = %s" % ("PASS" if all(s[1] in ("ok", "exists") for s in summary) else "FAIL"))
w("结束内存: 可用 %.0f MB (占用 %d%%)" % (av2, load2))
flush()
