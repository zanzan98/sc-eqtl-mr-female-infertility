# -*- coding: utf-8 -*-
"""s43c_run_smr.py — 任务三 3.3 步骤 2/3：SMR + HEIDI 主分析 + 敏感性，并汇总。

输入（均已通过闸门验证）：
  LD 参考   : 00_data_raw/smr_3p3/ld_ref/chr1_20_25mb.{bed,bim,fam}  (980 人 / 9,623 变异, GRCh37 ID)
  GWAS      : 00_data_raw/smr_3p3/gwas_GCST90483463.cojo.txt        (4,232 SNP, GCTA-COJI 7 列)
  eQTL BESD : 00_data_raw/smr_3p3/besd/<cell>.{besd,esi,epi}         (14 细胞类型 × 2 probe)

分析：
  [main] 默认参数（--peqtl-smr 5e-8, --peqtl-heidi 1.57e-3, --ld-upper-limit 0.9,
         --ld-lower-limit 0.05, --heidi-min-m 3, --heidi-max-m 20, --heidi-mtd 1, --cis-wind 2000）
  [sens] 仅把 --peqtl-smr 放宽到 1e-5（覆盖主分析因 cis-eQTL 未达 5e-8 而被跳过的 probe）

输出：
  00_data_raw/smr_3p3/smr_out/<cell>.smr        [main]
  00_data_raw/smr_3p3/smr_out_sens/<cell>.smr   [sens]
  tables/45_smr_3p3_results.csv                 [main 汇总]
  tables/45b_smr_3p3_results_sens.csv           [sens 汇总]
  日志 scripts/_s43c_run_smr.log
"""
import os, io, sys, glob, subprocess, ctypes, csv

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

PROJ = r"D:\endometriosis_project\11_sc_eqtl_mr_project"
SMR3 = os.path.join(PROJ, "00_data_raw", "smr_3p3")
BESD = os.path.join(SMR3, "besd")
LD = os.path.join(SMR3, "ld_ref", "chr1_20_25mb")
COJO = os.path.join(SMR3, "gwas_GCST90483463.cojo.txt")
OUT_MAIN = os.path.join(SMR3, "smr_out")
OUT_SENS = os.path.join(SMR3, "smr_out_sens")
TAB = os.path.join(PROJ, "tables")
SMR_EXE = r"D:\endometriosis_project\_tools\smr-1.3.1-win-x86_64\smr-1.3.1-win-x86_64\smr-1.3.1-win.exe"
LOG = os.path.join(PROJ, "scripts", "_s43c_run_smr.log")

buf = []
def w(s=""):
    buf.append(str(s))
    print(s)

def flush():
    os.makedirs(os.path.dirname(LOG), exist_ok=True)
    with io.open(LOG, "w", encoding="utf-8") as f:
        f.write("\n".join(buf) + "\n")

class MS(ctypes.Structure):
    _fields_ = [("dwLength", ctypes.c_ulong), ("dwMemoryLoad", ctypes.c_ulong),
                ("ullTotalPhys", ctypes.c_ulonglong), ("ullAvailPhys", ctypes.c_ulonglong),
                ("ullTotalPageFile", ctypes.c_ulonglong), ("ullAvailPageFile", ctypes.c_ulonglong),
                ("ullTotalVirtual", ctypes.c_ulonglong), ("ullAvailVirtual", ctypes.c_ulonglong),
                ("ullAvailExtendedVirtual", ctypes.c_ulonglong)]

def mem():
    st = MS(); st.dwLength = ctypes.sizeof(MS)
    ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(st))
    return (st.ullAvailPhys / 1048576.0, st.dwMemoryLoad)

def run(tag, outdir, extra, cells):
    os.makedirs(outdir, exist_ok=True)
    res = []
    for i, cell in enumerate(cells, 1):
        o = os.path.join(outdir, cell)
        cmd = [SMR_EXE, "--bfile", LD, "--gwas-summary", COJO,
               "--beqtl-summary", os.path.join(BESD, cell), "--out", o,
               "--thread-num", "4"] + extra
        p = subprocess.run(cmd, capture_output=True, text=True,
                           encoding="utf-8", errors="replace")
        so = (p.stdout or "") + "\n" + (p.stderr or "")
        io.open(os.path.join(outdir, cell + ".console.txt"), "w", encoding="utf-8").write(so)
        # 抽取关键告警行
        warns = [ln.strip() for ln in so.splitlines()
                 if ("no cis-eQTL" in ln.lower() or "skipped" in ln.lower()
                     or "error" in ln.lower() or "warning" in ln.lower()
                     or "excluded from the analysis" in ln.lower()
                     or "probes (with at least" in ln.lower())]
        nrow = 0
        sp = o + ".smr"
        if os.path.exists(sp):
            with io.open(sp, "r", encoding="utf-8", errors="replace") as f:
                nrow = max(0, len([l for l in f.read().splitlines() if l.strip()]) - 1)
        av, load = mem()
        w("[%s %2d/%d] %-12s rc=%d rows=%d  内存可用=%.0fMB(占%d%%)"
          % (tag, i, len(cells), cell, p.returncode, nrow, av, load))
        for x in warns:
            w("        ! %s" % x[:200])
        res.append((cell, p.returncode, nrow))
    return res

def parse(outdir, tag, cells):
    rows = []
    for cell in cells:
        sp = os.path.join(outdir, cell + ".smr")
        if not os.path.exists(sp):
            continue
        with io.open(sp, "r", encoding="utf-8", errors="replace") as f:
            hdr = f.readline().rstrip("\n").split("\t")
            for ln in f:
                if not ln.strip():
                    continue
                p = ln.rstrip("\n").split("\t")
                if len(p) != len(hdr):
                    w("  !! 列数不符 %s %s: %d vs %d" % (tag, cell, len(p), len(hdr)))
                    continue
                d = dict(zip(hdr, p))
                d["cell_type"] = cell
                d["analysis"] = tag
                rows.append(d)
    return rows

cells = sorted(os.path.basename(p).replace(".query.txt", "")
               for p in glob.glob(os.path.join(SMR3, "eqtl", "*.query.txt")))
av, load = mem()
w("=" * 90)
w("s43c_run_smr.py — SMR + HEIDI（任务三 3.3）")
w("cells (%d): %s" % (len(cells), ", ".join(cells)))
w("起始内存: 可用 %.0f MB (占 %d%%)" % (av, load))
w("=" * 90)

w()
w("### [main] 默认参数（--peqtl-smr 5e-8）")
resM = run("main", OUT_MAIN, [], cells)
w()
w("### [sens] 敏感性（--peqtl-smr 1e-5）")
resS = run("sens", OUT_SENS, ["--peqtl-smr", "1e-5"], cells)

rowsM = parse(OUT_MAIN, "main", cells)
rowsS = parse(OUT_SENS, "sens", cells)

os.makedirs(TAB, exist_ok=True)
COLS = ["cell_type", "analysis", "probeID", "ProbeChr", "Gene", "Probe_bp", "topSNP",
        "topSNP_chr", "topSNP_bp", "A1", "A2", "Freq", "b_GWAS", "se_GWAS", "p_GWAS",
        "b_eQTL", "se_eQTL", "p_eQTL", "b_SMR", "se_SMR", "p_SMR", "p_HEIDI", "nsnp_HEIDI"]

def dump(path, rows):
    with io.open(path, "w", encoding="utf-8", newline="") as f:
        wr = csv.writer(f)
        wr.writerow(COLS)
        for r in rows:
            wr.writerow([r.get(c, "NA") for c in COLS])
    w("  写出 %s (%d 行)" % (path, len(rows)))

w()
dump(os.path.join(TAB, "45_smr_3p3_results.csv"), rowsM)
dump(os.path.join(TAB, "45b_smr_3p3_results_sens.csv"), rowsS)
dump(os.path.join(TAB, "45c_smr_3p3_combined.csv"), rowsM + rowsS)

# ---- 核验：主分析每个细胞是否两 probe 都在；缺失者列出
w()
w("### 主分析（5e-8）逐细胞 probe 覆盖")
for cell in cells:
    got = sorted(r["Gene"] for r in rowsM if r["cell_type"] == cell)
    miss = sorted(set(["CDC42", "LINC00339"]) - set(got))
    w("  %-12s 分析到 %s%s" % (cell, ",".join(got) if got else "(无)",
                              ("  缺失: " + ",".join(miss)) if miss else ""))

w()
w("### 主分析结果（P_SMR / P_HEIDI）")
w("%-12s %-10s %-14s %10s %10s %10s %6s" % ("cell", "Gene", "topSNP", "b_SMR", "p_SMR", "p_HEIDI", "m"))
for cell in cells:
    for r in rowsM:
        if r["cell_type"] != cell:
            continue
        w("%-12s %-10s %-14s %10s %10s %10s %6s"
          % (cell, r["Gene"], r["topSNP"], r["b_SMR"], r["p_SMR"], r["p_HEIDI"], r["nsnp_HEIDI"]))

w()
w("### 敏感性（1e-5）结果")
w("%-12s %-10s %-14s %10s %10s %10s %6s" % ("cell", "Gene", "topSNP", "b_SMR", "p_SMR", "p_HEIDI", "m"))
for cell in cells:
    for r in rowsS:
        if r["cell_type"] != cell:
            continue
        w("%-12s %-10s %-14s %10s %10s %10s %6s"
          % (cell, r["Gene"], r["topSNP"], r["b_SMR"], r["p_SMR"], r["p_HEIDI"], r["nsnp_HEIDI"]))

av, load = mem()
w()
w("VERDICT = %s" % ("PASS" if len(rowsM) >= 1 and len(rowsS) >= 1 else "FAIL"))
w("结束内存: 可用 %.0f MB (占 %d%%)" % (av, load))
flush()
