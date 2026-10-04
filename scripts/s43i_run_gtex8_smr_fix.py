# -*- coding: utf-8 -*-
"""s43i_run_gtex8_smr_fix.py — 选项 A 正式运行：GTEx v8 全血（bulk）3 基因 SMR + HEIDI
（使用 s43h 修正 Freq 取向后的 query_fix）

暴露：GTEx v8 全血（eQTL Catalogue imported，n=670，bulk 全血 —— 非单细胞！）
结局：GCST90483463（女性不孕症，EUR，40,024 病例 / 665,658 对照）
LD 参考：OneK1K 980 供者子集 chr1:20–25 Mb（hg19；与 GTEx 供者不重叠）

参数：
  [A] default  --peqtl-smr 5e-8   （预期：WNT4 无 SNP 过阈 → 不检验；CDC42/LINC00339 可检验）
  [B] relaxed  --peqtl-smr 1e-4   （WNT4 唯一可行设置：minP = 5.526e-05）

输出：00_data_raw/smr_3p3/gtex8/{besd_fix,smr_fix_A,smr_fix_B}/...
      tables/47_gtex8_wb_smr_results.csv
"""
import os, io, sys, subprocess, ctypes, csv
import numpy as np
import pandas as pd

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

PROJ = r"D:\endometriosis_project"
SMR3 = os.path.join(PROJ, "11_sc_eqtl_mr_project", "00_data_raw", "smr_3p3")
G8 = os.path.join(SMR3, "gtex8")
QDIR = os.path.join(G8, "query_fix")
BESD = os.path.join(G8, "besd_fix")
LD = os.path.join(SMR3, "ld_ref", "chr1_20_25mb")
COJO = os.path.join(SMR3, "gwas_GCST90483463.cojo.txt")
TAB = os.path.join(PROJ, "11_sc_eqtl_mr_project", "tables")
EXE = os.path.join(PROJ, "_tools", "smr-1.3.1-win-x86_64",
                   "smr-1.3.1-win-x86_64", "smr-1.3.1-win.exe")
LOG = os.path.join(PROJ, "11_sc_eqtl_mr_project", "scripts", "_s43i_run_gtex8_smr_fix.log")

buf = []
def w(s=""):
    buf.append(str(s)); print(s)
def flush():
    io.open(LOG, "w", encoding="utf-8").write("\n".join(buf) + "\n")

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
GENES = ["WNT4", "CDC42", "LINC00339"]

w("=" * 100)
w("s43i — GTEx v8 全血（bulk）SMR + HEIDI（任务三 3.3 · 选项 A · 修正 Freq 取向后正式运行）")
w("=" * 100)
av, load = mem(); w("起始内存: 可用 %.0f MB（占 %d%%）" % (av, load))

# ---- N
src = os.path.join(PROJ, "01_data_raw", "eqtlcat_gtex_v8",
                   "Whole_Blood.chr1_21_23.3Mb_3genes.tsv")
an_vals = []
with io.open(src, "r", encoding="utf-8", errors="replace") as f:
    ix = {c: i for i, c in enumerate(f.readline().rstrip("\n").split("\t"))}
    for ln in f:
        p = ln.rstrip("\n").split("\t")
        if len(p) > ix["an"]:
            try:
                an_vals.append(int(p[ix["an"]]))
            except ValueError:
                pass
an_arr = np.array(an_vals)
N_GTEX = int(round(np.median(an_arr) / 2.0))
w("[0] an 中位数 = %d → GTEx 全血样本量 N = an/2 = %d" % (int(np.median(an_arr)), N_GTEX))

# ---- BESD
w()
w("[1] 建 BESD（--make-besd --add-n %d）" % N_GTEX)
for g in GENES:
    o = os.path.join(BESD, g)
    p = subprocess.run([EXE, "--qfile", os.path.join(QDIR, "%s.query.txt" % g),
                        "--make-besd", "--out", o, "--add-n", str(N_GTEX)],
                       capture_output=True, text=True, encoding="utf-8", errors="replace", cwd=BESD)
    ok = all(os.path.exists(o + e) for e in (".besd", ".epi", ".esi"))
    w("    %-10s rc=%d  生成=%s" % (g, p.returncode, ok))
    io.open(os.path.join(BESD, g + ".console.txt"), "w", encoding="utf-8").write(
        (p.stdout or "") + "\n" + (p.stderr or ""))
    sp = o + ".summary"
    if os.path.exists(sp):
        txt = io.open(sp, "r", encoding="utf-8", errors="replace").read().splitlines()
        w("        .summary: %s" % " || ".join([l for l in txt if l.startswith("{")]))

# ---- SMR
def run_smr(tag, extra, outname):
    outdir = os.path.join(G8, outname)
    os.makedirs(outdir, exist_ok=True)
    res = []
    for g in GENES:
        o = os.path.join(outdir, g)
        p = subprocess.run([EXE, "--bfile", LD, "--gwas-summary", COJO,
                            "--beqtl-summary", os.path.join(BESD, g), "--out", o,
                            "--thread-num", "4"] + extra,
                           capture_output=True, text=True, encoding="utf-8", errors="replace")
        so = (p.stdout or "") + "\n" + (p.stderr or "")
        io.open(os.path.join(outdir, "%s.console.txt" % g), "w", encoding="utf-8").write(so)
        key = [l.strip() for l in so.splitlines()
               if ("to be included after allele checking" in l
                   or "with at least a cis-eQTL" in l
                   or "allele frequency differences" in l
                   or "SNPs are excluded" in l
                   or "ERROR" in l
                   or "no SNP" in l.lower())]
        nrow = 0
        sp = o + ".smr"
        if os.path.exists(sp):
            nrow = max(0, len([l for l in io.open(sp, "r", encoding="utf-8",
                                                   errors="replace").read().splitlines()
                               if l.strip()]) - 1)
        av2, l2 = mem()
        w("[%s] %-10s rc=%d rows=%d 内存可用=%.0fMB(占%d%%)" % (tag, g, p.returncode, nrow, av2, l2))
        for k in key:
            w("        %s" % k[:200])
        res.append((g, p.returncode, nrow))
    return res

w()
w("=" * 100); w("[A] 默认参数（--peqtl-smr 5e-8）"); w("=" * 100)
run_smr("A", [], "smr_fix_A")
w()
w("=" * 100); w("[B] 放宽 target 阈值（--peqtl-smr 1e-4）—— WNT4 唯一可行设置"); w("=" * 100)
run_smr("B", ["--peqtl-smr", "1e-4"], "smr_fix_B")

# ---- 汇总
COLS = ["Gene", "analysis", "probeID", "Probe_bp", "topSNP", "topSNP_bp", "A1", "A2", "Freq",
        "b_GWAS", "se_GWAS", "p_GWAS", "b_eQTL", "se_eQTL", "p_eQTL",
        "b_SMR", "se_SMR", "p_SMR", "p_HEIDI", "nsnp_HEIDI"]
rows = []
for tag, outname in (("A(5e-8)", "smr_fix_A"), ("B(1e-4)", "smr_fix_B")):
    for g in GENES:
        sp = os.path.join(G8, outname, g + ".smr")
        if not os.path.exists(sp):
            continue
        with io.open(sp, "r", encoding="utf-8", errors="replace") as f:
            hdr = f.readline().rstrip("\n").split("\t")
            for ln in f:
                if not ln.strip():
                    continue
                p = ln.rstrip("\n").split("\t")
                if len(p) != len(hdr):
                    continue
                d = dict(zip(hdr, p)); d["analysis"] = tag; d["Gene"] = g
                rows.append(d)
os.makedirs(TAB, exist_ok=True)
with io.open(os.path.join(TAB, "47_gtex8_wb_smr_results.csv"), "w",
             encoding="utf-8", newline="") as f:
    wr = csv.writer(f); wr.writerow(COLS)
    for r in rows:
        wr.writerow([r.get(c, "NA") for c in COLS])

w()
w("=" * 100)
w("结果（GTEx v8 全血 bulk n=%d；结局 GCST90483463 女性不孕症）" % N_GTEX)
w("=" * 100)
w("%-9s %-10s %-14s %10s %11s %11s %4s" % ("analysis", "Gene", "topSNP", "b_SMR", "p_SMR", "p_HEIDI", "m"))
for r in rows:
    try:
        bs = "%.5g" % float(r["b_SMR"])
    except (ValueError, TypeError, KeyError):
        bs = "NA"
    w("%-9s %-10s %-14s %10s %11s %11s %4s"
      % (r["analysis"], r["Gene"], r.get("topSNP", "NA"), bs,
         r.get("p_SMR", "NA"), r.get("p_HEIDI", "NA"), r.get("nsnp_HEIDI", "NA")))

av3, l3 = mem()
w()
w("写出 tables/47_gtex8_wb_smr_results.csv（%d 行）" % len(rows))
w("结束内存: 可用 %.0f MB（占 %d%%）" % (av3, l3))
w("VERDICT = %s" % ("PASS" if len(rows) > 0 else "FAIL"))
flush()
