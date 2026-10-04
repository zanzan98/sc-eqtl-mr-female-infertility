# -*- coding: utf-8 -*-
"""s55c_run_gtex8_smr_L2L3.py -- 裁定1：GTEx v8 全血（bulk, n=670）L2/L3 规范 cis-SMR + HEIDI

严格复刻 s43i_run_gtex8_smr_fix.py 的参数口径：
  [A] default   --peqtl-smr 5e-8      （主分析）
  [B] relaxed   --peqtl-smr 1e-5      （敏感性；两基因主分析均已可行，此步仅留痕）
其余取 SMR 1.3.1 默认：--peqtl-heidi 1.57e-3 / --ld-upper-limit 0.9 / --ld-lower-limit 0.05
                        --heidi-min-m 3 / --heidi-max-m 20

LD 参考 = 同一位点的 OneK1K 980 供者 PLINK 子集（与 GTEx v8 供者不重叠 → 欧洲血统近似 LD，已知局限）
结局 GWAS = GCST90483463（与既有 chr1 GTEx 分析同一结局）

输出：00_data_raw/smr_L2L3/gtex8_besd/<TAG>.{besd,esi,epi}
      00_data_raw/smr_L2L3/gtex8_smr_out{,_sens}/<TAG>.{smr,console.txt}
      tables/59_gtex8_L2L3_smr_results.csv
"""
import ctypes
import csv
import io
import os
import subprocess
import sys

import numpy as np
import pandas as pd

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

PROJ = r"D:\endometriosis_project"
P11 = os.path.join(PROJ, "11_sc_eqtl_mr_project")
SMRD = os.path.join(P11, "00_data_raw", "smr_L2L3")
QDIR = os.path.join(SMRD, "gtex8_query")
BESD = os.path.join(SMRD, "gtex8_besd")
G8RAW = os.path.join(PROJ, "01_data_raw", "eqtlcat_gtex_v8")
TAB = os.path.join(P11, "tables")
EXE = os.path.join(PROJ, "_tools", "smr-1.3.1-win-x86_64",
                   "smr-1.3.1-win-x86_64", "smr-1.3.1-win.exe")
LOG = os.path.join(P11, "scripts", "_s55c_run_gtex8_smr_L2L3.log")
os.makedirs(BESD, exist_ok=True)

DS = [
    dict(tag="L2_YME1L1", gene="YME1L1", chrom="10",
         bfile=os.path.join(SMRD, "ld_ref", "L2_YME1L1_CD4_NC_cis"),
         cojo=os.path.join(SMRD, "gwas_L2_YME1L1_CD4_NC.cojo.txt"),
         src=os.path.join(G8RAW, "Whole_Blood.L2_YME1L1.tsv")),
    dict(tag="L3_ANXA4", gene="ANXA4", chrom="2",
         bfile=os.path.join(SMRD, "ld_ref", "L3_ANXA4_Mono_NC_cis"),
         cojo=os.path.join(SMRD, "gwas_L3_ANXA4_Mono_NC.cojo.txt"),
         src=os.path.join(G8RAW, "Whole_Blood.L3_ANXA4.tsv")),
]

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

w("=" * 104)
w("s55c — GTEx v8 全血（bulk）L2/L3 cis-SMR + HEIDI（裁定1）")
w("=" * 104)
av, load = mem(); w("起始内存: 可用 %.0f MB（占 %d%%）" % (av, load))

# ---- N（由源 TSV 的 an 中位数推） ----
NMAP = {}
for D in DS:
    an = []
    with io.open(D["src"], "r", encoding="utf-8", errors="replace") as f:
        ix = {c: i for i, c in enumerate(f.readline().rstrip("\n").split("\t"))}
        for ln in f:
            p = ln.rstrip("\n").split("\t")
            if len(p) <= ix["an"]:
                continue
            try:
                an.append(int(p[ix["an"]]))
            except ValueError:
                pass
    an = np.array(an)
    N = int(round(np.median(an) / 2.0))
    NMAP[D["tag"]] = N
    w("[0] %-11s an 中位=%d（唯一值 %d，范围 %d..%d）→ N = an/2 = %d"
      % (D["tag"], int(np.median(an)), len(set(an.tolist())), an.min(), an.max(), N))

# ---- BESD ----
w("")
w("[1] 建 BESD（--make-besd --add-n <N>）")
for D in DS:
    o = os.path.join(BESD, D["tag"])
    p = subprocess.run([EXE, "--qfile", os.path.join(QDIR, "%s.query.txt" % D["tag"]),
                        "--make-besd", "--out", o, "--add-n", str(NMAP[D["tag"]])],
                       capture_output=True, text=True, encoding="utf-8",
                       errors="replace", cwd=BESD)
    ok = all(os.path.exists(o + e) for e in (".besd", ".epi", ".esi"))
    io.open(o + ".console.txt", "w", encoding="utf-8").write(
        (p.stdout or "") + "\n" + (p.stderr or ""))
    w("    %-11s rc=%d 生成=%s" % (D["tag"], p.returncode, ok))
    sp = o + ".summary"
    if os.path.exists(sp):
        txt = io.open(sp, "r", encoding="utf-8", errors="replace").read().splitlines()
        w("        .summary: %s" % " || ".join([l for l in txt if l.startswith('{')]))
    epi = o + ".epi"
    if os.path.exists(epi):
        line = io.open(epi, "r", encoding="utf-8", errors="replace").read().splitlines()
        w("        .epi: %s" % (line[1] if len(line) > 1 else line))

# ---- SMR ----
def run_smr(tag, extra, outname):
    outdir = os.path.join(SMRD, outname)
    os.makedirs(outdir, exist_ok=True)
    for D in DS:
        o = os.path.join(outdir, D["tag"])
        cmd = [EXE, "--bfile", D["bfile"], "--gwas-summary", D["cojo"],
               "--beqtl-summary", os.path.join(BESD, D["tag"]), "--out", o,
               "--thread-num", "4"] + extra
        p = subprocess.run(cmd, capture_output=True, text=True,
                           encoding="utf-8", errors="replace")
        so = (p.stdout or "") + "\n" + (p.stderr or "")
        io.open(os.path.join(outdir, "%s.console.txt" % D["tag"]), "w",
                encoding="utf-8").write(so)
        key = [l.strip() for l in so.splitlines()
               if ("to be included after allele checking" in l
                   or "with at least a cis-eQTL" in l
                   or "allele frequency differences" in l
                   or "excluded" in l or "ERROR" in l or "Error" in l
                   or "no SNP" in l.lower() or "no probe" in l.lower())]
        nrow = 0
        sp = o + ".smr"
        if os.path.exists(sp):
            nrow = max(0, len([l for l in io.open(sp, "r", encoding="utf-8",
                                                   errors="replace").read().splitlines()
                               if l.strip()]) - 1)
        av2, l2 = mem()
        w("[%s] %-11s rc=%d rows=%d 内存可用=%.0fMB(占%d%%)"
          % (tag, D["tag"], p.returncode, nrow, av2, l2))
        for k in key:
            w("        %s" % k[:210])

w("")
w("=" * 104); w("[A] 主分析：--peqtl-smr 5e-8"); w("=" * 104)
run_smr("A", [], "gtex8_smr_out")
w("")
w("=" * 104); w("[B] 敏感性：--peqtl-smr 1e-5"); w("=" * 104)
run_smr("B", ["--peqtl-smr", "1e-5"], "gtex8_smr_out_sens")

# ---- 汇总 ----
COLS = ["Gene", "analysis", "probeID", "Probe_bp", "topSNP", "topSNP_bp", "A1", "A2", "Freq",
        "b_GWAS", "se_GWAS", "p_GWAS", "b_eQTL", "se_eQTL", "p_eQTL",
        "b_SMR", "se_SMR", "p_SMR", "p_HEIDI", "nsnp_HEIDI"]
rows = []
for tag, outname in (("A(5e-8)", "gtex8_smr_out"), ("B(1e-5)", "gtex8_smr_out_sens")):
    for D in DS:
        sp = os.path.join(SMRD, outname, D["tag"] + ".smr")
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
                d = dict(zip(hdr, p))
                d["analysis"] = tag
                d["Gene"] = D["gene"]
                d["_tag"] = D["tag"]
                rows.append(d)
os.makedirs(TAB, exist_ok=True)
with io.open(os.path.join(TAB, "59_gtex8_L2L3_smr_results.csv"), "w",
             encoding="utf-8", newline="") as f:
    wr = csv.writer(f); wr.writerow(COLS)
    for r in rows:
        wr.writerow([r.get(c, "NA") for c in COLS])

def fnum(x, fmt="%.5g"):
    try:
        return fmt % float(x)
    except (TypeError, ValueError):
        return "NA"

def verdict(ph):
    try:
        v = float(ph)
    except (TypeError, ValueError):
        return "不可评估"
    if np.isnan(v):
        return "不可评估"
    return "consistent（无显著异质性）" if v >= 0.05 else "heterogeneity（存在显著异质性）"

w("")
w("=" * 104)
w("结果：GTEx v8 全血 bulk（n=%s）；结局 GCST90483463 女性不孕症；LD=OneK1K 980 供者同区域子集"
  % "/".join(str(NMAP[d["tag"]]) for d in DS))
w("=" * 104)
w("%-9s %-9s %-14s %-9s %11s %11s %11s %6s  %s"
  % ("analysis", "Gene", "topSNP", "b_SMR", "p_SMR", "p_HEIDI", "b_GWAS", "m", "HEIDI 判读"))
for r in rows:
    w("%-9s %-9s %-14s %-9s %11s %11s %11s %6s  %s"
      % (r["analysis"], r["Gene"], r.get("topSNP", "NA"), fnum(r.get("b_SMR")),
         fnum(r.get("p_SMR"), "%.4e"), fnum(r.get("p_HEIDI"), "%.4e"),
         fnum(r.get("b_GWAS")), r.get("nsnp_HEIDI", "NA"), verdict(r.get("p_HEIDI"))))

w("")
w("已写 tables/59_gtex8_L2L3_smr_results.csv（%d 行）" % len(rows))
av3, l3 = mem()
w("结束内存: 可用 %.0f MB（占 %d%%）" % (av3, l3))
w("VERDICT = %s" % ("PASS" if len(rows) > 0 else "FAIL"))
flush()
print("WROTE " + LOG)
