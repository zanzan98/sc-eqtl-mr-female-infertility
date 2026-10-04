# -*- coding: utf-8 -*-
"""s43g_run_gtex8_smr.py — 选项 A 收尾：GTEx v8 全血（bulk）3 基因 SMR + HEIDI

暴露：GTEx v8 全血（eQTL Catalogue harmonised，n≈670，bulk 全血 —— 非单细胞！）
结局：GCST90483463（女性不孕症）
LD 参考：OneK1K 980 供者子集 chr1:20–25 Mb（hg19；与 GTEx 供者不重叠，为欧洲血统近似 LD）

参数化：
  [A] default  --peqtl-smr 5e-8   （预期：WNT4 无 SNP 过阈 → 不检验；CDC42/LINC00339 可检验）
  [B] relaxed  --peqtl-smr 1e-4   （WNT4 唯一可行设置：minP = 5.53e-05）

输出：00_data_raw/smr_3p3/gtex8/{besd,smr_out_A,smr_out_B}/...；tables/47_*.csv
"""
import os, io, sys, glob, subprocess, ctypes, csv
import numpy as np
import pandas as pd

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

PROJ = r"D:\endometriosis_project\11_sc_eqtl_mr_project"
SMR3 = os.path.join(PROJ, "00_data_raw", "smr_3p3")
G8 = os.path.join(SMR3, "gtex8")
QDIR = os.path.join(G8, "query")
BESD = os.path.join(G8, "besd")
LD = os.path.join(SMR3, "ld_ref", "chr1_20_25mb")
COJO = os.path.join(SMR3, "gwas_GCST90483463.cojo.txt")
TAB = os.path.join(PROJ, "tables")
EXE = r"D:\endometriosis_project\_tools\smr-1.3.1-win-x86_64\smr-1.3.1-win-x86_64\smr-1.3.1-win.exe"
LOG = os.path.join(PROJ, "scripts", "_s43g_run_gtex8_smr.log")

buf = []
def w(s=""):
    buf.append(str(s)); print(s)
def flush():
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
    return st.ullAvailPhys / 1048576.0, st.dwMemoryLoad

os.makedirs(BESD, exist_ok=True)
GENES = ["WNT4", "CDC42", "LINC00339"]

w("=" * 96)
w("s43g_run_gtex8_smr.py — GTEx v8 全血（bulk）SMR + HEIDI（任务三 3.3 · 选项 A）")
w("=" * 96)
av, load = mem(); w("起始内存: 可用 %.0f MB (占 %d%%)" % (av, load))

# ---- N（由 an 列推样本量）
w()
w("[0] 由 query 的 Freq/等位数列推 GTEx 全血样本量")
n_obs = []
for g in GENES:
    an = []
    with io.open(os.path.join(QDIR, "%s.query.txt" % g), "r", encoding="utf-8") as f:
        hdr = f.readline()
        for ln in f:
            p = ln.rstrip("\n").split("\t")
            if len(p) >= 14:
                pass
    # an 未写入 query（query 只有 14 列）；改为从源 TSV 读
w("    （an 列未入 query；改从源 TSV 读取）")
src = r"D:\endometriosis_project\01_data_raw\eqtlcat_gtex_v8\Whole_Blood.chr1_21_23.3Mb_3genes.tsv"
an_vals = []
with io.open(src, "r", encoding="utf-8", errors="replace") as f:
    ix = {c: i for i, c in enumerate(f.readline().rstrip("\n").split("\t"))}
    for ln in f:
        p = ln.rstrip("\n").split("\t")
        if len(p) > ix["an"]:
            an_vals.append(int(p[ix["an"]]))
an_arr = np.array(an_vals)
N_GTEX = int(round(np.median(an_arr) / 2.0))
w("    an 中位数 = %d → 样本量 N = an/2 = %d （an 唯一值数=%d，范围 %d..%d）"
  % (int(np.median(an_arr)), N_GTEX, len(set(an_arr.tolist())), an_arr.min(), an_arr.max()))

# ---- BESD
w()
w("[1] 建 BESD（--make-besd --add-n %d）" % N_GTEX)
for g in GENES:
    o = os.path.join(BESD, g)
    if os.path.exists(o + ".besd") and os.path.exists(o + ".epi") and os.path.exists(o + ".esi"):
        w("    %-10s 已存在，跳过" % g); continue
    p = subprocess.run([EXE, "--qfile", os.path.join(QDIR, "%s.query.txt" % g),
                        "--make-besd", "--out", o, "--add-n", str(N_GTEX)],
                       capture_output=True, text=True, encoding="utf-8", errors="replace", cwd=BESD)
    ok = all(os.path.exists(o + e) for e in (".besd", ".epi", ".esi"))
    w("    %-10s rc=%d  生成=%s" % (g, p.returncode, ok))
    io.open(os.path.join(BESD, g + ".console.txt"), "w", encoding="utf-8").write(
        (p.stdout or "") + "\n" + (p.stderr or ""))

# ---- 检查 .summary（cis 窗口与 SNP 数）
w()
for g in GENES:
    sp = os.path.join(BESD, g + ".summary")
    if os.path.exists(sp):
        txt = io.open(sp, "r", encoding="utf-8", errors="replace").read().splitlines()
        w("    %s.summary: %s" % (g, " || ".join([l for l in txt if l.startswith("{")])))

# ---- SMR 运行
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
        io.open(outdir + "\\%s.console.txt" % g, "w", encoding="utf-8").write(so)
        key = [l.strip() for l in so.splitlines()
               if ("to be included after allele checking" in l
                   or "with at least a cis-eQTL" in l
                   or "allele frequency differences" in l
                   or "Probes to be included" in l)]
        nrow = 0
        sp = o + ".smr"
        if os.path.exists(sp):
            nrow = max(0, len([l for l in io.open(sp, "r", encoding="utf-8",
                                                   errors="replace").read().splitlines() if l.strip()]) - 1)
        av2, l2 = mem()
        w("[%s] %-10s rc=%d rows=%d 内存可用=%.0fMB(占%d%%)" % (tag, g, p.returncode, nrow, av2, l2))
        for k in key:
            w("        %s" % k[:190])
        res.append((g, p.returncode, nrow))
    return res

w()
w("=" * 96)
w("[A] 默认参数（--peqtl-smr 5e-8）")
w("=" * 96)
resA = run_smr("A", [], "smr_out_A")
w()
w("=" * 96)
w("[B] 放宽 target 阈值（--peqtl-smr 1e-4）—— WNT4 唯一可行设置")
w("=" * 96)
resB = run_smr("B", ["--peqtl-smr", "1e-4"], "smr_out_B")

# ---- 汇总
COLS = ["Gene", "analysis", "probeID", "Probe_bp", "topSNP", "topSNP_bp", "A1", "A2", "Freq",
        "b_GWAS", "se_GWAS", "p_GWAS", "b_eQTL", "se_eQTL", "p_eQTL",
        "b_SMR", "se_SMR", "p_SMR", "p_HEIDI", "nsnp_HEIDI"]
rows = []
for tag, outname in (("A(5e-8)", "smr_out_A"), ("B(1e-4)", "smr_out_B")):
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
                    w("  !! %s %s 列数不符 %d vs %d" % (tag, g, len(p), len(hdr))); continue
                d = dict(zip(hdr, p)); d["analysis"] = tag; d["Gene"] = g
                rows.append(d)
os.makedirs(TAB, exist_ok=True)
with io.open(os.path.join(TAB, "47_gtex8_wb_smr_results.csv"), "w", encoding="utf-8", newline="") as f:
    wr = csv.writer(f); wr.writerow(COLS)
    for r in rows:
        wr.writerow([r.get(c, "NA") for c in COLS])
w()
w("写出 tables/47_gtex8_wb_smr_results.csv（%d 行）" % len(rows))

w()
w("=" * 96)
w("结果（GTEx v8 全血 bulk；结局 GCST90483463 女性不孕症）")
w("=" * 96)
w("%-9s %-10s %-14s %9s %11s %11s %5s" % ("analysis", "Gene", "topSNP", "b_SMR", "p_SMR", "p_HEIDI", "m"))
for r in rows:
    w("%-9s %-10s %-14s %9s %11s %11s %5s"
      % (r["analysis"], r["Gene"], r.get("topSNP", "NA"),
         ("%.5g" % float(r["b_SMR"])) if r.get("b_SMR") not in (None, "NA") else "NA",
         r.get("p_SMR", "NA"), r.get("p_HEIDI", "NA"), r.get("nsnp_HEIDI", "NA")))

av3, l3 = mem()
w()
w("结束内存: 可用 %.0f MB (占 %d%%)" % (av3, l3))
w("VERDICT = %s" % ("PASS" if len(rows) > 0 else "FAIL"))
flush()
