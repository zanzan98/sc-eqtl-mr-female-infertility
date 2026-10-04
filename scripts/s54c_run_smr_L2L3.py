# -*- coding: utf-8 -*-
"""s54c_run_smr_L2L3.py — 任务A：L2/L3 的 cis-SMR + HEIDI（复刻 s43c）

每个数据集独立 LD 参考（各自 cis 窗子集）：
  LD  = 00_data_raw/smr_L2L3/ld_ref/<TAG>_cis.{bed,bim,fam}
  GWAS= 00_data_raw/smr_L2L3/gwas_<TAG>.cojo.txt
  BESD= 00_data_raw/smr_L2L3/besd/<TAG>.{besd,esi,epi}

main: --peqtl-smr 5e-8 ; sens: --peqtl-smr 1e-5
其余参数与 chr1 链一致。
"""
import os, io, sys, glob, subprocess, ctypes, csv

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

PROJ = r"D:\endometriosis_project\11_sc_eqtl_mr_project"
SMRD = os.path.join(PROJ, "00_data_raw", "smr_L2L3")
BESD = os.path.join(SMRD, "besd")
LDD = os.path.join(SMRD, "ld_ref")
OUT_MAIN = os.path.join(SMRD, "smr_out")
OUT_SENS = os.path.join(SMRD, "smr_out_sens")
TAB = os.path.join(PROJ, "tables")
SMR_EXE = r"D:\endometriosis_project\_tools\smr-1.3.1-win-x86_64\smr-1.3.1-win-x86_64\smr-1.3.1-win.exe"
LOG = r"D:\endometriosis_project\_s54c_run_smr_L2L3.log"

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

def run(tag_label, outdir, extra, tags):
    os.makedirs(outdir, exist_ok=True)
    rows = []
    for i, tag in enumerate(tags, 1):
        ld = os.path.join(LDD, "%s_cis" % tag)
        cojo = os.path.join(SMRD, "gwas_%s.cojo.txt" % tag)
        be = os.path.join(BESD, tag)
        o = os.path.join(outdir, tag)
        cmd = [SMR_EXE, "--bfile", ld, "--gwas-summary", cojo, "--beqtl-summary", be,
               "--out", o, "--thread-num", "4"] + extra
        w("CMD: %s" % " ".join(cmd))
        p = subprocess.run(cmd, capture_output=True, text=True, encoding="utf-8", errors="replace")
        so = (p.stdout or "") + "\n" + (p.stderr or "")
        io.open(os.path.join(outdir, tag + ".console.txt"), "w", encoding="utf-8").write(so)
        sp = o + ".smr"
        n = 0
        if os.path.exists(sp):
            ls = [l for l in io.open(sp, encoding="utf-8", errors="replace").read().splitlines() if l.strip()]
            if len(ls) > 1:
                hdr = ls[0].split("\t")
                for ln in ls[1:]:
                    pp = ln.split("\t")
                    if len(pp) != len(hdr):
                        continue
                    d = dict(zip(hdr, pp)); d["dataset"] = tag; d["analysis"] = tag_label
                    rows.append(d); n += 1
        av, load = mem()
        w("[%s %d/%d] %-20s rc=%d rows=%d  内存可用=%.0fMB(占%d%%)"
          % (tag_label, i, len(tags), tag, p.returncode, n, av, load))
        for ln in so.splitlines():
            low = ln.lower()
            if any(k in low for k in ("error", "warning", "no cis-eqtl", "skipped",
                                      "excluded", "not enough")):
                w("        ! %s" % ln.strip()[:220])
    return rows

tags = sorted(os.path.basename(p).replace(".query.txt", "")
              for p in glob.glob(os.path.join(SMRD, "eqtl", "*.query.txt")))
av, load = mem()
w("=" * 90)
w("s54c — L2/L3 cis-SMR + HEIDI")
w("数据集 (%d): %s" % (len(tags), ", ".join(tags)))
w("起始内存: 可用 %.0f MB (占 %d%%)" % (av, load))
w("=" * 90)

w(); w("### [main] --peqtl-smr 5e-8")
rowsM = run("main(5e-8)", OUT_MAIN, [], tags)
w(); w("### [sens] --peqtl-smr 1e-5")
rowsS = run("sens(1e-5)", OUT_SENS, ["--peqtl-smr", "1e-5"], tags)

COLS = ["dataset", "analysis", "probeID", "ProbeChr", "Gene", "Probe_bp", "topSNP",
        "topSNP_chr", "topSNP_bp", "A1", "A2", "Freq", "b_GWAS", "se_GWAS", "p_GWAS",
        "b_eQTL", "se_eQTL", "p_eQTL", "b_SMR", "se_SMR", "p_SMR", "p_HEIDI", "nsnp_HEIDI"]

def dump(path, rows):
    with io.open(path, "w", encoding="utf-8", newline="") as f:
        wr = csv.writer(f); wr.writerow(COLS)
        for r in rows:
            wr.writerow([r.get(c, "NA") for c in COLS])
    w("  写出 %s (%d 行)" % (path, len(rows)))

w()
dump(os.path.join(TAB, "56_smr_L2L3_main.csv"), rowsM)
dump(os.path.join(TAB, "56b_smr_L2L3_sens.csv"), rowsS)

w(); w("### 结果明细")
for r in rowsM + rowsS:
    w("  %-20s %-9s %-10s topSNP=%-14s b_SMR=%10s p_SMR=%10s p_HEIDI=%10s m=%s"
      % (r["dataset"], r["analysis"], r.get("Gene", "?"), r.get("topSNP", "?"),
         r.get("b_SMR", "?"), r.get("p_SMR", "?"), r.get("p_HEIDI", "?"), r.get("nsnp_HEIDI", "?")))
av, load = mem()
w(); w("VERDICT = %s" % ("PASS" if len(rowsM) >= 1 or len(rowsS) >= 1 else "FAIL"))
w("结束内存: 可用 %.0f MB (占 %d%%)" % (av, load))
flush()
