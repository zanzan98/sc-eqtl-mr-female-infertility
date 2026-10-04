# -*- coding: utf-8 -*-
"""s31_sensitivity_outcome_mr.py —— Task A：P3 敏感性结局（GCST90483469）MR

顺序（严格）：
  1. 校验 GCST90483469.h.tsv.gz：字节数 == 官方 Content-Length、md5 == EBI md5sum.txt、
     gzip 流式读到 EOF 不抛 EOFError（三者全过才继续）
  2. 用与主分析**完全相同**的 s06_mr_extract.py 参数跑一遍，只换结局文件与 label
  3. 与主结局结果（10_discovery_MR_main.csv）逐（gene, cell_type, variant_id）对比

官方基准：EBI harmonised/md5sum.txt -> 322ed203bcc7799a6d7e0e8201a9ec80
          Content-Length -> 1,128,940,620
"""
import gzip, hashlib, os, subprocess, sys, time, datetime
import numpy as np
import pandas as pd

ROOT = r"D:\endometriosis_project\11_sc_eqtl_mr_project"
GW   = os.path.join(ROOT, "00_data_raw", "gwas")
TAB  = os.path.join(ROOT, "tables")
PY   = r"C:\Users\28144\.workbuddy\binaries\python\envs\default\Scripts\python.exe"

TARGET = os.path.join(GW, "GCST90483469.h.tsv.gz")
EXPECT_SIZE = 1128940620
EXPECT_MD5  = "322ed203bcc7799a6d7e0e8201a9ec80"

INSTR = os.path.join(TAB, "01b_discovery_instruments_grch38.csv")
BIMZ  = os.path.join(ROOT, "00_data_raw", "onek1k", "plink_genotype_980.zip")
OUT_SENS = os.path.join(TAB, "38_discovery_MR_sensitivity_GCST90483469.csv")
MAIN  = os.path.join(TAB, "10_discovery_MR_main.csv")
CMP   = os.path.join(TAB, "39_sensitivity_vs_main_GCST90483469.csv")
LOG   = r"D:\endometriosis_project\_s31_sensitivity.log"

L = []
def w(s=""):
    line = str(s)
    L.append(line)
    print(line)

def flush():
    with open(LOG, "w", encoding="utf-8") as fh:
        fh.write("\n".join(L) + "\n")

w("=== s31 sensitivity outcome MR (Task A) @ %s ===" % datetime.datetime.now().isoformat())

# ---------------------------------------------------------------- 1) verify
w("\n[1] 完整性校验 %s" % TARGET)
if not os.path.exists(TARGET):
    w("  FATAL: 文件不存在"); flush(); sys.exit(2)

with open(TARGET, "rb") as fh:
    fh.seek(0, os.SEEK_END)
    sz = fh.tell()
w("  size = %d (expect %d) -> %s" % (sz, EXPECT_SIZE, "OK" if sz == EXPECT_SIZE else "MISMATCH"))
if sz != EXPECT_SIZE:
    w("  FATAL: 字节数不符，拒绝用残缺文件运行。"); flush(); sys.exit(2)

t0 = time.time()
h = hashlib.md5()
with open(TARGET, "rb") as fh:
    for chunk in iter(lambda: fh.read(1 << 22), b""):
        h.update(chunk)
got = h.hexdigest()
w("  md5  = %s (expect %s) -> %s   [%.1fs]"
  % (got, EXPECT_MD5, "OK" if got == EXPECT_MD5 else "MISMATCH", time.time() - t0))
if got != EXPECT_MD5:
    w("  FATAL: md5 不符。"); flush(); sys.exit(2)

# gzip stream to EOF (a truncated file throws EOFError here)
t0 = time.time()
nbytes = 0
nline = 0
try:
    with gzip.open(TARGET, "rb") as fh:
        while True:
            b = fh.read(1 << 24)
            if not b:
                break
            nbytes += len(b)
            nline += b.count(b"\n")
    w("  gzip = OK, 解压 %.2f GB / %d 行  [%.1fs]"
      % (nbytes / 1e9, nline, time.time() - t0))
except EOFError as e:
    w("  FATAL: gzip 流截断 -> EOFError: %s" % e); flush(); sys.exit(2)
except Exception as e:
    w("  FATAL: gzip 读取异常 %s: %s" % (type(e).__name__, e)); flush(); sys.exit(2)

w("  === 三项校验全部通过 ===")
flush()

# ---------------------------------------------------------------- 2) run s06
w("\n[2] 运行 s06_mr_extract.py（与主分析同参数，仅换结局）")
cmd = [PY, os.path.join(ROOT, "scripts", "s06_mr_extract.py"),
       "--instruments", INSTR,
       "--outcome-gz", TARGET,
       "--bim-zip", BIMZ,
       "--out", OUT_SENS,
       "--label", "GCST90483469"]
w("  cmd = " + " ".join('"%s"' % c if " " in c else c for c in cmd))
t0 = time.time()
p = subprocess.run(cmd, capture_output=True, text=True,
                   encoding="utf-8", errors="replace", cwd=ROOT)
w("  exit = %d  [%.1f min]" % (p.returncode, (time.time() - t0) / 60))
for ln in (p.stdout or "").splitlines()[-25:]:
    w("  | " + ln)
for ln in (p.stderr or "").splitlines()[-10:]:
    w("  ! " + ln)
if not os.path.exists(OUT_SENS):
    w("  FATAL: 未生成 %s" % OUT_SENS); flush(); sys.exit(3)
flush()

# ---------------------------------------------------------------- 3) compare
w("\n[3] 与主结局逐对比较")
sens = pd.read_csv(OUT_SENS, dtype={"variant_id": str}, encoding="utf-8-sig")
main = pd.read_csv(MAIN, encoding="utf-8-sig")
w("  敏感性表行 = %d ; 主表行 = %d" % (len(sens), len(main)))

key = ["gene", "cell_type"]
sv = sens[sens["p_MR"].notna()].copy()
sv["variant_id"] = sv["variant_id"].astype(str)
mv = main.copy()
if "variant_id_grch38" in mv.columns:
    mv = mv.rename(columns={"variant_id_grch38": "variant_id"})
mv["variant_id"] = mv["variant_id"].astype(str)

j = sv.merge(mv, on=key, how="inner", suffixes=("_sens", "_main"))
w("  可连接对 = %d" % len(j))
# variant consistency within pair
j["same_variant"] = j["variant_id_sens"] == j["variant_id_main"]
w("  同变异一致 = %d / %d" % (int(j["same_variant"].sum()), len(j)))

j = j[j["same_variant"]].copy()
if len(j) == 0:
    w("  !! 无同变异配对，无法比较"); flush(); sys.exit(4)

j["b_main"] = j["b"]
j["b_sens"] = j["beta_MR"]
j["p_main"] = j["p"]
j["p_sens"] = j["p_MR"]
j["sign_consistent"] = np.sign(j["b_main"]) == np.sign(j["b_sens"])
j["sens_FDR05"] = j["p_FDR"] < 0.05
j["main_FDR05"] = j["qval"] < 0.05
j["both_FDR05"] = j["sens_FDR05"] & j["main_FDR05"]

w("")
w("  方向一致率            = %d / %d = %.4f"
  % (int(j["sign_consistent"].sum()), len(j), j["sign_consistent"].mean()))
w("  beta 相关系数         = %.4f" % j[["b_main", "b_sens"]].corr().iloc[0, 1])
w("  主结局 FDR<0.05       = %d" % int(j["main_FDR05"].sum()))
w("  敏感性 FDR<0.05       = %d" % int(j["sens_FDR05"].sum()))
w("  两者同时 FDR<0.05     = %d" % int(j["both_FDR05"].sum()))
w("  主显著、敏感不显著     = %d" % int((j["main_FDR05"] & ~j["sens_FDR05"]).sum()))
w("  敏感显著、主不显著     = %d" % int((~j["main_FDR05"] & j["sens_FDR05"]).sum()))

# L1 三对（本任务 D 关注）
w("\n  --- L1 座（CDC42 / LINC00339）逐对 ---")
l1 = j[j["gene"].isin(["CDC42", "LINC00339"]) & j["variant_id_sens"].isin(
    ["1:22088292", "1:22135618", "1:22132301", "1:22096228", "1:22113027", "1:22031964"])]
show = ["gene", "cell_type", "variant_id_sens", "b_main", "b_sens", "p_main", "p_sens",
        "sign_consistent", "main_FDR05", "sens_FDR05"]
w(l1[show].to_string(index=False))

keep = ["gene", "cell_type", "variant_id_sens", "beta_exp", "se_exp", "F",
        "beta_out_sens", "se_out_sens", "p_out_sens", "b_sens", "se_MR", "p_sens", "z_MR",
        "p_FDR", "status_sens", "allele_source_sens",
        "b_main", "se_main", "p_main", "qval", "F_stat",
        "sign_consistent", "main_FDR05", "sens_FDR05", "both_FDR05"]
keep = [c for c in keep if c in j.columns]
j[keep].sort_values("p_sens").to_csv(CMP, index=False, encoding="utf-8-sig")
w("\n  写出 %s (rows=%d)" % (CMP, len(j)))

w("\ntime %s" % datetime.datetime.now().isoformat())
flush()
