# -*- coding: utf-8 -*-
"""s40d_prep_negctl_ld.py -- 任务 3.1 补丁：为「程序性阴性对照」计算**配对正确的带符号 r**

背景（本轮自查发现的实现缺陷）：
  s40c 的负对照臂把 LD 列硬绑在 `r_to_cond`（即与 rs56318008 的 r），
  却用负对照变异自身的 z_k 去条件 —— 即 r 与 k 不配对。
  后果：条件 z 被 1/sqrt(1-r_{j,rs56318008}^2) 无谓放大（实测 max|z| 27.14 / 18.52），
  负对照的 PP.H4 ~ 1 不能作为「机制不会无差别摧毁共定位」的证据。
  （主口径 cond = rs56318008 的 LD 列本身配对正确，故主结果不受影响。）

本脚本：从既有 ld_geno.raw（OneK1K 980 供者, GRCh37, A1 剂量）重算
  r(·, 1:22535811 = rs2473247) 与 r(·, 1:22885644 = rs12048511)，
  并将 r(·, 1:22470407) 一并重算，与既有 r_to_rs56318008.csv **逐位比对**（管线自验证）。

产出：00_data_raw/onek1k/cond_coloc/r_to_negctl.csv
  列：snp_grch37,pos37,a1,a2,af_ld,sd_a1,
      r_to_rs2473247,r2_to_rs2473247,r_to_rs12048511,r2_to_rs12048511,
      r_to_rs56318008_check,r2_to_rs56318008_check,d_r_check
日志：D:\\endometriosis_project\\_s40d_prep_negctl_ld.log
纪律：只读既有产物；只新增文件；结果写文件再读（不靠 stdout / 退出码）。
"""
import io
import os
import subprocess
import time

import numpy as np
import pandas as pd

PROJ = r"D:\endometriosis_project\11_sc_eqtl_mr_project"
COND = os.path.join(PROJ, "00_data_raw", "onek1k", "cond_coloc")
BFILE = r"D:\endometriosis_project\_onek1k_plink\plink_merged_980_donors"
LOG = r"D:\endometriosis_project\_s40d_prep_negctl_ld.log"

SID_MAIN = "1:22470407"     # rs56318008
SID_CTRL = {"rs2473247": "1:22535811", "rs12048511": "1:22885644"}
LO37, HI37 = 22000000, 22900000

L = []


def p(s=""):
    L.append(str(s))
    print(s)


def free_gb():
    try:
        out = subprocess.check_output(
            ["powershell", "-NoProfile", "-Command",
             "(Get-CimInstance Win32_OperatingSystem).FreePhysicalMemory"],
            text=True, stderr=subprocess.DEVNULL, timeout=60).strip()
        return int(out) / 1024.0 / 1024.0
    except Exception:
        return float("nan")


T0 = time.time()
p("=== s40d 负对照配对 LD 补丁 ===")
p("启动 %s   可用内存 %.2f GB" % (time.strftime("%Y-%m-%d %H:%M:%S"), free_gb()))

raw = os.path.join(COND, "ld_geno.raw")
if not os.path.exists(raw):
    p("!! ld_geno.raw 不存在，终止（需先跑 s40b）")
    io.open(LOG, "w", encoding="utf-8", newline="\n").write("\n".join(L) + "\n")
    raise SystemExit(2)

rw = pd.read_csv(raw, sep=r"\s+")
hdr = list(rw.columns)
snp_cols = [c for c in hdr if c not in ("FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE")]
cnt_allele = {c.rsplit("_", 1)[0]: c.rsplit("_", 1)[1] for c in snp_cols}
ids = [c.rsplit("_", 1)[0] for c in snp_cols]
p("raw shape = %s   剂量列 = %d   示例 %s -> 计数等位 %s"
  % (rw.shape, len(snp_cols), ids[0], cnt_allele[ids[0]]))

bim = pd.read_csv(BFILE + ".bim", sep="\t", header=None, dtype={0: str},
                  names=["chr", "snp", "cm", "pos", "a1", "a2"])
bim = bim[(bim["chr"] == "1") & (bim["pos"] >= LO37) & (bim["pos"] <= HI37)].copy()
bim["a1"] = bim["a1"].str.upper()
bim["a2"] = bim["a2"].str.upper()
map_a1 = dict(zip(bim["snp"], bim["a1"]))
map_a2 = dict(zip(bim["snp"], bim["a2"]))
pos_of = dict(zip(bim["snp"], bim["pos"]))

D = rw[snp_cols].astype(np.float32)
flip_to_a1 = []
for sid in ids:
    a1, a2, counted = map_a1.get(sid), map_a2.get(sid), cnt_allele[sid]
    flip_to_a1.append(0 if (counted == a1) else (1 if (counted == a2) else -1))
flip_to_a1 = np.array(flip_to_a1, dtype=np.int8)
p("计数等位 == bim A1 的列 = %d / %d ；== A2 的列 = %d ；无法判定 = %d"
  % (int((flip_to_a1 == 0).sum()), len(ids), int((flip_to_a1 == 1).sum()),
     int((flip_to_a1 < 0).sum())))
if (flip_to_a1 < 0).any():
    p("!! 存在无法对齐的等位，终止")
    io.open(LOG, "w", encoding="utf-8", newline="\n").write("\n".join(L) + "\n")
    raise SystemExit(3)

X = D.values.copy()
idx_a2 = np.where(flip_to_a1 == 1)[0]
X[:, idx_a2] = 2.0 - X[:, idx_a2]
X = X.astype(np.float64)
p("缺失率 = %.6f" % np.isnan(X).mean())

sd = np.nanstd(X, axis=0, ddof=1)
af_ld = np.nanmean(X, axis=0) / 2.0
Xc = X - np.nanmean(X, axis=0)


def r_vs(sid):
    if sid not in ids:
        raise SystemExit("条件变异 %s 不在 LD 矩阵中" % sid)
    j = ids.index(sid)
    xk = X[:, j] - np.nanmean(X[:, j])
    num = np.nansum(Xc * xk[:, None], axis=0)
    den = (len(xk) - 1) * sd * sd[j]
    with np.errstate(invalid="ignore", divide="ignore"):
        r = num / den
    r[j] = 1.0
    return r


p("")
p("[1] 计算配对 r")
out = pd.DataFrame(dict(snp_grch37=ids, pos37=[pos_of.get(s) for s in ids],
                        a1=[map_a1.get(s) for s in ids], a2=[map_a2.get(s) for s in ids],
                        af_ld=af_ld, sd_a1=sd))
for rs, sid in SID_CTRL.items():
    rr = r_vs(sid)
    out["r_to_%s" % rs] = rr
    out["r2_to_%s" % rs] = rr ** 2
    p("   r(·, %s = %s) 完成" % (sid, rs))
r_main = r_vs(SID_MAIN)
out["r_to_rs56318008_check"] = r_main
out["r2_to_rs56318008_check"] = r_main ** 2

# ---- 管线自验证：与既有 r_to_rs56318008.csv 逐位比对 ----
p("")
p("[2] 管线自验证（vs 既有 r_to_rs56318008.csv）")
old = pd.read_csv(os.path.join(COND, "r_to_rs56318008.csv"))
old["snp_grch37"] = old["snp_grch37"].astype(str)
mg = out.merge(old[["snp_grch37", "r_to_cond", "r2_to_cond"]], on="snp_grch37", how="left")
mg["d_r_check"] = (mg["r_to_rs56318008_check"] - mg["r_to_cond"]).abs()
overlap = mg[mg["r_to_cond"].notna()]
p("   重叠变异 = %d / %d" % (len(overlap), len(out)))
p("   |Δ r| max = %.3e   median = %.3e   (>1e-9 的个数 = %d)"
  % (overlap["d_r_check"].max(), overlap["d_r_check"].median(),
     int((overlap["d_r_check"] > 1e-9).sum())))
ok_pipe = bool(overlap["d_r_check"].max() < 1e-9)
p("   -> 管线一致：%s" % ("是（新 r 与既有 r 逐位一致）" if ok_pipe else "否 !!"))

# r2 已知值核验（沿用 s40b 口径）
p("   已知 r2 核验：")
for rsid, want in [("rs12037376", 0.899359), ("rs10917151", 0.84732)]:
    mm0 = pd.read_csv(os.path.join(COND, "m_B_MEM.csv"))
    q = mm0[mm0["rsid"].astype(str) == rsid]
    if len(q) == 0:
        p("     %s 不在交集内，跳过" % rsid)
        continue
    sid = str(q.iloc[0]["snp_grch37"])
    rr = out[out["snp_grch37"] == sid]
    if len(rr) == 0:
        p("     %s 不在 LD 矩阵内，跳过" % rsid)
        continue
    got = float(rr.iloc[0]["r2_to_rs56318008_check"])
    p("     r2(%s, rs56318008) = %.6f  (口径 %.6f, 差 %.6f) %s"
      % (rsid, got, want, abs(got - want), "OK" if abs(got - want) < 1e-4 else "!!"))

# 负对照自身 LD 概览（决定 keep 集的膨胀范围）
p("")
p("[3] 负对照的 LD 结构（决定各自 keep 集与最大膨胀因子）")
for rs, sid in SID_CTRL.items():
    rr = out["r_to_%s" % rs].values
    jj = ids.index(sid)
    mask = np.arange(len(rr)) != jj
    for kr in (1e-9, 0.01, 0.1):
        keep = rr ** 2 <= 1 - kr
        infl = 1.0 / np.sqrt(np.clip(1 - rr[keep] ** 2, 1e-300, None))
        p("   %s  1-r2>=%-6g -> 保留 %d（剔除 %d）  max_infl=%.3f（|r|=%.4f）"
          % (rs, kr, int(keep.sum()), int(len(rr) - keep.sum()), infl.max(),
             np.nanmax(np.abs(rr[keep]))))
    p("      |r| max（含自身）=%.4f   >0.99（不含自身）=%d   >0.9=%d"
      % (np.nanmax(np.abs(rr)), int((np.abs(rr[mask]) > 0.99).sum()),
         int((np.abs(rr) > 0.9).sum())))

out.to_csv(os.path.join(COND, "r_to_negctl.csv"), index=False, encoding="utf-8-sig")
p("")
p("写入 r_to_negctl.csv  rows=%d  cols=%d" % (len(out), len(out.columns)))
p("总耗时 %.1f s   末次可用内存 %.2f GB" % (time.time() - T0, free_gb()))
p("VERDICT = %s" % ("PASS" if ok_pipe else "FAIL(管线不一致)"))
io.open(LOG, "w", encoding="utf-8", newline="\n").write("\n".join(L) + "\n")
print("WROTE " + LOG)
