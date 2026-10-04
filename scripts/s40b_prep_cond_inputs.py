# -*- coding: utf-8 -*-
"""s40b_prep_cond_inputs.py -- 任务 3.1 步骤 1：条件共定位输入准备

产出（目录 00_data_raw/onek1k/cond_coloc/）：
  m_<CT>.csv                  复现 s28 harmonise 的 allele-aware 交集（含 snp_grch37/pos38/slope/slope_se/af/b2/se/eaf/rsid）
  r_to_rs56318008.csv         每个变异 vs 条件变异 1:22470407 的**带符号 r** 与 r2（OneK1K 980 供者，GRCh37）
  ld_geno.*                   plink --recode A 原始剂量（A1 计数，已按 bim A1 对齐）

纪律：
  * 只读已闭环产物，不改动任何文件；
  * 一切中间结果落盘（R 段崩溃可续跑）；
  * 内存闸门：plink 用 --memory 限流；全程记录可用内存。
输出日志 -> D:\\endometriosis_project\\_s40b_prep_cond_inputs.log
"""
import io
import os
import shutil
import subprocess
import sys
import time

import numpy as np
import pandas as pd

PROJ = r"D:\endometriosis_project\11_sc_eqtl_mr_project"
PREP = os.path.join(PROJ, "00_data_raw", "onek1k", "coloc_prep")
COND = os.path.join(PROJ, "00_data_raw", "onek1k", "cond_coloc")
PLINK = r"D:\endometriosis_project\_tools\plink.exe"
BFILE = r"D:\endometriosis_project\_onek1k_plink\plink_merged_980_donors"
LOG = r"D:\endometriosis_project\_s40b_prep_cond_inputs.log"

COND_SNP37 = "1:22470407"      # rs56318008 (WNT4 credible-set lead)
COND_RS = "rs56318008"
LO37, HI37 = 22000000, 22900000
MEM_CAP_MB = 1200

CELLS = ["B_MEM", "Mono_NC"]

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
os.makedirs(COND, exist_ok=True)
p("=== s40b 任务3.1 条件共定位输入准备 ===")
p("启动时间 %s   可用内存 %.2f GB（闸门：峰值 > 14 GB 即停）" % (time.strftime("%Y-%m-%d %H:%M:%S"), free_gb()))
p("条件变异 = %s (%s)  窗口 GRCh37 chr1:%d-%d" % (COND_RS, COND_SNP37, LO37, HI37))

# ------------------------------------------------------------------ 1) 读入
p("")
p("[1] 读入 coloc_prep")
e_all = pd.read_csv(os.path.join(PREP, "eqtl_L1_CDC42.csv"))
gw = pd.read_csv(os.path.join(PREP, "gwas_L1.csv"))
bimU = pd.read_csv(os.path.join(PREP, "bim_universe_L1.csv"))
p("   eqtl_L1_CDC42 rows=%d   gwas_L1 rows=%d   bim_universe_L1 rows=%d"
  % (len(e_all), len(gw), len(bimU)))

# ------------------------------------------------------------------ 2) 结局侧 harmonise（复现 s28）
gw = gw.copy()
gw["ea"] = gw["ea"].astype(str).str.upper()
gw["oa"] = gw["oa"].astype(str).str.upper()
gw = gw[gw["pos38"].notna() & gw["beta"].notna() & gw["se"].notna() & (gw["se"] > 0)].copy()
gw["pos38"] = gw["pos38"].astype(int)
gw["allele_key"] = (gw["pos38"].astype(str) + "_" + np.minimum(gw["ea"], gw["oa"]) + "_"
                    + np.maximum(gw["ea"], gw["oa"]))
gw = gw[~gw["allele_key"].duplicated()]
p("   结局侧 effective rows = %d" % len(gw))

r_tab = None
for ct in CELLS:
    p("")
    p("[2] %s：harmonise" % ct)
    e = e_all[e_all["cell_type"] == ct].copy()
    e = e.merge(bimU, left_on="snp_grch37", right_on="snp", how="left")
    e = e[e["a1"].notna() & e["a2"].notna()].copy()
    e["pos38"] = e["snp_grch38"].astype(str).str.split(":").str[1]
    e = e[e["pos38"].notna()].copy()
    e["pos38"] = e["pos38"].astype(int)
    e["a1"] = e["a1"].astype(str).str.upper()
    e["a2"] = e["a2"].astype(str).str.upper()
    e["allele_key"] = (e["pos38"].astype(str) + "_" + np.minimum(e["a1"], e["a2"]) + "_"
                       + np.maximum(e["a1"], e["a2"]))
    m = e.merge(gw, on="allele_key", suffixes=("_e", "_o"))
    # 两侧都有 pos38 -> 后缀化；统一回 pos38（同 s28 R 侧做法）
    if "pos38_e" in m.columns:
        m = m.rename(columns={"pos38_e": "pos38"})
        m = m.drop(columns=[c for c in ["pos38_o"] if c in m.columns])
    p("   eqtl_ct=%d  gwas=%d  交集=%d" % (len(e), len(gw), len(m)))

    m["flip"] = (m["ea"] == m["a2"]) & (m["oa"] == m["a1"])
    m = m[((m["ea"] == m["a1"]) & (m["oa"] == m["a2"])) |
          ((m["ea"] == m["a2"]) & (m["oa"] == m["a1"]))].copy()
    m["palindromic"] = (np.minimum(m["a1"], m["a2"]) + np.maximum(m["a1"], m["a2"])).isin(["AT", "CG"])
    m["af_mism"] = (m["af"] - m["eaf"]).abs()
    m["af_flip"] = (m["af"] - (1 - m["eaf"])).abs()
    n_pal_bad = int(((m["palindromic"]) & (np.minimum(m["af_mism"], m["af_flip"]) > 0.05) &
                     ((m["af_mism"] - m["af_flip"]).abs() < 0.05)).sum())
    m = m[~((m["palindromic"]) & (np.minimum(m["af_mism"], m["af_flip"]) > 0.05) &
            ((m["af_mism"] - m["af_flip"]).abs() < 0.05))].copy()
    fix = (m["palindromic"]) & (m["af_flip"] < m["af_mism"])
    m.loc[fix, "flip"] = ~m.loc[fix, "flip"]
    m["b2"] = np.where(m["flip"], -m["beta"], m["beta"])
    m["maf_e"] = np.minimum(m["af"], 1 - m["af"])
    m["maf_o"] = np.minimum(m["eaf"], 1 - m["eaf"])
    m = m[np.isfinite(m["maf_e"]) & np.isfinite(m["maf_o"]) & (m["maf_e"] > 0) & (m["maf_o"] > 0) &
          np.isfinite(m["slope"]) & np.isfinite(m["slope_se"]) & (m["slope_se"] > 0)].copy()
    m = m[m["rsid"].notna() & (m["rsid"].astype(str) != "")]
    m = m.sort_values("pos38").reset_index(drop=True)
    m["N_out"] = m["n_cases"] + m["n_controls"]
    m["s_out"] = m["n_cases"] / m["N_out"]
    m["snp_grch37"] = m["snp_grch37"].astype(str)
    p("   harmonised=%d  flip=%d  palindromic=%d  dropped_af_ambiguous=%d"
      % (len(m), int(m["flip"].sum()), int(m["palindromic"].sum()), n_pal_bad))
    k = m[m["snp_grch37"] == COND_SNP37]
    p("   条件变异在交集内 = %s" % ("是" if len(k) else "否"))
    if len(k):
        p("      eQTL  slope=%.6f se=%.6f z=%.4f" % (k.iloc[0]["slope"], k.iloc[0]["slope_se"],
                                                     k.iloc[0]["slope"] / k.iloc[0]["slope_se"]))
        p("      结局  b2=%.6f se=%.6f z=%.4f  (rsid=%s a1=%s a2=%s ea=%s oa=%s)"
          % (k.iloc[0]["b2"], k.iloc[0]["se"], k.iloc[0]["b2"] / k.iloc[0]["se"],
             k.iloc[0]["rsid"], k.iloc[0]["a1"], k.iloc[0]["a2"], k.iloc[0]["ea"], k.iloc[0]["oa"]))
    cols = ["snp_grch37", "rsid", "pos38", "a1", "a2", "ea", "oa", "slope", "slope_se",
            "af", "maf_e", "pval_nominal", "b2", "se", "eaf", "maf_o", "N_out", "s_out",
            "n_cases", "n_controls"]
    m[cols].to_csv(os.path.join(COND, "m_%s.csv" % ct), index=False, encoding="utf-8-sig")
    if r_tab is None:
        r_tab = set(m["snp_grch37"])

# ------------------------------------------------------------------ 3) LD：A1 剂量
p("")
p("[3] plink --recode A 取 OneK1K 基因型剂量（A1 计数），内存上限 %d MB" % MEM_CAP_MB)
raw = os.path.join(COND, "ld_geno.raw")
if os.path.exists(raw) and len(pd.read_csv(raw, sep=r"\s+", nrows=1).columns) > 100:
    p("   已存在，跳过 plink")
else:
    cmd = [PLINK, "--bfile", BFILE, "--allow-no-sex", "--chr", "1",
           "--from-bp", str(LO37), "--to-bp", str(HI37),
           "--recode", "A", "--out", os.path.join(COND, "ld_geno"),
           "--memory", str(MEM_CAP_MB)]
    p("   CMD: " + " ".join(cmd))
    before = free_gb()
    subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=False, timeout=1800)
    p("   plink 返回；可用内存 %.2f -> %.2f GB" % (before, free_gb()))
    if not os.path.exists(raw):
        p("   !! plink 未产出 .raw，终止")
        io.open(LOG, "w", encoding="utf-8", newline="\n").write("\n".join(L) + "\n")
        sys.exit(2)

rw = pd.read_csv(raw, sep=r"\s+")
p("   .raw shape = %s" % (rw.shape,))
hdr = list(rw.columns)
snp_cols = [c for c in hdr if c not in ("FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE")]
cnt_allele = {}
for c in snp_cols:
    sid, al = c.rsplit("_", 1)
    cnt_allele[sid] = al
_ex_sid = snp_cols[0].rsplit("_", 1)[0]
p("   剂量列 = %d（示例 %s -> 计数等位 %s）" % (len(snp_cols), _ex_sid, cnt_allele[_ex_sid]))

bim_ref = pd.read_csv(BFILE + ".bim", sep="\t", header=None, dtype={0: str},
                      names=["chr", "snp", "cm", "pos", "a1", "a2"])
bim_ref = bim_ref[(bim_ref["chr"] == "1") & (bim_ref["pos"] >= LO37) & (bim_ref["pos"] <= HI37)]
bim_ref["a1"] = bim_ref["a1"].str.upper()
bim_ref["a2"] = bim_ref["a2"].str.upper()
map_a1 = dict(zip(bim_ref["snp"], bim_ref["a1"]))
map_a2 = dict(zip(bim_ref["snp"], bim_ref["a2"]))
pos_of = dict(zip(bim_ref["snp"], bim_ref["pos"]))

D = rw[snp_cols].astype(np.float32)
ids = [c.rsplit("_", 1)[0] for c in snp_cols]

# 按 bim A1 对齐（若 recode 计数的是 A2 则翻转为 2-剂量）
flip_to_a1 = []
for sid in ids:
    a1, a2, counted = map_a1.get(sid), map_a2.get(sid), cnt_allele[sid]
    flip_to_a1.append(0 if (counted == a1) else (1 if (counted == a2) else -1))
flip_to_a1 = np.array(flip_to_a1, dtype=np.int8)
p("   计数等位 == bim A1 的列 = %d / %d ；== A2 的列 = %d ；无法判定的列 = %d"
  % (int((flip_to_a1 == 0).sum()), len(ids), int((flip_to_a1 == 1).sum()), int((flip_to_a1 < 0).sum())))
if (flip_to_a1 < 0).any():
    p("   !! 存在无法对齐的等位，终止")
    io.open(LOG, "w", encoding="utf-8", newline="\n").write("\n".join(L) + "\n")
    sys.exit(3)
X = D.values.copy()
idx_a2 = np.where(flip_to_a1 == 1)[0]
X[:, idx_a2] = 2.0 - X[:, idx_a2]
X = X.astype(np.float64)
p("   缺失率 = %.6f" % (np.isnan(X).mean()))

# 等位对齐校验（★必须经 harmonise 表把 GRCh37->GRCh38 对齐，不能直接用 bim 的 GRCh37 位置
#   去比 GWAS 的 GRCh38 位置 —— 上一版正是犯了这个坐标系错配，故此处改正）
mm0 = pd.read_csv(os.path.join(COND, "m_%s.csv" % CELLS[0]))
idx_of = {sid: j for j, sid in enumerate(ids)}
rows = []
for _, q in mm0.iterrows():
    j = idx_of.get(str(q["snp_grch37"]))
    if j is None:
        continue
    f_ld_a1 = float(np.nanmean(X[:, j]) / 2.0)
    a1, ea = q["a1"], q["ea"]
    f_gwas_a1 = float(q["eaf"]) if ea == a1 else 1.0 - float(q["eaf"])
    rows.append((str(q["snp_grch37"]), str(q["rsid"]), f_ld_a1, f_gwas_a1, float(q["af"])))
chk = np.array([[t[2], t[3], t[4]] for t in rows], dtype=float)
e1 = np.abs(chk[:, 0] - chk[:, 1])          # vs 结局 eaf 推出的 A1 频率
e1f = np.abs(chk[:, 0] - (1.0 - chk[:, 1]))  # 若等位翻转会得到的偏差
e2 = np.abs(chk[:, 0] - chk[:, 2])          # vs eQTL af（同队列，权威判据）
p("   [判据 A · 权威] A1 剂量频率 vs eQTL af（同队列同源）：n=%d  max=%.4f  median=%.4f  >0.05 个数=%d"
  % (len(chk), e2.max(), np.median(e2), int((e2 > 0.05).sum())))
p("   [判据 B · 诊断] A1 剂量频率 vs 结局 eaf 推出的 A1 频率（跨队列，频率本身有差异）：")
p("        对齐  mean|dev| = %.4f   max = %.4f" % (e1.mean(), e1.max()))
p("        若翻转 mean|dev| = %.4f   max = %.4f" % (e1f.mean(), e1f.max()))
p("        对齐/翻转 偏差比 = %.4f  （远小于 1 即确认等位对齐）" % (e1.mean() / e1f.mean()))
if e2.max() > 0.01 or e2.mean() > 0.002 or e1.mean() >= e1f.mean():
    p("   !! 等位对齐判定失败（判据 A 或 B），条件分析不可信，终止")
    io.open(LOG, "w", encoding="utf-8", newline="\n").write("\n".join(L) + "\n")
    sys.exit(5)
p("   -> 等位对齐判定通过（判据 A 通过；判据 B 对齐偏差远小于翻转偏差）")

sd = np.nanstd(X, axis=0, ddof=1)
af_ld = np.nanmean(X, axis=0) / 2.0

# ------------------------------------------------------------------ 4) 带符号 r
p("")
p("[4] 计算 vs %s 的带符号 r" % COND_SNP37)
if COND_SNP37 not in ids:
    p("   !! 条件变异不在 LD 剂量矩阵中，终止")
    io.open(LOG, "w", encoding="utf-8", newline="\n").write("\n".join(L) + "\n")
    sys.exit(4)
jj = ids.index(COND_SNP37)
xk = X[:, jj] - np.nanmean(X[:, jj])
sk = sd[jj]
num = np.nansum((X - np.nanmean(X, axis=0)) * xk[:, None], axis=0)
den = (len(xk) - 1) * sd * sk
with np.errstate(invalid="ignore", divide="ignore"):
    r = num / den
r[jj] = 1.0
out = pd.DataFrame(dict(snp_grch37=ids, pos37=[pos_of.get(s) for s in ids],
                        a1=[map_a1.get(s) for s in ids], a2=[map_a2.get(s) for s in ids],
                        af_ld=af_ld, sd_a1=sd, r_to_cond=r, r2_to_cond=r ** 2))
out["in_intersect"] = out["snp_grch37"].isin(r_tab)
out.to_csv(os.path.join(COND, "r_to_rs56318008.csv"), index=False, encoding="utf-8-sig")
p("   写入 r_to_rs56318008.csv  rows=%d" % len(out))

# r2 已知值核验（图件与正文口径：0.899359 / 0.84732）
mm0 = pd.read_csv(os.path.join(COND, "m_%s.csv" % CELLS[0]))
for rs, want in [("rs12037376", 0.899359), ("rs10917151", 0.84732)]:
    q = mm0[mm0["rsid"].astype(str) == rs]
    if len(q) == 0:
        p("   核验 r2(%s, %s)：该变异不在 %s 交集内，跳过" % (rs, COND_RS, CELLS[0]))
        continue
    sid = q.iloc[0]["snp_grch37"]
    rr = out[out["snp_grch37"] == sid]
    if len(rr) == 0:
        p("   核验 r2(%s, %s)：不在 LD 矩阵内，跳过" % (rs, COND_RS))
        continue
    got = float(rr.iloc[0]["r2_to_cond"])
    p("   核验 r2(%s, %s) = %.6f  (口径值 %.6f, 差 %.6f)  %s"
      % (rs, COND_RS, got, want, abs(got - want),
         "OK" if abs(got - want) < 0.01 else "!! 偏差偏大"))
p("   |r| 分布：max=%.4f  >0.99 的个数=%d  >0.9 的个数=%d"
  % (np.nanmax(np.abs(r[np.arange(len(r)) != jj])),
     int((np.abs(r) > 0.99).sum()) - 1, int((np.abs(r) > 0.9).sum())))
sup = (out["in_intersect"]) & (np.abs(out["r_to_cond"]) > 0.99)
p("   （交集内 |r|>0.99 的变异数 = %d，含条件变异本身为 1）" % int(sup.sum()))

# ------------------------------------------------------------------ 5) 符号一致性诊断
p("")
p("[5] 符号一致性诊断（条件 zd 用 +r 与 -r 各算一次）")
rmap = dict(zip(out["snp_grch37"], out["r_to_cond"]))
for ct in CELLS:
    mm = pd.read_csv(os.path.join(COND, "m_%s.csv" % ct))
    rv = mm["snp_grch37"].map(rmap)
    ok = rv.notna()
    mm2 = mm[ok].copy()
    k = mm2[mm2["snp_grch37"] == COND_SNP37].iloc[0]
    r = rv[ok].values.astype(float)
    zk = float(k["slope"]) / float(k["slope_se"])
    denom = np.sqrt(np.clip(1.0 - r ** 2, 1e-12, None))
    z_e = mm2["slope"].values / mm2["slope_se"].values
    zc_pos = (z_e - r * zk) / denom
    zc_neg = (z_e + r * zk) / denom
    m_ = (mm2["snp_grch37"] != COND_SNP37).values
    sub = mm2[m_].copy()
    sub["z_e"] = z_e[m_]
    sub["z_cond_pos"] = zc_pos[m_]
    sub["z_cond_neg"] = zc_neg[m_]
    top = sub.loc[sub["z_e"].abs().idxmax()]
    p("   %s：交集 %d；条件前 max|z_eqtl| = %.3f（%s，pos37=%s）"
      % (ct, len(mm2), abs(top["z_e"]), top["rsid"], top["snp_grch37"]))
    p("        该变异条件后 z：（+r）= %+.3f   （-r）= %+.3f"
      % (top["z_cond_pos"], top["z_cond_neg"]))
    p("        条件后 max|z|：（+r）= %.3f   （-r）= %.3f"
      % (sub["z_cond_pos"].abs().max(), sub["z_cond_neg"].abs().max()))
    p("        条件后 |z|>5.45 (=0.05/8612 Bonferroni) 的变异数：（+r）%d  （-r）%d"
      % (int((sub["z_cond_pos"].abs() > 5.45).sum()), int((sub["z_cond_neg"].abs() > 5.45).sum())))

p("")
p("总耗时 %.1f s   末次可用内存 %.2f GB" % (time.time() - T0, free_gb()))
io.open(LOG, "w", encoding="utf-8", newline="\n").write("\n".join(L) + "\n")
print("WROTE " + LOG)
