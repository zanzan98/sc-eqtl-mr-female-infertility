# -*- coding: utf-8 -*-
"""chr1 基因座 LD 结构数据构建
1) 探查 OneK1K bim 的 SNP ID 格式
2) 用 plink 在 chr1:22.25-22.57 Mb (GRCh37) 内计算两两 r2
3) 把 OneK1K GRCh37 位置 liftover 到 GRCh38（复用已缓存 hg19ToHg38 chain）
4) 识别结局 top 变异（GRCh38 chr1:22,139,327）在 OneK1K 中的对应变异
5) 扫描结局文件 chr1 区域（文件按 chr+pos 排序，可从前端截读）取 p 值
6) 输出绘图数据 _chr1_plotdata.csv 与 _chr1_matrix.csv
"""
import gzip, os, subprocess, sys, datetime
import numpy as np
import pandas as pd
from pyliftover import LiftOver

ROOT = r"D:\endometriosis_project\11_sc_eqtl_mr_project"
PLINK = r"D:\endometriosis_project\_tools\plink.exe"
BFILE = r"D:\endometriosis_project\_onek1k_plink\plink_merged_980_donors"
OUTCOME = os.path.join(ROOT, "00_data_raw", "gwas", "GCST90483463.h.tsv.gz")
WORK = r"D:\endometriosis_project\_chr1_ld"
os.makedirs(WORK, exist_ok=True)

INDEX_GRCH38 = 22139327          # 结局 top 变异（GWAS Catalog harmonised, GRCh38）
WIN37 = (22250000, 22570000)     # GRCh37 窗口
INSTRUMENTS_GRCH37 = ["1:22358457", "1:22439520", "1:22414785",
                      "1:22462111", "1:22458794", "1:22422721"]

LOG = []


def w(s):
    LOG.append(str(s))


# ---------- 1) bim 探查 ----------
w("=== 1) OneK1K bim 探查 @ %s ===" % datetime.datetime.now().strftime("%H:%M:%S"))
with open(BFILE + ".bim", "r", encoding="utf-8", errors="replace") as f:
    head = [next(f).rstrip("\n") for _ in range(3)]
for h in head:
    w("   " + h)

# ---------- 2) plink 两两 r2 ----------
w("")
w("=== 2) plink 两两 r2 (chr1 %d-%d, GRCh37) ===" % WIN37)
cmd = [PLINK, "--bfile", BFILE,
       "--chr", "1", "--from-bp", str(WIN37[0]), "--to-bp", str(WIN37[1]),
       "--r2", "--ld-window-kb", "2000", "--ld-window", "99999", "--ld-window-r2", "0",
       "--out", os.path.join(WORK, "chr1_r2")]
try:
    p = subprocess.run(cmd, capture_output=True, text=True, timeout=1800)
    w("   exit = %d" % p.returncode)
    for ln in (p.stdout or "").splitlines()[-25:]:
        w("   | " + ln)
    for ln in (p.stderr or "").splitlines()[-10:]:
        w("   ! " + ln)
except Exception as ex:
    w("   plink 调用失败: %s: %s" % (type(ex).__name__, ex))

r2f = os.path.join(WORK, "chr1_r2.ld")
if not os.path.exists(r2f):
    w("   !! 未生成 %s" % r2f)
    open(r"D:\endometriosis_project\_chr1_ld.txt", "w", encoding="utf-8").write("\n".join(LOG))
    sys.exit(1)

r2 = pd.read_csv(r2f, sep=r"\s+")
w("   r2 对数 = %d ; 列 = %s" % (len(r2), list(r2.columns)))

# ---------- 3) liftover GRCh37 -> GRCh38 ----------
w("")
w("=== 3) liftover GRCh37 -> GRCh38 ===")
lo = LiftOver('hg19', 'hg38')
allpos = sorted(set(r2["BP_A"]).union(set(r2["BP_B"])))
p37_to_38 = {}
failed = []
for pp in allpos:
    res = lo.convert_coordinate('chr1', int(pp))
    if res:
        p37_to_38[int(pp)] = int(res[0][1])
    else:
        failed.append(pp)
w("   唯一位置 = %d ; liftover 成功 = %d (%.2f%%) ; 失败 = %d"
  % (len(allpos), len(p37_to_38), 100.0 * len(p37_to_38) / max(len(allpos), 1), len(failed)))

r2["BP_A_38"] = r2["BP_A"].map(p37_to_38)
r2["BP_B_38"] = r2["BP_B"].map(p37_to_38)

# ---------- 4) 索引变异 ----------
w("")
w("=== 4) 索引变异（结局 top, GRCh38 chr1:%d）===" % INDEX_GRCH38)
cand = pd.DataFrame({"p37": list(p37_to_38.keys()), "p38": list(p37_to_38.values())})
cand["dist"] = (cand["p38"] - INDEX_GRCH38).abs()
cand = cand.sort_values("dist")
w(cand.head(5).to_string(index=False))
best = cand.iloc[0]
idx_p37 = int(best["p37"])
idx_p38 = int(best["p38"])
w("   选定 OneK1K 变异: GRCh37:%d <-> GRCh38:%d (距离 %d bp)" % (idx_p37, idx_p38, int(best["dist"])))

# 索引 vs 区域内全部
sel = r2[(r2["SNP_A"] == "1:%d" % idx_p37) | (r2["SNP_B"] == "1:%d" % idx_p37)].copy()
w("   索引相关对数 = %d" % len(sel))
if len(sel):
    sel["other_38"] = np.where(sel["SNP_A"] == "1:%d" % idx_p37, sel["BP_B_38"], sel["BP_A_38"])
    sel["other_37"] = np.where(sel["SNP_A"] == "1:%d" % idx_p37, sel["BP_B"], sel["BP_A"])
    sel = sel[["other_37", "other_38", "R2"]].rename(columns={"R2": "r2_index"})

# ---------- 5) 结局区域 p 值 ----------
w("")
w("=== 5) 结局区域扫描 (GRCh38 chr1 %.2f-%.2f Mb) ===" % (21.8, 22.6))
ALIAS = {"chrom": ["chromosome", "chr"], "pos": ["base_pair_location"], "ea": ["effect_allele"],
         "oa": ["other_allele"], "beta": ["beta"], "se": ["standard_error"],
         "eaf": ["effect_allele_frequency"], "p": ["p_value"], "rsid": ["rsid", "variant_id"]}
LO, HI = 21_800_000, 22_600_000
rows = []
nread = 0
with gzip.open(OUTCOME, "rt", encoding="utf-8", errors="replace") as f:
    hdr = f.readline().rstrip("\n").split("\t")
    h = [c.strip().lower().lstrip("#") for c in hdr]
    ix = {}
    for k, alts in ALIAS.items():
        for a in alts:
            if a in h:
                ix[k] = h.index(a)
                break
    w("   列映射 = %s" % ix)
    for line in f:
        nread += 1
        f2 = line.rstrip("\n").split("\t")
        if len(f2) < 6:
            continue
        try:
            chrom = f2[ix["chrom"]]
            pos = int(f2[ix["pos"]])
        except Exception:
            continue
        if chrom not in ("1", "chr1"):
            if nread > 1000:
                break
            continue
        if pos > HI:
            break
        if pos < LO:
            continue
        def fnum(x):
            try:
                return float(x)
            except Exception:
                return np.nan
        rows.append(dict(pos_38=pos, rsid=f2[ix["rsid"]] if "rsid" in ix else "",
                         p=fnum(f2[ix["p"]]), beta=fnum(f2[ix["beta"]]),
                         ea=f2[ix["ea"]], oa=f2[ix["oa"]],
                         eaf=fnum(f2[ix["eaf"]]) if "eaf" in ix else np.nan))
region = pd.DataFrame(rows)
region["neglogp"] = -np.log10(region["p"].clip(lower=1e-300))
w("   读取行数 = %d ; 区域内变异 = %d ; p 范围 %.3g ~ %.3g"
  % (nread, len(region), region["p"].min(), region["p"].max()))

region = region.merge(sel, left_on="pos_38", right_on="other_38", how="left")
w("   有 LD 值的变异 = %d / %d" % (region["r2_index"].notna().sum(), len(region)))

# ---------- 6) 7x7 矩阵 ----------
w("")
w("=== 6) 索引+6 工具的 r2 矩阵 ===")
ids37 = INSTRUMENTS_GRCH37 + ["1:%d" % idx_p37]
lab38 = {v: p37_to_38.get(int(v.split(":")[1]), np.nan) for v in ids37}
mat = pd.DataFrame(np.nan, index=ids37, columns=ids37)
for v in ids37:
    mat.loc[v, v] = 1.0
for _, r in r2.iterrows():
    if r["SNP_A"] in mat.index and r["SNP_B"] in mat.index:
        mat.loc[r["SNP_A"], r["SNP_B"]] = r["R2"]
        mat.loc[r["SNP_B"], r["SNP_A"]] = r["R2"]
mat.to_csv(os.path.join(WORK, "chr1_matrix.csv"), encoding="utf-8-sig")
w(mat.round(3).to_string())
w("")
w("标签映射 (GRCh37 id -> GRCh38 pos): %s" % lab38)
w("矩阵缺失值数 = %d（缺失 = 该变异不在 OneK1K 或超出窗口）" % int(mat.isna().sum().sum()))

region.to_csv(os.path.join(WORK, "chr1_plotdata.csv"), index=False, encoding="utf-8-sig")
with open(r"D:\endometriosis_project\_chr1_ld.txt", "w", encoding="utf-8") as fh:
    fh.write("\n".join(LOG))
print("done")
