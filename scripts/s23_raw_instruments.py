# -*- coding: utf-8 -*-
"""s23_raw_instruments.py —— 从 OneK1K raw parquet 提取三座工具变量（P<5e-8）

裁定 (a)：多 SNP 工具筛选统一用主分析 **P<5e-8**（raw 无 qval）+ clumping r²<0.001 / 10 Mb。
本脚本只做「提取候选工具 + 记录每细胞类型的 cis 变异总数与最小 p」，clumping 由 PLINK 单独做。

输出：
  tables/27_raw_instruments_preclump.csv   （gene, cell_type, variant_id_grch37(chr:pos), ...）
  tables/27b_raw_cis_totals.csv            （gene, cell_type, n_cis_variants, min_pval, n_P5e8）
"""
import os, io, re, sys, time
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from pyliftover import LiftOver

PROJ = r"D:\endometriosis_project\11_sc_eqtl_mr_project"
SUB = os.path.join(PROJ, "00_data_raw", "onek1k", "parquet_subset")
OUT_IV = os.path.join(PROJ, "tables", "27_raw_instruments_preclump.csv")
OUT_TOT = os.path.join(PROJ, "tables", "27b_raw_cis_totals.csv")
LOG = r"D:\endometriosis_project\_s23_ivw.log"
PTHR = 5e-8

# 座 -> 染色体 -> 目标基因
LOCI = {
    "1": ["CDC42", "LINC00339"],
    "2": ["ANXA4"],
    "10": ["YME1L1"],
}

L = []
def p(s=""):
    L.append(str(s)); print(s)

T0 = time.time()
p("=== s23 从 raw 提取三座工具变量（P<5e-8） ===")
p("时间 %s" % time.strftime("%Y-%m-%d %H:%M:%S"))

lo = LiftOver("hg19", "hg38")   # ★ 只用正向链（裁定：拒绝反向链）

rows = []
tot_rows = []
files = sorted(os.listdir(SUB))
for fn in files:
    if not (fn.endswith(".parquet") and "cis_qtl_pairs" in fn):
        continue
    m = re.search(r"chr([0-9XYM]+)\.parquet$", fn)
    if not m or m.group(1) not in LOCI:
        continue
    chrm = m.group(1)
    ct = fn.split(".cis_qtl_pairs")[0].replace("OneK1K_", "")
    target = set(LOCI[chrm])
    fp = os.path.join(SUB, fn)
    t = pq.read_table(fp, columns=["phenotype_id", "variant_id", "tss_distance", "af",
                                   "slope", "slope_se", "pval_nominal", "ma_samples", "ma_count"])
    df = t.to_pandas()
    df = df[df["phenotype_id"].isin(target)]
    for gene, g in df.groupby("phenotype_id"):
        n_cis = len(g)
        minp = float(g["pval_nominal"].min())
        gg = g[g["pval_nominal"] < PTHR].copy()
        tot_rows.append(dict(gene=gene, cell_type=ct, chrom=chrm,
                             n_cis_variants=n_cis, min_pval=minp, n_P5e8=len(gg)))
        for r in gg.itertuples(index=False):
            rows.append(dict(gene=gene, cell_type=ct, chrom=chrm, chrom_pos_grch37=r.variant_id,
                             tss_distance=r.tss_distance, af=r.af, beta=r.slope, se=r.slope_se,
                             pval_nominal=r.pval_nominal, ma_samples=r.ma_samples, ma_count=r.ma_count))
    p("  %-46s 目标基因行=%d" % (fn, len(df)))

raw_iv = pd.DataFrame(rows)
p("")
p("P<5e-8 候选工具行数 = %d" % len(raw_iv))
p("按基因：")
p(raw_iv.groupby("gene").size().to_string())

# ---- 正向 liftover GRCh37 -> GRCh38（拒绝反向链）----
p("")
p("正向 liftover hg19->hg38 …")
uniq = raw_iv[["chrom", "chrom_pos_grch37"]].drop_duplicates().copy()
uniq["pos37"] = uniq["chrom_pos_grch37"].str.split(":").str[1].astype(int)
mapped = {}
for r in uniq.itertuples(index=False):
    res = lo.convert_coordinate("chr" + str(r.chrom), int(r.pos37) - 1)
    mapped[r.chrom_pos_grch37] = ("%s:%d" % (r.chrom, res[0][1] + 1)) if res else None
raw_iv["variant_id_grch38"] = raw_iv["chrom_pos_grch37"].map(mapped)
n_ok = raw_iv["variant_id_grch38"].notna().sum()
p("  唯一位点 = %d ; 映射成功 = %d ; 失败 = %d"
  % (len(uniq), n_ok, len(uniq) - n_ok))

raw_iv = raw_iv.rename(columns={"chrom_pos_grch37": "variant_id_grch37"})
raw_iv["F_approx"] = np.nan  # 见下：用 af/beta 近似，公式单独写
tau = 2 * raw_iv["af"] * (1 - raw_iv["af"])
R2 = (tau * raw_iv["beta"] ** 2) / (1.0 + tau * raw_iv["beta"] ** 2)
raw_iv["F_approx"] = (980 - 2) * R2 / (1 - R2)

raw_iv = raw_iv[["gene", "cell_type", "chrom", "variant_id_grch37", "variant_id_grch38",
                 "tss_distance", "af", "beta", "se", "pval_nominal", "F_approx",
                 "ma_samples", "ma_count"]].sort_values(
    ["chrom", "gene", "cell_type", "pval_nominal"]).reset_index(drop=True)

os.makedirs(os.path.dirname(OUT_IV), exist_ok=True)
raw_iv.to_csv(OUT_IV, index=False, encoding="utf-8-sig")
pd.DataFrame(tot_rows).sort_values(["chrom", "gene", "cell_type"]).to_csv(
    OUT_TOT, index=False, encoding="utf-8-sig")

p("")
p("已写出 %s (%d 行)" % (OUT_IV, len(raw_iv)))
p("已写出 %s (%d 行)" % (OUT_TOT, len(tot_rows)))
p("")
p("=== 每 gene×cell_type 的候选工具数分布 ===")
cnt = raw_iv.groupby(["gene", "cell_type"]).size().rename("n_P5e8").reset_index()
p(cnt.to_string(index=False))
p("")
p("总耗时 %.1f min" % ((time.time() - T0) / 60))
open(LOG, "w", encoding="utf-8").write("\n".join(L) + "\n")
