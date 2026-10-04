# -*- coding: utf-8 -*-
"""s51a2_coord_forensics.py -- 判定 parquet variant_id 的坐标体系，并定位真正的工具变量

背景：31_/15_ 表把 instrument 记为 `10:27154694` / `2:69729589`（GRCh38），
      而 parquet_subset 的 variant_id 疑似 GRCh37（`1:<pos37>`，与 _onek1k_plink bim 对齐）。
      本脚本用**已知答案**的 chr1 座（L1）做标定：
        CDC42/CD4_NC 的权威记录 =  GRCh38 `1:22088292` ↔ GRCh37 `1:22414785`
        LINC00339/NK  = GRCh38 `1:22031964` ↔ GRCh37 `1:22358457`
      哪个字符串能在 parquet 里以正确统计量命中，就说明 parquet 是哪个体系。

同时输出 L2/L3 的真·工具变量（parquet 坐标 + liftover 到另一体系）。
只读，不写任何已交付产物。
"""
import io
import os

import numpy as np
import pandas as pd
from pyliftover import LiftOver

PROJ = r"D:\endometriosis_project\11_sc_eqtl_mr_project"
PQ = os.path.join(PROJ, "00_data_raw", "onek1k", "parquet_subset")
TAB = os.path.join(PROJ, "tables")
BIM = r"D:\endometriosis_project\_onek1k_plink\plink_merged_980_donors.bim"
LOG = r"D:\endometriosis_project\_s51a2_forensics.txt"

L = []


def p(s=""):
    L.append(str(s))
    print(s)


lo37_38 = LiftOver("hg19", "hg38")
lo38_37 = LiftOver("hg38", "hg19")


def read_pq(ct, chrom, gene, cols=None):
    f = os.path.join(PQ, "OneK1K_%s.cis_qtl_pairs.chr%s.parquet" % (ct, chrom))
    df = pd.read_parquet(f)
    e = df[df["phenotype_id"].astype(str) == gene].copy()
    e["vid"] = e["variant_id"].astype(str)
    return e


p("=== s51a2 坐标体系取证 ===")
p("time-free (no network)")

# ---------------------------------------------------------------- 标定：chr1
p("")
p("### [标定] chr1 座：parquet 里 1:22088292 (GRCh38) 还是 1:22414785 (GRCh37) ?")
e1 = read_pq("CD4_NC", "1", "CDC42")
p("  CDC42/CD4_NC parquet 变异数 = %d  vid 前缀 = %s" % (len(e1), e1["vid"].iloc[0].rsplit(":", 1)[0]))
for v in ["1:22088292", "1:22414785"]:
    q = e1[e1["vid"] == v]
    if len(q):
        r = q.iloc[0]
        p("  %-14s 命中  p=%.4e  slope=%+.6f se=%.6f  |z|=%.3f" %
          (v, r["pval_nominal"], r["slope"], r["slope_se"], abs(r["slope"] / r["slope_se"])))
    else:
        p("  %-14s 不存在" % v)
p("  权威（36_finemap_pair.csv CDC42_CD4_NC）: GRCh38=1:22088292  GRCh37=1:22414785  "
  "01b eQTL beta/se ?")

e1b = read_pq("NK", "1", "LINC00339")
p("  LINC00339/NK parquet 变异数 = %d" % len(e1b))
for v in ["1:22031964", "1:22358457"]:
    q = e1b[e1b["vid"] == v]
    if len(q):
        r = q.iloc[0]
        p("  %-14s 命中  p=%.4e  slope=%+.6f se=%.6f" % (v, r["pval_nominal"], r["slope"], r["slope_se"]))
    else:
        p("  %-14s 不存在" % v)

# 检查 01b 里 CDC42/CD4_NC 的一致行
p("")
p("### [标定] 01b / 15_ / 31_ 三表对 CDC42·CD4_NC 的记录")
b01 = pd.read_csv(os.path.join(TAB, "01b_discovery_instruments_grch38.csv"))
sel = b01[(b01["gene"] == "CDC42") & (b01["cell_type"] == "CD4_NC")]
p(b01.columns.tolist())
p(sel.to_string(index=False))
t15 = pd.read_csv(os.path.join(TAB, "15_discovery_significant_with_locus.csv"))
p(t15[(t15["gene"] == "CDC42") & (t15["cell_type"] == "CD4_NC")][
    ["locus_id", "gene", "cell_type", "b", "se", "F_stat", "eQTL_pval_nominal",
     "variant_id_grch38", "tss_pos_grch38", "tss_distance"]].to_string(index=False))

# ---------------------------------------------------------------- L2 / L3
p("")
p("### [L2] YME1L1 / CD4_NC / chr10")
e2 = read_pq("CD4_NC", "10", "YME1L1")
imin = int(e2["pval_nominal"].idxmin())
r = e2.loc[imin]
p("  parquet 最小 p 变异 = %s  p=%.6e slope=%+.6f se=%.6f F=%.4f"
  % (r["vid"], r["pval_nominal"], r["slope"], r["slope_se"], (r["slope"] / r["slope_se"]) ** 2))
p("  01b 记录:")
s01 = b01[(b01["gene"] == "YME1L1") & (b01["cell_type"] == "CD4_NC")]
p("    %s" % (s01[["variant_id", "beta", "se", "pval_nominal", "F", "pos37", "pos38",
                   "variant_id_grch37"]].to_string(index=False) if len(s01) else "    (无)"))
p("  15_ 记录:")
s15 = t15[(t15["gene"] == "YME1L1") & (t15["cell_type"] == "CD4_NC")]
p("    %s" % (s15[["variant_id_grch38", "tss_pos_grch38", "tss_distance", "F_stat",
                   "eQTL_pval_nominal"]].to_string(index=False) if len(s15) else "    (无)"))
for v in ["10:27443623", "10:27154694"]:
    q = e2[e2["vid"] == v]
    p("    检查 %-14s -> %s" % (v, ("命中 p=%.4e slope=%+.6f" % (q.iloc[0]["pval_nominal"], q.iloc[0]["slope"])) if len(q) else "不存在"))
lt = lo37_38.convert_coordinate("chr10", 27443623)
p("    liftover chr10:27443623 GRCh37 -> GRCh38 = %s" % (lt[0][1] if lt else "FAIL"))
lt2 = lo38_37.convert_coordinate("chr10", 27154694)
p("    liftover chr10:27154694 GRCh38 -> GRCh37 = %s" % (lt2[0][1] if lt2 else "FAIL"))

p("")
p("### [L3] ANXA4 / Mono_NC / chr2")
e3 = read_pq("Mono_NC", "2", "ANXA4")
imin = int(e3["pval_nominal"].idxmin())
r = e3.loc[imin]
p("  parquet 最小 p 变异 = %s  p=%.6e slope=%+.6f se=%.6f F=%.4f"
  % (r["vid"], r["pval_nominal"], r["slope"], r["slope_se"], (r["slope"] / r["slope_se"]) ** 2))
s01 = b01[(b01["gene"] == "ANXA4") & (b01["cell_type"] == "Mono_NC")]
p("  01b 记录:")
p("    %s" % (s01[["variant_id", "beta", "se", "pval_nominal", "F", "pos37", "pos38",
                   "variant_id_grch37"]].to_string(index=False) if len(s01) else "    (无)"))
s15 = t15[(t15["gene"] == "ANXA4") & (t15["cell_type"] == "Mono_NC")]
p("  15_ 记录:")
p("    %s" % (s15[["variant_id_grch38", "tss_pos_grch38", "tss_distance", "F_stat",
                   "eQTL_pval_nominal"]].to_string(index=False) if len(s15) else "    (无)"))
for v in ["2:69956721", "2:69729589", "2:69946670", "2:69719538"]:
    q = e3[e3["vid"] == v]
    p("    检查 %-14s -> %s" % (v, ("命中 p=%.4e slope=%+.6f" % (q.iloc[0]["pval_nominal"], q.iloc[0]["slope"])) if len(q) else "不存在"))
lt = lo37_38.convert_coordinate("chr2", 69956721)
p("    liftover chr2:69956721 GRCh37 -> GRCh38 = %s" % (lt[0][1] if lt else "FAIL"))
lt2 = lo38_37.convert_coordinate("chr2", 69729589)
p("    liftover chr2:69729589 GRCh38 -> GRCh37 = %s" % (lt2[0][1] if lt2 else "FAIL"))
lt3 = lo37_38.convert_coordinate("chr2", 69946670)
p("    liftover chr2:69946670 GRCh37 -> GRCh38 = %s" % (lt3[0][1] if lt3 else "FAIL"))

# ---------------------------------------------------------------- bim 侧
p("")
p("### [bim] _onek1k_plink 是否含这些 SNP（bim SNP ID 体系）")
bim = pd.read_csv(BIM, sep="\t", header=None, dtype={0: str},
                  names=["chr", "snp", "cm", "pos", "a1", "a2"])
p("  bim 总行 = %d   chr10 行 = %d   chr2 行 = %d"
  % (len(bim), (bim["chr"] == "10").sum(), (bim["chr"] == "2").sum()))
for v in ["10:27443623", "10:27154694", "2:69956721", "2:69729589", "2:69946670"]:
    q = bim[bim["snp"] == v]
    p("    bim 含 %-14s -> %s" % (v, ("是 (pos=%d, a1=%s, a2=%s)" % (q.iloc[0]["pos"], q.iloc[0]["a1"], q.iloc[0]["a2"])) if len(q) else "否"))

io.open(LOG, "w", encoding="utf-8", newline="\n").write("\n".join(L) + "\n")
print("WROTE " + LOG)
