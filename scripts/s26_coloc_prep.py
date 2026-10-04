# -*- coding: utf-8 -*-
"""s26_coloc_prep.py —— 共定位数据准备（Python 侧）

裁定与红线：
  * 坐标系统一用 **GRCh37（eQTL 侧为准）**，liftover **只用正向链 hg19->hg38**（拒绝反向链恢复 GRCh37）。
  * LD 窗口从 OneK1K **bim** 按 `chr:pos` 取（GRCh37，与 eQTL 同源同 build）。
  * 与结局连接走 **rsID**（结局侧提供 rsid + GRCh38 位置）；本脚本用「位置直连」仅作候选，
    最终以 rsID 一致为判据落表。

产出（目录 00_data_raw/onek1k/coloc_prep/）：
  gwas_L1.csv / gwas_L2.csv / gwas_L3.csv        结局窗口（GRCh38）
  eqtl_<locus>_<gene>.csv                         长表：14 细胞类型的 eQTL cis 变异
  bim_universe_L<k>.csv                           窗口内 bim 变异（chr:pos, A1, A2）
  window.json                                     各座窗口（GRCh37/GRCh38 双向留痕）
"""
import gzip, json, os, re, time
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from pyliftover import LiftOver

PROJ = r"D:\endometriosis_project\11_sc_eqtl_mr_project"
SUB = os.path.join(PROJ, "00_data_raw", "onek1k", "parquet_subset")
BIM = os.path.join(SUB, "plink_merged.bim")
GWAS = os.path.join(PROJ, "00_data_raw", "gwas", "GCST90483463.h.tsv.gz")
OUT = os.path.join(PROJ, "00_data_raw", "onek1k", "coloc_prep")
LOG = r"D:\endometriosis_project\_s26_prep.log"
os.makedirs(OUT, exist_ok=True)

# 座定义：GRCh37 区域（含基因体 + 两侧余量）+ GRCh38 结局抽取窗口（宽松）
LOCI = {
    "L1": dict(chrom="1", lo37=22_000_000, hi37=22_900_000,
               lo38=21_550_000, hi38=22_700_000, genes=["CDC42", "LINC00339"]),
    "L2": dict(chrom="10", lo37=26_900_000, hi37=27_900_000,
               lo38=26_500_000, hi38=27_750_000, genes=["YME1L1"]),
    "L3": dict(chrom="2", lo37=69_300_000, hi37=70_200_000,
               lo38=68_950_000, hi38=70_100_000, genes=["ANXA4"]),
}

L = []
def p(s=""):
    L.append(str(s)); print(s)

T0 = time.time()
p("=== s26 共定位数据准备 ===")
p("时间 %s" % time.strftime("%Y-%m-%d %H:%M:%S"))
lo = LiftOver("hg19", "hg38")   # ★ 仅正向链

# ---------- 1) bim 窗口宇宙 ----------
p("")
p("[1] 从 bim 取各座窗口内变异（chr:pos, A1, A2；GRCh37）")
bim = pd.read_csv(BIM, sep="\t", header=None, dtype={0: str},
                  names=["chr", "snp", "cm", "pos", "a1", "a2"])
bim["chr"] = bim["chr"].astype(str).str.replace("chr", "", regex=False)
universe = {}
for k, v in LOCI.items():
    s = bim[(bim["chr"] == v["chrom"]) & (bim["pos"] >= v["lo37"]) & (bim["pos"] <= v["hi37"])].copy()
    s = s[s["snp"].astype(str).str.startswith(str(v["chrom"]) + ":")]
    universe[k] = s
    sp = os.path.join(OUT, "bim_universe_%s.csv" % k)
    s[["snp", "a1", "a2", "pos"]].to_csv(sp, index=False, encoding="utf-8-sig")
    p("   %s chr%s:%d-%d(GRCh37)  变异=%d -> %s" % (k, v["chrom"], v["lo37"], v["hi37"], len(s), os.path.basename(sp)))

# ---------- 2) eQTL：逐基因 × 14 细胞类型 ----------
p("")
p("[2] eQTL cis 变异（raw parquet）")
gene2locus = {}
for k, v in LOCI.items():
    for g in v["genes"]:
        gene2locus[g] = k

for gene, k in gene2locus.items():
    v = LOCI[k]
    rows = []
    files = [f for f in sorted(os.listdir(SUB))
             if f.endswith(".parquet") and re.search(r"chr%s\.parquet$" % v["chrom"], f)]
    for fn in files:
        ct = fn.split(".cis_qtl_pairs")[0].replace("OneK1K_", "")
        df = pq.read_table(os.path.join(SUB, fn),
                           columns=["phenotype_id", "variant_id", "tss_distance", "af",
                                    "slope", "slope_se", "pval_nominal"]).to_pandas()
        g = df[df["phenotype_id"] == gene]
        if g.empty:
            continue
        g = g[(g["variant_id"].str.split(":").str[1].astype(int) >= v["lo37"]) &
              (g["variant_id"].str.split(":").str[1].astype(int) <= v["hi37"])]
        g = g.assign(cell_type=ct)
        rows.append(g)
    e = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()
    e = e.rename(columns={"phenotype_id": "gene", "variant_id": "snp_grch37"})
    # 正向 liftover（每个唯一位置只算一次）
    upos = e["snp_grch37"].drop_duplicates().tolist()
    mp = {}
    for s in upos:
        pp = int(s.split(":")[1])
        r = lo.convert_coordinate("chr" + v["chrom"], pp - 1)
        mp[s] = ("%s:%d" % (v["chrom"], r[0][1] + 1)) if r else None
    e["snp_grch38"] = e["snp_grch37"].map(mp)
    tss37 = float((e["snp_grch37"].str.split(":").str[1].astype(int) - e["tss_distance"]).mode().iloc[0])
    e["tss_grch37"] = tss37
    sp = os.path.join(OUT, "eqtl_%s_%s.csv" % (k, gene))
    e.to_csv(sp, index=False, encoding="utf-8-sig")
    p("   %-10s(%s) 行=%7d  细胞类型=%2d  cis唯一位置=%5d  liftover成功=%5d  TSS(GRCh37)=%d"
      % (gene, k, len(e), e["cell_type"].nunique(), len(upos),
         int(e["snp_grch38"].notna().sum() and e.loc[e["snp_grch38"].notna(), "snp_grch37"].nunique()), tss37))
    p("      -> %s" % os.path.basename(sp))

# ---------- 3) 结局窗口（GRCh38） ----------
p("")
p("[3] 结局 GCST90483463 三座窗口抽取（一次流式扫描）")
want = {}
for k, v in LOCI.items():
    want[v["chrom"]] = (v["lo38"], v["hi38"])
hdr = None
recs = {k: [] for k in LOCI}
with gzip.open(GWAS, "rt", errors="replace") as fh:
    first = fh.readline()
    hdr = first.lstrip("#").rstrip("\n").split("\t")
    idx = {c: i for i, c in enumerate(hdr)}
    need = ["chromosome", "base_pair_location", "effect_allele", "other_allele",
            "beta", "standard_error", "effect_allele_frequency", "p_value", "rsid",
            "n_cases", "n_controls"]
    miss = [c for c in need if c not in idx]
    if miss:
        raise SystemExit("FATAL 结局缺列 %s；实际=%s" % (miss, hdr))
    n_line = 0
    for ln in fh:
        n_line += 1
        f = ln.rstrip("\n").split("\t")
        if len(f) <= idx["n_controls"]:
            continue
        c = f[idx["chromosome"]].replace("chr", "")
        if c not in want:
            continue
        try:
            pos = int(float(f[idx["base_pair_location"]]))
        except ValueError:
            continue
        lo38, hi38 = want[c]
        if not (lo38 <= pos <= hi38):
            continue
        recs[[k for k in LOCI if LOCI[k]["chrom"] == c][0]].append(
            dict(chrom=c, pos38=pos, ea=f[idx["effect_allele"]].upper(), oa=f[idx["other_allele"]].upper(),
                 beta=f[idx["beta"]], se=f[idx["standard_error"]],
                 eaf=f[idx["effect_allele_frequency"]], p=f[idx["p_value"]],
                 rsid=f[idx["rsid"]], n_cases=f[idx["n_cases"]], n_controls=f[idx["n_controls"]]))
    p("   扫描行数 = %d" % n_line)

for k in LOCI:
    d = pd.DataFrame(recs[k])
    for c in ["beta", "se", "eaf", "p"]:
        d[c] = pd.to_numeric(d[c], errors="coerce")
    d["n_cases"] = pd.to_numeric(d["n_cases"], errors="coerce")
    d["n_controls"] = pd.to_numeric(d["n_controls"], errors="coerce")
    sp = os.path.join(OUT, "gwas_%s.csv" % k)
    d.to_csv(sp, index=False, encoding="utf-8-sig")
    rs = d["rsid"].astype(str).str.startswith("rs").sum()
    p("   %s 窗口行=%6d  含 rsid 行=%6d  beta/se 非空=%6d -> %s"
      % (k, len(d), int(rs), int(d[["beta", "se"]].notna().all(axis=1).sum()), os.path.basename(sp)))

win = {k: dict(chrom=v["chrom"], lo37=v["lo37"], hi37=v["hi37"], lo38=v["lo38"], hi38=v["hi38"],
               genes=v["genes"]) for k, v in LOCI.items()}
with open(os.path.join(OUT, "window.json"), "w", encoding="utf-8") as fh:
    json.dump(win, fh, ensure_ascii=False, indent=1)

p("")
p("总耗时 %.1f min" % ((time.time() - T0) / 60))
open(LOG, "w", encoding="utf-8").write("\n".join(L) + "\n")
