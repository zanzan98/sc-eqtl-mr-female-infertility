# -*- coding: utf-8 -*-
"""s33_coloc_sensitivity_prep.py —— 为 coloc 敏感性分析准备 3 窗口 × 输入数据

针对 H4_verdict=strong 的 5 个位点，围绕 lead variant 准备 3 种窗口（500Kb/1Mb/2Mb）的：
  - eQTL 数据（从 OneK1K parquet 提取，含 af/slope/slope_se/pval_nominal）
  - GWAS 数据（从结局缓存 GCST90483463.h.tsv.gz 提取）
  - 等位基因映射（从 plink bim 提取）
输出到 00_data_raw/onek1k/coloc_sens_prep/，供 R 端 coloc.abf 调用。

窗口定义（GRCh37，围绕 lead）：
  500Kb = lead ± 250000
  1Mb   = lead ± 500000
  2Mb   = lead ± 1000000

红线遵守：仅 chr1/chr2 同染色体内扩窗口，不扩展 parquet_subset 的染色体范围。
"""
import os, gzip
import numpy as np
import pandas as pd
from pyliftover import LiftOver

PROJ = r"D:\endometriosis_project\11_sc_eqtl_mr_project"
PREP = os.path.join(PROJ, "00_data_raw", "onek1k", "coloc_prep")
SENS = os.path.join(PROJ, "00_data_raw", "onek1k", "coloc_sens_prep")
PQ   = os.path.join(PROJ, "00_data_raw", "onek1k", "parquet_subset")
GW   = os.path.join(PROJ, "00_data_raw", "gwas", "GCST90483463.h.tsv.gz")
BIM  = r"D:\endometriosis_project\_onek1k_plink\plink_merged_980_donors.bim"
os.makedirs(SENS, exist_ok=True)

lo37_38 = LiftOver('hg19', 'hg38')
lo38_37 = LiftOver('hg38', 'hg19')

# 5 strong 位点：gene, cell_type, lead_pos38, chrom
LOCI = [
    ("CDC42", "B_IN",   "1", 22132301),
    ("CDC42", "B_MEM",  "1", 22135618),
    ("CDC42", "Mono_C", "1", 22096228),
    ("CDC42", "Mono_NC","1", 22096228),
    ("ANXA4", "Mono_NC","2", 69729589),
]

WINDOWS = {"500Kb": 250000, "1Mb": 500000, "2Mb": 1000000}

# ---------- 1) 读 plink bim（全 chr1/chr2 等位基因映射）----------
bim = {}
print("读取 plink bim ...")
with open(BIM) as f:
    for L in f:
        c = L.rstrip().split("\t")
        ch, snp, _, pos, a1, a2 = c[0], c[1], c[2], int(c[3]), c[4], c[5]
        if ch in ("1", "2"):
            bim[(ch, pos)] = (snp, a1, a2)
print(f"bim 载入: chr1+chr2 共 {len(bim)} SNP")

# ---------- 2) 结局 GWAS 缓存（一次性预加载 chr1/chr2 目标大区间，避免 15 次全扫）----------
# 目标区间（GRCh38）：chr1 21.0-23.6Mb（覆盖 CDC42 所有窗口），chr2 68.4-71.1Mb（覆盖 ANXA4）
GW_PRELOAD = {"1": (20900000, 23700000), "2": (68300000, 71200000)}
_gwas_cache = {"1": [], "2": []}
print("预加载结局 GWAS 缓存（chr1+chr2 目标区间，单次流式扫）...")
with gzip.open(GW, "rt") as f:
    hdr = f.readline().rstrip("\n").split("\t")
    ix = {c: i for i, c in enumerate(hdr)}
    for line in f:
        c = line.rstrip("\n").split("\t")
        ch = c[ix["chromosome"]]
        if ch not in GW_PRELOAD:
            continue
        p = int(c[ix["base_pair_location"]])
        lo, hi = GW_PRELOAD[ch]
        if p < lo or p > hi:
            continue
        _gwas_cache[ch].append(dict(
            chrom=int(ch), pos38=p,
            ea=c[ix["effect_allele"]], oa=c[ix["other_allele"]],
            beta=float(c[ix["beta"]]), se=float(c[ix["standard_error"]]),
            eaf=float(c[ix["effect_allele_frequency"]]), p=float(c[ix["p_value"]]),
            rsid=c[ix["rsid"]], n_cases=int(c[ix["n_cases"]]), n_controls=int(c[ix["n_controls"]]),
        ))
_gwas_cache = {k: pd.DataFrame(v) for k, v in _gwas_cache.items()}
print(f"GWAS 预载: chr1={len(_gwas_cache['1'])} SNP, chr2={len(_gwas_cache['2'])} SNP")

def load_gwas(chrom, lo, hi):
    """从预载缓存切片提取 chrom 上 [lo,hi] 范围（pos 是 GRCh38）。"""
    df = _gwas_cache[str(chrom)]
    return df[(df["pos38"] >= lo) & (df["pos38"] <= hi)].copy()

# ---------- 3) eQTL 从 parquet 提取 ----------
def load_eqtl(cell_type, gene, chrom):
    f = os.path.join(PQ, f"OneK1K_{cell_type}.cis_qtl_pairs.chr{chrom}.parquet")
    df = pd.read_parquet(f, columns=["phenotype_id", "variant_id", "af", "slope", "slope_se", "pval_nominal"])
    df = df[df["phenotype_id"] == gene].copy()
    df["pos37"] = df["variant_id"].astype(str).str.split(":").str[1].astype(int)
    return df

# ---------- 4) liftover 一个 pos37 -> pos38 ----------
def lift37_38(chrom, pos37):
    r = lo37_38.convert_coordinate(f"chr{chrom}", pos37)
    return r[0][1] if r else None

summary = []
for gene, ct, chrom, lead38 in LOCI:
    chrom_i = int(chrom)
    # lead pos37
    r = lo38_37.convert_coordinate(f"chr{chrom}", lead38)
    lead37 = r[0][1] if r else None
    # 载入 eQTL 完整 ±1Mb
    eq = load_eqtl(ct, gene, chrom_i)
    for wname, half in WINDOWS.items():
        lo37 = lead37 - half; hi37 = lead37 + half
        # eQTL 按 pos37 过滤
        e = eq[(eq["pos37"] >= lo37) & (eq["pos37"] <= hi37)].copy()
        # liftover eQTL pos37 -> pos38
        e["pos38"] = e["pos37"].apply(lambda p: lift37_38(chrom_i, p))
        e = e.dropna(subset=["pos38"]).copy()
        e["pos38"] = e["pos38"].astype(int)
        # 等位基因映射
        e["snp_grch37"] = e["pos37"].map(lambda p: f"{chrom_i}:{p}")
        def get_alleles(p37):
            key = (chrom, p37)
            return bim.get(key, (None, None, None))
        e[["bim_snp", "a1", "a2"]] = e["pos37"].apply(lambda p: pd.Series(get_alleles(p)))
        e = e[e["a1"].notna()].copy()
        # GWAS 窗口（pos38）
        lo38 = lead38 - half; hi38 = lead38 + half
        g = load_gwas(chrom_i, lo38, hi38)
        # 保存
        tag = f"{gene}_{ct}_{wname}"
        e_out = e[["snp_grch37", "pos37", "pos38", "af", "slope", "slope_se", "pval_nominal", "a1", "a2"]].copy()
        e_out.insert(0, "gene", gene); e_out.insert(1, "cell_type", ct)
        e_out.to_csv(os.path.join(SENS, f"eqtl_{tag}.csv"), index=False, encoding="utf-8-sig")
        g_out = g.copy()
        g_out.to_csv(os.path.join(SENS, f"gwas_{tag}.csv"), index=False, encoding="utf-8-sig")
        summary.append(dict(gene=gene, cell_type=ct, window=wname,
                            lead37=lead37, n_eqtl=len(e), n_gwas=len(g)))
        print(f"{tag}: lead37={lead37} 窗口[{lo37},{hi37}] eqtl={len(e)} gwas={len(g)}")

sdf = pd.DataFrame(summary)
sdf.to_csv(os.path.join(SENS, "_summary.csv"), index=False, encoding="utf-8-sig")
print("\n=== 提取汇总 ===")
print(sdf.to_string(index=False))
print("\nDONE")
