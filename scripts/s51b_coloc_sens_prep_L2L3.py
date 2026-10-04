# -*- coding: utf-8 -*-
"""s51b_coloc_sens_prep_L2L3.py -- 任务A：为 L2 (YME1L1/CD4_NC) 生成 coloc 敏感性输入

严格复用 s33_coloc_sensitivity_prep.py 的口径与代码路径，仅两点改动：
  1) LOCI 增加 L2 位点 (YME1L1, CD4_NC, chr10, lead38=27154694)；
     并保留 L3 (ANXA4, Mono_NC, chr2, lead38=69729589) 以便同口径重生成（写入独立文件名，
     与原 s33 产物**逐字节应为同一内容**，从而验证管线可复现）。
  2) GW_PRELOAD 增加 chr10；输出文件名加 `_L2` 后缀的汇总，**不覆盖** `_summary.csv`。

窗口（GRCh37，围绕 lead）：500Kb = ±250kb，1Mb = ±500kb，2Mb = ±1000kb
只读已有数据，不下载。
"""
import gzip
import io
import os

import numpy as np
import pandas as pd
from pyliftover import LiftOver

PROJ = r"D:\endometriosis_project\11_sc_eqtl_mr_project"
SENS = os.path.join(PROJ, "00_data_raw", "onek1k", "coloc_sens_prep")
SENS_RECHECK = os.path.join(PROJ, "00_data_raw", "onek1k", "coloc_sens_prep_recheck")
PQ = os.path.join(PROJ, "00_data_raw", "onek1k", "parquet_subset")
GW = os.path.join(PROJ, "00_data_raw", "gwas", "GCST90483463.h.tsv.gz")
BIM = r"D:\endometriosis_project\_onek1k_plink\plink_merged_980_donors.bim"
LOG = r"D:\endometriosis_project\_s51b_prep_L2L3.log"
os.makedirs(SENS, exist_ok=True)
os.makedirs(SENS_RECHECK, exist_ok=True)

L = []


def p(s=""):
    L.append(str(s))
    print(s)


lo37_38 = LiftOver("hg19", "hg38")
lo38_37 = LiftOver("hg38", "hg19")

# (gene, cell_type, chrom, lead_pos38, outdir)
#   L2 -> 写入 coloc_sens_prep（新文件，不覆盖任何既有产物）
#   L3 -> 写入 coloc_sens_prep_recheck（**不碰** s33 已产出的原文件；用于同口径可复现核验）
#   lead_pos38 取自 31_final_MR_with_method.csv
LOCI = [
    ("YME1L1", "CD4_NC", "10", 27154694, SENS),
    ("ANXA4", "Mono_NC", "2", 69729589, SENS_RECHECK),
]
WINDOWS = {"500Kb": 250000, "1Mb": 500000, "2Mb": 1000000}

p("=== s51b L2/L3 敏感性输入准备 ===")
p("时间 %s" % pd.Timestamp.now())

# ---------- 1) plink bim（chr2/chr10 等位基因映射）----------
bim = {}
p("[1] 读取 plink bim ...")
with open(BIM) as f:
    for line in f:
        c = line.rstrip().split("\t")
        ch, snp, _, pos, a1, a2 = c[0], c[1], c[2], int(c[3]), c[4], c[5]
        if ch in ("2", "10"):
            bim[(ch, pos)] = (snp, a1, a2)
p("    bim 载入 chr2+chr10 共 %d SNP" % len(bim))

# ---------- 2) 结局 GWAS 预载（GRCh38 目标区间）----------
# chr10 lead38=27154694，±1Mb 窗口 -> [26154694, 28154694]，留余量
# chr2  lead38=69729589，±1Mb 窗口 -> [68729589, 70729589]，留余量
GW_PRELOAD = {"10": (26050000, 28250000), "2": (68650000, 70800000)}
_cache = {k: [] for k in GW_PRELOAD}
p("[2] 预载结局 GWAS（GCST90483463.h.tsv.gz）chr2/chr10 目标区间，单次流式扫 ...")
with gzip.open(GW, "rt") as f:
    hdr = f.readline().rstrip("\n").split("\t")
    ix = {c: i for i, c in enumerate(hdr)}
    p("    表头: %s" % hdr)
    for line in f:
        c = line.rstrip("\n").split("\t")
        ch = c[ix["chromosome"]]
        if ch not in GW_PRELOAD:
            continue
        pos = int(c[ix["base_pair_location"]])
        lo, hi = GW_PRELOAD[ch]
        if pos < lo or pos > hi:
            continue
        _cache[ch].append(dict(
            chrom=int(ch), pos38=pos,
            ea=c[ix["effect_allele"]], oa=c[ix["other_allele"]],
            beta=float(c[ix["beta"]]), se=float(c[ix["standard_error"]]),
            eaf=float(c[ix["effect_allele_frequency"]]), p=float(c[ix["p_value"]]),
            rsid=c[ix["rsid"]], n_cases=int(c[ix["n_cases"]]), n_controls=int(c[ix["n_controls"]]),
        ))
_cache = {k: pd.DataFrame(v) for k, v in _cache.items()}
for k, v in _cache.items():
    p("    已预载 chr%s: %d SNP" % (k, len(v)))


def load_gwas(chrom, lo, hi):
    df = _cache[str(chrom)]
    return df[(df["pos38"] >= lo) & (df["pos38"] <= hi)].copy()


# ---------- 3) eQTL ----------
def load_eqtl(cell_type, gene, chrom):
    f = os.path.join(PQ, "OneK1K_%s.cis_qtl_pairs.chr%s.parquet" % (cell_type, chrom))
    df = pd.read_parquet(f, columns=["phenotype_id", "variant_id", "af", "slope",
                                     "slope_se", "pval_nominal"])
    df = df[df["phenotype_id"].astype(str) == gene].copy()
    df["pos37"] = df["variant_id"].astype(str).str.split(":").str[1].astype(int)
    return df


def lift37_38(chrom, pos37):
    r = lo37_38.convert_coordinate("chr%s" % chrom, int(pos37))
    return r[0][1] if r else None


summary = []
p("")
p("[3] 逐位点生成窗口输入")
for gene, ct, chrom, lead38, outdir in LOCI:
    chrom_i = int(chrom)
    r = lo38_37.convert_coordinate("chr%s" % chrom, lead38)
    lead37 = int(r[0][1]) if r else None
    p("")
    p("--- %s / %s / chr%s   lead38=%d -> lead37=%s   输出目录=%s"
      % (gene, ct, chrom, lead38, lead37, outdir))
    eq = load_eqtl(ct, gene, chrom_i)
    p("    eQTL 全量 = %d （pos37 %d-%d）" % (len(eq), eq["pos37"].min(), eq["pos37"].max()))
    for wname, half in WINDOWS.items():
        lo37, hi37 = lead37 - half, lead37 + half
        e = eq[(eq["pos37"] >= lo37) & (eq["pos37"] <= hi37)].copy()
        e["pos38"] = e["pos37"].apply(lambda x: lift37_38(chrom_i, x))
        e = e.dropna(subset=["pos38"]).copy()
        e["pos38"] = e["pos38"].astype(int)
        e["snp_grch37"] = e["pos37"].map(lambda x: "%d:%d" % (chrom_i, x))

        def get_alleles(x, _ch=chrom):
            return bim.get((_ch, int(x)), (None, None, None))
        e[["bim_snp", "a1", "a2"]] = e["pos37"].apply(lambda x: pd.Series(get_alleles(x)))
        n0 = len(e)
        e = e[e["a1"].notna()].copy()
        lo38, hi38 = lead38 - half, lead38 + half
        g = load_gwas(chrom_i, lo38, hi38)
        tag = "%s_%s_%s" % (gene, ct, wname)
        e_out = e[["snp_grch37", "pos37", "pos38", "af", "slope", "slope_se",
                   "pval_nominal", "a1", "a2"]].copy()
        e_out.insert(0, "gene", gene)
        e_out.insert(1, "cell_type", ct)
        e_out.to_csv(os.path.join(outdir, "eqtl_%s.csv" % tag), index=False, encoding="utf-8-sig")
        g.to_csv(os.path.join(outdir, "gwas_%s.csv" % tag), index=False, encoding="utf-8-sig")
        summary.append(dict(gene=gene, cell_type=ct, window=wname, lead37=lead37,
                            n_eqtl_raw=n0, n_eqtl=len(e), n_gwas=len(g)))
        p("    %-24s lead37=%s 窗口[%d,%d] eqtl=%d/%d gwas=%d"
          % (tag, lead37, lo37, hi37, len(e), n0, len(g)))

sdf = pd.DataFrame(summary)
sdf.to_csv(os.path.join(SENS, "_summary_L2L3.csv"), index=False, encoding="utf-8-sig")
sdf.to_csv(os.path.join(SENS_RECHECK, "_summary.csv"), index=False, encoding="utf-8-sig")
p("")
p("=== 汇总 ===")
p(sdf.to_string(index=False))
io.open(LOG, "w", encoding="utf-8", newline="\n").write("\n".join(L) + "\n")
print("WROTE " + LOG)
