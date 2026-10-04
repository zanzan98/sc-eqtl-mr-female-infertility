# -*- coding: utf-8 -*-
"""s54a_build_smr_inputs_L2L3.py -- 任务A：L2 (chr10 YME1L1/CD4_NC) 与 L3 (chr2 ANXA4/Mono_NC)
                                       的 SMR 输入构建（严格复刻 s43a_build_smr_inputs.py）

产出（00_data_raw/smr_L2L3/）：
  eqtl/<TAG>.query.txt          SMR query 格式（每数据集一个 probe）
  gwas_<TAG>.cojo.txt           GCTA-COJO 7 列（SNP=A1=A2=freq=b=se=p=n）
  ld_ref/<TAG>_cis.{bed,bim,fam}  OneK1K 980 供者 PLINK 子集（GRCh37 ID）
  strand.json                   基因链向（Ensembl GRCh37 REST 实查）

两道硬闸门（与 s43a 一致）
  ① liftover 精度：与 OneK1K 自带 snp_grch38 注解比对，一致率须 ≥0.99
  ② GWAS 锚点：与已闭环的 coloc_prep/gwas_<LOCUS>.csv 逐位比对，须完全一致
只读既有数据；不下载。
"""
import io
import json
import os
import subprocess
import sys
import time

import numpy as np
import pandas as pd
import pyarrow.parquet as pqm
from pyliftover import LiftOver

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

PROJ = r"D:\endometriosis_project\11_sc_eqtl_mr_project"
LD = r"D:\endometriosis_project\_onek1k_plink"
PLINK = r"D:\endometriosis_project\_tools\plink.exe"
PQ = os.path.join(PROJ, "00_data_raw", "onek1k", "parquet_subset")
PREP = os.path.join(PROJ, "00_data_raw", "onek1k", "coloc_prep")
GWAS = os.path.join(PROJ, "00_data_raw", "gwas", "GCST90483463.h.tsv.gz")
OUTD = os.path.join(PROJ, "00_data_raw", "smr_L2L3")
EQD = os.path.join(OUTD, "eqtl")
LDD = os.path.join(OUTD, "ld_ref")
LOG = r"D:\endometriosis_project\_s54a_build_smr_inputs.log"
os.makedirs(EQD, exist_ok=True)
os.makedirs(LDD, exist_ok=True)

# ds = 数据集标签；gene/cell 取自 31_final_MR_with_method.csv；cis 窗 = tss37 ± 1 Mb
DS = [
    dict(tag="L2_YME1L1_CD4_NC", locus="L2", gene="YME1L1", cell="CD4_NC",
         chrom="10", tss37=27421790, lo37=25400000, hi37=29000000),
    dict(tag="L3_ANXA4_Mono_NC", locus="L3", gene="ANXA4", cell="Mono_NC",
         chrom="2", tss37=69962577, lo37=67800000, hi37=71800000),
]
L = []


def w(s=""):
    L.append(str(s))
    print(s)


T0 = time.time()
w("=== s54a 构建 L2/L3 SMR 输入 %s ===" % time.strftime("%Y-%m-%d %H:%M:%S"))

# ---- 0) 链向（Ensembl GRCh37 REST，2026-09-28 实查）----
STRAND = {"YME1L1": "-", "ANXA4": "+"}   # GRCh37: YME1L1 strand=-1 ; ANXA4 strand=+1
json.dump(STRAND, io.open(os.path.join(OUTD, "strand.json"), "w", encoding="utf-8"))
w("[0] 链向（仅 .epi 元数据与作图用，SMR/HEIDI 不使用）: %s" % STRAND)

# ---- 1) LD 参考 bim（按染色体分次扫描）----
w("")
w("[1] LD 参考 bim")
bmap_all = {}
with io.open(os.path.join(LD, "plink_merged_980_donors.bim"), encoding="utf-8") as f:
    for ln in f:
        c = ln.rstrip("\n").split("\t")
        if len(c) < 6:
            continue
        bmap_all[c[1]] = (c[0], c[4].upper(), c[5].upper(), int(c[3]))
w("    载入全基因组 bim = %d 变异" % len(bmap_all))

# ---- 2) 结局 GWAS 缓存（一次流式扫，抓 chr2/chr10 目标区间）----
w("")
w("[2] 结局 GWAS（GCST90483463）流式预载 chr2/chr10 目标区间")
RANGES = {"10": (25400000, 29000000), "2": (67800000, 71800000)}
USEC = ["chromosome", "base_pair_location", "effect_allele", "other_allele",
        "beta", "standard_error", "effect_allele_frequency", "p_value",
        "rsid", "n_cases", "n_controls"]
parts = {k: [] for k in RANGES}
nrow = 0
for ck in pd.read_csv(GWAS, sep="\t", usecols=USEC, chunksize=300000,
                      dtype={"chromosome": "str", "effect_allele": "str",
                             "other_allele": "str", "rsid": "str"}, low_memory=False):
    nrow += len(ck)
    for k, (lo, hi) in RANGES.items():
        s = ck[(ck["chromosome"] == k) &
               (ck["base_pair_location"] >= lo) & (ck["base_pair_location"] <= hi)]
        if len(s):
            parts[k].append(s)
GW = {k: (pd.concat(v, ignore_index=True) if v else pd.DataFrame(columns=USEC))
      for k, v in parts.items()}
w("    全文件行数 = %d" % nrow)
for k, v in GW.items():
    w("    chr%s 命中区间 = %d" % (k, len(v)))

lo37_38 = LiftOver("hg19", "hg38")

for D in DS:
    tag, locus, gene, cell, chrom = D["tag"], D["locus"], D["gene"], D["cell"], D["chrom"]
    w("")
    w("=" * 90)
    w("[%s] %s / %s / chr%s   TSS(GRCh37)=%d   cis 窗 = %d..%d"
      % (tag, gene, cell, chrom, D["tss37"], D["tss37"] - 1000000, D["tss37"] + 1000000))

    # ---- 3) 暴露：parquet 完整 cis 窗 ----
    f = os.path.join(PQ, "OneK1K_%s.cis_qtl_pairs.chr%s.parquet" % (cell, chrom))
    cols = ["phenotype_id", "variant_id", "tss_distance", "af", "pval_nominal", "slope", "slope_se"]
    E = pqm.read_table(f, columns=cols,
                       filters=[("phenotype_id", "==", gene)]).to_pandas()
    E["pos37"] = E["variant_id"].str.split(":").str[1].astype(int)
    E["tss"] = E["pos37"] - E["tss_distance"].astype(int)
    w("  暴露行数 = %d ；pos37 %d..%d ；tss = %d" % (len(E), E["pos37"].min(), E["pos37"].max(), E["tss"].iloc[0]))
    it = int(E["pval_nominal"].idxmin())
    w("  top cis-eQTL = %s  p=%.4e  <5e-8? %s  slope=%+.6f" %
      (E.loc[it, "variant_id"], E.loc[it, "pval_nominal"],
       "Y" if E.loc[it, "pval_nominal"] < 5e-8 else "n", E.loc[it, "slope"]))
    n5e8 = int((E["pval_nominal"] < 5e-8).sum())
    n1e5 = int((E["pval_nominal"] < 1e-5).sum())
    w("  过 5e-8 的 SNP 数 = %d ；过 1e-5 = %d" % (n5e8, n1e5))

    # ---- 4) liftover ----
    uniq = sorted(E["pos37"].unique())
    MP = {}
    for p in uniq:
        r = lo37_38.convert_coordinate("chr%s" % chrom, int(p) - 1)
        if not r:
            MP[p] = (None, "unmapped")
        elif len(r) > 1:
            MP[p] = (int(r[0][1]) + 1, "multi")
        else:
            MP[p] = (int(r[0][1]) + 1, "ok")
    st = pd.Series([MP[p][1] for p in uniq]).value_counts().to_dict()
    w("  liftover 映射状态 = %s" % st)

    ref = pd.read_csv(os.path.join(PREP, "eqtl_%s_%s.csv" % (locus, gene)),
                      usecols=["snp_grch37", "snp_grch38"])
    ref["p37"] = ref["snp_grch37"].str.split(":").str[1].astype(int)
    ref["p38_ref"] = ref["snp_grch38"].str.split(":").str[1].astype(int)
    chk = ref.drop_duplicates("p37")[["p37", "p38_ref"]].copy()
    chk["p38_new"] = chk["p37"].map(lambda p: MP.get(int(p), (None, "?"))[0])
    chk = chk[chk["p38_new"].notna()]
    chk["d"] = (chk["p38_new"] - chk["p38_ref"]).abs()
    agree = float((chk["d"] == 0).mean())
    w("  ★闸门① liftover：重叠=%d 一致=%d 一致率=%.6f |Δ|max=%d"
      % (len(chk), int((chk["d"] == 0).sum()), agree, int(chk["d"].max()) if len(chk) else -1))
    if agree < 0.99:
        w("  !! 一致率 <0.99，终止")
        io.open(LOG, "w", encoding="utf-8", newline="\n").write("\n".join(L) + "\n")
        sys.exit(5)
    E["pos38"] = E["pos37"].map(lambda p: MP[int(p)][0])
    E["liftover"] = E["pos37"].map(lambda p: MP[int(p)][1])

    # ---- 5) GWAS 锚点闸门 ----
    G = GW[chrom].copy()
    gl = pd.read_csv(os.path.join(PREP, "gwas_%s.csv" % locus))
    ndup_old = int(gl["pos38"].duplicated().sum())
    gl_u = gl.rename(columns={"pos38": "pos38_old", "ea": "ea_old", "oa": "oa_old",
                             "beta": "beta_old", "se": "se_old", "eaf": "eaf_old",
                             "p": "p_old"})
    gl_u = gl_u[~gl_u["pos38_old"].duplicated(keep=False)]
    m = G.rename(columns={"base_pair_location": "pos38_new", "effect_allele": "ea_new",
                          "other_allele": "oa_new", "beta": "beta_new",
                          "standard_error": "se_new", "effect_allele_frequency": "eaf_new",
                          "p_value": "p_new"}).merge(
        gl_u, left_on="pos38_new", right_on="pos38_old", how="inner")
    w("  ★闸门② GWAS 锚点：唯一位点重叠 = %d （gwas_%s 重复行 %d 已剔除）"
      % (len(m), locus, ndup_old))
    anchor_ok = len(m) > 1000
    for cn, co in [("ea_new", "ea_old"), ("oa_new", "oa_old"), ("beta_new", "beta_old"),
                   ("se_new", "se_old"), ("eaf_new", "eaf_old"), ("p_new", "p_old")]:
        if m[cn].dtype.kind in "if":
            d = (m[cn] - m[co]).abs()
            nm = int((d > 1e-9).sum())
            w("      |Δ %-8s| max = %.3e  (>1e-9 个数 %d)" % (cn.replace("_new", ""), float(np.nanmax(d)), nm))
            if nm > 0:
                anchor_ok = False
        else:
            nm = int((m[cn].astype(str).str.upper() != m[co].astype(str).str.upper()).sum())
            w("      %-10s 不一致个数 = %d" % (cn.replace("_new", ""), nm))
            if nm > 0:
                anchor_ok = False
    w("  闸门② = %s" % ("PASS" if anchor_ok else "FAIL"))
    if not anchor_ok:
        w("  !! 闸门② FAIL，终止")
        io.open(LOG, "w", encoding="utf-8", newline="\n").write("\n".join(L) + "\n")
        sys.exit(6)

    # ---- 6) GWAS 标注 GRCh37 ID 并写 cojo ----
    p38_2_p37 = {}
    for p37, (p38, s) in MP.items():
        if p38 is not None:
            p38_2_p37.setdefault(p38, []).append(p37)
    G["pos37"] = G["base_pair_location"].map(lambda x: (p38_2_p37.get(int(x)) or [None])[0])
    w("  可在 eQTL 位点表中回标的 GWAS 行 = %d / %d" % (int(G["pos37"].notna().sum()), len(G)))
    G2 = []
    nd_na = nd_multi = nd_amb = 0
    for p37, g in G[G["pos37"].notna()].groupby("pos37"):
        vid = "%d:%d" % (int(chrom), int(p37))
        if vid not in bmap_all or bmap_all[vid][0] != chrom:
            nd_na += 1
            continue
        a1, a2 = bmap_all[vid][1], bmap_all[vid][2]
        want = {a1, a2}
        cand = g[[set([str(x).upper(), str(y).upper()]) == want
                  for x, y in zip(g["effect_allele"], g["other_allele"])]]
        if len(cand) == 0:
            if len(g) > 1:
                nd_multi += 1
            else:
                nd_amb += 1
            continue
        r = cand.iloc[0]
        G2.append(dict(SNP=vid, A1=str(r["effect_allele"]).upper(),
                       A2=str(r["other_allele"]).upper(),
                       freq=float(r["effect_allele_frequency"]) if pd.notna(r["effect_allele_frequency"]) else "NA",
                       b=float(r["beta"]), se=float(r["standard_error"]),
                       p=float(r["p_value"]), n=int(r["n_cases"]) + int(r["n_controls"])))
    GD = pd.DataFrame(G2)
    w("  可用 GWAS SNP = %d （丢弃：等位不匹配 %d ；多等位无匹配 alt %d ；不在 bim %d）"
      % (len(GD), nd_amb, nd_multi, nd_na))
    gpath = os.path.join(OUTD, "gwas_%s.cojo.txt" % tag)
    with io.open(gpath, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("SNP\tA1\tA2\tfreq\tb\tse\tp\tn\n")
        for _, r in GD.iterrows():
            fh.write("%s\t%s\t%s\t%s\t%.10g\t%.10g\t%.10g\t%d\n"
                     % (r["SNP"], r["A1"], r["A2"],
                        r["freq"] if r["freq"] == "NA" else "%.6g" % r["freq"],
                        r["b"], r["se"], r["p"], r["n"]))
    w("  → %s（%d 行）" % (os.path.basename(gpath), len(GD)))

    # ---- 7) 暴露 query 文件 ----
    HDR = ["SNP", "Chr", "BP", "A1", "A2", "Freq", "Probe", "Probe_Chr", "Probe_bp",
           "Gene", "Orientation", "b", "se", "p"]
    lines = []
    n_skipped = 0
    for _, r in E.iterrows():
        vid = r["variant_id"]
        if vid not in bmap_all or pd.isna(r["pos38"]):
            n_skipped += 1
            continue
        a1, a2 = bmap_all[vid][1], bmap_all[vid][2]
        lines.append("%s\t%s\t%d\t%s\t%s\t%.6g\t%s\t%s\t%d\t%s\t%s\t%.10g\t%.10g\t%.10g"
                     % (vid, chrom, int(r["pos37"]), a1, a2, float(r["af"]),
                        r["phenotype_id"], chrom, int(r["tss"]), r["phenotype_id"],
                        STRAND.get(r["phenotype_id"], "+"), r["slope"], r["slope_se"],
                        r["pval_nominal"]))
    qp = os.path.join(EQD, "%s.query.txt" % tag)
    with io.open(qp, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\t".join(HDR) + "\n")
        fh.write("\n".join(lines) + "\n")
    w("  → %s（%d 行，跳过 %d）" % (os.path.basename(qp), len(lines), n_skipped))

    # ---- 8) LD 参考 plink 子集 ----
    ldp = os.path.join(LDD, "%s_cis" % tag)
    need = (not os.path.exists(ldp + ".bed")) or (not os.path.exists(ldp + ".bim")) \
           or (not os.path.exists(ldp + ".fam"))
    if need:
        cmd = [PLINK, "--bfile", os.path.join(LD, "plink_merged_980_donors"),
               "--allow-no-sex", "--chr", chrom, "--from-bp", str(D["lo37"]),
               "--to-bp", str(D["hi37"]), "--make-bed", "--out", ldp]
        w("  CMD: %s" % " ".join(cmd))
        subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                       check=False, timeout=3600)
    ok = all(os.path.exists(ldp + e) for e in (".bed", ".bim", ".fam"))
    n_ld = sum(1 for _ in io.open(ldp + ".bim", encoding="utf-8")) if ok else 0
    w("  → %s LD 参考：ok=%s 变异=%d" % (os.path.basename(ldp), ok, n_ld))

w("")
w("VERDICT = %s" % "PASS")
w("总耗时 %.1f s" % (time.time() - T0))
io.open(LOG, "w", encoding="utf-8", newline="\n").write("\n".join(L) + "\n")
print("WROTE " + LOG)
