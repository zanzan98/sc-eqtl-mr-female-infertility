# -*- coding: utf-8 -*-
"""s43a_build_smr_inputs.py -- 任务 3.3 SMR 输入构建

暴露（OneK1K sc-eQTL，逐细胞类型，**直接取原始 parquet 的完整 ±1 Mb cis 窗**）
  → 00_data_raw/smr_3p3/eqtl/<cell>.query.txt（SMR query 输出格式，含 2 个 probe）
结局（GCST90483463，女性不孕症主结局）
  → 00_data_raw/smr_3p3/gwas_GCST90483463.cojo.txt（GCTA-COJO 格式）

★ 为什么不用既有的 `coloc_prep/eqtl_L1_*.csv`：该表是 GRCh37 22.0-22.9 Mb 的**截断**，
  而原始 cis 窗是 TSS ±1 Mb（21.4-23.4 Mb）。为不损失 HEIDI 可用位点，改从 parquet 重建。

坐标系：LD 参考 bim 的 SNP ID 为 `1:<pos37>`（GRCh37）→ 全链统一用 GRCh37 ID 空间；
  GWAS 原生 GRCh38 → 用 pyliftover（本项目既有实现）把 eQTL 的 GRCh37 位点抬到 GRCh38，
  以此把 GWAS 行标注回 `1:<pos37>`。

两道硬闸门
  ① liftover 精度：与 OneK1K 自带 `snp_grch38` 注解在 2,068 个重叠位点上比对，一致率须 ≥0.99；
  ② GWAS 回归锚点：与已闭环的 `coloc_prep/gwas_L1.csv` 在重叠位点上逐位比对，须完全一致。
"""
import csv
import io
import os
import sys
import time

import numpy as np
import pandas as pd
import pyarrow.parquet as pqm
from pyliftover import LiftOver

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

PROJ = r"D:\endometriosis_project\11_sc_eqtl_mr_project"
LD = r"D:\endometriosis_project\_onek1k_plink"
PQ = os.path.join(PROJ, "00_data_raw", "onek1k", "parquet_subset")
PREP = os.path.join(PROJ, "00_data_raw", "onek1k", "coloc_prep")
GWAS = os.path.join(PROJ, "00_data_raw", "gwas", "GCST90483463.h.tsv.gz")
OUTD = os.path.join(PROJ, "00_data_raw", "smr_3p3")
EQD = os.path.join(OUTD, "eqtl")
LOG = r"D:\endometriosis_project\_s43a_build_smr_inputs.log"

GENES = ["CDC42", "LINC00339", "WNT4"]
# ★ 由 Ensembl REST 实查（GRCh38，见 00_data_raw/smr_3p3/strand.json）：
#   CDC42 (ENSG00000070831) strand=+ ; LINC00339 (ENSG00000218510) strand=+ ; WNT4 (ENSG00000162552) strand=-
#   注：Orientation 仅用于 .epi 元数据与 locus plot，SMR/HEIDI 检验不使用。
STRAND = {"CDC42": "+", "LINC00339": "+", "WNT4": "-"}
L = []


def w(s=""):
    L.append(str(s))
    print(s)


T0 = time.time()
w("=== s43a 构建 SMR 输入 %s ===" % time.strftime("%Y-%m-%d %H:%M:%S"))
os.makedirs(EQD, exist_ok=True)

# ================= 1) LD 参考 bim（只取 chr1） =================
w("")
w("[1] LD 参考 bim（chr1）")
bim = os.path.join(LD, "plink_merged_980_donors.bim")
bmap = {}
n_scan = 0
with io.open(bim, encoding="utf-8") as f:
    for ln in f:
        n_scan += 1
        p = ln.rstrip("\n").split("\t")
        if len(p) < 6:
            continue
        if p[0] != "1":
            if bmap:
                break                     # bim 按染色体排序，chr1 读完即止
            continue
        bmap[p[1]] = (p[4].upper(), p[5].upper(), int(p[3]))
w("   扫描行数=%d ；chr1 变异=%d" % (n_scan, len(bmap)))
w("   样例 = %s" % list(bmap.items())[:2])

# ================= 2) 暴露：从 parquet 取完整 cis 窗 =================
w("")
w("[2] 暴露：OneK1K 原始 parquet（phenotype ∈ %s）" % GENES)
cols = ["phenotype_id", "variant_id", "tss_distance", "af", "pval_nominal", "slope", "slope_se"]
frames = []
wnt4_n = 0
pqs = sorted([os.path.join(PQ, f) for f in os.listdir(PQ) if f.endswith(".chr1.parquet")])
for fp in pqs:
    cell = os.path.basename(fp).split(".")[0].replace("OneK1K_", "")
    t = pqm.read_table(fp, columns=cols, filters=[("phenotype_id", "in", GENES)])
    d = t.to_pandas()
    wnt4_n += int((d["phenotype_id"] == "WNT4").sum())
    d = d[d["phenotype_id"].isin(["CDC42", "LINC00339"])].copy()
    d["cell"] = cell
    frames.append(d)
E = pd.concat(frames, ignore_index=True)
w("   行数=%d ；细胞=%d ；基因=%s" % (len(E), E["cell"].nunique(), sorted(E["phenotype_id"].unique())))
w("   ★ WNT4 在 14 个细胞类型中的 cis 对数 = %d（0 = 该基因在 OneK1K 血液免疫细胞无 cis-eQTL）" % wnt4_n)
E["pos37"] = E["variant_id"].str.split(":").str[1].astype(int)
E["tss"] = E["pos37"] - E["tss_distance"].astype(int)
top = E.loc[E.groupby(["phenotype_id", "cell"])["pval_nominal"].idxmin()][
    ["phenotype_id", "cell", "variant_id", "pval_nominal", "slope", "slope_se", "af", "tss"]]
w("   每 (基因,细胞) 的 top cis-eQTL：")
for _, r in top.sort_values(["phenotype_id", "pval_nominal"]).iterrows():
    w("     %-10s %-10s %-13s p=%.3e  <5e-8? %s  tss=%d"
      % (r["phenotype_id"], r["cell"], r["variant_id"], r["pval_nominal"],
         "Y" if r["pval_nominal"] < 5e-8 else "n", r["tss"]))
n_pass = int((top["pval_nominal"] < 5e-8).sum())
w("   ★ top cis-eQTL p<5e-8 的 (基因,细胞) 组合数 = %d / %d" % (n_pass, len(top)))

# ================= 3) liftover GRCh37 -> GRCh38 =================
w("")
w("[3] liftover GRCh37 -> GRCh38（pyliftover hg19ToHg38，本项目既有实现）")
lo = LiftOver("hg19", "hg38")
uniq = sorted(E["pos37"].unique())
w("   唯一 GRCh37 位点 = %d ；范围 %d..%d" % (len(uniq), uniq[0], uniq[-1]))
MP = {}
for p in uniq:
    r = lo.convert_coordinate("chr1", int(p) - 1)
    if not r:
        MP[p] = (None, "unmapped")
    elif len(r) > 1:
        MP[p] = (int(r[0][1]) + 1, "multi")
    else:
        MP[p] = (int(r[0][1]) + 1, "ok")
st = pd.Series([MP[p][1] for p in uniq]).value_counts().to_dict()
w("   映射状态 = %s" % st)

# ★ 闸门 ①：与 OneK1K 自带 snp_grch38 注解比对
ref = pd.read_csv(os.path.join(PREP, "eqtl_L1_CDC42.csv"), usecols=["snp_grch37", "snp_grch38"])
ref["p37"] = ref["snp_grch37"].str.split(":").str[1].astype(int)
ref["p38_ref"] = ref["snp_grch38"].str.split(":").str[1].astype(int)
chk = ref.drop_duplicates("p37")[["p37", "p38_ref"]].copy()
chk["p38_new"] = chk["p37"].map(lambda p: MP.get(p, (None, "?"))[0])
chk = chk[chk["p38_new"].notna()]
chk["d"] = (chk["p38_new"] - chk["p38_ref"]).abs()
agree = float((chk["d"] == 0).mean())
w("   ★ 闸门① liftover 精度：重叠位点=%d ；与 OneK1K 注解一致=%d ；一致率=%.6f ；|Δ| max=%d"
  % (len(chk), int((chk["d"] == 0).sum()), agree, int(chk["d"].max()) if len(chk) else -1))
if agree < 0.99:
    w("   !! 一致率 <0.99，终止")
    io.open(LOG, "w", encoding="utf-8", newline="\n").write("\n".join(L) + "\n")
    sys.exit(5)

E["pos38"] = E["pos37"].map(lambda p: MP[p][0])
E["liftover"] = E["pos37"].map(lambda p: MP[p][1])
w("   暴露表：pos38 缺失 = %d" % int(E["pos38"].isna().sum()))

# ================= 4) 结局 GWAS =================
w("")
w("[4] 结局 GWAS：GCST90483463（女性不孕症主结局）")
lo38 = int(E["pos38"].min()) - 100000
hi38 = int(E["pos38"].max()) + 100000
w("   提取范围（GRCh38）= %d .. %d（由暴露 pos38 范围 ±100 kb）" % (lo38, hi38))
USEC = ["chromosome", "base_pair_location", "effect_allele", "other_allele",
        "beta", "standard_error", "effect_allele_frequency", "p_value",
        "rsid", "n_cases", "n_controls"]
parts, nrow = [], 0
for ck in pd.read_csv(GWAS, sep="\t", usecols=USEC, chunksize=300000,
                      dtype={"chromosome": "str", "effect_allele": "str", "other_allele": "str",
                             "rsid": "str"}, low_memory=False):
    nrow += len(ck)
    ck = ck[ck["chromosome"] == "1"]
    ck = ck[(ck["base_pair_location"] >= lo38) & (ck["base_pair_location"] <= hi38)]
    if len(ck):
        parts.append(ck)
G = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame(columns=USEC)
w("   全文件行数=%d ；命中范围行数=%d" % (nrow, len(G)))

# ★ 闸门 ②：与已闭环的 gwas_L1.csv 逐位比对
#   注意 gwas_L1.csv 自身含 4 个**多等位**位点（8 行，同一 pos38 两个 alt）→ 多对多连接会造假差异，
#   故锚点只在**唯一位点**上比对，多等位位点单独登记。
gl = pd.read_csv(os.path.join(PREP, "gwas_L1.csv"))
ndup_old = int(gl["pos38"].duplicated().sum())
gl_u = gl.rename(columns={"pos38": "pos38_old", "ea": "ea_old", "oa": "oa_old",
                          "beta": "beta_old", "se": "se_old", "eaf": "eaf_old",
                          "p": "p_old", "n_cases": "nc_old", "n_controls": "nct_old"})
gl_u = gl_u[~gl_u["pos38_old"].duplicated(keep=False)]
m = G.rename(columns={"base_pair_location": "pos38_new", "effect_allele": "ea_new",
                      "other_allele": "oa_new", "beta": "beta_new", "standard_error": "se_new",
                      "effect_allele_frequency": "eaf_new", "p_value": "p_new"}).merge(
    gl_u, left_on="pos38_new", right_on="pos38_old", how="inner")
w("   ★ 闸门② GWAS 锚点：唯一位点重叠 = %d（gwas_L1 重复行 %d 已剔除）" % (len(m), ndup_old))
anchor_ok = len(m) > 12000
for cn, co in [("ea_new", "ea_old"), ("oa_new", "oa_old"), ("beta_new", "beta_old"),
               ("se_new", "se_old"), ("eaf_new", "eaf_old"), ("p_new", "p_old")]:
    if m[cn].dtype.kind in "if":
        d = (m[cn] - m[co]).abs()
        nm = int((d > 1e-9).sum())
        w("      |Δ %-10s| max = %.3e  (>1e-9 的个数 %d)" % (cn.replace("_new", ""), float(np.nanmax(d)), nm))
        if nm > 0:
            anchor_ok = False
    else:
        nm = int((m[cn].astype(str).str.upper() != m[co].astype(str).str.upper()).sum())
        w("      %-13s 不一致个数 = %d" % (cn.replace("_new", ""), nm))
        if nm > 0:
            anchor_ok = False
w("   闸门② 结果 = %s" % ("PASS" if anchor_ok else "FAIL"))

# ================= 5) GWAS 按 eQTL 位点标注 GRCh37 ID =================
w("")
w("[5] 把 GWAS 行标注回 `1:<pos37>`（SMR 三个输入必须共用同一 SNP ID）")
p38_2_p37 = {}
for p37, (p38, s) in MP.items():
    if p38 is not None:
        p38_2_p37.setdefault(p38, []).append(p37)
G["pos37"] = G["base_pair_location"].map(lambda x: (p38_2_p37.get(int(x)) or [None])[0])
G["n_ambi"] = G["base_pair_location"].map(lambda x: len(p38_2_p37.get(int(x), [])))
w("   位置可在 eQTL 位点表中找到的 GWAS 行 = %d / %d" % (int(G["pos37"].notna().sum()), len(G)))
w("   一对多的位置（GRCh37 合并）= %d 行" % int((G["n_ambi"] > 1).sum()))

Ekey = E[["variant_id", "phenotype_id", "cell"]].drop_duplicates()
al = E.drop_duplicates("variant_id").set_index("variant_id")
G2 = []
ndrop_amb = ndrop_na = ndrop_multi = 0
for p37, g in G[G["pos37"].notna()].groupby("pos37"):
    vid = "1:%d" % int(p37)
    if vid not in bmap:
        ndrop_na += 1
        continue
    a1, a2 = bmap[vid][0], bmap[vid][1]
    want = {a1, a2}
    # ★ 只接受**等位集合完全一致**的行（多等位位点必须选到同一个变异，否则宁可丢弃）
    cand = g[[set([str(x).upper(), str(y).upper()]) == want
              for x, y in zip(g["effect_allele"], g["other_allele"])]]
    if len(cand) == 0:
        if len(g) > 1:
            ndrop_multi += 1
        else:
            ndrop_amb += 1
        continue
    r = cand.iloc[0]
    G2.append(dict(SNP=vid, A1=str(r["effect_allele"]).upper(), A2=str(r["other_allele"]).upper(),
                   freq=float(r["effect_allele_frequency"]) if pd.notna(r["effect_allele_frequency"]) else "NA",
                   b=float(r["beta"]), se=float(r["standard_error"]), p=float(r["p_value"]),
                   n=int(r["n_cases"]) + int(r["n_controls"])))
GD = pd.DataFrame(G2)
w("   可用 GWAS SNP = %d" % len(GD))
w("     丢弃：等位集合不匹配（单位点）%d ；多等位位点无匹配 alt %d ；位点不在参考 bim %d"
  % (ndrop_amb, ndrop_multi, ndrop_na))
gpath = os.path.join(OUTD, "gwas_GCST90483463.cojo.txt")
with io.open(gpath, "w", encoding="utf-8", newline="\n") as f:
    f.write("SNP\tA1\tA2\tfreq\tb\tse\tp\tn\n")
    for _, r in GD.iterrows():
        f.write("%s\t%s\t%s\t%s\t%.10g\t%.10g\t%.10g\t%d\n"
                % (r["SNP"], r["A1"], r["A2"], r["freq"] if r["freq"] == "NA" else "%.6g" % r["freq"],
                   r["b"], r["se"], r["p"], r["n"]))
w("   → %s（%d 行 + 表头）" % (gpath, len(GD)))

# ================= 6) 暴露：SMR query 格式，逐细胞类型 =================
w("")
w("[6] 暴露文件（SMR query 格式，每细胞类型一个，含 2 个 probe）")
HDR = ["SNP", "Chr", "BP", "A1", "A2", "Freq", "Probe", "Probe_Chr", "Probe_bp", "Gene", "Orientation", "b", "se", "p"]
rows_all = []
out_files = []
for cell, g in E.groupby("cell"):
    lines = []
    for _, r in g.iterrows():
        vid = r["variant_id"]
        if vid not in bmap or pd.isna(r["pos38"]):
            continue
        a1, a2 = bmap[vid][0], bmap[vid][1]
        lines.append("%s\t1\t%d\t%s\t%s\t%.6g\t%s\t1\t%d\t%s\t%s\t%.10g\t%.10g\t%.10g"
                     % (vid, int(r["pos37"]), a1, a2, float(r["af"]), r["phenotype_id"], int(r["tss"]),
                        r["phenotype_id"], STRAND.get(r["phenotype_id"], "+"),
                        r["slope"], r["slope_se"], r["pval_nominal"]))
        rows_all.append(dict(cell=cell, gene=r["phenotype_id"], snp=vid, p=r["pval_nominal"]))
    fp = os.path.join(EQD, "%s.query.txt" % cell)
    with io.open(fp, "w", encoding="utf-8", newline="\n") as f:
        f.write("\t".join(HDR) + "\n")
        f.write("\n".join(lines) + "\n")
    out_files.append(fp)
    w("   %-10s 行=%5d  cells=%s  probes=%s" % (cell, len(lines), cell,
                                                sorted(g["phenotype_id"].unique())))
RA = pd.DataFrame(rows_all)
w("   暴露总行数 = %d ；文件数 = %d" % (len(RA), len(out_files)))

# ================= 7) 落盘小结 =================
w("")
w("[7] 小结")
w("   liftover 一致率 = %.6f（闸门①）" % agree)
w("   GWAS 锚点 = %s（闸门②）" % ("PASS" if anchor_ok else "FAIL"))
w("   WNT4 cis 对数 = %d → %s" % (wnt4_n, "WNT4 无法运行血液 SMR（缺暴露工具）" if wnt4_n == 0 else "可运行"))
w("   总耗时 %.1f s" % (time.time() - T0))
verdict = "PASS" if (agree >= 0.99 and anchor_ok and wnt4_n == 0) else "WARN"
w("VERDICT = %s" % verdict)
io.open(LOG, "w", encoding="utf-8", newline="\n").write("\n".join(L) + "\n")
print("WROTE " + LOG)
