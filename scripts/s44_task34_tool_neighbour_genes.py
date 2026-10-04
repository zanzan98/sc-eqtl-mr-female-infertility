# -*- coding: utf-8 -*-
"""s44_task34_tool_neighbour_genes.py — 任务 3.4：工具邻近基因 cis 效应检查

问题：rs12037376 与 rs10917151 各自是否只调控 CDC42，还是同时调控 CDC42 / LINC00339 / WNT4？
      若是后者 ⇒ 它们是「区域标签（region tag）」而非「基因标签（gene tag）」。

数据源（两层，互相独立）：
  ① OneK1K 单细胞 cis-eQTL（14 细胞类型，GRCh37，纯免疫细胞）
     `00_data_raw/onek1k/parquet_subset/OneK1K_<cell>.cis_qtl_pairs.chr1.parquet`
     列：phenotype_id / variant_id / tss_distance / af / ma_samples / ma_count / pval_nominal / slope / slope_se
  ② GTEx v8 全血 bulk cis-eQTL（eQTL Catalogue imported，n=670，GRCh38 源）
     `01_data_raw/eqtlcat_gtex_v8/Whole_Blood.chr1_21_23.3Mb_3genes.tsv`

坐标系（★两套必须换算后并列，不得混用）：
  rs12037376 : GRCh37 hg19 `1:22462111`  ↔  GRCh38 `chr1:22135618`  (δ = pos38 − pos37 = −326,493)
  rs10917151 : GRCh37 hg19 `1:22422721`  ↔  GRCh38 `chr1:22096228`
  OneK1K 原生 GRCh37；GTEx v8 原生 GRCh38；LD 参考 ID 空间 = `1:<pos37>`

效应等位判定：
  · GTEx v8 / eQTL Catalogue：`alt` 恒为效应等位（官方 FAQ）
  · OneK1K：`slope` 与 `af` 同为「剂量等位」；用 LD 参考实测频率反查该等位是哪个碱基
    （避免凭空假设 alt，本项目既有纪律：第三方文件列语义必须实证）

输出：
  tables/49_task34_tool_neighbour_cis_eqtl.csv     长表（逐 SNP×基因×细胞/组织）
  tables/49b_task34_region_tag_summary.csv         宽表（逐 SNP×基因 汇总 + 区域标签判定）
  tables/49c_task34_gwas_context.csv               工具在结局 GWAS 中的效应（区域标签论证的补充）
"""
import os, io, sys, csv
import numpy as np
import pandas as pd
import pyarrow.parquet as pq

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

PROJ = r"D:\endometriosis_project\11_sc_eqtl_mr_project"
PQ = os.path.join(PROJ, "00_data_raw", "onek1k", "parquet_subset")
SMR3 = os.path.join(PROJ, "00_data_raw", "smr_3p3")
LD = os.path.join(SMR3, "ld_ref", "chr1_20_25mb")
COJO = os.path.join(SMR3, "gwas_GCST90483463.cojo.txt")
G8SRC = os.path.join(PROJ, "..", "01_data_raw", "eqtlcat_gtex_v8",
                     "Whole_Blood.chr1_21_23.3Mb_3genes.tsv")
G8SRC = os.path.abspath(G8SRC)
TAB = os.path.join(PROJ, "tables")
LOG = os.path.join(PROJ, "scripts", "_s44_task34.log")
os.makedirs(TAB, exist_ok=True)

buf = []
def w(s=""):
    buf.append(str(s)); print(s)
def flush():
    io.open(LOG, "w", encoding="utf-8").write("\n".join(buf) + "\n")

DELTA = -326493                     # pos38 - pos37
GENES = ["CDC42", "LINC00339", "WNT4"]
TOOLS = {
    "rs12037376": {"pos37": 22462111, "pos38": 22135618},
    "rs10917151": {"pos37": 22422721, "pos38": 22096228},
}
V37 = {v["pos37"]: k for k, v in TOOLS.items()}

w("=" * 104)
w("s44 — 任务 3.4：工具邻近基因 cis 效应检查（rs12037376 / rs10917151）")
w("=" * 104)

# ---------- 0. LD 参考（用于判定 OneK1K 效应等位）----------
bim = pd.read_csv(LD + ".bim", sep=r"\s+", header=None,
                  names=["CHR", "SNP", "CM", "BP", "A1", "A2"]).set_index("SNP")
frq = pd.read_csv(os.path.join(SMR3, "ld_ref", "_ldfreq.frq"), sep=r"\s+").set_index("SNP")
w("[0] LD 参考 OneK1K 980 人频率（用于反查 OneK1K 剂量等位碱基）")
ldfreq = {}
for rs, t in TOOLS.items():
    snp = "1:%d" % t["pos37"]
    a1, a2 = bim.at[snp, "A1"], bim.at[snp, "A2"]
    maf = float(frq.at[snp, "MAF"])
    ldfreq[snp] = {a1: maf, a2: 1.0 - maf}
    w("    %-12s %-12s bim A1=%s(%.4f) A2=%s(%.4f)" % (rs, snp, a1, maf, a2, 1.0 - maf))

# ---------- 1. OneK1K 14 细胞类型 ----------
w()
w("[1] OneK1K：14 细胞类型 parquet 过滤提取（variant_id ∈ 两个工具 × phenotype_id ∈ 3 基因）")
rows = []
cells = sorted(f.split(".")[0].replace("OneK1K_", "")
               for f in os.listdir(PQ) if f.endswith(".chr1.parquet"))
targets = ["1:%d" % t["pos37"] for t in TOOLS.values()]
for cell in cells:
    fp = os.path.join(PQ, "OneK1K_%s.cis_qtl_pairs.chr1.parquet" % cell)
    t = pq.read_table(fp,
                      columns=["phenotype_id", "variant_id", "tss_distance", "af",
                               "pval_nominal", "slope", "slope_se"],
                      filters=[("variant_id", "in", targets),
                               ("phenotype_id", "in", GENES)])
    d = t.to_pandas()
    n = 0
    for r in d.itertuples():
        rs = V37.get(int(r.variant_id.split(":")[1]), r.variant_id)
        # 反查剂量等位
        fmap = ldfreq.get(r.variant_id, {})
        a1, a2 = ("?", "?")
        if fmap:
            cands = sorted(fmap.items(), key=lambda kv: abs(kv[1] - r.af))
            a1 = cands[0][0]          # 与 af 最接近的碱基 = 剂量等位
            a2 = cands[1][0]
        rows.append({"SNP_rsID": rs, "pos37_hg19": int(r.variant_id.split(":")[1]),
                     "pos38_GRCh38": int(r.variant_id.split(":")[1]) + DELTA,
                     "dataset": "OneK1K_scRNA", "cell_or_tissue": cell,
                     "gene": r.phenotype_id, "effect_allele": a1, "other_allele": a2,
                     "af": float(r.af), "b_eQTL": float(r.slope), "se_eQTL": float(r.slope_se),
                     "p_eQTL": float(r.pval_nominal), "tss_distance": int(r.tss_distance),
                     "n_obs": np.nan})
        n += 1
    w("    %-10s 命中 %d 行" % (cell, n))
one = pd.DataFrame(rows)
w("    OneK1K 合计 %d 行；基因分布 = %s" % (len(one), one.gene.value_counts().to_dict()))

# ---------- 2. GTEx v8 全血 ----------
w()
w("[2] GTEx v8 全血 bulk：源 TSV 取两个 b38 位置")
w("    源 = %s" % G8SRC)
trow = []
with io.open(G8SRC, "r", encoding="utf-8", errors="replace") as f:
    ix = {c: i for i, c in enumerate(f.readline().rstrip("\n").split("\t"))}
    for ln in f:
        q = ln.rstrip("\n").split("\t")
        if len(q) <= ix["gene_symbol"]:
            continue
        try:
            pos38 = int(q[ix["position"]])
        except ValueError:
            continue
        if pos38 not in [t["pos38"] for t in TOOLS.values()]:
            continue
        rs = [k for k, v in TOOLS.items() if v["pos38"] == pos38][0]
        trow.append({"SNP_rsID": rs, "pos37_hg19": pos38 - DELTA, "pos38_GRCh38": pos38,
                     "dataset": "GTEx_v8_WholeBlood_bulk", "cell_or_tissue": "WholeBlood",
                     "gene": q[ix["gene_symbol"]],
                     "effect_allele": q[ix["alt"]], "other_allele": q[ix["ref"]],
                     "af": float(q[ix["ac"]]) / float(q[ix["an"]]),
                     "b_eQTL": float(q[ix["beta"]]), "se_eQTL": float(q[ix["se"]]),
                     "p_eQTL": float(q[ix["pvalue"]]),
                     "tss_distance": np.nan,
                     "n_obs": int(round(float(q[ix["an"]]) / 2.0))})
g8 = pd.DataFrame(trow)
w("    GTEx v8 命中 %d 行" % len(g8))

long = pd.concat([one, g8], ignore_index=True)
long = long.sort_values(["SNP_rsID", "gene", "dataset", "cell_or_tissue"])
long.to_csv(os.path.join(TAB, "49_task34_tool_neighbour_cis_eqtl.csv"),
            index=False, encoding="utf-8-sig")
w("    已写 tables/49_task34_tool_neighbour_cis_eqtl.csv（%d 行）" % len(long))

# ---------- 3. 宽表汇总 + 区域标签判定 ----------
w()
w("=" * 104)
w("[3] 逐 SNP × 基因 汇总")
w("=" * 104)
summ = []
for rs in TOOLS:
    for g in GENES:
        a = one[(one.SNP_rsID == rs) & (one.gene == g)]
        b = g8[(g8.SNP_rsID == rs) & (g8.gene == g)]
        n_ok = int((a.p_eQTL < 0.05).sum()) if len(a) else 0
        sig = a[a.p_eQTL < 0.05] if len(a) else a
        if len(sig) and (sig.b_eQTL > 0).all():
            sc = "all_pos"
        elif len(sig) and (sig.b_eQTL < 0).all():
            sc = "all_neg"
        elif len(sig):
            sc = "mixed"
        else:
            sc = "NA"
        summ.append({
            "SNP_rsID": rs,
            "pos37_hg19": TOOLS[rs]["pos37"], "pos38_GRCh38": TOOLS[rs]["pos38"],
            "gene": g,
            "OneK1K_n_celltypes": int(len(a)),
            "OneK1K_n_sig": n_ok,
            "OneK1K_n_available": int(len(a)),
            "OneK1K_min_p": float(a.p_eQTL.min()) if len(a) else np.nan,
            "OneK1K_sig_celltypes": ",".join(
                sig.sort_values("p_eQTL").cell_or_tissue.tolist()) if len(sig) else "",
            # ★ 符号一致性只在 p<0.05 子集内判定（否则会把 CD4_SOX4 等零效应算进去）
            "OneK1K_sign_consistency_sig_only": sc,
            "OneK1K_b_range_all": ("%.4f..%.4f" % (a.b_eQTL.min(), a.b_eQTL.max())) if len(a) else "NA",
            "OneK1K_b_range_sig": ("%.4f..%.4f" % (sig.b_eQTL.min(), sig.b_eQTL.max())) if len(sig) else "NA",
            "GTEx_v8_b_eQTL": float(b.b_eQTL.iloc[0]) if len(b) else np.nan,
            "GTEx_v8_se": float(b.se_eQTL.iloc[0]) if len(b) else np.nan,
            "GTEx_v8_p_eQTL": float(b.p_eQTL.iloc[0]) if len(b) else np.nan,
            "GTEx_v8_effect_allele": b.effect_allele.iloc[0] if len(b) else "NA",
            "GTEx_v8_af": float(b.af.iloc[0]) if len(b) else np.nan,
            "GTEx_v8_signif": bool(len(b) and b.p_eQTL.iloc[0] < 0.05),
        })
S = pd.DataFrame(summ)
S.to_csv(os.path.join(TAB, "49b_task34_region_tag_summary.csv"),
         index=False, encoding="utf-8-sig")

hdr = "%-12s %-10s %6s %8s %12s %18s %14s %12s %12s %10s"
w(hdr % ("SNP", "gene", "n_cell", "n_p<.05", "OneK1K_minP", "OneK1K_sig方向", "OneK1K_b范围(p<.05)",
         "GTEx_b", "GTEx_p", "GTEx_p<.05"))
for r in summ:
    w(hdr % (r["SNP_rsID"], r["gene"], r["OneK1K_n_celltypes"], r["OneK1K_n_sig"],
             ("%.4g" % r["OneK1K_min_p"]) if not np.isnan(r["OneK1K_min_p"]) else "NA",
             r["OneK1K_sign_consistency_sig_only"], r["OneK1K_b_range_sig"],
             ("%+.4f" % r["GTEx_v8_b_eQTL"]) if not np.isnan(r["GTEx_v8_b_eQTL"]) else "NA",
             ("%.4g" % r["GTEx_v8_p_eQTL"]) if not np.isnan(r["GTEx_v8_p_eQTL"]) else "NA",
             r["GTEx_v8_signif"]))

# ---------- 3b. 效应等位一致性 ----------
w()
w("=" * 104)
w("[3b] 效应等位一致性（OneK1K 剂量等位 vs GTEx v8 alt vs 结局 GWAS A1）")
w("=" * 104)
gwas_a1 = {}
with io.open(COJO, "r", encoding="utf-8") as f:
    h = f.readline().rstrip("\n").split("\t")
    ic = {c: i for i, c in enumerate(h)}
    for ln in f:
        p = ln.rstrip("\n").split("\t")
        if len(p) == len(h) and p[ic["SNP"]] in targets:
            gwas_a1[V37[int(p[ic["SNP"]].split(":")[1])]] = p[ic["A1"]]
for rs in TOOLS:
    a = one[(one.SNP_rsID == rs)]
    ea_one = a.effect_allele.mode().iloc[0] if len(a) else "NA"
    b = g8[g8.SNP_rsID == rs]
    ea_gt = b.effect_allele.mode().iloc[0] if len(b) else "NA"
    g1 = gwas_a1.get(rs, "NA")
    same = len({ea_one, ea_gt, g1} - {"NA"}) == 1
    w("    %-12s OneK1K 剂量等位=%s ｜ GTEx v8 alt(效应等位)=%s ｜ GWAS A1=%s  → %s"
      % (rs, ea_one, ea_gt, g1, "三层一致 ✅" if same else "⚠不一致，需换算符号"))

# ---------- 4. 区域标签判定 ----------
w()
w("=" * 104)
w("[4] 区域标签 vs 基因标签 判定（阈值 p < 0.05）")
w("=" * 104)
for rs in TOOLS:
    s = S[S.SNP_rsID == rs]
    hit = []
    for g in GENES:
        r = s[s.gene == g].iloc[0]
        layers = []
        if r.OneK1K_n_sig > 0:
            layers.append("OneK1K(%d/14)" % r.OneK1K_n_sig)
        if r.GTEx_v8_signif:
            layers.append("GTEx_v8")
        if layers:
            hit.append("%s[%s]" % (g, "+".join(layers)))
    n_genes = len(hit)
    w("    %-12s hg19=%d / GRCh38=%d" % (rs, TOOLS[rs]["pos37"], TOOLS[rs]["pos38"]))
    w("        显著调控基因数 = %d / 3   → %s" % (n_genes, " ｜ ".join(hit) if hit else "无"))
    w("        判定 = %s"
      % ("**区域标签（region tag）**：同时调控 >=2 个基因"
         if n_genes >= 2 else
         ("单基因标签（gene tag）" if n_genes == 1 else "两层均未达 p<0.05")))

# ---------- 5. 结局 GWAS 上下文 ----------
w()
w("=" * 104)
w("[5] 工具在结局 GWAS（GCST90483463）中的效应（区域标签论证的补充）")
w("=" * 104)
grow = []
with io.open(COJO, "r", encoding="utf-8") as f:
    h = f.readline().rstrip("\n").split("\t")
    ic = {c: i for i, c in enumerate(h)}
    for ln in f:
        p = ln.rstrip("\n").split("\t")
        if len(p) != len(h):
            continue
        snp = p[ic["SNP"]]
        if snp not in targets:
            continue
        rs = V37[int(snp.split(":")[1])]
        grow.append({"SNP_rsID": rs, "pos37_hg19": int(snp.split(":")[1]),
                     "GWAS_A1": p[ic["A1"]], "GWAS_A2": p[ic["A2"]],
                     "GWAS_freq_A1": float(p[ic["freq"]]),
                     "b_GWAS": float(p[ic["b"]]), "se_GWAS": float(p[ic["se"]]),
                     "p_GWAS": float(p[ic["p"]]), "n_GWAS": int(p[ic["n"]])})
        w("    %-12s %-12s A1=%s A2=%s freq(A1)=%.4f  b=%+.5f se=%.5f p=%.4g n=%s"
          % (rs, snp, p[ic["A1"]], p[ic["A2"]], float(p[ic["freq"]]),
             float(p[ic["b"]]), float(p[ic["se"]]), float(p[ic["p"]]), p[ic["n"]]))
# 逐基因的 Wald ratio（仅 GTEx 层工具可算；OneK1K 层已在 3.3 主分析给出）
w()
w("    [GTEx v8 层] 逐基因 单 SNP Wald ratio = b_GWAS / b_eQTL（仅用于展示「同一工具、同一 b_GWAS，")
w("     不同基因 |b_MR| 只由 b_eQTL 决定」这一区域标签逻辑）")
for rs in TOOLS:
    gr = [x for x in grow if x["SNP_rsID"] == rs]
    if not gr:
        continue
    gr = gr[0]
    w("      %s（GWAS b=%+.5f, A1=%s）" % (rs, gr["b_GWAS"], gr["GWAS_A1"]))
    for g in GENES:
        b = g8[(g8.SNP_rsID == rs) & (g8.gene == g)]
        if len(b) == 0:
            continue
        be = float(b.b_eQTL.iloc[0])
        se = float(b.se_eQTL.iloc[0])
        mre = gr["b_GWAS"] / be
        mre_se = float(np.sqrt(gr["se_GWAS"] ** 2 / be ** 2
                               + gr["b_GWAS"] ** 2 * se ** 2 / be ** 4))
        star = ""
        if gr["GWAS_A1"] != b.effect_allele.iloc[0]:
            star = "  ⚠等位不一致（GWAS A1=%s vs 效应等位=%s）" % (
                gr["GWAS_A1"], b.effect_allele.iloc[0])
        w("          %-10s b_eQTL=%+.5f → b_MR=%+.5f (se=%.5f)%s" % (g, be, mre, mre_se, star))
pd.DataFrame(grow).to_csv(os.path.join(TAB, "49c_task34_gwas_context.csv"),
                          index=False, encoding="utf-8-sig")
w()
w("已写 tables/49c_task34_gwas_context.csv（%d 行）" % len(grow))
w("VERDICT = PASS")
flush()
