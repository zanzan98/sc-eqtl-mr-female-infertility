# -*- coding: utf-8 -*-
"""s52b1_diag_L3_N.py -- 决定性检验：parquet af 的分母是否为 980 供者

parquet schema 含 ma_samples / ma_count（minor allele 计数）。
若 ma_count/(2*N) 与 PLINK A1 剂量频率一致，而 af 不一致，
则 af 的分母 N 与 PLINK 面板不同 => LD 参考与 eQTL 面板非同源样本。
"""
import io, os
import numpy as np
import pandas as pd

PROJ = r"D:\endometriosis_project\11_sc_eqtl_mr_project"
PQ = os.path.join(PROJ, "00_data_raw", "onek1k", "parquet_subset")
COND = os.path.join(PROJ, "00_data_raw", "onek1k", "cond_coloc_L2L3")
BFILE = r"D:\endometriosis_project\_onek1k_plink\plink_merged_980_donors"
LOG = r"D:\endometriosis_project\_s52b1_diag_L3N.log"
L = []

def p(s=""):
    s = str(s); L.append(s)
    try: print(s)
    except UnicodeEncodeError: print(s.encode("ascii", "backslashreplace").decode("ascii"))

COLS = ["phenotype_id", "variant_id", "af", "ma_samples", "ma_count", "slope", "slope_se", "pval_nominal"]

for locus, chrom, ct, lo, hi in [("L2", "10", "CD4_NC", 26900000, 27900000),
                                 ("L3", "2", "Mono_NC", 69300000, 70200000)]:
    p("")
    p("=" * 78)
    fp = os.path.join(PQ, "OneK1K_%s.cis_qtl_pairs.chr%s.parquet" % (ct, chrom))
    df = pd.read_parquet(fp, columns=COLS)
    df["pos37"] = df["variant_id"].astype(str).str.split(":").str[1].astype(int)
    sub = df[(df["pos37"] >= lo) & (df["pos37"] <= hi)].copy()
    p("[%s] parquet 全域 rows=%d ；窗口 rows=%d ；unique variant=%d ；dtype af=%s ma_count=%s"
      % (locus, len(df), len(sub), sub["variant_id"].nunique(), sub["af"].dtype, sub["ma_count"].dtype))

    # 用 ma_count 反推总等位计数 N*2
    sub["mac_tot"] = sub["ma_count"] / np.maximum(sub["af"], 1e-12)
    sub["n_hat_af"] = sub["mac_tot"] / 2.0
    sub["n_hat_ma"] = sub["ma_count"] / np.maximum(np.minimum(sub["af"], 1 - sub["af"]), 1e-12) / 2.0
    p("    反推 N（ma_count/af/2）: median=%.2f  q10=%.2f q90=%.2f"
      % (sub["n_hat_af"].median(), sub["n_hat_af"].quantile(.1), sub["n_hat_af"].quantile(.9)))

    # 与 PLINK 剂量比较
    raw = os.path.join(COND, "ld_geno_%s.raw" % locus)
    rw = pd.read_csv(raw, sep=r"\s+")
    snp_cols = [c for c in rw.columns if c not in ("FID","IID","PAT","MAT","SEX","PHENOTYPE")]
    ids = [c.rsplit("_", 1)[0] for c in snp_cols]
    cnt = {c.rsplit("_", 1)[0]: c.rsplit("_", 1)[1] for c in snp_cols}
    X = rw[snp_cols].astype(np.float64).values
    bim = pd.read_csv(BFILE + ".bim", sep="\t", header=None, dtype={0: str},
                      names=["chr","snp","cm","pos","a1","a2"])
    bim = bim[(bim["chr"] == chrom) & (bim["pos"] >= lo) & (bim["pos"] <= hi)].copy()
    bim["a1"] = bim["a1"].str.upper(); bim["a2"] = bim["a2"].str.upper()
    ma1 = dict(zip(bim["snp"], bim["a1"])); ma2 = dict(zip(bim["snp"], bim["a2"]))
    flip = np.array([0 if cnt[s] == ma1.get(s) else (1 if cnt[s] == ma2.get(s) else -1) for s in ids],
                    dtype=np.int8)
    X2 = X.copy(); X2[:, np.where(flip == 1)[0]] = 2.0 - X2[:, np.where(flip == 1)[0]]
    fdose = np.nanmean(X2, axis=0) / 2.0
    p("    PLINK 面板：rows=%d  SNP 列=%d  N=980 ；A1剂量频率粒度=1/1960=%.6f" % (len(rw), len(ids), 1/1960))

    # 聚合到 variant 层
    g = sub.groupby("variant_id").agg(af=("af", "first"), ma_count=("ma_count", "first"),
                                      n_hat=("n_hat_af", "median")).reset_index()
    d = pd.DataFrame(dict(variant_id=ids, fdose=fdose))
    j = g.merge(d, on="variant_id", how="inner")
    j["af_x1960"] = j["af"] * 1960
    j["ma_count_on_af"] = j["ma_count"]
    j["dose_x1960"] = j["fdose"] * 1960
    j["round_dev_af"] = (j["af_x1960"] - j["af_x1960"].round()).abs()
    j["round_dev_dose"] = (j["dose_x1960"] - j["dose_x1960"].round()).abs()
    p("    交集 variant = %d" % len(j))
    p("    af   *1960 距最近整数:  median=%.4f  max=%.4f   （若为 980 供者网格应 ~0）"
      % (j["round_dev_af"].median(), j["round_dev_af"].max()))
    p("    dose *1960 距最近整数:  median=%.4f  max=%.4f"
      % (j["round_dev_dose"].median(), j["round_dev_dose"].max()))
    p("    |dose - af|          :  max=%.5f  mean=%.5f" % ((j["fdose"]-j["af"]).abs().max(),
                                                           (j["fdose"]-j["af"]).abs().mean()))
    # 关键：ma_count 是否等于 dose*1960（即 PLINK A1 计数）
    j["ma_count_x2"] = j["ma_count"] * 2
    j["dose_cnt"] = (j["fdose"] * 1960)
    p("    |ma_count*2 - dose*1960|: median=%.2f  max=%.2f   （若同源应为 0）"
      % ((j["ma_count_x2"] - j["dose_cnt"]).abs().median(),
         (j["ma_count_x2"] - j["dose_cnt"]).abs().max()))
    p("    ma_count 分布样本（前 8）：")
    for _, r in j.head(8).iterrows():
        p("       %-14s af=%.6f  af*1960=%.2f  ma_count=%.0f  dose=%.6f  dose*1960=%.2f  n_hat=%.2f"
          % (r["variant_id"], r["af"], r["af_x1960"], r["ma_count"], r["fdose"],
             r["dose_x1960"], r["n_hat"]))
    # 反推 L3 的 af 分母：用 af 与 dose 的频率差是否集中在整数网格
    j["n_grid"] = (j["af"] * 2 * 980)
    p("    af 相对 dose 的偏置（af-dose）: mean=%+.5f  median=%+.5f  >0 占比=%.3f"
      % ((j["af"]-j["fdose"]).mean(), (j["af"]-j["fdose"]).median(),
         float((j["af"] > j["fdose"]).mean())))

io.open(LOG, "w", encoding="utf-8", newline="\n").write("\n".join(L) + "\n")
print("WROTE " + LOG)
