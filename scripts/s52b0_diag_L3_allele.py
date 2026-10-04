# -*- coding: utf-8 -*-
"""s52b0_diag_L3_allele.py -- 诊断 L3 判据A 偏差来源（只读）

背景：s52a 阶段B 中 L2 判据A max=0.0000（完美），L3 max=0.0199 > 0.01 阈值被拦。
本脚本定位：哪些变异偏差大、是否集中、与 MAF / 等位 / 位置的关系。
"""
import io, os, sys
import numpy as np
import pandas as pd

PROJ = r"D:\endometriosis_project\11_sc_eqtl_mr_project"
COND = os.path.join(PROJ, "00_data_raw", "onek1k", "cond_coloc_L2L3")
BFILE = r"D:\endometriosis_project\_onek1k_plink\plink_merged_980_donors"
LOG = r"D:\endometriosis_project\_s52b0_diag_L3.log"
L = []

def p(s=""):
    s = str(s); L.append(s)
    try: print(s)
    except UnicodeEncodeError: print(s.encode("ascii", "backslashreplace").decode("ascii"))

def dose_freq(locus, chrom, lo, hi):
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
    map_a1 = dict(zip(bim["snp"], bim["a1"])); map_a2 = dict(zip(bim["snp"], bim["a2"]))
    pos_of = dict(zip(bim["snp"], bim["pos"]))
    # 与 s52a 完全相同的翻转逻辑
    flip = np.array([0 if cnt[s] == map_a1.get(s) else (1 if cnt[s] == map_a2.get(s) else -1)
                     for s in ids], dtype=np.int8)
    X2 = X.copy()
    X2[:, np.where(flip == 1)[0]] = 2.0 - X2[:, np.where(flip == 1)[0]]
    f = np.nanmean(X2, axis=0) / 2.0
    n_missing = int(np.isnan(X2).sum())
    return pd.DataFrame(dict(snp_grch37=ids, pos37=[pos_of.get(s) for s in ids],
                             a1_bim=[map_a1.get(s) for s in ids],
                             a2_bim=[map_a2.get(s) for s in ids],
                             recode_cnt_allele=[cnt[s] for s in ids],
                             flip=flip, dose_freq=f)), n_missing, bim

p("=== s52b0 诊断：L3 判据A 偏差定位 ===")
for locus, chrom, lo, hi in [("L2", "10", 26900000, 27900000), ("L3", "2", 69300000, 70200000)]:
    p("")
    p("=" * 78)
    m = pd.read_csv(os.path.join(COND, "m_%s_%s.csv" % (locus, "CD4_NC" if locus == "L2" else "Mono_NC")))
    d, nmiss, bim = dose_freq(locus, chrom, lo, hi)
    p("[%s] m rows=%d  ld rows=%d  missing_dosage=%d" % (locus, len(m), len(d), nmiss))
    j = m.merge(d.rename(columns={"flip": "flip_ld"}), on="snp_grch37", how="inner")
    p("     交集 = %d" % len(j))
    dev = np.abs(j["dose_freq"] - j["af"])
    p("     dev vs eQTL af : max=%.5f mean=%.5f median=%.5f  n(>0.01)=%d n(>0.005)=%d n(>0.002)=%d"
      % (dev.max(), dev.mean(), np.median(dev), int((dev > 0.01).sum()),
         int((dev > 0.005).sum()), int((dev > 0.002).sum())))
    # 对照：dev vs 1-af（是否等价于等位翻转）
    dev_f = np.abs(j["dose_freq"] - (1 - j["af"]))
    p("     dev vs 1-eQTL af: max=%.5f mean=%.5f" % (dev_f.max(), dev_f.mean()))
    j["dev"] = dev; j["dev_flip"] = dev_f
    j["maf"] = np.minimum(j["af"], 1 - j["af"])
    bad = j[j["dev"] > 0.01].sort_values("dev", ascending=False)
    p("")
    p("     >>> dev>0.01 的变异（共 %d）：" % len(bad))
    cols = ["snp_grch37", "pos37", "a1_bim", "a2_bim", "recode_cnt_allele", "flip_ld",
            "af", "dose_freq", "dev", "dev_flip", "maf", "eaf"]
    for _, r in bad.head(30).iterrows():
        p("       %-14s pos=%d  bim a1=%s a2=%s  recode=%s flip=%d | eQTL af=%.4f  dose=%.4f "
          "dev=%.4f dev_flip=%.4f  MAF_e=%.3f  GWAS eaf=%.4f"
          % (r["snp_grch37"], int(r["pos37"]), r["a1_bim"], r["a2_bim"], r["recode_cnt_allele"],
             int(r["flip_ld"]), r["af"], r["dose_freq"], r["dev"], r["dev_flip"], r["maf"], r["eaf"]))
    if len(bad):
        p("")
        p("     偏差变异的位置分布：min_pos=%d max_pos=%d ；中位 pos=%d"
          % (bad["pos37"].min(), bad["pos37"].max(), int(bad["pos37"].median())))
        # 是否与 MAF 相关：分 MAF 箱统计 dev
        p("")
        p("     按 eQTL MAF 分箱的 dev（看偏差是否随 MAF 增大）：")
        for lo_m, hi_m in [(0.0, 0.05), (0.05, 0.1), (0.1, 0.2), (0.2, 0.3), (0.3, 0.5)]:
            sel = j[(j["maf"] >= lo_m) & (j["maf"] < hi_m)]
            if len(sel) == 0:
                continue
            p("       MAF [%.2f,%.2f) n=%-5d  dev max=%.5f mean=%.5f  dev_flip mean=%.5f"
              % (lo_m, hi_m, len(sel), sel["dev"].max(), sel["dev"].mean(), sel["dev_flip"].mean()))

io.open(LOG, "w", encoding="utf-8", newline="\n").write("\n".join(L) + "\n")
print("WROTE " + LOG)
