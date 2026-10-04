# -*- coding: utf-8 -*-
"""Discovery 小任务交付物：
19_locus_summary.csv            —— 15 个显著对压缩为 3 个独立信号座
20_celltype_power.csv           —— 14 种细胞类型的工具数/检验数/零显著标记
21_celltype_mapping.csv         —— OneK1K 14 精细亚型 -> scBloodNL 9 粗分群 映射
22_replication_test_plan.csv    —— 复现检验清单（按 Tier 分层 + 粗分群折叠）
"""
import os
import numpy as np
import pandas as pd

ROOT = r"D:\endometriosis_project\11_sc_eqtl_mr_project"
T = os.path.join(ROOT, "tables")
LOG = []


def w(s):
    LOG.append(str(s))


sig = pd.read_csv(os.path.join(T, "15_discovery_significant_with_locus.csv"))
allres = pd.read_csv(os.path.join(T, "10_discovery_MR_main.csv"))
iv = pd.read_csv(os.path.join(T, "01b_discovery_instruments_grch38.csv"),
                 dtype={"variant_id": str, "variant_id_grch37": str})
qa = pd.read_csv(os.path.join(T, "13_discovery_sensitivity_qval05_significant.csv"))

w("=== 输入 ===")
w("15 显著对 = %d 行 ; 10 全量 = %d 行 ; 01b 工具 = %d 行 ; 13 qval敏感 = %d 行"
  % (len(sig), len(allres), len(iv), len(qa)))

# ---------------- 19. 基因座汇总 ----------------
rows = []
for lid, g in sig.groupby("locus_id", sort=True):
    pos = g["tss_pos_grch38"].dropna()
    genes = sorted(g["gene"].unique())
    rows.append(dict(
        locus=lid,
        chr=str(g["chr"].iloc[0]),
        position_Mb="%.2f-%.2f" % (pos.min() / 1e6, pos.max() / 1e6),
        genes=";".join(genes),
        n_genes=len(genes),
        cell_types=";".join(sorted(g["cell_type"].unique())),
        n_cell_types=g["cell_type"].nunique(),
        n_pairs=len(g),
        max_b=round(g["b"].max(), 4),
        max_abs_b=round(g["b"].abs().max(), 4),
        b_direction_of_gene=";" .join("%s=%s" % (gn, "+" if gg["b"].mean() > 0 else "-")
                                      for gn, gg in g.groupby("gene")),
        min_p="%.3g" % g["p"].min(),
        max_F=round(g["F_stat"].max(), 1),
        p_outcome_range="%.3g ~ %.3g" % (g["p_outcome"].min(), g["p_outcome"].max()),
        n_pairs_in_outcome_GWS=int(g["outcome_GWS"].sum()),
        in_qval_sensitivity=int(g["in_qval_sensitivity"].sum()),
        is_MHC=bool(g["is_MHC"].any()),
        distinct_instruments=g["variant_id_grch38"].nunique(),
        instruments=";".join(sorted(g["variant_id_grch38"].unique())),
    ))
locus = pd.DataFrame(rows)
locus.to_csv(os.path.join(T, "19_locus_summary.csv"), index=False, encoding="utf-8-sig")
w("")
w("=== 19. 基因座汇总 ===")
w(locus[["locus", "position_Mb", "n_genes", "genes", "n_pairs",
         "max_abs_b", "min_p", "p_outcome_range", "is_MHC"]].to_string(index=False))
w("独立基因座数 = %d（显著对 %d 个 -> 压缩比 %.1f:1）" % (len(locus), len(sig), len(sig) / len(locus)))
w("唯一工具变异总数 = %d" % sig["variant_id_grch38"].nunique())

# ---------------- 20. 细胞类型功效表 ----------------
cnt = iv.groupby("cell_type").agg(n_instruments=("variant_id", "size"),
                                  n_genes=("gene", "nunique")).reset_index()
tst = allres.groupby("cell_type").agg(n_tested=("p", "size"),
                                      n_nominal=("p", lambda s: int((s < 0.05).sum())),
                                      n_FDR05=("qval", lambda s: int((s < 0.05).sum()))).reset_index()
pw = cnt.merge(tst, on="cell_type", how="outer").fillna(0)
pw["n_instruments"] = pw["n_instruments"].astype(int)
pw["pct_nominal"] = (100.0 * pw["n_nominal"] / pw["n_tested"].replace(0, np.nan)).round(2)
pw["power_flag"] = np.where(pw["n_FDR05"] > 0, "signal",
                            np.where(pw["n_instruments"] < 200, "LOW_POWER_DO_NOT_READ_AS_NULL", "no_signal"))
pw = pw.sort_values("n_instruments", ascending=False)
pw.to_csv(os.path.join(T, "20_celltype_power.csv"), index=False, encoding="utf-8-sig")
w("")
w("=== 20. 细胞类型功效 ===")
w(pw.to_string(index=False))
zero = pw[pw["n_FDR05"] == 0]
w("")
w("零显著细胞类型 = %d : %s" % (len(zero), "; ".join(
    "%s(%d 工具)" % (r["cell_type"], r["n_instruments"]) for _, r in zero.iterrows())))

# ---------------- 21. 细胞类型映射 ----------------
MAP = [
    ("CD4_NC", "CD4T", "sc-eQTL 发布层为 lowerres；naive/central-memory CD4T 合并"),
    ("CD4_ET", "CD4T", "effector/Th1/Th2/Th17 CD4T 合并入 CD4T"),
    ("CD4_SOX4", "CD4T", "SOX4+ CD4T 并入 CD4T"),
    ("CD8_NC", "CD8T", "naive/central-memory CD8T"),
    ("CD8_ET", "CD8T", "effector CD8T"),
    ("CD8_S100B", "CD8T", "S100B+ CD8T 并入 CD8T"),
    ("NK", "NK", "NK 主群"),
    ("NK_R", "NK", "NKdim/NKbright 并入 NK"),
    ("B_MEM", "B", "scBloodNL 发布层未分 naive/memory B"),
    ("B_IN", "B", "同上"),
    ("Mono_NC", "monocyte", "mono 1/2/4 合并入 monocyte"),
    ("Mono_C", "monocyte", "同上"),
    ("DC", "DC", "mDC+pDC 合并入 DC"),
    ("Plasma", "plasma B", "★ scBloodNL 有 plasma B 注释（710 细胞）但未发布独立 eQTL 文件 -> 无法复现"),
]
mp = pd.DataFrame(MAP, columns=["cell_type_OneK1K", "cell_type_scBloodNL", "note"])
mp["available_in_scbloodnl"] = mp["cell_type_scBloodNL"].isin(
    ["CD4T", "CD8T", "NK", "B", "monocyte", "DC", "megakaryocyte"])
mp.to_csv(os.path.join(T, "21_celltype_mapping.csv"), index=False, encoding="utf-8-sig")
w("")
w("=== 21. 细胞类型映射（多对一）===")
w(mp.to_string(index=False))

# ---------------- 22. 复现检验清单 ----------------
TIER1 = {("YME1L1", "CD4_NC"), ("ANXA4", "Mono_NC")}
pool = pd.read_csv(os.path.join(T, "18_replication_candidate_pool.csv"),
                   dtype={"variant_id": str, "variant_id_grch37": str,
                          "variant_id_grch38": str})
t = sig.copy()
t["tier"] = np.where([(a, b) in TIER1 for a, b in zip(t["gene"], t["cell_type"])],
                     "Tier1", "Tier2_chr1")
t = t.merge(pool[["gene", "cell_type", "variant_id_grch37", "effect_allele",
                  "other_allele", "af"]], on=["gene", "cell_type"], how="left")
t = t.merge(mp[["cell_type_OneK1K", "cell_type_scBloodNL", "available_in_scbloodnl"]],
            left_on="cell_type", right_on="cell_type_OneK1K", how="left")
plan = t[["tier", "locus_id", "gene", "cell_type", "cell_type_scBloodNL",
          "available_in_scbloodnl", "variant_id_grch38", "variant_id_grch37",
          "b", "se", "p", "qval", "F_stat",
          "beta_outcome", "se_outcome", "p_outcome", "outcome_GWS"]].rename(columns={
              "b": "discovery_b", "se": "discovery_se", "p": "discovery_p", "qval": "discovery_qval",
              "F_stat": "discovery_F"})
plan = plan.sort_values(["tier", "gene", "cell_type_scBloodNL", "gene"])
plan.to_csv(os.path.join(T, "22_replication_test_plan.csv"), index=False, encoding="utf-8-sig")
w("")
w("=== 22. 复现检验清单 ===")
w(plan[["tier", "gene", "cell_type", "cell_type_scBloodNL", "variant_id_grch37",
        "discovery_b", "discovery_p", "p_outcome"]].to_string(index=False))
w("")
coll = (plan.groupby(["gene", "cell_type_scBloodNL"])
        .agg(n_discovery_pairs=("discovery_p", "size"),
             n_instruments=("variant_id_grch37", "nunique")).reset_index().sort_values(["gene", "cell_type_scBloodNL"]))
w("折叠后复现检验数（按 基因 x 粗分群 去重）= %d" % len(coll))
w(coll.to_string(index=False))

with open(r"D:\endometriosis_project\_small_tasks.txt", "w", encoding="utf-8") as fh:
    fh.write("\n".join(LOG))
print("done")
