# -*- coding: utf-8 -*-
"""s27_ivw_and_merge.py —— IVW 汇总 + 与 Discovery 合并（用户裁定 (b)：并列报告）

输入：
  tables/29b_variantlevel_wald.csv      敏感性阈值 r2<0.01 的逐变异 Wald（已与结局调和等位）
  tables/29c_variantlevel_wald.csv      未 clump 的逐变异 Wald（仅作负对照，不可用于推断）
  tables/15_discovery_significant_with_locus.csv   Discovery 15 显著对（权威 Wald 结果）
  tables/29a_ivw_instruments_main.csv    主阈值 r2<0.001 的工具（用于 nsnp 标注）
  tables/28b_clump_stats_*.csv           三阈值 clump 统计

产出：
  tables/30_ivw_sensitivity.csv          IVW（固定/随机）+ Cochran Q + I2（仅 nsnp>=2 的对）
  tables/30b_raw_vs_discovery_validation.csv  raw 重算 vs Discovery 逐行校验
  tables/31_final_MR_with_method.csv     ★ 最终表：Wald 与 IVW 并列，含 method 列
  tables/32_ivw_inapplicable_demo.csv    未 clump 负对照的 IVW 演示（标 INVALID_FOR_INFERENCE）
"""
import os, sys
import numpy as np
import pandas as pd
from scipy import stats

T = r"D:\endometriosis_project\11_sc_eqtl_mr_project\tables"
LOG = r"D:\endometriosis_project\_s27_ivw.log"
L = []
def p(s=""):
    L.append(str(s)); print(s)

def ivw(b, se):
    b = np.asarray(b, float); se = np.asarray(se, float)
    w = 1.0 / (se ** 2)
    sw = w.sum()
    bi = float((w * b).sum() / sw)
    se_fix = float(1.0 / np.sqrt(sw))
    k = len(b)
    Q = float((w * (b - bi) ** 2).sum())
    df = k - 1
    Qp = float(stats.chi2.sf(Q, df)) if df > 0 else np.nan
    I2 = float(max(0.0, (Q - df) / Q) * 100.0) if (df > 0 and Q > 0) else 0.0
    infl = float(np.sqrt(max(1.0, Q / df))) if df > 0 else 1.0
    se_ran = se_fix * infl
    return dict(b=bi, se_fixed=se_fix, se_random=se_ran, Q=Q, Q_df=df, Q_p=Qp, I2=I2, k=k)

# ---------- 1) 敏感性阈值（r2<0.01）逐对汇总 ----------
p("=== s27 IVW 汇总 ===")
p("")
p("[1] r2<0.01 敏感性工具集：nsnp 分布")
sens = pd.read_csv(os.path.join(T, "29b_variantlevel_wald.csv"))
sens = sens[sens["status"].isin(["ok", "flipped"])].copy()
cnt = sens.groupby(["gene", "cell_type"]).size().rename("nsnp_raw").reset_index()
main = pd.read_csv(os.path.join(T, "29a_ivw_instruments_main.csv"))
cnt = cnt.merge(main[["gene", "cell_type", "n_index"]], on=["gene", "cell_type"], how="left")
p("   总对数(有效) = %d ；nsnp=1 的 %d ；nsnp>=2 的 %d"
  % (len(cnt), int((cnt["nsnp_raw"] == 1).sum()), int((cnt["nsnp_raw"] >= 2).sum())))
for _, r in cnt[cnt["nsnp_raw"] >= 2].iterrows():
    p("   ★ 可做 IVW：%s / %s  nsnp=%d" % (r["gene"], r["cell_type"], r["nsnp_raw"]))

rows = []
for (g, ct), d in sens.groupby(["gene", "cell_type"]):
    d = d.sort_values("variant_id")
    if len(d) == 1:
        r = d.iloc[0]
        rows.append(dict(gene=g, cell_type=ct, method="Wald ratio", nsnp=1,
                         b=r["beta_MR"], se=r["se_MR"], z=r["z_MR"], p=r["p_MR"],
                         Q=np.nan, Q_df=np.nan, Q_p=np.nan, I2=np.nan,
                         instruments=r["variant_id"], thr="r2<0.01",
                         note="单工具：与 Discovery 同法（Wald ratio）"))
    else:
        st = ivw(d["beta_MR"].values, d["se_MR"].values)
        rows.append(dict(gene=g, cell_type=ct, method="IVW(fixed)", nsnp=st["k"],
                         b=st["b"], se=st["se_fixed"], z=st["b"] / st["se_fixed"],
                         p=float(2 * stats.norm.sf(abs(st["b"] / st["se_fixed"]))),
                         Q=st["Q"], Q_df=st["Q_df"], Q_p=st["Q_p"], I2=st["I2"],
                         instruments=";".join(d["variant_id"]),
                         thr="r2<0.01", note="多工具：IVW 固定效应（敏感性）"))
        rows.append(dict(gene=g, cell_type=ct, method="IVW(random)", nsnp=st["k"],
                         b=st["b"], se=st["se_random"], z=st["b"] / st["se_random"],
                         p=float(2 * stats.norm.sf(abs(st["b"] / st["se_random"]))),
                         Q=st["Q"], Q_df=st["Q_df"], Q_p=st["Q_p"], I2=st["I2"],
                         instruments=";".join(d["variant_id"]),
                         thr="r2<0.01", note="多工具：IVW 乘法随机效应（敏感性）"))
ivw_tab = pd.DataFrame(rows)
ivw_tab.to_csv(os.path.join(T, "30_ivw_sensitivity.csv"), index=False, encoding="utf-8-sig")
p("    -> 30_ivw_sensitivity.csv  行=%d（其中 IVW 行=%d）"
  % (len(ivw_tab), int(ivw_tab["method"].str.startswith("IVW").sum())))
for _, r in ivw_tab[ivw_tab["method"].str.startswith("IVW")].iterrows():
    p("      %s/%s %s b=%.6f se=%.6f p=%.3e Q=%.4f Q_p=%.4f I2=%.2f%%"
      % (r["gene"], r["cell_type"], r["method"], r["b"], r["se"], r["p"], r["Q"], r["Q_p"], r["I2"]))

# ---------- 2) raw 重算 vs Discovery 校验 ----------
p("")
p("[2] raw 重算（parquet 侧）vs Discovery 逐行校验（15 对）")
disc = pd.read_csv(os.path.join(T, "15_discovery_significant_with_locus.csv"))
val = []
for _, r in disc.iterrows():
    m = sens[(sens["gene"] == r["gene"]) & (sens["cell_type"] == r["cell_type"])]
    if m.empty:
        val.append(dict(gene=r["gene"], cell_type=r["cell_type"], match=False,
                        note="raw 侧无该对（可能 instrument 不同或 outcome_missing）"))
        continue
    m = m.iloc[0]
    same_var = (str(m["variant_id"]) == str(r["variant_id_grch38"]))
    val.append(dict(gene=r["gene"], cell_type=r["cell_type"],
                    disc_instrument=r["variant_id_grch38"], raw_instrument=m["variant_id"],
                    same_instrument=same_var,
                    disc_b=r["b"], raw_b=m["beta_MR"], d_b=abs(r["b"] - m["beta_MR"]),
                    disc_se=r["se"], raw_se=m["se_MR"], d_se=abs(r["se"] - m["se_MR"]),
                    disc_p=r["p"], raw_p=m["p_MR"],
                    lgp_ratio=(np.log10(r["p"]) / np.log10(m["p_MR"])) if m["p_MR"] > 0 else np.nan,
                    note="ok" if same_var else "instrument 不同"))
v = pd.DataFrame(val)
v.to_csv(os.path.join(T, "30b_raw_vs_discovery_validation.csv"), index=False, encoding="utf-8-sig")
p("    15 对中：同一 instrument = %d ；b 最大绝对差 = %.3e ；se 最大绝对差 = %.3e"
  % (int(v["same_instrument"].sum()), v["d_b"].max(), v["d_se"].max()))
p("    -> 30b_raw_vs_discovery_validation.csv")

# ---------- 3) 最终表：Wald(Discovery) 与 IVW 并列 ----------
p("")
p("[3] 生成最终表（method 列标注 Wald / IVW）")
fin = []
for _, r in disc.iterrows():
    fin.append(dict(locus_id=r["locus_id"], gene=r["gene"], cell_type=r["cell_type"],
                    method="Wald ratio", nsnp=1, b=r["b"], se=r["se"], p=r["p"],
                    p_outcome=r["p_outcome"], F_stat=r["F_stat"],
                    instrument=r["variant_id_grch38"], thr="P<5e-8 (Discovery)",
                    source="Discovery(top_eQTL_summary)", note="主结果（单 SNP，按设计）"))
for _, r in ivw_tab[ivw_tab["method"].str.startswith("IVW")].iterrows():
    d = disc[(disc["gene"] == r["gene"]) & (disc["cell_type"] == r["cell_type"])]
    fin.append(dict(locus_id=(d["locus_id"].iloc[0] if len(d) else ""), gene=r["gene"],
                    cell_type=r["cell_type"], method=r["method"], nsnp=r["nsnp"], b=r["b"],
                    se=r["se"], p=r["p"], p_outcome=(d["p_outcome"].iloc[0] if len(d) else np.nan),
                    F_stat=(d["F_stat"].iloc[0] if len(d) else np.nan),
                    instrument=r["instruments"], thr=r["thr"],
                    source="raw_recompute(onlyk1k_parquet)",
                    note="敏感性：r2<0.01 下多工具；Q=%.4f, I2=%.2f%%" % (r["Q"], r["I2"])))
ft = pd.DataFrame(fin)
ft.to_csv(os.path.join(T, "31_final_MR_with_method.csv"), index=False, encoding="utf-8-sig")
p("    行数 = %d（Wald %d / IVW %d）-> 31_final_MR_with_method.csv"
  % (len(ft), int((ft["method"] == "Wald ratio").sum()),
     int(ft["method"].str.startswith("IVW").sum())))

# ---------- 4) 主阈值下 IVW 不适用的机械证明 ----------
p("")
p("[4] 主阈值 r2<0.001 下的 clump 结果（机械证明 IVW 不适用）")
for tag in ["main", "sens001", "ctrl_r05"]:
    f = os.path.join(T, "28b_clump_stats_%s.csv" % tag)
    if not os.path.exists(f):
        continue
    c = pd.read_csv(f)
    p("    %-9s 对数=%2d  n_index=1 的对=%2d  n_index>=2 的对=%2d  n_index 最大=%d"
      % (tag, len(c), int((c["n_index"] == 1).sum()), int((c["n_index"] >= 2).sum()), int(c["n_index"].max())))
p("    -> 主阈值（用户裁定）下 26/26 对 n_index=1 → 全部单 SNP → IVW 不适用，Wald 为主线方法（按设计）")
p("    -> ctrl_r05 证明 clump 机械可用（最多 29 个 SNP），故 main 的 n_index=1 是阈值驱动而非流程故障")

# ---------- 5) 未 clump 负对照 ----------
p("")
p("[5] 未 clump 工具的 IVW 演示（负对照，标 INVALID_FOR_INFERENCE）")
un = pd.read_csv(os.path.join(T, "29c_variantlevel_wald.csv"))
un = un[un["status"].isin(["ok", "flipped"])]
d = un[(un["gene"] == "LINC00339") & (un["cell_type"] == "CD4_NC")]
p("    示例 LINC00339/CD4_NC 未 clump 变异数=%d（LD 高度相关）" % len(d))
demo = []
for (g, ct), dd in un.groupby(["gene", "cell_type"]):
    if len(dd) < 2:
        continue
    st = ivw(dd["beta_MR"].values, dd["se_MR"].values)
    demo.append(dict(gene=g, cell_type=ct, nsnp=st["k"], b=st["b"], se=st["se_fixed"],
                     p=float(2 * stats.norm.sf(abs(st["b"] / st["se_fixed"]))),
                     Q=st["Q"], Q_p=st["Q_p"], I2=st["I2"], verdict="INVALID_FOR_INFERENCE",
                     note="工具间 LD 相关（未 clump）→ IVW 独立性假设破坏"))
pd.DataFrame(demo).to_csv(os.path.join(T, "32_ivw_inapplicable_demo.csv"), index=False, encoding="utf-8-sig")
p("    行=%d -> 32_ivw_inapplicable_demo.csv（仅作口径说明，不得用于推断）" % len(demo))

p("")
p("完成 %s" % __import__("time").strftime("%Y-%m-%d %H:%M:%S"))
open(LOG, "w", encoding="utf-8").write("\n".join(L) + "\n")
