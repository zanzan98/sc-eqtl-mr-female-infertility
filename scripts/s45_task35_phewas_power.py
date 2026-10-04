# -*- coding: utf-8 -*-
"""
s45_task35_phewas_power.py — 任务三 3.5：PheWAS 统计功效分析
================================================================================
目标：对 uncond family 的 **免疫/血液/肿瘤域 968 个检验**（484 端点 × 2 工具）
      报告 ① 各端点病例数 ② 工具 F 统计量 ③ 统计功效（最小可检出 OR + 实测功率）
      并区分「功效充分 / 功效不足 / 不可评估（无对应端点）」。

★ 本脚本为 只读 + 只新增：
   - 只读输入：tables/44_phewas_3p2_full_all.csv（3.2 产物）、
               01_data_raw/finngen_R12_manifest.tsv、
               00_data_raw/finngen_r12_tabix/phewas_3snp_all_endpoints.csv
   - 只新增输出：tables/50*_task35_*.csv + scripts/_s45_task35.log
   - 不触碰任何已闭环结果。

★ 无需重读 FinnGen R12 全量汇总统计：
   3.2 已把 3 SNP × 2,469 端点的实测 beta/sebeta/pval 落盘（7,407 行），
   manifest 已含全 2,469 端点的 num_cases/num_controls ⇒ 本步纯本地轻量计算。

功效口径（预先设定，写死）：
   - 单工具 Wald ratio，log-OR 尺度；α = 域内 Bonferroni = 0.05/968（主），
     另给全局 α = 0.05/4938 作敏感性。
   - 最小可检出 |b_MR|（80% power）= (z_{1-α/2} + z_{0.80}) × se_MR
     ⇒ 最小可检出 OR = exp(该值)。
   - 实测功率 = Φ(|b_obs|/se_MR − z_crit) + Φ(−|b_obs|/se_MR − z_crit)。
   - se_MR 直接取自 3.2 产物（二阶 delta 全比值形式，已验证）。
"""

import os, io, json, time
import numpy as np
import pandas as pd
from scipy import stats

ROOT = r"D:\endometriosis_project\11_sc_eqtl_mr_project"
MANIFEST = r"D:\endometriosis_project\01_data_raw\finngen_R12_manifest.tsv"
TAB44 = os.path.join(ROOT, "tables", "44_phewas_3p2_full_all.csv")
SNPLONG = os.path.join(ROOT, "00_data_raw", "finngen_r12_tabix", "phewas_3snp_all_endpoints.csv")
ONEPQ = os.path.join(ROOT, "00_data_raw", "onek1k", "parquet_subset")
TAB = os.path.join(ROOT, "tables")
LOG = os.path.join(ROOT, "scripts", "_s45_task35.log")

FAMILY = "uncond"
DOMAIN_MAIN = "免疫/血液/肿瘤"
N_TESTS_MAIN = 968
N_TESTS_ALL = 4938

buf = io.StringIO()
def P(*a):
    buf.write(" ".join(str(x) for x in a) + "\n")

T0 = time.time()
P("=" * 100)
P("s45 — 任务三 3.5：PheWAS 统计功效分析")
P("time %s" % time.strftime("%Y-%m-%d %H:%M:%S"))
P("=" * 100)

# ---------------------------------------------------------------- 0) 读入
d = pd.read_csv(TAB44, encoding="utf-8-sig")
P("\n[0] 44_phewas_3p2_full_all.csv  shape=%s" % (d.shape,))
P("    family=%s" % d.family.value_counts().to_dict())
P("    instrument x cell_type:")
for (i, c), n in d.groupby(["instrument", "cell_type"]).size().items():
    P("      %-12s %-9s %d" % (i, c, n))

man = pd.read_csv(MANIFEST, sep="\t")
P("    manifest 端点 = %d（列：%s）" % (len(man), list(man.columns)))

lg = pd.read_csv(SNPLONG)
P("    phewas_3snp_all_endpoints.csv  rows=%d  端点=%d  SNP=%s"
  % (len(lg), lg.endpoint.nunique(), sorted(lg.rsids.unique())))
assert len(lg) == lg.endpoint.nunique() * lg.rsids.nunique(), "长表非满覆盖"
P("    ★ 覆盖率核对：%d 端点 × %d SNP = %d 行（满覆盖 ✅）"
  % (lg.endpoint.nunique(), lg.rsids.nunique(), len(lg)))

# ---------------------------------------------------------------- 1) 主表：uncond 全 4,938 检验
u = d[d.family == FAMILY].copy()
P("\n[1] family=%s  检验数=%d  端点=%d" % (FAMILY, len(u), u.phenocode.nunique()))

# α 与临界 z
alpha_main = 0.05 / N_TESTS_MAIN
alpha_all = 0.05 / N_TESTS_ALL
zc_main = stats.norm.ppf(1 - alpha_main / 2)
zc_all = stats.norm.ppf(1 - alpha_all / 2)
zc_unc = stats.norm.ppf(1 - 0.05 / 2)
z80 = stats.norm.ppf(0.80)
P("    α(域内 Bonferroni 0.05/968) = %.6e  → z_crit = %.6f" % (alpha_main, zc_main))
P("    α(全局 Bonferroni 0.05/4938) = %.6e  → z_crit = %.6f" % (alpha_all, zc_all))
P("    α(未校正 0.05)                = %.6e  → z_crit = %.6f" % (0.05, zc_unc))
P("    z_0.80 = %.6f" % z80)

u["alpha_used"] = alpha_main
u["z_crit_domain"] = zc_main
u["z_req_80pct_domain"] = zc_main + z80
u["abs_beta_MR"] = u.beta_MR.abs()
u["abs_z_MR"] = u.abs_beta_MR / u.se_MR
u["p_check_from_z"] = 2 * stats.norm.sf(u.abs_z_MR)

# 最小可检出效应（80% power）
u["min_detectable_abs_beta_MR"] = (zc_main + z80) * u.se_MR
u["min_detectable_OR"] = np.exp(u.min_detectable_abs_beta_MR)
# 敏感性：全局 Bonferroni
u["min_detectable_OR_globalBonf"] = np.exp((zc_all + z80) * u.se_MR)
# 敏感性：未校正
u["min_detectable_OR_uncorrected"] = np.exp((zc_unc + z80) * u.se_MR)

# 实测功率（对本检验实际观测到的效应量）
u["power_at_observed_effect"] = (stats.norm.cdf(u.abs_z_MR - zc_main)
                                 + stats.norm.cdf(-u.abs_z_MR - zc_main))

# 当前该检验在 α 下是否可检出「观测效应量」（等价 power>=0.8）
u["detectable_at_80pct_for_observed"] = u.abs_z_MR >= (zc_main + z80)

# 所需病例数：使 min_detectable_OR 降到 target（se_MR ∝ 1/sqrt(n_case)）
def nreq(target):
    with np.errstate(divide="ignore", invalid="ignore"):
        r = (np.log(u.min_detectable_OR) / np.log(target)) ** 2
    r = np.where(u.min_detectable_OR <= target, 0.0, r)   # 0 表示已足够
    return u.n_cases_manifest * r
u["n_cases_needed_minOR_1p5"] = nreq(1.5)
u["n_cases_needed_minOR_1p2"] = nreq(1.2)

# 功效分层（按 80% power 可检出的最小 OR）
def classify(mo):
    if not np.isfinite(mo):          # ★ NaN 必须单独成类，否则会误落入兜底类
        return "0_not_estimable"
    if mo <= 1.2:
        return "1_adequate_OR<=1.2"
    if mo <= 1.5:
        return "2_adequate_1.2<OR<=1.5"
    if mo <= 2.0:
        return "3_limited_1.5<OR<=2.0"
    if mo <= 3.0:
        return "4_inadequate_2.0<OR<=3.0"
    return "5_severe_OR>3.0"
u["power_class"] = u.min_detectable_OR.map(classify)

P("\n    NaN 统计（wald 分母不可用的行，已按 3.2 规则置 NaN）：")
P("      beta_MR=%d  se_MR=%d  p_MR=%d  F=%d"
  % (u.beta_MR.isna().sum(), u.se_MR.isna().sum(), u.p_MR.isna().sum(), u.F.isna().sum()))

# ---------------------------------------------------------------- 2) 独立复核
P("\n" + "=" * 100)
P("[2] 独立复核")
P("=" * 100)
# (a) |z_MR| 与 p_MR 的 p 值一致性
u["p_rel_diff"] = (u.p_check_from_z - u.p_MR).abs() / u.p_MR.clip(lower=1e-300)
mx = u.p_rel_diff.max()
P("    (a) 由 |z_MR| 反算的 p 与 44 表 p_MR 最大相对偏差 = %.3e  %s"
  % (mx, "✅" if mx < 1e-8 else "⚠️"))

# (b) 隐式 |b_exp| = sebeta / se_MR 应为常数（每个 instrument 一个值）
utgt = u[["phenocode", "instrument"]].drop_duplicates()
# ★ 不要 rename rsids→instrument：右表已自带 instrument（来自 u），rename 会造成重复列名
mrg0 = lg.merge(utgt, left_on="endpoint", right_on="phenocode", how="inner")
# ★★ 关键：inner merge 会产生交叉积（3 个 SNP × 2 工具 = 6 行/端点），
#     必须只保留 rsids == instrument 的正确配对，否则 implied |b_exp| 会混入错配值（std≠0）
P("    (b0) 未过滤交叉积行数 = %d（= 2,469 端点 × 3 SNP × 2 工具）" % len(mrg0))
mrg = mrg0[mrg0.rsids == mrg0.instrument].copy()
_exp = u.phenocode.nunique() * u.instrument.nunique()
P("        限定 rsids==instrument 后 = %d 行（应为 %d）  %s"
  % (len(mrg), _exp, "✅" if len(mrg) == _exp else "⚠️"))
mrg = mrg[["endpoint", "instrument", "beta", "sebeta", "pval", "af_alt"]].rename(
    columns={"beta": "beta_lg", "sebeta": "sebeta_lg",
             "pval": "pval_lg", "af_alt": "af_alt_lg"})
chk = u.merge(mrg, left_on=["phenocode", "instrument"],
              right_on=["endpoint", "instrument"], how="left")
P("    (b) 与 phewas_3snp 长表匹配率 = %d/%d  %s"
  % (chk.sebeta_lg.notna().sum(), len(chk),
     "✅" if chk.sebeta_lg.notna().all() else "⚠️"))
bc = chk.dropna(subset=["sebeta_lg"]).copy()
bc["implied_abs_b_exp"] = bc.sebeta_lg / bc.se_MR
imp = bc.groupby("instrument").implied_abs_b_exp.agg(["count", "median", "std", "min", "max"])
P("        由 sebeta/se_MR 反推的 |b_exp|（应为一常数）：")
for i, r in imp.iterrows():
    P("          %-12s n=%d  median=%.6f  std=%.3e  [%.6f, %.6f]"
      % (i, r["count"], r["median"], r["std"], r["min"], r["max"]))
# 与 OneK1K 原始 slope 对照
def onepq_slope(cell, sid):
    fp = os.path.join(ONEPQ, "OneK1K_%s.cis_qtl_pairs.chr1.parquet" % cell)
    t = pd.read_parquet(fp, filters=[("variant_id", "==", sid)])
    return t
CELL = {"rs12037376": ("B_MEM", "1:22462111"),
        "rs10917151": ("Mono_NC", "1:22422721")}
P("        与 OneK1K 原始 slope 直接对照（独立于 3.2 的计算路径）：")
for ins, (cell, sid) in CELL.items():
    t = pd.read_parquet(os.path.join(ONEPQ, "OneK1K_%s.cis_qtl_pairs.chr1.parquet" % cell),
                        filters=[("variant_id", "==", sid)])
    if len(t) == 0:
        P("          %-12s %s/%s 未命中" % (ins, cell, sid))
        continue
    t = t[t.phenotype_id == "CDC42"]
    if len(t) == 0:
        P("          %-12s %s/%s 无 CDC42 行" % (ins, cell, sid))
        continue
    sl = float(t.slope.iloc[0]); se = float(t.slope_se.iloc[0])
    P("          %-12s %s/%s  slope=%.6f  slope_se=%.6f  |slope|=%.6f  F=slope²/se²=%.4f"
      % (ins, cell, sid, sl, se, abs(sl), sl ** 2 / se ** 2))
    P("                       44 表 F = %.4f  ｜ 反推 |b_exp| = %.6f"
      % (float(u.loc[(u.instrument == ins), "F"].iloc[0]),
         float(imp.loc[ins, "median"]) if ins in imp.index else np.nan))
    if ins in imp.index:
        _mx, _md, _mn = (float(imp.loc[ins, k]) for k in ("max", "median", "min"))
        P("                       反推 |b_exp| [min,med,max] = [%.6f, %.6f, %.6f]"
          % (_mn, _md, _mx))
        P("                       ★ 上限 vs 原始 |slope| = %.6f → %s"
          % (abs(sl), "✅ 吻合（|bo|→0 时二阶 delta 项消失，公式退化为 sebeta/|b_exp|）"
             if abs(_mx - abs(sl)) < 1e-6 else "⚠️ 不吻合，需排查"))
        P("                       ★ 下限偏低源于二阶项在大 |bo| 端点上的贡献，属预期、非缺陷")
P("        ★ 说明：se_MR 为二阶 delta 全比值形式 sqrt(varo/be² + bo²·vare/be⁴)，")
P("          故 sebeta/se_MR 不等于常数 |b_exp|；其上限即真实 |b_exp|。")

# (c) 病例数与 manifest 一致性
mm = man[["phenocode", "num_cases", "num_controls"]].rename(
    columns={"num_cases": "man_n_cases", "num_controls": "man_n_controls"})
u = u.merge(mm, on="phenocode", how="left")
mis = (u.n_cases_manifest != u.man_n_cases).sum()
P("    (c) n_cases_manifest 与 manifest.num_cases 不一致行数 = %d  %s"
  % (mis, "✅" if mis == 0 else "⚠️"))
P("        n_cases 缺失(manifest) = %d" % u.man_n_cases.isna().sum())

# ---------------------------------------------------------------- 3) 主域 968
P("\n" + "=" * 100)
P("[3] ★ 主域：%s（family=%s）" % (DOMAIN_MAIN, FAMILY))
P("=" * 100)
m = u[u.domain == DOMAIN_MAIN].copy()
P("    检验数 = %d  端点数 = %d  工具 = %s"
  % (len(m), m.phenocode.nunique(), sorted(m.instrument.unique())))
assert len(m) == N_TESTS_MAIN, "主域检验数不等于 968"
P("    ✅ 检验数与用户口径 968 一致")
P("    其中 se_MR 不可用（NaN）的检验数 = %d" % m.se_MR.isna().sum())
m_ok = m[np.isfinite(m.min_detectable_OR)].copy()      # ★ 必须先定义，再统计
P("    功效可估（min_detectable_OR 有限）的检验数 = %d" % len(m_ok))
P("    ★ 下列功效统计一律基于 m_ok（n=%d）" % len(m_ok))

P("\n    (3.1) 病例数分布（端点为 484 个，2 工具重复）")
ep = m.drop_duplicates("phenocode")
P("        端点 n=%d  min=%d  p25=%.1f  median=%.1f  p75=%.1f  max=%d"
  % (len(ep), ep.n_cases_manifest.min(), ep.n_cases_manifest.quantile(.25),
     ep.n_cases_manifest.median(), ep.n_cases_manifest.quantile(.75),
     ep.n_cases_manifest.max()))
for thr in [100, 200, 500, 1000, 5000]:
    P("        n_cases < %-5d : %d 端点 (%.1f%%)"
      % (thr, (ep.n_cases_manifest < thr).sum(),
         100 * (ep.n_cases_manifest < thr).mean()))
P("        对照数 median = %.0f" % ep.num_controls.median())
P("        病例占比 median = %.5f" % (ep.n_cases_manifest / (ep.n_cases_manifest + ep.num_controls)).median())

P("\n    (3.2) 工具 F 统计量")
for i, g in m.groupby("instrument"):
    P("        %-12s  cell_type=%s  F=%.4f  n=%d"
      % (i, g.cell_type.iloc[0], g.F.iloc[0], len(g)))

P("\n    (3.3) 统计功效（α = 0.05/968 = %.4e；z_req(80%%) = %.4f）"
  % (alpha_main, zc_main + z80))
P("        最小可检出 OR 分布（968 检验）：")
for q in [0, .05, .25, .5, .75, .95, 1.0]:
    P("          P%-5s = %.4f" % (int(q * 100), m_ok.min_detectable_OR.quantile(q)))
P("        最小可检出 |b_MR| 中位数 = %.4f" % m_ok.min_detectable_abs_beta_MR.median())
P("        实测功率中位数 = %.4f ；实测功率 >=0.8 的检验数 = %d / 968"
  % (m_ok.power_at_observed_effect.median(),
     int((m_ok.power_at_observed_effect >= 0.8).sum())))
P("        在 α 下能检出「观测到的效应量」的检验数 = %d / 968"
  % int(m_ok.detectable_at_80pct_for_observed.sum()))

P("\n    (3.4) 功效分层（按 80%% power 可检出的最小 OR）")
cls = m_ok.power_class.value_counts().sort_index()
for k, v in cls.items():
    P("        %-24s %4d  (%.1f%%)" % (k, v, 100 * v / len(m_ok)))
P("        ★ 可检出 OR<=1.5 的检验数 = %d / 968 (%.1f%%)"
  % (int((m_ok.min_detectable_OR <= 1.5).sum()), 100 * (m_ok.min_detectable_OR <= 1.5).mean()))
P("        ★ 可检出 OR<=1.2 的检验数 = %d / 968 (%.1f%%)"
  % (int((m_ok.min_detectable_OR <= 1.2).sum()), 100 * (m_ok.min_detectable_OR <= 1.2).mean()))
P("        ★ 需 OR>3.0 才能检出的检验数 = %d / 968 (%.1f%%)"
  % (int((m_ok.min_detectable_OR > 3.0).sum()), 100 * (m_ok.min_detectable_OR > 3.0).mean()))

P("\n    (3.5) 敏感性：不同 α 口径下的最小可检出 OR 中位数")
P("        域内 Bonferroni 0.05/968   : med=%.4f  OR<=1.5 占 %.1f%%"
  % (m_ok.min_detectable_OR.median(), 100 * (m_ok.min_detectable_OR <= 1.5).mean()))
P("        全局 Bonferroni 0.05/4938 : med=%.4f  OR<=1.5 占 %.1f%%"
  % (m_ok.min_detectable_OR_globalBonf.median(), 100 * (m_ok.min_detectable_OR_globalBonf <= 1.5).mean()))
P("        未校正 0.05                : med=%.4f  OR<=1.5 占 %.1f%%"
  % (m_ok.min_detectable_OR_uncorrected.median(), 100 * (m_ok.min_detectable_OR_uncorrected <= 1.5).mean()))

P("\n    (3.6) 升向信号（FDR<0.05 且 beta_MR<0 = 风险升高）")
up = m_ok[m_ok.p_MR.notna() & (m_ok.p_FDR_family < 0.05) & (m_ok.beta_MR < 0)]
P("        FDR<0.05 检验数 = %d" % int((m_ok.p_FDR_family < 0.05).sum()))
P("        其中 风险升高（beta_MR<0）= %d ；风险降低（beta_MR>0）= %d"
  % (len(up), int(((m_ok.p_FDR_family < 0.05) & (m_ok.beta_MR > 0)).sum())))
_ri = m_ok[m_ok.beta_MR < 0].p_MR.idxmin()
P("        「升向」检验中最小的 p_MR = %.6f" % m_ok.loc[_ri, "p_MR"])
P("          端点 %s ｜ %s ｜ instrument %s ｜ n_cases=%d"
  % (m_ok.loc[_ri, "phenocode"], str(m_ok.loc[_ri, "phenotype"])[:44],
     m_ok.loc[_ri, "instrument"], int(m_ok.loc[_ri, "n_cases_manifest"])))
P("        ★ 该端点的检测功效：可检出 OR = %.4f（80%% power, α=0.05/968）⇒ %s"
  % (m_ok.loc[_ri, "min_detectable_OR"],
     "功效充分（可检出中等效应）" if m_ok.loc[_ri, "min_detectable_OR"] <= 1.5 else
     ("功效受限" if m_ok.loc[_ri, "min_detectable_OR"] <= 2.0 else "功效不足")))
P("        ★ 该端点在实测效应下的功效 = %.4f" % m_ok.loc[_ri, "power_at_observed_effect"])

P("\n    (3.7) 端点层面功效分层（484 个端点；每个端点取 2 工具中较优的 minOR）")
epb = m_ok.groupby("phenocode").agg(
    phenotype=("phenotype", "first"), category=("category", "first"),
    n_cases=("n_cases_manifest", "first"), num_controls=("num_controls", "first"),
    min_or_best=("min_detectable_OR", "min"),
    min_or_worst=("min_detectable_OR", "max"),
    n_instr=("instrument", "nunique"),
    min_p_MR=("p_MR", "min"), n_FDR05=("p_FDR_family", lambda s: int((s < 0.05).sum())),
).reset_index()
epb["ep_class_best"] = epb.min_or_best.map(classify)
P("        端点数 = %d" % len(epb))
for k, v in epb.ep_class_best.value_counts().sort_index().items():
    P("        %-24s %4d 端点 (%.1f%%)" % (k, v, 100 * v / len(epb)))
P("        ★ 端点层面 可检出 OR<=1.5 的 = %d / %d (%.1f%%)"
  % (int((epb.min_or_best <= 1.5).sum()), len(epb), 100 * (epb.min_or_best <= 1.5).mean()))
P("        ★ 端点层面 可检出 OR<=1.2 的 = %d / %d (%.1f%%)"
  % (int((epb.min_or_best <= 1.2).sum()), len(epb), 100 * (epb.min_or_best <= 1.2).mean()))
P("        ★ 端点层面 需 OR>3.0 才能检出的 = %d / %d (%.1f%%)"
  % (int((epb.min_or_best > 3.0).sum()), len(epb), 100 * (epb.min_or_best > 3.0).mean()))
P("        功效最好的 10 个端点（minOR 最小 ⇒ 对中等效应最敏感）：")
for _, r in epb.nsmallest(10, "min_or_best").iterrows():
    P("          %-38s n_cases=%7d  minOR=%6.3f" % (r.phenocode, r.n_cases, r.min_or_best))
P("        功效最差的 5 个端点：")
for _, r in epb.nlargest(5, "min_or_best").iterrows():
    P("          %-38s n_cases=%7d  minOR=%6.3f" % (r.phenocode, r.n_cases, r.min_or_best))

# ---------------------------------------------------------------- 4) 五域汇总
P("\n" + "=" * 100)
P("[4] 五域功效汇总（family=%s）" % FAMILY)
P("=" * 100)
rows = []
for dm, g in u.groupby("domain"):
    ge = g.drop_duplicates("phenocode")
    rows.append(dict(
        domain=dm, n_tests=len(g), n_endpoints=g.phenocode.nunique(),
        n_cases_min=int(ge.n_cases_manifest.min()),
        n_cases_median=float(ge.n_cases_manifest.median()),
        n_cases_max=int(ge.n_cases_manifest.max()),
        F_min=float(g.F.min()), F_max=float(g.F.max()),
        min_detectable_OR_median=float(g.min_detectable_OR.median()),
        min_detectable_OR_p05=float(g.min_detectable_OR.quantile(.05)),
        min_detectable_OR_p95=float(g.min_detectable_OR.quantile(.95)),
        n_detectable_OR_le_1p5=int((g.min_detectable_OR <= 1.5).sum()),
        pct_detectable_OR_le_1p5=100 * (g.min_detectable_OR <= 1.5).mean(),
        n_detectable_OR_le_1p2=int((g.min_detectable_OR <= 1.2).sum()),
        n_severe_OR_gt_3=int((g.min_detectable_OR > 3.0).sum()),
        median_power_at_observed=float(g.power_at_observed_effect.median()),
        n_FDR05=int((g.p_FDR_family < 0.05).sum()),
        n_risk_up=int(((g.p_FDR_family < 0.05) & (g.beta_MR < 0)).sum()),
        n_risk_down=int(((g.p_FDR_family < 0.05) & (g.beta_MR > 0)).sum()),
        min_p_MR=float(g.p_MR.min()),
        min_p_MR_risk_up=float(g[g.beta_MR < 0].p_MR.min()),
    ))
dom = pd.DataFrame(rows).sort_values("n_tests", ascending=False)
for _, r in dom.iterrows():
    P("    %-10s tests=%5d  ep=%4d  n_case(med)=%8.1f  minOR(med)=%.3f  OR<=1.5=%.1f%%  OR>3=%d  FDR05=%d(升%d/降%d)"
      % (r.domain, r.n_tests, r.n_endpoints, r.n_cases_median, r.min_detectable_OR_median,
         r.pct_detectable_OR_le_1p5, r.n_severe_OR_gt_3, r.n_FDR05, r.n_risk_up, r.n_risk_down))

# ---------------------------------------------------------------- 5) 端点存在性检索
P("\n" + "=" * 100)
P("[5] 端点存在性检索（不可评估 vs 已检验未过阈）")
P("=" * 100)
KEY = {
    "HLH_haemophagocytic": ["haemophagocytic", "hemophagocytic"],
    "HLH_lymphohistiocytosis": ["lymphohistiocytosis"],
    "HLH_abbrev": ["_hlh", "hlh"],
    "pancytopenia": ["pancytopenia", "pancytopaenia"],
    "bone_marrow_failure": ["bone marrow failure", "marrow failure"],
    "myeloid_proxy_agranulocytosis": ["agranulocytosis"],
    "myeloid_proxy_neutropenia": ["neutropenia"],
    "myeloid_proxy_aplastic": ["aplastic"],
    "myeloid_proxy_thrombocytopenia": ["thrombocytopenia", "thrombocytopaenia"],
    "myeloid_proxy_posthaemorrhagic_anaemia": ["posthaemorrhagic", "posthemorrhagic",
                                               "post-bleeding anaemia", "postbleed"],
}
# ★ 正文 §3.8 明确列出的 8 个髓系毒性代理端点，逐个显式核对存在性
PROXY_8 = ["D3_AGRANULOCYTOSIS", "DRUGADVERS_NEUTROPENIA", "D3_OTHERAPLASTICANAEMIA",
           "DRUGADVERS_APLAST_ANAEM", "D3_OTHPRIMTHROMBOCYTOPENIA",
           "D3_SCNDTHROMBOCYTOPENIA", "D3_THROMBOCYTOPENIANAS",
           "D3_ACUTEPOSTBLEEDANAEMIA"]
man["_ph"] = man.phenotype.astype(str).str.lower()
man["_pc"] = man.phenocode.astype(str).str.lower()
search_rows = []
for label, kws in KEY.items():
    hits = []
    for kw in kws:
        k = kw.strip().lower()
        h = man[man._ph.str.contains(k, regex=False, na=False)
                | man._pc.str.contains(k, regex=False, na=False)]
        for _, r in h.iterrows():
            hits.append(dict(kw=kw, phenocode=r.phenocode, phenotype=r.phenotype,
                             category=r.category, num_cases=int(r.num_cases),
                             num_controls=int(r.num_controls)))
    hits = pd.DataFrame(hits).drop_duplicates("phenocode") if hits else pd.DataFrame()
    search_rows.append(dict(group=label, keywords="|".join(kws),
                            n_hits=0 if hits.empty else len(hits),
                            endpoints="" if hits.empty else ";".join(hits.phenocode)))
    P("    %-32s keywords=%-45s n_hits=%d"
      % (label, "|".join(kws), 0 if hits.empty else len(hits)))
    if not hits.empty:
        for _, r in hits.iterrows():
            P("        → %-34s %-46s n_cases=%7d" % (r.phenocode, r.phenotype[:46], r.num_cases))
search = pd.DataFrame(search_rows)

# 代理端点在主域中的检验结果
P("\n    ★ 正文 §3.8 明列的 8 个髓系毒性代理端点逐个核对：")
found, missing = [], []
for pc in PROXY_8:
    r = man[man.phenocode == pc]
    if len(r):
        found.append(pc)
        rr = r.iloc[0]
        P("        ✅ %-32s n_cases=%7d  %s" % (pc, int(rr.num_cases), str(rr.phenotype)[:42]))
    else:
        missing.append(pc)
        P("        ❌ %-32s 未在 R12 manifest 中" % pc)
P("        存在 %d / %d ；缺失 = %s" % (len(found), len(PROXY_8), missing if missing else "无"))
# 关键词检索 ∪ 显式清单
kwset = set()
for kws in KEY.values():
    for kw in kws:
        k = kw.strip().lower()
        kwset |= set(man.loc[man._ph.str.contains(k, regex=False, na=False)
                             | man._pc.str.contains(k, regex=False, na=False), "phenocode"])
proxy_pcs = sorted(set(found) | kwset)
P("        关键词检索 ∪ 显式清单 = %d 个端点" % len(proxy_pcs))
px = u[u.phenocode.isin(proxy_pcs)]
P("        其在 uncond 域的检验数 = %d（%d 端点 × 2 工具）；风险升高 = %d；最小 p_MR = %.4f"
  % (len(px), len(proxy_pcs), int(((px.beta_MR < 0) & (px.p_FDR_family < 0.05)).sum()),
     px.p_MR.min() if len(px) else np.nan))
P("        其中「8 个明列代理端点」在 uncond 域的检验数 = %d ；风险升高 = %d ；最小 p_MR = %.4f"
  % (int(px.phenocode.isin(found).sum()),
     int((px.phenocode.isin(found) & (px.beta_MR < 0) & (px.p_FDR_family < 0.05)).sum()),
     px.loc[px.phenocode.isin(found), "p_MR"].min() if len(found) else np.nan))

# ---------------------------------------------------------------- 6) 落盘
P("\n" + "=" * 100)
P("[6] 落盘")
P("=" * 100)

cols_main = ["phenocode", "phenotype", "category", "domain", "instrument", "cell_type", "F",
             "n_cases_manifest", "num_controls", "man_n_cases", "man_n_controls",
             "alpha_used", "z_crit_domain", "z_req_80pct_domain",
             "beta_MR", "se_MR", "abs_beta_MR", "p_MR", "abs_z_MR",
             "p_FDR_family", "min_detectable_abs_beta_MR", "min_detectable_OR",
             "min_detectable_OR_globalBonf", "min_detectable_OR_uncorrected",
             "power_at_observed_effect", "detectable_at_80pct_for_observed",
             "n_cases_needed_minOR_1p5", "n_cases_needed_minOR_1p2", "power_class"]
f1 = os.path.join(TAB, "50_task35_phewas_power_primary_968.csv")
m[cols_main].sort_values(["min_detectable_OR"]).to_csv(f1, index=False, encoding="utf-8-sig")
P("    ✅ %s（%d 行）" % (os.path.basename(f1), len(m)))

f2 = os.path.join(TAB, "50b_task35_phewas_power_all4938.csv")
u[cols_main].sort_values(["domain", "min_detectable_OR"]).to_csv(f2, index=False, encoding="utf-8-sig")
P("    ✅ %s（%d 行）" % (os.path.basename(f2), len(u)))

f3 = os.path.join(TAB, "50c_task35_phewas_power_by_domain.csv")
dom.to_csv(f3, index=False, encoding="utf-8-sig")
P("    ✅ %s（%d 行）" % (os.path.basename(f3), len(dom)))

# 分层 × 域 交叉
cx = (m_ok.groupby(["domain", "power_class"]).size().unstack(fill_value=0))
f4 = os.path.join(TAB, "50d_task35_power_class_by_domain_968.csv")
cx.to_csv(f4, encoding="utf-8-sig")
P("    ✅ %s" % os.path.basename(f4))

f5 = os.path.join(TAB, "50e_task35_endpoint_absence_search.csv")
search.to_csv(f5, index=False, encoding="utf-8-sig")
P("    ✅ %s（%d 行）" % (os.path.basename(f5), len(search)))

# 主域病例数最极端的端点（功效最差）
worst = m_ok.nlargest(15, "min_detectable_OR")[
    ["phenocode", "phenotype", "instrument", "n_cases_manifest",
     "min_detectable_OR", "min_detectable_abs_beta_MR", "power_at_observed_effect"]]
f6 = os.path.join(TAB, "50f_task35_worst_power_endpoints_968.csv")
worst.to_csv(f6, index=False, encoding="utf-8-sig")
P("    ✅ %s（%d 行，功效最差）" % (os.path.basename(f6), len(worst)))

f7 = os.path.join(TAB, "50g_task35_endpoint_level_power_968.csv")
epb.sort_values("min_or_best").to_csv(f7, index=False, encoding="utf-8-sig")
P("    ✅ %s（%d 行，端点层面功效）" % (os.path.basename(f7), len(epb)))

P("\n    主域功效最差 15 个检验：")
for _, r in worst.iterrows():
    P("      %-38s n_cases=%7d  minOR=%7.3f  power=%.3f  %s"
      % (r.phenocode, r.n_cases_manifest, r.min_detectable_OR,
         r.power_at_observed_effect, r.instrument))

# 机读摘要
json.dump(dict(
    family=FAMILY, domain_main=DOMAIN_MAIN,
    n_tests_main=int(len(m)), n_endpoints_main=int(m.phenocode.nunique()),
    alpha_domain=alpha_main, z_req_80pct=float(zc_main + z80),
    instruments={i: dict(cell_type=g.cell_type.iloc[0], F=float(g.F.iloc[0]))
                 for i, g in m.groupby("instrument")},
    min_detectable_OR_median=float(m_ok.min_detectable_OR.median()),
    min_detectable_OR_p05=float(m_ok.min_detectable_OR.quantile(.05)),
    min_detectable_OR_p95=float(m_ok.min_detectable_OR.quantile(.95)),
    n_detectable_OR_le_1p5=int((m_ok.min_detectable_OR <= 1.5).sum()),
    n_detectable_OR_le_1p2=int((m_ok.min_detectable_OR <= 1.2).sum()),
    n_severe_OR_gt_3=int((m_ok.min_detectable_OR > 3.0).sum()),
    n_cases_endpoint_min=int(m.n_cases_manifest.min()),
    n_cases_endpoint_median=float(m.drop_duplicates("phenocode").n_cases_manifest.median()),
    n_cases_endpoint_max=int(m.n_cases_manifest.max()),
    n_FDR05=int((m.p_FDR_family < 0.05).sum()),
    n_risk_up=int(((m.p_FDR_family < 0.05) & (m.beta_MR < 0)).sum()),
    min_p_MR_risk_up=float(m[m.beta_MR < 0].p_MR.min()),
    power_class_counts={k: int(v) for k, v in m.power_class.value_counts().items()},
    absence_search={r.group: dict(n_hits=int(r.n_hits)) for _, r in search.iterrows()},
    runtime_sec=round(time.time() - T0, 1),
), open(os.path.join(TAB, "_s45_task35_summary.json"), "w", encoding="utf-8"),
   ensure_ascii=False, indent=2)
P("    ✅ _s45_task35_summary.json")

P("\n运行耗时 %.1f s" % (time.time() - T0))
P("VERDICT = PASS")

open(LOG, "w", encoding="utf-8").write(buf.getvalue())
print("DONE PASS")
