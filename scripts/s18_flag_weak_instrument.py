# -*- coding: utf-8 -*-
"""s18 —— 复现表的「弱工具效应量」标记（Yan 2026-09-23 裁定）

背景
----
复现队列（1M-scBloodNL, N≈119）中，当工具变异在该队列的 eQTL **未被检出**
（`instrument_qc == "not_detected"`，即 scb_pval 不达名义显著）时，Wald ratio
`replication_b = beta_outcome / scb_beta` 的**分母趋零**：

  * 分母趋零 -> 比值的抽样分布重尾/无界，`replication_b` 失去效应量语义；
  * 二阶 delta SE 虽能放大区间（例：CDC42/B_IN 的 se=3912），但**点估计本身不可用**。

故按 Yan 裁定（2026-09-23 二次裁定：采用**广义口径**）：
  ★ 这些行**不得作为效应量报告**，一律标记 `INVALID_DUE_TO_WEAK_INSTRUMENT`。
  ★ Yan 明确点名 CDC42/B_IN 与 LINC00339/CD4T；本脚本按**缺陷类别**（而非个案）
    把该标记施加于**全部** `not_detected` 行——因为它们的 Wald 分母缺陷同源。
    点名行是此集合的子集。
  ★★ Yan 2026-09-23 追加要求：**不要收窄为点名的 2 个**，按统一客观标准
    （`instrument_qc == "not_detected"`）处理 → 实际命中 **4 个唯一检验**；
    且 `replication_effect_note` 必须**保留逐行差异说明**，不得一刀切。
    故本脚本按 F 量级把缺陷分为两档写入 note：
      · A 档｜严格数值无界：F < 1e-2（Wald 分母≈0，点估计=纯噪声伪影）
      · B 档｜有限但不可靠：F ≥ 1e-2，点估计数值有限但区间无信息价值
    A/B 只写在 note 文本里；`replication_effect_status` **仍是单一客观标记**，
    以保证「按统一标准处理，逻辑自证」，避免审稿人质疑主观挑选。

对 `weak_nominal` 行（工具名义可检出但未校正显著）：效应量可报，但**必须带工具质量限定语**，
标为 `REPORTABLE_WITH_INSTRUMENT_CAVEAT`。

另：Yan 要求复现计数按「唯一复现检验 12 个」口径。本脚本同时产出
`24b_replication_unique_tests.csv`（按 gene × cell_type_replication × rsid 折叠）。

用法：
  python s18_flag_weak_instrument.py
"""
import os
import pandas as pd

ROOT = r"D:\endometriosis_project\11_sc_eqtl_mr_project\tables"
SRC = os.path.join(ROOT, "24_replication_MR_UT.csv")
SRC_T1 = os.path.join(ROOT, "24_replication_MR_UT_tier1.csv")
BAK = os.path.join(ROOT, "24_replication_MR_UT_preflag_backup.csv")
OUT_UNIQ = os.path.join(ROOT, "24b_replication_unique_tests.csv")
LOG = r"D:\endometriosis_project\_s18_flag.log"

msgs = []


def log(s):
    msgs.append(str(s))


log("=== s18 弱工具效应量标记 :: %s ===" % SRC)

# 幂等：若备份存在，一律从**备份**（未添加标记列的原始表）重新派生，
# 避免反复在已标记表上叠加导致口径漂移。
if os.path.exists(BAK):
    src = BAK
    log("检测到备份，从原始表重新派生 -> %s" % BAK)
else:
    src = SRC
    log("无备份，首次运行，读入 -> %s" % SRC)

df = pd.read_csv(src, encoding="utf-8-sig", dtype={"chr": str})
log("读入 %d 行 × %d 列" % (df.shape[0], df.shape[1]))

if not os.path.exists(BAK):
    df.to_csv(BAK, index=False, encoding="utf-8-sig")
    log("原表已备份 -> %s" % BAK)
else:
    log("备份已存在，跳过备份（%s）" % BAK)

# ---- 1) 效应量可报性标记 ----
bad = df["instrument_qc"].astype(str).eq("not_detected")
df["effect_reportable"] = (~bad).map({True: "TRUE", False: "FALSE"})
df["replication_effect_status"] = [
    "INVALID_DUE_TO_WEAK_INSTRUMENT" if b else "REPORTABLE_WITH_INSTRUMENT_CAVEAT"
    for b in bad
]

# 逐行说明（便于审稿人/合作者复核）—— 必须保留 A/B 档差异，不得一刀切
UNBOUNDED_F_CUT = 1e-2  # F 小于此值视为 Wald 分母≈0，点估计严格无界

notes = []
for _, r in df.iterrows():
    if str(r["instrument_qc"]) == "not_detected":
        f = float(r["scb_F"])
        b = float(r["replication_b"])
        se = abs(float(r["replication_se"])) if pd.notna(r["replication_se"]) else float("nan")
        ratio = (se / abs(b)) if b else float("inf")
        if f < UNBOUNDED_F_CUT:
            # A 档：严格数值无界
            notes.append(
                "[A/严格数值无界] 工具(%s)在复现队列未检出(β=%.4g, P=%.3g, F=%.3g)："
                "Wald 分母 = scb_beta 趋零(F=%.3g)，比值**严格数值无界**；"
                "实际报出的 replication_b=%.4g 与 se=%.4g（se/|b|=%.3g）均为分母噪声伪影，"
                "任何量级都不可解读 —— 禁止作为效应量报告，仅留存原始数值备核"
                % (r["rsid"], r["scb_beta"], r["scb_pval"], f, f, b, se, ratio)
            )
        else:
            # B 档：数值有限但不可靠
            notes.append(
                "[B/有限但不可靠] 工具(%s)在复现队列未检出(β=%.4g, P=%.3g, F=%.3g)："
                "比值数值有限(replication_b=%.4g, se=%.4g, se/|b|=%.3g)，"
                "但点估计由未检出工具驱动、区间无信息价值 —— 禁止作为效应量报告，"
                "仅留存原始数值备核"
                % (r["rsid"], r["scb_beta"], r["scb_pval"], f, b, se, ratio)
            )
    else:
        notes.append(
            "工具(%s)在复现队列名义可检出(β=%.4g, P=%.3g, F=%.3g) -> 效应量可报, "
            "但工具质量偏弱, 须带 weak_nominal 限定语"
            % (r["rsid"], r["scb_beta"], r["scb_pval"], r["scb_F"])
        )
df["replication_effect_note"] = notes

df.to_csv(SRC, index=False, encoding="utf-8-sig")
log("已写回 %s（新增 effect_reportable / replication_effect_status / replication_effect_note）" % SRC)

# ---- 2) Tier1 单表同步 ----
if os.path.exists(SRC_T1):
    t1 = df[df["tier"].astype(str).str.upper().eq("TIER1")].copy()
    t1.to_csv(SRC_T1, index=False, encoding="utf-8-sig")
    log("Tier1 单表已按同一口径重写（%d 行）：%s" % (t1.shape[0], SRC_T1))
else:
    log("!! 未找到 Tier1 单表，跳过")

# ---- 3) 唯一复现检验折叠（12 个）----
key = ["gene", "cell_type_replication", "rsid"]
uniq = df.drop_duplicates(subset=key).copy()
# 记录被折叠进来的 OneK1K 精细亚型
grp = df.groupby(key)["cell_type_discovery"].apply(lambda s: ";".join(sorted(set(map(str, s)))))
uniq["n_discovery_pairs"] = df.groupby(key)["cell_type_discovery"].size().values
uniq["collapsed_discovery_cell_types"] = [grp.get(k, "") for k in
                                          zip(uniq["gene"], uniq["cell_type_replication"], uniq["rsid"])]
uniq = uniq.sort_values(["tier", "gene", "cell_type_replication", "rsid"])

cols = ["tier", "locus_id", "gene", "cell_type_replication", "rsid",
        "variant_id_grch37", "variant_id_grch38",
        "discovery_b", "discovery_p", "discovery_qval", "discovery_F",
        "replication_b", "replication_p", "direction_consistent",
        "replication_status", "instrument_qc",
        "effect_reportable", "replication_effect_status", "replication_effect_note",
        "scb_beta", "scb_se", "scb_F", "scb_pval", "scb_fdr",
        "n_discovery_pairs", "collapsed_discovery_cell_types",
        "allele_check", "rsid_consistent", "method"]
uniq = uniq[[c for c in cols if c in uniq.columns]]
uniq.to_csv(OUT_UNIQ, index=False, encoding="utf-8-sig")
log("唯一复现检验表 -> %s（%d 个检验，来自 %d 个 discovery 对）"
    % (OUT_UNIQ, uniq.shape[0], df.shape[0]))

# ---- 4) 核对计数 ----
log("")
log("=== 计数核对（Yan 要求：按 12 个唯一检验报）===")
log("行级（15 行）:")
for k, v in df["replication_status"].value_counts().items():
    log("   %-20s %d" % (k, v))
log("唯一检验级（%d 个）:" % uniq.shape[0])
for k, v in uniq["replication_status"].value_counts().items():
    log("   %-20s %d" % (k, v))

log("")
log("=== 被标记 INVALID_DUE_TO_WEAK_INSTRUMENT 的检验（不得作效应量）===")
log("    客观标准: instrument_qc == not_detected（广义口径，非点名）")
inv = uniq[uniq["replication_effect_status"].eq("INVALID_DUE_TO_WEAK_INSTRUMENT")]
for _, r in inv.iterrows():
    sev = "A/严格无界" if float(r["scb_F"]) < UNBOUNDED_F_CUT else "B/有限不可靠"
    log("   %-10s %-10s %-12s rep_b=%12.4g  rep_p=%.4g  scb_F=%.4g  [%s]  (%s)"
        % (r["gene"], r["cell_type_replication"], r["rsid"],
           r["replication_b"], r["replication_p"], r["scb_F"], sev, r["variant_id_grch37"]))
nA = 0
nB = 0
for _, r in inv.iterrows():
    if float(r["scb_F"]) < UNBOUNDED_F_CUT:
        nA += 1
    else:
        nB += 1
log("    >> 唯一检验标记数 = %d（A 档严格无界 %d / B 档有限不可靠 %d）" % (inv.shape[0], nA, nB))
log("    >> 行级标记数 = %d（含 OneK1K 精细亚型重复行；CD4_ET/CD4_NC 折叠为 1 个唯一检验）"
    % int(df["replication_effect_status"].eq("INVALID_DUE_TO_WEAK_INSTRUMENT").sum()))
log("    >> replication_effect_note 已保留 A/B 逐行差异说明，未作一刀切表述")

log("")
log("=== 可报但需限定语的检验 ===")
ok = uniq[uniq["replication_effect_status"].eq("REPORTABLE_WITH_INSTRUMENT_CAVEAT")]
for _, r in ok.iterrows():
    log("   %-10s %-10s %-12s rep_b=%8.4f  rep_p=%.4g  %s"
        % (r["gene"], r["cell_type_replication"], r["rsid"],
           r["replication_b"], r["replication_p"], r["replication_status"]))

log("")
log("★ 红线 17-① 复述：CDC42 与 LINC00339 在同一座（chr1）内**必须分别陈述**，")
log("  严禁合并为「chr1 座共定位阳性」一句话。二者在本表中始终各占独立行，未做任何合并。")
log("★ 共定位尚未运行；本脚本不触碰 coloc / IVW。")

with open(LOG, "w", encoding="utf-8") as fh:
    fh.write("\n".join(msgs) + "\n")

print("done")
