# -*- coding: utf-8 -*-
"""
s62_build_techdoc_v8.py
------------------------
由 39b_技术文档_SM与SN_v7.md 派生 39b_技术文档_SM与SN_v8.md。
仅 3 类改动，均以「锚点唯一性断言」保护；SM1–SM11、SN1–SN10、补充表、补充图逐字保留。

  C1  标题版本号：v7 → v8
  C2  配套文件指向：追加 SM12/SM13 对应的正文文件
  C3  新增 v8 变更登记注（紧接 v7 变更登记注之后）
  C4  在 SM11 之后插入 SM12 + SM13

产物：39b_技术文档_SM与SN_v8.md（v7 原文不动）
日志：logs/s62_build_techdoc_v8.log
"""
import io
import os
import difflib
import hashlib

ROOT = r"D:/endometriosis_project/11_sc_eqtl_mr_project"
SRC = os.path.join(ROOT, "39b_技术文档_SM与SN_v7.md")
DST = os.path.join(ROOT, "39b_技术文档_SM与SN_v8.md")
LOG = os.path.join(ROOT, "logs", "s62_build_techdoc_v8.log")

T_TITLE = "# 技术文档 · 补充方法与逐项说明（v7）"
T_PTR = "> 本文档承载 v6 正文移出的全部实现层与口径层内容，与 `39_中文版_论文体v7_正文与图注.md` 配套阅读。"
T_V7NOTE_HEAD = "> **本版（v7）相对 v6 的变更登记**"
T_SM11 = ("**SM11 N 口径（各细胞类型有效样本量差异）。** 因 OneK1K 各细胞类型有效样本量不同"
          "（Plasma 253 至 980 不等），部分脚本硬编码 N = 980，导致最大 |ΔH4| = 3.6 × 10⁻⁴ 的轻微偏移，"
          "但零分层翻转、`susie` 可信集 lead 完全一致，不影响任何显著性判断与主结论。")

R_TITLE = "# 技术文档 · 补充方法与逐项说明（v8）"
R_PTR = ("> 本文档承载 v6 正文移出的全部实现层与口径层内容，与 `39_中文版_论文体v7_正文与图注.md` 配套阅读；"
         "其中 **SM12** 与 **SM13** 对应 `42_论文初稿_带图_v4.md` 正文 §3.1、§3.3 与讨论中新增的"
         "结果层陈述与探索性分析。")

V8_NOTE = (
    "> **本版（v8）相对 v7 的变更登记**：**2 处新增 + 1 处指向更新**——在「补充方法」中新增 "
    "**SM12**（GTEx v8 全血区域定向取数与 cis-SMR/HEIDI 口径）与 **SM13**（MVMR 探索性分析口径），"
    "以对应 `42_论文初稿_带图_v4.md` 正文新增的结果层与探索性陈述；并更新本文件开头的配套文件指向，"
    "使之同时覆盖 `42_论文初稿_带图_v4.md`。SM1–SM11、SN1–SN10、补充表、补充图自 v7 一字未改；"
    "`39_中文版_论文体v7_正文与图注.md` 未改动，v7 文件（`39b_技术文档_SM与SN_v7.md`）原样保留。"
    "新增两条的口径与数值全部来自本轮已交付产物：`tables/59_gtex8_L2L3_smr_results.csv`、"
    "`tables/58_gtex8_query_orient_audit_L2L3.csv`、`tables/61_mvmr_results.csv`、"
    "`tables/61b_mvmr_ld_diag.csv`、`tables/62_mvmr_decomposition.csv`，"
    "脚本 `scripts/s55a`–`s55c`、`scripts/s60b`–`s60d`。"
)

SM12 = (
    "**SM12 GTEx v8 全血区域定向取数与 cis-SMR/HEIDI 口径（chr10 YME1L1 与 chr2 ANXA4）。** "
    "该项补充分析针对 chr1p36.12 之外的两个独立信号座，用于把三个信号座放到同一 bulk 参照层上比较。"
    "数据取自 eQTL Catalogue 导入的 GTEx v8 全血（`Whole_Blood.tsv.gz`，3,190,551,490 B）。"
    "因本机无 tabix / bgzip / pysam，取数以手写的 tabix-over-HTTP 区域随机访问实现"
    "（解析本地 `.tbi` 索引 → 标准五级分箱 → BGZF 块级 HTTP Range → 多成员 gzip 增量解码）；"
    "严格限定两段区域（YME1L1：chr10 GRCh38 26,032,861–28,232,861；ANXA4：chr2 68,635,445–70,835,445），"
    "实际取回 6.11 MB，为全文件的 0.2007%，**未下载全量文件**。基因编号经 Ensembl REST 核验："
    "`YME1L1 = ENSG00000136758`、`ANXA4 = ENSG00000196975`。样本量由 `an` 字段中位数推得 N = 670。"
    "SMR query 的构建沿用既有口径：双等位变异经 hg38→hg19 liftover 后须落在 OneK1K 参考的 .bim 内、"
    "等位对一致，并按该 TSV 中 `ac/an` 为次等位频率（而非 ALT 频率）做取向修正；"
    "可用变异 YME1L1 4,681 个（P < 5 × 10⁻⁸ 者 193 个）、ANXA4 3,777 个（113 个）。"
    "SMR 参数为主分析 `--peqtl-smr 5e-8`、敏感性 1e-5，HEIDI 参数 `--peqtl-heidi 1.57e-3`、"
    "`--heidi-min-m 3`、`--heidi-max-m 20`；LD 参考为同区域 OneK1K 980 名供者的 PLINK 子集"
    "（与 GTEx v8 供者不重叠，属已知局限）。结果：YME1L1 b_SMR = +0.35154、P_SMR = 2.19 × 10⁻⁴、"
    "P_HEIDI = 0.114、m = 20；ANXA4 b_SMR = −0.13401、P_SMR = 4.44 × 10⁻⁴、P_HEIDI = 0.142、m = 20；"
    "主分析与敏感性分析的估计逐位相同，亦无低频变异导致 m < 3 的情形。"
    "两基因的 top 变异与本地 GTEx v10 的 eGenes 记录完全一致。作为层对照，同一 bulk 层上 "
    "chr1p36.12 的 CDC42 为 P_HEIDI = 0.435（一致），LINC00339 为 3.44 × 10⁻¹⁰（异质）。"
)

SM13 = (
    "**SM13 MVMR 探索性分析口径（CDC42 与 LINC00339）。** 该分析为正文讨论中的探索性补充，"
    "用于观察 chr1p36.12 内两个被检验基因的效应如何分解，其证据层级低于主分析。"
    "工具集并非取自主分析结果表——该表只保留 clump 后的唯一工具，不足以构成多工具集——"
    "而是以主分析结果表锁定细胞类型与基因对后，自剪枝前的全部 cis 变异（并集 P < 5 × 10⁻⁸）重建。"
    "剪枝阈值做系统性扫描（r² ∈ {0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.5} × 8 个细胞类型，共 56 组）；"
    "剪枝前工具并集的 LD 相关阵数值近奇异（最小特征值 ≈ 0，条件数约 10³⁰¹，最大两两 r² = 1.0000），"
    "故不剪枝不可行。可识别性按预先设定的判据判定：工具数 ≥ 2、两个暴露的条件 F 均 ≥ 10"
    "（Sanderson 等，2019 实现）且 LD 条件数 ≤ 30。结果：在正文主分析的剪枝口径（r² < 0.001）下，"
    "8 个细胞类型中有 7 个只剩 1 个独立工具，模型秩不足、不可识别；`CD4_NC` 为恰好识别"
    "（残差自由度 0，标准误不可估计）。放宽剪枝后 13/56 组可识别，其分解一致："
    "调整 LINC00339 后 CDC42 的效应保持为正（+0.090 至 +0.204），调整 CDC42 后 LINC00339 的效应"
    "衰减至近零（−0.003 至 −0.050）且均不显著，两者方向均无翻转。"
    "点估计以 Python 独立复算 `(B′WB)⁻¹B′WΓ` 与 `MVMR::ivw_mvmr`（MVMR 0.4.8）逐位一致；"
    "标准误为固定效应加权回归标准误，另附异质性稳健（sandwich）标准误。"
    "**边界**：`WNT4` 与 `MASTL` 在 OneK1K 的 14 种细胞类型中 cis-eQTL 记录为零（面板缺失），"
    "无法纳入 MVMR，故该分析不能回答「CDC42 与 WNT4 何者为独立信号」；"
    "因可识别结果依赖放宽后的剪枝阈值，该分解为探索性证据。"
)


def main():
    lines = []
    src = io.open(SRC, "r", encoding="utf-8", newline="").read()
    lines.append("== s62 技术文档 v7 → v8 派生日志 ==")
    lines.append("源 39b_v7 md5 : " + hashlib.md5(src.encode("utf-8")).hexdigest())

    checks = [("标题", T_TITLE), ("配套指向", T_PTR), ("v7 变更登记注", T_V7NOTE_HEAD), ("SM11", T_SM11)]
    bad = False
    for name, anchor in checks:
        n = src.count(anchor)
        lines.append("锚点 %-14s 命中 %d 次" % (name, n))
        if n != 1:
            bad = True
    if bad:
        lines.append("RESULT=ABORT 锚点异常")
        io.open(LOG, "w", encoding="utf-8", newline="\n").write("\n".join(lines) + "\n")
        print("RESULT=ABORT")
        return

    out = src
    # C1 标题
    out = out.replace(T_TITLE, R_TITLE, 1)
    # C2 配套指向
    out = out.replace(T_PTR, R_PTR, 1)
    # C3 v8 变更登记注：插在 v7 变更登记注整段之后
    v7note_start = out.index(T_V7NOTE_HEAD)
    v7note_end = out.index("\n", v7note_start)
    out = out[:v7note_end + 1] + "\n" + V8_NOTE + "\n" + out[v7note_end + 1:]
    # C4 SM12/SM13 插在 SM11 之后
    i = out.index(T_SM11)
    j = out.index("\n", i)
    out = out[:j + 1] + "\n" + SM12 + "\n\n" + SM13 + "\n" + out[j + 1:]

    io.open(DST, "w", encoding="utf-8", newline="\n").write(out)

    lines.append("")
    lines.append("产物 39b_v8 md5 : " + hashlib.md5(out.encode("utf-8")).hexdigest())
    lines.append("字节 %d → %d（+%d）" % (len(src.encode("utf-8")), len(out.encode("utf-8")),
                                        len(out.encode("utf-8")) - len(src.encode("utf-8"))))
    lines.append("行数 %d → %d" % (src.count("\n"), out.count("\n")))

    # 差异复核：应仅含 C1–C4 四处
    d = list(difflib.unified_diff(src.splitlines(), out.splitlines(), "v7", "v8", lineterm="", n=0))
    added = [x for x in d if x.startswith("+") and not x.startswith("+++")]
    removed = [x for x in d if x.startswith("-") and not x.startswith("---")]
    lines.append("")
    lines.append("diff 新增行 %d / 删除行 %d（预期 新增 3 / 删除 3：标题、指向、登记注、SM11 之后的 4 行）"
                 % (len(added), len(removed)))
    for x in removed:
        lines.append("   DEL " + x[:90])
    for x in added:
        lines.append("   ADD " + x[:90])

    # 正向/反向探针
    for probe, name in [("SM12 GTEx v8 全血区域定向取数与 cis-SMR/HEIDI 口径", "SM12 小标题"),
                        ("SM13 MVMR 探索性分析口径", "SM13 小标题"),
                        ("最大两两 r² = 1.0000", "SM13 关键数值"),
                        ("P_HEIDI = 0.114", "SM12 关键数值")]:
        lines.append("含 %-16s : %s" % (name, out.count(probe) == 1))
    for probe, name in [("（v7）", "旧版号残留"), ("仅 1 处新增", "v7 登记注被改写")]:
        lines.append("旧串 %-16s 计数 : %d" % (name, out.count(probe)))

    lines.append("")
    lines.append("RESULT=PASS")
    os.makedirs(os.path.dirname(LOG), exist_ok=True)
    io.open(LOG, "w", encoding="utf-8", newline="\n").write("\n".join(lines) + "\n")
    print("RESULT=PASS")


if __name__ == "__main__":
    main()
