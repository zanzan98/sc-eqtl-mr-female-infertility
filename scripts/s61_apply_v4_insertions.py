# -*- coding: utf-8 -*-
"""
s61_apply_v4_insertions.py
--------------------------
按用户裁定 1/3/4 在 41_论文初稿_带图_v3.md 的副本（42_论文初稿_带图_v4.md）上
落实 4 处插入。每处插入均先断言锚点在文件中唯一出现（count == 1），
否则中止写入并报错，避免误插。

插入清单（用户裁定 → 41 稿落位映射）：
  I1  裁定"§3.1 末尾加三位点比较说明"     → 41 稿 §3.1 末尾（【插入图 2】之前）
  I2  裁定"§3.9 末尾加 GTEx 补充陈述"     → 41 稿 §3.3 末尾（【插入图 4】之前；
                                           41 稿 §3.3 即 SMR/HEIDI/GTEx 段所在）
  I3  裁定"§4.1 末尾加方法学讨论"         → 41 稿 §4「方法学贡献」段末尾
  I4  裁定"§4.3 末尾加 MVMR 探索性一段"   → 41 稿 §4「WNT4 归属」段之后
                                           （对应 39_v7 §4.2 → §4.3 的相邻序）
输出：
  - 覆盖写 42_论文初稿_带图_v4.md
  - 日志 logs/s61_apply_v4_insertions.log（UTF-8）
"""
import io
import os
import hashlib

ROOT = r"D:/endometriosis_project/11_sc_eqtl_mr_project"
TARGET = os.path.join(ROOT, "42_论文初稿_带图_v4.md")
BASELINE = os.path.join(ROOT, "41_论文初稿_带图_v3.md")
LOG = os.path.join(ROOT, "logs", "s61_apply_v4_insertions.log")

# ---------------------------------------------------------------- 锚点与插入内容
A1 = "记忆 B 细胞的复核状态因而有待更大规模的单细胞队列确认。\n\n**【插入图 2】**"
I1 = (
    "三个信号座在结构上并不对称。另两个座——chr10 的 YME1L1 与 chr2 的 ANXA4——各只承载一个检验对，"
    "而 chr1p36.12 承载 13 个，并且是唯一在同一区域内同时出现两个基因、且二者效应方向相反的座。"
    "把三个座放到同一参照层上比较可以更清楚地看出这一点：在 GTEx v8 全血中，YME1L1（P_HEIDI = 0.114）、"
    "ANXA4（P_HEIDI = 0.142）与 CDC42（P_HEIDI = 0.435）的 cis-SMR 均通过 HEIDI 的一致性判读，"
    "而同一层上的 LINC00339 判为异质（P_HEIDI = 3.4 × 10⁻¹⁰）。"
    "因此 chr1p36.12 的特殊性不止在于它同时携带两个方向相反的信号，还在于这种判读分层只出现在该区域内部——"
    "这构成本研究把该区域单独展开的理由。"
)

A2 = ("综合精细定位、条件分析与多组织 eQTL 注释，归属证据一致指向 WNT4，"
      "而 CDC42 更可能是这一区域效应在特定免疫细胞中的“标签”。\n\n**【插入图 4】**")
I2 = (
    "在 GTEx v8 全血（n = 670）中，YME1L1 与 ANXA4 均通过 HEIDI，而 LINC00339 判为异质；"
    "在 OneK1K 单细胞层，YME1L1 的 HEIDI 判定为临界（P = 0.049）。"
    "三位点在两层暴露上的判读分层，进一步支持 chr1p36.12 在三个信号座中的特殊性。"
)

A3 = ("整体 eQTL 与单细胞 eQTL 不是替代关系，而是互补关系：整体层确立可检出性，"
      "单细胞层确立细胞归属，两层共同解释了 WNT4 在全血可检出却在单细胞血液面板中缺席这一格局。")
I3 = (
    "三个信号座的并列比较为这一分层提供了内部校准：同一套流程在 chr10 的 YME1L1 与 chr2 的 ANXA4 上"
    "给出跨层一致的判读——两个基因在 GTEx v8 全血中均通过 HEIDI——而在 chr1p36.12 上给出座内分层的判读"
    "（CDC42 为一致、LINC00339 为异质）。"
    "判读差异因此更可能来自位点结构本身（单一强 LD 区块、多基因共享同一单倍型），而不是流程的系统性偏倚。"
    "这也说明，在强 LD 区域报告区域层级结论、在单基因信号座报告基因层级结论，"
    "是同一条推断链上的两种分辨率，而不是两套标准。"
)

A4 = ("这一解读与条件共定位的结果一致：以 WNT4 的 lead 为条件后，CDC42 的共享信号从 0.996 降至 0.058，"
      "而改用位点自身的 lead 作为条件时信号保持不变。"
      "共享成分由区域 lead 标记的单倍型承载，而非被检验基因自身的 eQTL。")
I4 = (
    "作为对该归属边界的补充探索，我们把 CDC42 与 LINC00339 的 cis-eQTL 工具同时纳入"
    "多变量孟德尔随机化（MVMR），以观察两者的效应如何分解。"
    "在本文主分析的剪枝口径（r² < 0.001）下，8 种细胞类型中有 7 种只剩一个独立工具，模型秩不足、无法识别；"
    "其原因是 chr1p36.12 构成单一强 LD 区块，工具并集内最大两两 r² 达 1.0000。"
    "将剪枝阈值放宽到 r² = 0.05–0.3 后，56 个「细胞类型 × 阈值」组合中有 13 个可识别，其分解结果高度一致："
    "调整 LINC00339 后 CDC42 的效应保持为正（+0.090 至 +0.204），"
    "而调整 CDC42 后 LINC00339 的效应衰减至近零（−0.003 至 −0.050）且均不显著。"
    "这一结果支持「在 CDC42 与 LINC00339 之间，前者是更接近独立信号载体的一方」，"
    "但它不构成基因归属的判决：WNT4 与 MASTL 在 OneK1K 的 14 种细胞类型中 cis-eQTL 记录为零（面板缺失），"
    "无法纳入 MVMR，因此该分析无法回答「CDC42 与 WNT4 何者为独立信号」。"
    "由于可识别结果依赖放宽后的剪枝阈值，上述分解应视为探索性证据，与主分析的证据层级不同。"
)

# (标签, 锚点, 替换文本, 位置说明)
JOBS = [
    ("I1", A1, "记忆 B 细胞的复核状态因而有待更大规模的单细胞队列确认。\n\n" + I1 + "\n\n**【插入图 2】**",
     "41 稿 §3.1 末尾（【插入图 2】之前）"),
    ("I2", A2, ("综合精细定位、条件分析与多组织 eQTL 注释，归属证据一致指向 WNT4，"
                "而 CDC42 更可能是这一区域效应在特定免疫细胞中的“标签”。\n\n" + I2 + "\n\n**【插入图 4】**"),
     "41 稿 §3.3 末尾（【插入图 4】之前）"),
    ("I3", A3, A3 + "\n\n" + I3, "41 稿 §4「方法学贡献」段末尾"),
    ("I4", A4, A4 + "\n\n" + I4, "41 稿 §4「WNT4 归属」段之后"),
]


def main():
    lines = []
    src = io.open(TARGET, "r", encoding="utf-8", newline="").read()
    base = io.open(BASELINE, "r", encoding="utf-8", newline="").read()

    lines.append("== s61 插入执行日志 ==")
    lines.append("底本 41_v3 md5 : " + hashlib.md5(base.encode("utf-8")).hexdigest())
    lines.append("改写前 42_v4 md5: " + hashlib.md5(src.encode("utf-8")).hexdigest())
    lines.append("")

    ok = True
    for tag, anchor, repl, where in JOBS:
        n = src.count(anchor)
        lines.append("%s  锚点命中 %d 次  → %s" % (tag, n, where))
        if n != 1:
            lines.append("     !! 锚点数 != 1，中止，不写入")
            ok = False
    if not ok:
        lines.append("")
        lines.append("RESULT=ABORT")
        os.makedirs(os.path.dirname(LOG), exist_ok=True)
        io.open(LOG, "w", encoding="utf-8", newline="\n").write("\n".join(lines) + "\n")
        print("RESULT=ABORT")
        return

    out = src
    for tag, anchor, repl, where in JOBS:
        out = out.replace(anchor, repl, 1)

    io.open(TARGET, "w", encoding="utf-8", newline="\n").write(out)

    lines.append("")
    lines.append("改写后 42_v4 md5: " + hashlib.md5(out.encode("utf-8")).hexdigest())
    lines.append("字节数 %d → %d（+%d）" % (len(src.encode("utf-8")), len(out.encode("utf-8")),
                                          len(out.encode("utf-8")) - len(src.encode("utf-8"))))
    lines.append("段落数 %d → %d" % (src.count("\n\n") + 1, out.count("\n\n") + 1))
    lines.append("")
    # 复核：4 处新段确实落位，且未破坏图注标记
    for probe, name in [(I1[:20], "I1"), (I2[:20], "I2"), (I3[:20], "I3"), (I4[:20], "I4")]:
        lines.append("含 %s 首句: %s" % (name, out.count(probe) == 1))
    for marker in ["**【插入图 1】**", "**【插入图 2】**", "**【插入图 3】**",
                   "**【插入图 4】**", "**【插入图 5】**"]:
        lines.append("图占位 %s : %d" % (marker, out.count(marker)))
    lines.append("")
    lines.append("RESULT=PASS 4 insertions applied")
    os.makedirs(os.path.dirname(LOG), exist_ok=True)
    io.open(LOG, "w", encoding="utf-8", newline="\n").write("\n".join(lines) + "\n")
    print("RESULT=PASS")


if __name__ == "__main__":
    main()
