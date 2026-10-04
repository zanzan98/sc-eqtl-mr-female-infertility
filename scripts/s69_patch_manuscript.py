# -*- coding: utf-8 -*-
"""s69_patch_manuscript.py — 按用户裁定（2026-09-29）修改当前主线稿件  [v2: 幂等 + 行尾自适应]

三类改动，全部为**字节级精确替换**（保留各行尾与 BOM 原样）：
  A. Fig. 4(a) 图注参考面板误述：1000 Genomes European reference -> OneK1K（用户裁定第 3 项）
  B. Results 末尾加一句集合式引用（用户裁定第 2 项）
  C. 补充图注由 S1-S2 扩为 S1-S4，并把图 S1 展开为 (a-d)（用户裁定第 1 项）

v2 修正（2026-09-29 首次运行失败）：
  * 39_英文正文_论文体v6.md 为 CRLF，锚点改为按文件自适应；
  * 加入幂等守卫：old 计数为 0 且 new 已存在 => 视为已应用，跳过（不再断言失败）。

范围 = 当前主线 6 个文件。历史草稿（含 36_* 已被 P3-3 交付清单登记字节数者）不回填，列入待裁清单。
"""
import io
import os
import shutil

ROOT = r"D:/endometriosis_project/11_sc_eqtl_mr_project"
BAK = r"D:/_transfer_logs/_ms_backup_20260929_s1s4"

F_EN = ["38_英文正文_论文体重写_v5.md", "39_英文正文_论文体v6.md"]
F_ZH = ["38_中文版_论文体重写_v5_正文与图注.md", "39_中文版_论文体v7_正文与图注.md",
        "43_论文初稿_带图_v5.md"]
F_LEG_EN = ["38_英文正文_论文体重写_v5.md"]
F_LEG_ZH = ["38_中文版_论文体重写_v5_正文与图注.md", "39b_技术文档_SM与SN_v8.md"]
ALL6 = ["38_英文正文_论文体重写_v5.md", "39_英文正文_论文体v6.md",
        "38_中文版_论文体重写_v5_正文与图注.md", "39_中文版_论文体v7_正文与图注.md",
        "43_论文初稿_带图_v5.md", "39b_技术文档_SM与SN_v8.md"]

CAP_EN_OLD = b"(r\xc2\xb2, 1000 Genomes European reference)"
CAP_EN_NEW = b"(r\xc2\xb2, OneK1K reference panel, 980 donors)"
CAP_ZH_OLD = u"（r²，1000 Genomes 欧洲参考）".encode("utf-8")
CAP_ZH_NEW = u"（r²，OneK1K 参考面板，980 名供者）".encode("utf-8")

CITE_EN = ("Additional quality-control and stratified views are provided in "
           "Supplementary Figs. S1\u2013S4.").encode("utf-8")
CITE_ZH = u"质量控制与分层视图另见补充图 S1–S4。".encode("utf-8")
ANCH_EN = b"## 4. Discussion"
ANCH_ZH = u"## 4 讨论".encode("utf-8")

LEG_EN_OLD = (b"Fig. S1 Instrument power and analytical gating by cell type (instrument counts, "
              b"nominally and FDR-significant test counts, and power flags). Fig. S2 Comparison "
              b"with the sensitivity outcome `GCST90483469` (15/15 concordant).")

LEG_EN_NEW = """Fig. S1 Instrument power and analytical gating by cell type. **(a)** Number of cis-eQTL instruments per cell type, coloured by the power flag (`signal`, cell types with at least one test passing FDR < 0.05; `low power`, fewer than 200 instruments and no test passing FDR < 0.05); the low-power cell types (`DC`, `NK_R`, `Plasma`, `CD4_SOX4`) are flagged as not interpretable as null results. **(b)** Number of nominally significant (P < 0.05) tests. **(c)** Number of tests passing FDR < 0.05. **(d)** -log10 of the smallest q value per cell type; dashed line, FDR = 0.05. Fig. S2 Comparison with the sensitivity outcome `GCST90483469`. **(a)** Main-outcome versus sensitivity-outcome MR effect sizes for the 15 pairs significant in both (dashed line, identity; Pearson r = 0.9999, slope = 1.001); **(b)** paired forest plot of the same 15 pairs (filled, main outcome; open, sensitivity outcome; bars, 95% CI). All 15 pairs are direction-consistent and FDR < 0.05 in both outcomes. Fig. S3 Cell-type-stratified Manhattan plots of the 8,612 discovery MR tests. **(a-n)** One panel per cell type (14 in total); x axis, chromosome with position normalised within each chromosome; y axis, -log10(P). Dashed line, Bonferroni threshold (p = 0.05/8,612), identical to Fig. 2a; diamonds, pairs passing FDR < 0.05. Fig. S4 Distribution of discovery MR effect sizes. MR effect size (beta per SD of expression, log-odds of female infertility) against -log10(P) for all 8,612 tests; points, P >= 0.05 and nominal P < 0.05; diamonds, FDR < 0.05 (n = 15), annotated once per gene with the corresponding cell type. Dashed line, Bonferroni threshold (p = 0.05/8,612); dotted line, Benjamini-Hochberg threshold (p = 6.43e-05) corresponding to FDR < 0.05.""".encode("utf-8")

LEG_ZH_OLD = u"图 S1 各细胞类型的工具变量效力与分析门控（各细胞类型的工具数、名义与 FDR 显著检验数，以及效力标记）。图 S2 与敏感性结局 `GCST90483469` 的对照（15/15 一致）。".encode("utf-8")

LEG_ZH_NEW = u"""图 S1 各细胞类型的工具变量效力与分析门控。**（a）** 各细胞类型的 cis-eQTL 工具数，按效力标记区分（`signal` 指该细胞类型至少有 1 个检验通过 FDR < 0.05；`low power` 指工具数 < 200 且无检验通过 FDR < 0.05）；低功效细胞类型（`DC`、`NK_R`、`Plasma`、`CD4_SOX4`）已标注为不可读作阴性。**（b）** 名义显著（P < 0.05）检验数。**（c）** 通过 FDR < 0.05 的检验数。**（d）** 各细胞类型最小 q 值的 -log10；虚线为 FDR = 0.05。图 S2 与敏感性结局 `GCST90483469` 的对照。**（a）** 15 对在两种结局下均显著的 MR 效应量一致性（虚线为等值线；Pearson r = 0.9999，斜率 1.001）；**（b）** 同一 15 对的主/敏感双臂森林图（实心为主结局，空心为敏感结局；误差棒为 95% CI）。15 对方向全部一致，且两种结局下 FDR 均 < 0.05。图 S3 各细胞类型分层的 Discovery 曼哈顿图（8,612 个检验）。**（a–n）** 每个细胞类型一个面板（共 14 个）；横轴为染色体（染色体内位置归一化），纵轴为 -log10(P)；虚线为 Bonferroni 阈值（p = 0.05/8,612，与图 2a 一致），菱形为通过 FDR < 0.05 的对。图 S4 Discovery MR 效应量分布。8,612 个检验的 MR 效应量（每 SD 表达量的 beta，女性不孕症 log-odds）对 -log10(P)；显著性分三级：未达名义显著（P ≥ 0.05）、名义显著（P < 0.05），以及以菱形标记的 FDR < 0.05（n = 15），并逐基因标注一次基因与对应细胞类型；虚线为 Bonferroni 阈值（p = 0.05/8,612），点线为对应 FDR < 0.05 的 Benjamini-Hochberg 阈值（p = 6.43e-05）。""".encode("utf-8")


def read(p):
    return open(p, "rb").read()


def nl_of(d):
    """返回该文件的主行尾（比对 CRLF 与 LF 计数）。"""
    c_crlf = d.count(b"\r\n")
    c_all = d.count(b"\n")
    return b"\r\n" if (c_crlf > 0 and c_crlf == c_all) else b"\n"


os.makedirs(BAK, exist_ok=True)
LOG = []


def backup(name):
    src = os.path.join(ROOT, name)
    dst = os.path.join(BAK, name)
    if not os.path.exists(dst):
        shutil.copy2(src, dst)


def patch(name, pairs):
    p = os.path.join(ROOT, name)
    data = read(p)
    for old, new, must, note in pairs:
        n = data.count(old)
        if n == 0:
            k = data.count(new)
            assert k > 0, "%s: NEITHER old nor new present (%s)" % (name, note)
            LOG.append("    [skip: already applied] %-10s %s" % (note, name))
            continue
        assert n == must, "%s: %r count=%d (expect %d)" % (name, note, n, must)
        data = data.replace(old, new)
        LOG.append("    [applied]                %-10s %s" % (note, name))
    with io.open(p, "wb") as fh:
        fh.write(data)
    return data


# ---------------- 预扫描 ----------------
LOG.append("[pre-scan] 目标 6 文件现状")
for f in ALL6:
    d = read(os.path.join(ROOT, f))
    LOG.append("  %-44s bytes=%-7d nl=%-4s capEN=%d capZH=%d legEN=%d legZH=%d anchEN=%d anchZH=%d cite=%d"
               % (f, len(d), ("CRLF" if nl_of(d) == b"\r\n" else "LF"),
                  d.count(CAP_EN_OLD), d.count(CAP_ZH_OLD),
                  d.count(LEG_EN_OLD), d.count(LEG_ZH_OLD),
                  d.count(ANCH_EN), d.count(ANCH_ZH), d.count(CITE_EN) + d.count(CITE_ZH)))

TOUCH = sorted(set(F_EN + F_ZH + F_LEG_EN + F_LEG_ZH))
for t in TOUCH:
    backup(t)
LOG.append("[backup] -> %s  (%d files)" % (BAK, len(TOUCH)))

# ---------------- A. 图注误述修正 ----------------
LOG.append("[A] caption source correction (1000 Genomes -> OneK1K)")
for f in F_EN:
    patch(f, [(CAP_EN_OLD, CAP_EN_NEW, 1, "caption-EN")])
for f in F_ZH:
    patch(f, [(CAP_ZH_OLD, CAP_ZH_NEW, 1, "caption-ZH")])

# ---------------- B. Results 末尾集合式引用 ----------------
# ★ v3 修正（2026-10-03）：本步原**不幂等** —— anchor 在插入后仍存在且 new 以 anchor 结尾，
#   故 patch() 的 `n == 0` 守卫失效，重跑会二次插入、静默重复一句。改为先查 CITE 是否已在。
LOG.append("[B] Results closing citation (S1-S4)")
for f in F_EN:
    d = read(os.path.join(ROOT, f))
    if CITE_EN in d:
        LOG.append("    [skip: already applied] %-10s %s" % ("cite-EN", f))
        continue
    nl = nl_of(d)
    anchor = b"---" + nl + nl + ANCH_EN
    new = CITE_EN + nl + nl + anchor
    patch(f, [(anchor, new, 1, "cite-EN")])
for f in F_ZH:
    d = read(os.path.join(ROOT, f))
    if CITE_ZH in d:
        LOG.append("    [skip: already applied] %-10s %s" % ("cite-ZH", f))
        continue
    nl = nl_of(d)
    anchor = b"---" + nl + nl + ANCH_ZH
    new = CITE_ZH + nl + nl + anchor
    patch(f, [(anchor, new, 1, "cite-ZH")])

# ---------------- C. 补充图注扩为 S1-S4 ----------------
LOG.append("[C] supplementary legends expanded to S1-S4")
for f in F_LEG_EN:
    patch(f, [(LEG_EN_OLD, LEG_EN_NEW, 1, "legend-EN")])
for f in F_LEG_ZH:
    patch(f, [(LEG_ZH_OLD, LEG_ZH_NEW, 1, "legend-ZH")])

# ---------------- 复核 ----------------
LOG.append("")
LOG.append("--- 复核（期望：1000Genomes=0 ; S1-S4≥1 ; 其余科学文本 0 改动）---")
for f in TOUCH:
    d = read(os.path.join(ROOT, f))
    LOG.append("  %-44s 1000Genomes=%-2d  S1-S4=%-2d  Fig.S4=%-2d  bytes=%d"
               % (f, d.count(b"1000 Genomes"), d.count(b"S1\xe2\x80\x93S4"),
                  d.count(b"Fig. S4") + d.count(u"图 S4".encode("utf-8")), len(d)))

p = os.path.join(ROOT, "logs", "_s69_patch_manuscript.log")
with io.open(p, "w", encoding="utf-8", newline="\n") as fh:
    fh.write("\n".join(LOG) + "\n")
print("VERDICT: s69 v2 done -> %s" % p)
