#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""s38_export_supplementary_data.py -- 补充数据导出（P3 版）

两个交付族
----------
A. **panel_source_data（图面板源数据）**：Fig1–Fig5 每个面板所绘元素的逐行数值，
   口径与 `s37_final_figures.py` 完全一致（复用同一 load()/筛选/排序/派生逻辑），
   保证 CSV 逐行对应图上所绘元素。P3 新增 Fig3d / Fig4d（口径修正）/ Fig4e / Fig4f /
   Fig5d / Fig5e / Fig5f 七个面板。
B. **supplementary_table（补充表 S9–S15）**：P1 阶段（3.1–3.5）新增分析的**完整**结果表，
   编号与稿件「补充表」小节一致（S1–S8 见 `tables/`，本脚本导出 S9–S15）。
   该族为源表**忠实转写**（除文件名与编号外不改动任何单元格），使稿件编号与文件一一对应。
   ★ S9 为 P3-2 验收时发现的**编号缺口**（只有描述、无实体文件），
   见 `tables_s09()` 的缺口说明。

设计原则
--------
1. **只读源表，不修改任何 tables/ 文件。**
2. 全部结果写入文件（不依赖 stdout），符合本机段错误环境下「写文件再读」的纪律。
3. manifest 必须排除自身（`_manifest.csv`），否则自指污染。
4. zip 写 UTF-8 文件名标志位（0x800），避免中文名乱码。

P3 关键修正
-----------
- **删除 `Fig4d_sumpip_nearest_gene.csv`**：P0 审稿修正（硬伤‑4）已废除「按最近基因把 PIP 相加」
  的归属排序，改为「13 个 chr1 检验对的可信集 lead 落点计数」。旧文件与稿件口径矛盾，
  由 `Fig4d_credible_set_lead_placement_counts.csv` + `..._pairs.csv` 取代。
- 落点计数取 **回文清理后** 的可信集 lead（`cs_lead_outcome_nopal_gene`）：清理前 13/13 落
  `CDC42`，清理后为 `WNT4` 11 / `LINC00339` 侧翼 2 / `CDC42` 0 —— 与 Fig4(d) 完全一致。
"""
import os, sys, csv, math, re, json, hashlib, zipfile, datetime
from collections import defaultdict, OrderedDict

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)                       # 11_sc_eqtl_mr_project
BASE = os.path.join(ROOT, "tables")
OUTDIR = os.path.join(ROOT, "supplementary_data")
os.makedirs(OUTDIR, exist_ok=True)
LOG = []
META = []          # manifest 累积：[file, family, label, content, source, bytes, sha256]


def log(*a):
    s = " ".join(str(x) for x in a)
    LOG.append(s)
    print(s)


def load(f):
    """与 s37 完全一致的读取（utf-8-sig，DictReader）。"""
    with open(os.path.join(BASE, f), encoding="utf-8-sig", newline="") as fh:
        return list(csv.DictReader(fh))


def num(x, d=float("nan")):
    try:
        return float(x)
    except Exception:
        return d


def write_csv(name, header, rows, family, label, content, source):
    p = os.path.join(OUTDIR, name)
    with open(p, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        w.writerows(rows)
    nb = os.path.getsize(p)
    META.append([name, family, label, content, source, nb, sha256(p)])
    log("  [CSV] %-52s %6d 行 %9d B  %s" % (name, len(rows), nb, label))
    return p


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def copy_table(src, dst, family, label, content):
    """源表忠实转写：列名与列序保持不变，仅改文件名与编号。"""
    rs = load(src)
    if not rs:
        raise RuntimeError("空表: %s" % src)
    hdr = list(rs[0].keys())
    rows = [[r.get(k, "") for k in hdr] for r in rs]
    return write_csv(dst, hdr, rows, family, label, content, src)


# ============================================================ Fig 1
def fig1_params():
    rows = [
        ["exposure", "OneK1K sc-eQTL (TensorQTL re-derived)", "reference"],
        ["exposure_cell_types", "14", "immune cell types"],
        ["exposure_donors", "980", "donors"],
        ["exposure_build", "GRCh37", "cis-eQTL coordinates"],
        ["replication", "1M-scBloodNL unstimulated (UT)", "direction-consistency only"],
        ["replication_cell_types", "9", "cell types"],
        ["replication_N", "119", "v2=79; v3=40"],
        ["instrument_rule", "cis-eQTL P<5e-8; PLINK clumping r2<0.001 within 10 Mb", "primary"],
        ["instrument_n_per_pair", "1", "26/26 pairs n_index=1"],
        ["outcome_primary", "GCST90483463", "female infertility (all causes); 40,024 cases / 665,658 controls; EUR"],
        ["outcome_primary_pmid", "40229599", "Venkatesh SS et al., Nat Genet 2025"],
        ["outcome_sensitivity", "GCST90483469", "multi-ancestry sensitivity outcome"],
        ["outcome_snp_count", "26550829", "variants in GCST90483463"],
        ["method_mr", "single-SNP Wald ratio (primary); IVW not applicable; Steiger direction OK; MR-Egger not runnable", "MR"],
        ["method_coloc", "coloc.abf (primary); coloc.susie (diagnostic); p12 x window sensitivity", "Colocalization"],
        ["method_conditional", "coloc.abf and PheWAS re-estimated after conditioning on the WNT4 credible-set lead rs56318008 or on locus-specific alternative leads", "Conditional (task 3.1-3.2)"],
        ["method_finemap", "susie_rss (external LD); palindrome cleanup; credible set -> gene attribution", "Fine-mapping"],
        ["method_smr", "cis-SMR with HEIDI (exposures: OneK1K single-cell and GTEx v8 whole blood)", "SMR-HEIDI (task 3.3)"],
        ["method_annotation", "GTEx v10 / eQTLGen; Open Targets; drug and safety annotation", "Annotation & safety"],
        ["safety_endpoints", "2469", "FinnGen R12 disease endpoints"],
        ["safety_instruments", "2", "cis-eQTL instruments of the two robust loci"],
        ["safety_tests", "4938", "2469 endpoints x 2 instruments"],
        ["safety_fdr05_total", "52", "12 risk-increasing + 40 risk-decreasing"],
        ["safety_fdr05_increasing", "12", "risk-increasing"],
        ["safety_fdr05_decreasing", "40", "risk-decreasing"],
        ["safety_immune_blood_tumour_increasing", "0", "0 risk-increasing in immune/blood/tumour domain"],
        ["safety_immune_blood_tumour_median_minOR", "3.15", "80% power, domain Bonferroni; weakest-powered of the five domains (task 3.5)"],
        ["safety_conditional_signals", "0", "after conditioning on the WNT4 credible-set lead, all 52 FDR<0.05 signals vanish (task 3.2)"],
        ["safety_adverse_lineages", "4", "genital prolapse; postpartum haemorrhage; preterm birth; knee osteoarthritis"],
    ]
    write_csv("Fig1_study_design_parameters.csv", ["parameter", "value", "note"], rows,
              "panel_source_data", "Fig1", "Design parameters annotated in the Fig1 schematic",
              "-")


# ============================================================ Fig 2
def fig2_all():
    disc = load("10_discovery_MR_main.csv")
    allf = load("31_final_MR_with_method.csv")
    fin = [r for r in allf if r["thr"].startswith("P<5e-8")]
    ivw = [r for r in allf if not r["thr"].startswith("P<5e-8")][:1]
    rep = load("24_replication_MR_UT.csv")

    byc = defaultdict(list)
    for r in disc:
        byc[int(r["chr"])].append(r)
    order = sorted(byc); off, cur, GAP = {}, 0.0, 6e7
    for c in order:
        off[c] = cur
        cur += max(int(x["tss_pos_grch38"]) for x in byc[c]) + GAP
    thr = -math.log10(0.05 / len(disc))
    rows = []
    for c in order:
        for r in sorted(byc[c], key=lambda x: int(x["tss_pos_grch38"])):
            p = num(r["p"], 1.0)
            rows.append([r["gene"], r["cell_type"], c,
                         int(r["tss_pos_grch38"]),
                         off[c] + int(r["tss_pos_grch38"]),
                         "%.10g" % p,
                         "%.6f" % (-math.log10(max(p, 1e-15))),
                         r.get("F_stat", ""), r.get("qval", ""), r.get("sig_P_FDR", ""),
                         "%.6f" % thr])
    write_csv("Fig2a_manhattan_source_8612tests.csv",
              ["gene", "cell_type", "chr", "tss_pos_grch38", "plot_x_offset",
               "p", "neglog10p", "F", "qval", "sig_P_FDR", "bonferroni_neglog10p"], rows,
              "panel_source_data", "Fig2a", "Gene-level Manhattan source data (8,612 discovery tests)",
              "10_discovery_MR_main.csv")

    ordg = {"CDC42": 0, "LINC00339": 1, "YME1L1": 2, "ANXA4": 3}
    fin = sorted(fin, key=lambda r: (r["locus_id"], ordg.get(r["gene"], 9), r["cell_type"]))
    rows = []
    for i, r in enumerate(fin):
        b, se = num(r["b"]), num(r["se"])
        rows.append(["discovery", i, r["gene"], r["cell_type"], r.get("locus_id", ""),
                     "%.10g" % b, "%.10g" % se, "%.10g" % (b - 1.96 * se), "%.10g" % (b + 1.96 * se),
                     r.get("method", ""), r.get("nsnp", ""), r.get("thr", "")])
    for r in ivw:
        b, se = num(r["b"]), num(r["se"])
        rows.append(["IVW (r2<0.01)", len(fin), r["gene"], r["cell_type"], r.get("locus_id", ""),
                     "%.10g" % b, "%.10g" % se, "%.10g" % (b - 1.96 * se), "%.10g" % (b + 1.96 * se),
                     r.get("method", ""), r.get("nsnp", ""), r.get("thr", "")])
    write_csv("Fig2b_forest_discovery_15pairs_plus_IVW.csv",
              ["series", "plot_y_index", "gene", "cell_type", "locus_id",
               "b", "se", "ci_lo", "ci_hi", "method", "nsnp", "thr"], rows,
              "panel_source_data", "Fig2b", "Forest plot of the 15 FDR-significant pairs plus the only runnable IVW",
              "31_final_MR_with_method.csv")

    rep = sorted(rep, key=lambda r: (ordg.get(r["gene"], 9), r["cell_type_discovery"]))
    rows = []
    for r in rep:
        rows.append([r["gene"], r["cell_type_discovery"],
                     r.get("cell_type_replication", ""),
                     "%.10g" % num(r["discovery_b"]),
                     ("%.10g" % num(r["replication_b"])) if r.get("effect_reportable") == "TRUE" else "",
                     r.get("replication_b", ""), r.get("replication_se", ""),
                     r.get("replication_p", ""), r.get("effect_reportable", ""),
                     r.get("instrument_qc", ""), r.get("effect_status", "")])
    write_csv("Fig2c_replication_comparison.csv",
              ["gene", "cell_type_discovery", "cell_type_replication",
               "discovery_b", "replication_b_plot", "replication_b_raw", "replication_se",
               "replication_p", "effect_reportable", "instrument_qc", "effect_status"], rows,
              "panel_source_data", "Fig2c", "Discovery vs 1M-scBloodNL direction comparison (merged cell-type mapping)",
              "24_replication_MR_UT.csv")
    log("  [check] Manhattan=%d; forest=%d+1IVW; replication=%d" % (len(disc), len(fin), len(rep)))


# ============================================================ Fig 3
def fig3_all():
    col = load("34_coloc_summary.csv")
    sen = load("40_coloc_sensitivity_matrix.csv")
    ROBUST = ["CDC42_B_MEM", "CDC42_Mono_NC"]

    ordg = {"CDC42": 0, "LINC00339": 1, "YME1L1": 2, "ANXA4": 3}
    col = sorted(col, key=lambda r: (ordg.get(r["gene"], 9), r["cell_type"]))
    rows = [[r["gene"], r["cell_type"], r.get("locus_id", ""),
             r.get("abf_H0", ""), r.get("abf_H1", ""), r.get("abf_H2", ""),
             r.get("abf_H3", ""), r.get("abf_H4", ""), r.get("abf_H4_over_H3H4", ""),
             r.get("abf_verdict", ""), r.get("susie_maxH4", ""), r.get("susie_verdict", ""),
             r.get("MR_method", ""), r.get("MR_b", ""), r.get("MR_p", ""), r.get("MR_instrument", "")]
            for r in col]
    write_csv("Fig3a_coloc_pph4_by_pair.csv",
              ["gene", "cell_type", "locus_id", "abf_H0", "abf_H1", "abf_H2", "abf_H3",
               "abf_H4", "abf_H4_over_H3H4", "abf_verdict", "susie_maxH4", "susie_verdict",
               "MR_method", "MR_b", "MR_p", "MR_instrument"], rows,
              "panel_source_data", "Fig3a", "coloc.abf PP.H0-H4 for the 15 pairs with two-tier verdict",
              "34_coloc_summary.csv")

    loci, wl, pl = [], ["500Kb", "1Mb", "2Mb"], ["1e-4", "1e-5", "1e-6"]
    for r in sen:
        if r["locus"] not in loci:
            loci.append(r["locus"])
    loci = [x for x in ROBUST if x in loci] + [x for x in loci if x not in ROBUST]
    M = {}
    for r in sen:
        M[(r["locus"], r["window"], r["p12"])] = num(r["PP.H4"])
    rows = []
    for i, loc in enumerate(loci):
        for j, w in enumerate(wl):
            for k, p in enumerate(pl):
                v = M.get((loc, w, p), float("nan"))
                rows.append([loc, "robust" if loc in ROBUST else "unstable",
                             i, w, p, j * len(pl) + k,
                             "%.6g" % v, "TRUE" if v < 0.8 else "FALSE"])
    write_csv("Fig3b_sensitivity_matrix_5loci_x9settings.csv",
              ["locus", "robust_class", "plot_row", "window", "p12", "plot_col",
               "PP_H4", "below_0.8"], rows,
              "panel_source_data", "Fig3b", "PP.H4 sensitivity heatmap (5 loci x 9 window/prior settings)",
              "40_coloc_sensitivity_matrix.csv")

    rows = []
    for loc in loci:
        for k, p in enumerate(pl):
            vs = [num(r["PP.H4"]) for r in sen if r["locus"] == loc and r["p12"] == p]
            rows.append([loc, "robust" if loc in ROBUST else "unstable", p, k,
                         "%.6g" % (sum(vs) / len(vs)), len(vs)])
    write_csv("Fig3c_pph4_vs_p12_mean_over_windows.csv",
              ["locus", "robust_class", "p12", "plot_x_index", "PP_H4_mean", "n_windows"], rows,
              "panel_source_data", "Fig3c", "PP.H4 vs prior p12 (mean over windows)",
              "40_coloc_sensitivity_matrix.csv")

    # ---- P3 新增 (d) 条件共定位（task 3.1）----
    Q1 = "rs56318008"
    Q2 = {"B_MEM": "rs2473247", "Mono_NC": "rs12048511"}
    cells = ["B_MEM", "Mono_NC"]
    series = {"uncond": "unconditional",
              "cond_eqtl_only": "conditioned on eQTL only (not plotted)",
              "cond_both": "conditioned on eQTL and the conditioning SNP"}
    cc = load("43_cond_coloc_CDC42.csv")
    rows = []
    for c in cells:
        for snp, tag in [(Q1, "uncond"), (Q1, "cond_eqtl_only"), (Q1, "cond_both"),
                         (Q2[c], "uncond"), (Q2[c], "cond_eqtl_only"), (Q2[c], "cond_both")]:
            for r in cc:
                if (r["scope"] == c and r["keep_rule"] == "0.01" and
                        r["cond_snp"] == snp and r["tag"] == tag):
                    plotted = ""
                    if snp == Q1 and tag in ("uncond", "cond_both"):
                        plotted = {"uncond": "unconditional",
                                   "cond_both": "+ WNT4 credible-set lead"}[tag]
                    elif snp == Q2[c] and tag == "cond_both":
                        plotted = "+ alternative lead"
                    rows.append([c, snp, "conditioning SNP" if snp == Q1 else "alternative lead",
                                 tag, series[tag], plotted,
                                 r["nsnps"], r["PP.H0"], r["PP.H1"], r["PP.H2"], r["PP.H3"], r["PP.H4"],
                                 r["max_abs_z_e_pre"], r["n_e_pre_bonf"],
                                 r["max_abs_z_e_post"], r["n_e_post_bonf"],
                                 r["max_abs_z_o_pre"], r["n_o_pre_bonf"],
                                 r["max_abs_z_o_post"], r["n_o_post_bonf"], r["zcheck"]])
    write_csv("Fig3d_conditional_colocalization.csv",
              ["cell_type", "cond_snp", "cond_snp_role", "tag", "family", "plotted_series",
               "nsnps", "PP_H0", "PP_H1", "PP_H2", "PP_H3", "PP_H4",
               "max_abs_z_eqtl_pre", "n_eqtl_pre_bonf", "max_abs_z_eqtl_post", "n_eqtl_post_bonf",
               "max_abs_z_outcome_pre", "n_outcome_pre_bonf",
               "max_abs_z_outcome_post", "n_outcome_post_bonf", "zcheck"], rows,
              "panel_source_data", "Fig3d",
              "Conditional colocalization at the two robust CDC42 loci (PP.H4 and max|z| before/after conditioning)",
              "43_cond_coloc_CDC42.csv")
    k = {(r[0], r[1], r[3]): r[11] for r in rows}
    log("  [check] Fig3d uncond/cond_both(WNT4 lead)=%s/%s ; B_MEM alt=%s ; Mono_NC alt=%s"
        % (k[("B_MEM", Q1, "uncond")], k[("B_MEM", Q1, "cond_both")],
           k[("B_MEM", Q2["B_MEM"], "cond_both")], k[("Mono_NC", Q2["Mono_NC"], "cond_both")]))
    log("  [check] coloc 15 对；sensitivity loci=%d x 9；Fig3d=%d 行" % (len(loci), len(rows)))


# ============================================================ Fig 4
def fig4_all():
    ld = load("23_chr1_LD_matrix.csv")
    V = load("35_finemap_variant.csv")
    o = [r for r in V if r["scope"] == ""]

    keys = list(ld[0].keys())
    labels = [r[keys[0]] for r in ld]
    short = [l.split(" / ")[0] for l in labels]
    rows = []
    for i, r in enumerate(ld):
        for j, k in enumerate(keys[1:]):
            rows.append([short[i], short[j], i, j, "%.6g" % num(r[k], 0.0)])
    write_csv("Fig4a_LD_matrix_7x7.csv",
              ["variant_row", "variant_col", "row_idx", "col_idx", "r2"], rows,
              "panel_source_data", "Fig4a", "Pairwise LD (r2) among the 7 variants analysed at the locus",
              "23_chr1_LD_matrix.csv")

    genes = [("CDC42", 22388297, 22397199),
             ("LINC00339", 22352108, 22355978),
             ("WNT4", 22444975, 22470451)]
    marks = [(22462111, "rs12037376", "B_MEM instrument", -2.05),
             (22422721, "rs10917151", "Mono_NC instrument", -1.25),
             (22470407, "rs56318008", "outcome credible-set lead", -2.95)]
    rows = [["gene", g, s, e, "gene_body"] for g, s, e in genes]
    rows += [["variant", rs, p, "", tag] for p, rs, tag, _my in marks]
    write_csv("Fig4b_gene_variant_track.csv",
              ["type", "name", "grch37_pos_start", "grch37_pos_end", "role"], rows,
              "panel_source_data", "Fig4b", "Gene and variant track for chr1p36.12 (GRCh37)",
              "coordinates hard-coded from GRCh37 annotation (as plotted)")

    top = sorted(o, key=lambda r: -num(r["PIP"]))[0]
    rows = [[r["rsid"], r["pos37"], r.get("snp37", ""), r.get("gene_annot", ""),
             "%.10g" % num(r["PIP"], 0.0),
             "TRUE" if r["rsid"] == top["rsid"] else "FALSE"] for r in o]
    write_csv("Fig4c_finemap_pip_outcome_only.csv",
              ["rsid", "pos37", "snp37", "gene_annot", "PIP", "is_top"], rows,
              "panel_source_data", "Fig4c", "Outcome-only fine-mapping PIP across the 900-kb window",
              "35_finemap_variant.csv")
    nn = len(rows)

    # ---- P0 口径修正 (d)：可信集 lead 落点计数（取代已废除的 ΣPIP 最近基因排序）----
    fp = load("36_finemap_pair.csv")
    pairs = [r for r in fp if r["scope"] != "outcome_only"]
    rows = [[r["scope"], r["nsnp"], r["n_cs_outcome"], r["n_cs_outcome_nopal"], r["n_palindromic"],
             r["cs_lead_outcome_rsid"], r["cs_lead_outcome_pos37"], r["cs_lead_outcome_gene"],
             r["cs_lead_outcome_nopal_pos37"], r["cs_lead_outcome_nopal_gene"],
             r["r2_instr_vs_cs_lead_outcome"], r["max_pip_outcome"]
             ] for r in pairs]
    write_csv("Fig4d_credible_set_lead_placement_pairs.csv",
              ["scope", "nsnp", "n_cs_outcome", "n_cs_outcome_nopal", "n_palindromic",
               "cs_lead_outcome_rsid", "cs_lead_outcome_pos37", "cs_lead_outcome_gene",
               "cs_lead_outcome_palindrome_cleaned_pos37", "cs_lead_outcome_palindrome_cleaned_gene",
               "r2_instrument_vs_cs_lead_outcome", "max_pip_outcome"], rows,
              "panel_source_data", "Fig4d (per pair)",
              "Credible-set lead placement for the 13 chr1 test pairs, before and after palindrome cleanup",
              "36_finemap_pair.csv")
    grp = OrderedDict([("WNT4", 0), ("LINC00339", 0), ("CDC42", 0)])
    for r in pairs:
        g = r["cs_lead_outcome_nopal_gene"]
        key = "LINC00339" if g.startswith("flank") and "LINC00339" in g else g
        if key in grp:
            grp[key] += 1
    rows = [[g, grp[g], i] for i, g in enumerate(grp)]
    write_csv("Fig4d_credible_set_lead_placement_counts.csv",
              ["lead_gene_interval", "n_pairs", "plot_row"], rows,
              "panel_source_data", "Fig4d (bars)",
              "Fig4(d) bars: placement of the palindrome-cleaned credible-set lead across 13 chr1 pairs",
              "36_finemap_pair.csv")
    log("  [check] Fig4d pairs=%d, counts=%s (合计 %d)；LD 矩阵 %d 行；outcome_only=%d"
        % (len(pairs), dict(grp), sum(grp.values()), len(ld) * len(ld), nn))

    # ---- P3 新增 (e) cis-SMR + HEIDI 森林图（task 3.3）----
    sm = load("45c_smr_3p3_combined.csv")
    rows = []
    for r in sm:
        b, se = num(r["b_SMR"]), num(r["se_SMR"])
        rows.append(["OneK1K single-cell", r["analysis"], r["Gene"], r["cell_type"],
                     r["topSNP"], b, se, b - 1.96 * se, b + 1.96 * se,
                     r["p_SMR"], r["fdr_smr"], r["direction"], r["p_HEIDI"], r["nsnp_HEIDI"],
                     r["heidi_verdict"], "TRUE" if r["heidi_verdict"] == "consistent" else "FALSE"])
    g8 = load("47_gtex8_wb_smr_results.csv")
    for r in g8:
        if r["Gene"] != "WNT4":
            continue
        b, se = num(r["b_SMR"]), num(r["se_SMR"])
        rows.append(["GTEx v8 whole blood", r["analysis"], r["Gene"], "whole blood",
                     r["topSNP"], b, se, b - 1.96 * se, b + 1.96 * se,
                     r["p_SMR"], "", ("negative" if b < 0 else "positive"),
                     r["p_HEIDI"], r["nsnp_HEIDI"], "n/a (not FDR-significant)", "FALSE"])
    write_csv("Fig4e_smr_heidi_forest.csv",
              ["source_dataset", "analysis", "gene", "cell_or_tissue", "topSNP",
               "b_SMR", "se_SMR", "ci_lo", "ci_hi", "p_SMR", "fdr_smr", "direction",
               "p_HEIDI", "nsnp_HEIDI", "heidi_verdict", "marker_filled_heidi_consistent"], rows,
              "panel_source_data", "Fig4e",
              "cis-SMR with HEIDI: b_SMR (95% CI) for CDC42 / LINC00339 in OneK1K and WNT4 in GTEx v8 whole blood",
              "45c_smr_3p3_combined.csv; 47_gtex8_wb_smr_results.csv")

    # ---- P3 新增 (f) 工具邻近基因检查（task 3.4）----
    rt = load("49b_task34_region_tag_summary.csv")
    gp = OrderedDict([("CDC42", 0), ("LINC00339", 1), ("WNT4", 2)])
    snps = sorted({r["SNP_rsID"] for r in rt})
    idx = {(r["SNP_rsID"], r["gene"]): r for r in rt}
    rows = []
    for s in snps:
        for g, gx in gp.items():
            r = idx.get((s, g), {})
            n = int(num(r.get("OneK1K_n_sig"), 0))
            rows.append([s, g, gx, r.get("OneK1K_n_celltypes", ""),
                         r.get("OneK1K_n_available", ""), n,
                         r.get("OneK1K_min_p", ""), r.get("OneK1K_sig_celltypes", ""),
                         r.get("OneK1K_sign_consistency_sig_only", ""),
                         r.get("GTEx_v8_p_eQTL", ""), r.get("GTEx_v8_signif", ""),
                         "FALSE" if g == "WNT4" else "TRUE"])
    n_e = len(sm) + 1          # 45c 全部行 + GTEx v8 的 WNT4 行
    n_f = len(rows)            # Fig4f：2 个工具 x 3 个基因
    write_csv("Fig4f_tool_neighbour_gene_check.csv",
              ["instrument", "gene", "plot_x_index", "n_celltypes_tested",
               "n_celltypes_in_panel", "n_significant_cis_eQTL",
               "min_p", "significant_celltypes", "sign_consistency",
               "GTEx_v8_whole_blood_p", "GTEx_v8_significant", "available_in_OneK1K_panel"], rows,
              "panel_source_data", "Fig4f",
              "Instrument neighbouring-gene cis-eQTL check: both instruments are region tags; WNT4 is absent from the OneK1K panel",
              "49b_task34_region_tag_summary.csv")
    log("  [check] Fig4e=%d 行（45c 的 %d + GTEx v8 WNT4 1）；Fig4f=%d 行"
        % (n_e, len(sm), n_f))


# ============================================================ Fig 5
def fig5_all():
    ot = load("41d_opentargets_diseases_full.csv")
    dom = load("42b_phewas_safety_domains.csv")
    up = load("42_phewas_safety_assessment.csv")
    pw = load("50c_task35_phewas_power_by_domain.csv")
    pc = load("50d_task35_power_class_by_domain_968.csv")
    cd = load("44c_phewas_3p2_full_domains.csv")

    ot = sorted(ot, key=lambda r: -num(r["overall_score"]))[:18][::-1]
    rows = []
    for i, r in enumerate(ot):
        ta = r["therapeutic_areas"].lower()
        star = any(k in ta for k in ["reproductive", "hematologic", "immune"])
        rows.append([i, r["gene"], r["disease_name"], "%.6g" % num(r["overall_score"]),
                     r.get("therapeutic_areas", ""), "TRUE" if star else "FALSE"])
    write_csv("Fig5a_opentargets_top18.csv",
              ["plot_row", "gene", "disease_name", "overall_score",
               "therapeutic_areas", "flagged_repro_heme_immune"], rows,
              "panel_source_data", "Fig5a", "Top 18 Open Targets gene-disease associations",
              "41d_opentargets_diseases_full.csv")

    LG = OrderedDict([
        ("A. 盆底支持结构/结缔组织", "A pelvic-floor / connective tissue"),
        ("B. 胎盘附着/剥离异常与产科出血", "B placental / obstetric haemorrhage"),
        ("C. 早产", "C preterm birth"),
        ("D. 关节软骨退变", "D knee osteoarthritis"),
        # ★★ 2026-10-02（M5 落地）：显示标签由未完成状态词 "E pending verification"
        #    改为真实谱系名。★ key 必须逐字符等于 42_phewas_safety_assessment.csv 的
        #    lineage_group 值（join key，由 s36 生成，s36 依赖外置盘 E:\ 本轮不在线）
        #    ⇒ key 保持原样，只改显示标签。
        ("E. 可疑/待复核", "E residual foreign body in soft tissue"),
        ("F. 非疾病生育表型(不计为风险)", "F non-disease phenotype"),
    ])
    inc = [r for r in up if r["risk_if_cdc42_inhibited"] == "increased"]
    inc = sorted(inc, key=lambda r: -(-num(r["beta_MR"])))
    rows = []
    for i, r in enumerate(inc):
        m = -num(r["beta_MR"]); e = 1.96 * num(r["se_MR"])
        rows.append([i, r.get("lineage_group", ""), LG.get(r.get("lineage_group", ""), ""),
                     r.get("phenotype", ""), r.get("cell_type", ""),
                     "%.6g" % num(r["beta_MR"]), "%.6g" % num(r["se_MR"]),
                     "%.6g" % m, "%.6g" % (m - e), "%.6g" % (m + e),
                     r.get("p_FDR", ""), r.get("risk_if_cdc42_inhibited", "")])
    write_csv("Fig5bc_phewas_risk_increasing_12.csv",
              ["plot_row", "lineage_group", "lineage_label_en", "phenotype", "cell_type",
               "beta_MR", "se_MR", "neg_beta_MR_plot", "ci_lo_plot", "ci_hi_plot",
               "p_FDR", "risk_if_cdc42_inhibited"], rows,
              "panel_source_data", "Fig5b", "The 12 region-related risk-increasing signals (plotted as -beta_MR)",
              "42_phewas_safety_assessment.csv")

    order = ["免疫/血液/肿瘤", "女性生殖/不孕", "妊娠/分娩/产褥", "肌肉骨骼/结缔组织", "其他"]
    dm = {r["domain"]: r for r in dom}
    rows = []
    for i, d in enumerate(order):
        r = dm.get(d, {})
        rows.append([d, r.get("n_tests", ""), r.get("n_sig", ""),
                     r.get("n_risk_up", ""), r.get("n_risk_down", ""), r.get("min_p_MR", ""), i])
    write_csv("Fig5bc_phewas_domain_rollup.csv",
              ["domain", "n_tests", "n_fdr05", "n_risk_up", "n_risk_down", "min_p", "plot_row"], rows,
              "panel_source_data", "Fig5c", "Domain-level rollup of FDR < 0.05 signals (risk-increasing above axis)",
              "42b_phewas_safety_domains.csv")

    # ---- P3 新增 (d) 五域中位最小可检出 OR（task 3.5）----
    dmP = {r["domain"]: r for r in pw}
    dorder = ["免疫/血液/肿瘤", "其他", "肌肉骨骼/结缔组织", "女性生殖/不孕", "妊娠/分娩/产褥"]
    rows = []
    for i, d in enumerate(dorder):
        r = dmP.get(d, {})
        rows.append([d, r.get("n_tests", ""), r.get("n_endpoints", ""),
                     r.get("n_cases_median", ""), r.get("min_detectable_OR_median", ""),
                     r.get("min_detectable_OR_p05", ""), r.get("min_detectable_OR_p95", ""),
                     r.get("n_detectable_OR_le_1p5", ""), r.get("pct_detectable_OR_le_1p5", ""),
                     r.get("n_severe_OR_gt_3", ""), r.get("n_FDR05", ""), r.get("n_risk_up", ""),
                     r.get("n_risk_down", ""), r.get("min_p_MR_risk_up", ""), i])
    write_csv("Fig5d_domain_minOR_power.csv",
              ["domain", "n_tests", "n_endpoints", "n_cases_median", "min_detectable_OR_median",
               "min_detectable_OR_p05", "min_detectable_OR_p95", "n_detectable_OR_le_1p5",
               "pct_detectable_OR_le_1p5", "n_severe_OR_gt_3", "n_FDR05",
               "n_risk_up", "n_risk_down", "min_p_MR_risk_up", "plot_row"], rows,
              "panel_source_data", "Fig5d",
              "Fig5(d) bars: median minimum detectable OR (80% power, domain Bonferroni) by disease domain",
              "50c_task35_phewas_power_by_domain.csv")

    # ---- P3 新增 (e) 968 检验功效分层 ----
    cls = pc[0]
    keys5 = [("1_adequate_OR<=1.2", "<=1.2"), ("2_adequate_1.2<OR<=1.5", "1.2-1.5"),
             ("3_limited_1.5<OR<=2.0", "1.5-2.0"), ("4_inadequate_2.0<OR<=3.0", "2.0-3.0"),
             ("5_severe_OR>3.0", ">3.0")]
    tot = sum(int(cls[k]) for k, _l in keys5)
    rows = []
    for i, (k, lab) in enumerate(keys5):
        v = int(cls[k])
        rows.append([lab, v, "%.4f" % (100.0 * v / tot), i])
    write_csv("Fig5e_power_class_968.csv",
              ["power_class_smallest_detectable_OR", "n_tests", "pct_of_968", "plot_x_index"], rows,
              "panel_source_data", "Fig5e",
              "Fig5(e) bars: distribution of the smallest detectable OR across the 968 immune/blood/tumour tests",
              "50d_task35_power_class_by_domain_968.csv")

    # ---- P3 新增 (f) 条件 PheWAS 按域（task 3.2）----
    fam = {}
    for r in cd:
        fam.setdefault(r["domain"], {})[r["family"]] = r
    fams = [("uncond", "unconditional", "TRUE"),
            ("cond_eqtl_only", "conditioned on eQTL", "TRUE"),
            ("cond_both", "conditioned on eQTL and region lead", "TRUE"),
            ("region", "region instrument arm (not plotted)", "FALSE")]
    dof = ["免疫/血液/肿瘤", "女性生殖/不孕", "妊娠/分娩/产褥", "肌肉骨骼/结缔组织", "其他"]
    rows = []
    for d in dof:
        for fk, lab, pl in fams:
            r = fam.get(d, {}).get(fk, {})
            rows.append([d, fk, lab, pl, r.get("n_tests", ""), r.get("n_sig", ""),
                         r.get("n_risk_up", ""), r.get("n_risk_down", ""), r.get("min_p_MR", "")])
    write_csv("Fig5f_conditional_phewas_by_domain.csv",
              ["domain", "family", "family_label", "plotted", "n_tests", "n_fdr05",
               "n_risk_up", "n_risk_down", "min_p_MR"], rows,
              "panel_source_data", "Fig5f",
              "Fig5(f) bars: conditional PheWAS by domain; conditioning removes all 52 FDR<0.05 signals",
              "44c_phewas_3p2_full_domains.csv")
    u = sum(int(fam[d]["uncond"]["n_sig"]) for d in dof)
    cb = sum(int(fam[d]["cond_both"]["n_sig"]) for d in dof)
    ce = sum(int(fam[d]["cond_eqtl_only"]["n_sig"]) for d in dof)
    log("  [check] OT top18; risk-increasing=%d; domains=%d" % (len(inc), len(order)))
    log("  [check] Fig5e 累计=%d；Fig5f uncond=%d cond_eqtl_only=%d cond_both=%d" % (tot, u, ce, cb))


# ============================================================ 补充表 S9
def tables_s09():
    """Table S9 · FinnGen R12 PheWAS 全量关联表（此前只有描述、从未导出实体文件）。

    缺口背景（P3-2 验收时暴露）：稿件「补充表」小节早在 P0 即已列出 Table S9
    「4,938 检验全表 + 52 条显著项 + 8 个髓系代理端点核查」，但该表组在补充包中
    **不存在任何实体 CSV**（S10–S15 已导出 22 件，S1–S8/S9 一直只有编号）。
    在 E 盘（迁出包 `生信7` 与项目旧副本 `生信2`）均检索不到早期文件名
    `20_phewas_all_all.csv` —— 那两张盘都是 2026-09-25 前后的快照，只到 `42*`。
    本机 `tables/44_phewas_3p2_full_all.csv`（5,079,481 B = 5.08 MB，与旧记录
    「`20_phewas_all_all.csv` 4,938 行 / ≈5.08 MB」字节数一致）即该表的**同源后继**
    （P3-2 经 tabix 重跑并扩为 4 臂），故 S9 缺口在本机闭环，无需外部盘。

    导出内容
    --------
    S9a  4 臂全量关联（uncond / cond_eqtl_only / cond_both / region）= 19,752 行
    S9b  主扫描（uncond 臂）= 4,938 行（与稿件「4,938 检验全表」字面一致）
    S9c  52 条 FDR<0.05 显著项（12 风险升高 + 40 风险降低）

    ★「8 个髓系代理端点核查」不在此重复导出：该表即 `TableS15d`（源 `42c`），
    稿件 S9 条目已改为交叉引用 S15d，避免同一文件在包内出现两次。
    """
    full = load("44_phewas_3p2_full_all.csv")
    hdr = list(full[0].keys())
    rows = [[r.get(k, "") for k in hdr] for r in full]
    write_csv("TableS9a_phewas_full_association_all_arms.csv", hdr, rows,
              "supplementary_table", "Table S9a",
              "FinnGen R12 PheWAS full association table: 4,938 tests x 4 analysis arms "
              "(uncond / cond_eqtl_only / cond_both / region) = 19,752 rows",
              "44_phewas_3p2_full_all.csv")

    un = [r for r in full if r.get("family") == "uncond"]
    rows = [[r.get(k, "") for k in hdr] for r in un]
    write_csv("TableS9b_phewas_primary_scan_uncond_4938.csv", hdr, rows,
              "supplementary_table", "Table S9b",
              "FinnGen R12 PheWAS primary scan (unconditional arm): all 4,938 tests",
              "44_phewas_3p2_full_all.csv")

    sig = load("42_phewas_safety_assessment.csv")
    copy_table("42_phewas_safety_assessment.csv", "TableS9c_phewas_fdr_significant_52.csv",
               "supplementary_table", "Table S9c",
               "The 52 FDR<0.05 PheWAS signals (12 risk-increasing + 40 risk-decreasing) "
               "with MR effect, FDR tier and harmony status")

    print("  [check] TableS9: 源表 %d 行 / %d 臂；S9a=%d 行；S9b(uncond)=%d 行；S9c=%d 行"
          % (len(full), len({r.get("family") for r in full}), len(full), len(un), len(sig)))
    LOG.append("  [check] TableS9: 源表 %d 行 / %d 臂；S9a=%d 行；S9b(uncond)=%d 行；S9c=%d 行"
               % (len(full), len({r.get("family") for r in full}), len(full), len(un), len(sig)))


# ============================================================ 补充表 S10–S15
def tables_s10_s15():
    """P1 阶段（3.1–3.5）新增分析的完整结果表，编号与稿件「补充表」小节一致。

    该族为源表忠实转写（列名与列序不变，仅改文件名与编号），
    再打包一行 provenance 记录进 manifest,保证「稿件 Table S10a = tables/43_*.csv」可追。
    """
    # ---- Table S10 · 条件共定位（3.1）----
    copy_table("43_cond_coloc_CDC42.csv", "TableS10a_cond_coloc_PPH4_all_settings.csv",
               "supplementary_table", "Table S10a",
               "Conditional colocalization PP.H0-H4 and max|z| before/after conditioning (all settings)")
    copy_table("43b_cond_coloc_variant_detail.csv", "TableS10b_cond_coloc_variant_detail.csv",
               "supplementary_table", "Table S10b",
               "Per-variant z-scores before/after conditioning at the conditioning SNPs")
    copy_table("43c_cond_coloc_negctl.csv", "TableS10c_cond_coloc_negative_control.csv",
               "supplementary_table", "Table S10c",
               "Conditional colocalization negative control")

    # ---- Table S11 · 条件 PheWAS（3.2）----
    copy_table("44b_phewas_3p2_full_summary.csv", "TableS11a_conditional_phewas_summary.csv",
               "supplementary_table", "Table S11a",
               "Conditional PheWAS: FDR counts per analysis family (uncond / cond_eqtl_only / cond_both / region)")
    copy_table("44c_phewas_3p2_full_domains.csv", "TableS11b_conditional_phewas_by_domain.csv",
               "supplementary_table", "Table S11b",
               "Conditional PheWAS: FDR counts per analysis family x disease domain")
    copy_table("44d_phewas_3p2_full_risk12.csv", "TableS11c_conditional_phewas_risk_increasing_12.csv",
               "supplementary_table", "Table S11c",
               "The 12 region-related risk-increasing signals with beta and FDR in all four analysis arms")
    copy_table("44e_phewas_3p2_full_signals.csv", "TableS11d_conditional_phewas_all_fdr_signals.csv",
               "supplementary_table", "Table S11d",
               "All 100 FDR-significant signals across the four analysis arms")
    copy_table("44f_phewas_3p2_region_variant_raw.csv", "TableS11e_conditional_phewas_region_arm_raw.csv",
               "supplementary_table", "Table S11e",
               "Region-arm raw association values for all 2,469 FinnGen R12 endpoints")

    # ---- Table S12 · SMR + HEIDI（3.3）----
    copy_table("45c_smr_3p3_combined.csv", "TableS12a_smr_heidi_combined.csv",
               "supplementary_table", "Table S12a",
               "cis-SMR with HEIDI for CDC42 and LINC00339 in OneK1K (main and sensitivity analysis)")
    copy_table("46_smr_3p3_gene_summary.csv", "TableS12b_smr_gene_summary.csv",
               "supplementary_table", "Table S12b",
               "Per-gene SMR/HEIDI summary (cell types tested, FDR-significant count, b_SMR range, HEIDI verdicts)")
    copy_table("47_gtex8_wb_smr_results.csv", "TableS12c_smr_WNT4_GTEx_v8_whole_blood.csv",
               "supplementary_table", "Table S12c",
               "GTEx v8 whole-blood cis-SMR substitute analysis (WNT4 not detected: p_SMR = 0.148, p_HEIDI = 0.292, m = 11)")

    # ---- Table S13 · 工具邻近基因 cis 效应检查（3.4）----
    copy_table("49_task34_tool_neighbour_cis_eqtl.csv", "TableS13a_tool_neighbour_cis_eqtl.csv",
               "supplementary_table", "Table S13a",
               "cis-eQTL of both instruments for CDC42 / LINC00339 / WNT4 across OneK1K and GTEx v8 whole blood")
    copy_table("49b_task34_region_tag_summary.csv", "TableS13b_region_tag_summary.csv",
               "supplementary_table", "Table S13b",
               "Region-tag summary per instrument x gene (OneK1K significant cell-type count; WNT4 absent from panel)")
    copy_table("49c_task34_gwas_context.csv", "TableS13c_instrument_gwas_context.csv",
               "supplementary_table", "Table S13c",
               "Outcome GWAS context for both instruments (alleles, frequency, beta, p, N)")

    # ---- Table S14 · PheWAS 功效分析（3.5）----
    copy_table("50_task35_phewas_power_primary_968.csv", "TableS14a_phewas_power_endpoint_level_968.csv",
               "supplementary_table", "Table S14a",
               "Power analysis of the 968 immune/blood/tumour tests: minOR, case counts, power class, cases needed")
    copy_table("50g_task35_endpoint_level_power_968.csv", "TableS14b_phewas_power_endpoint_collapsed_484.csv",
               "supplementary_table", "Table S14b",
               "Endpoint-level (484) collapsed power summary across the two instruments")
    copy_table("50f_task35_worst_power_endpoints_968.csv", "TableS14c_phewas_power_worst15_endpoints.csv",
               "supplementary_table", "Table S14c",
               "The 15 worst-powered endpoints in the 968-test domain")
    copy_table("50b_task35_phewas_power_all4938.csv", "TableS14d_phewas_power_all4938.csv",
               "supplementary_table", "Table S14d",
               "Power analysis for all 4,938 tests (full coverage, mirrors the scope of Table S9)")

    # ---- Table S15 · 域级汇总与安全性核查（对应主表新增结果）----
    copy_table("50d_task35_power_class_by_domain_968.csv", "TableS15a_power_class_by_domain_968.csv",
               "supplementary_table", "Table S15a",
               "Power-class distribution of the 968 immune/blood/tumour tests (the five-domain comparison)")
    copy_table("50c_task35_phewas_power_by_domain.csv", "TableS15b_phewas_power_by_domain.csv",
               "supplementary_table", "Table S15b",
               "Domain-level power and signal summary (median minOR, case-count distribution, FDR counts)")
    copy_table("50e_task35_endpoint_absence_search.csv", "TableS15c_endpoint_absence_search.csv",
               "supplementary_table", "Table S15c",
               "Existence search for named safety endpoints in FinnGen R12 (not assessable vs not tested)")
    copy_table("42c_phewas_myeloid_proxy_check.csv", "TableS15d_myeloid_proxy_endpoint_check.csv",
               "supplementary_table", "Table S15d",
               "Myeloid toxicity proxy endpoints: 8/8 present, 16 tests, 0 risk-increasing")


# ============================================================ Fig S1 / S2（2026-09-29 新增）
def fig_s1_s2():
    """补充图 Fig. S1 / Fig. S2 的 panel 源数据。

    源 = scripts/s67_make_supp_figs.py 从**已闭环表**派生的 64a / 64b：
      64a <- tables/20_celltype_power.csv
      64b <- tables/39_sensitivity_vs_main_GCST90483469.csv
             + tables/15_discovery_significant_with_locus.csv（主分析 se / locus_id）
    忠实转写（列名与列序不变），仅改文件名，沿用 copy_table。
    """
    copy_table("64a_FigS1_celltype_power_gating.csv", "FigS1_celltype_power_gating.csv",
               "panel_source_data", "FigS1",
               "Fig. S1(a-d): per-cell-type instrument counts, nominal and FDR-significant test "
               "counts, minimum q value and power flags (14 cell types; 182.4 x 70.6 mm)")
    copy_table("64b_FigS2_sensitivity_vs_main_15pairs.csv", "FigS2_sensitivity_vs_main_15pairs.csv",
               "panel_source_data", "FigS2",
               "Fig. S2(a-b): main vs sensitivity outcome (GCST90483469) betas for the 15 "
               "FDR<0.05 concordant pairs, with main-outcome SE and locus (182.5 x 77.7 mm)")
    copy_table("64d_FigS3S4_discovery_thresholds.csv", "FigS3S4_discovery_thresholds.csv",
               "panel_source_data", "FigS3;FigS4",
               "Fig. S3 / Fig. S4: per-cell-type discovery test counts and the two significance "
               "thresholds used (Bonferroni p = 0.05/8612; BH p = 6.4277e-05 for FDR < 0.05). "
               "Point-level source data are the same 8,612 discovery tests as Fig2a, exported as "
               "Fig2a_manhattan_source_8612tests.csv and not duplicated here")


# ============================================================ 主流程
def main():
    log("=" * 78)
    log("s38 补充数据导出（P3 版）  %s" % datetime.datetime.now().isoformat(timespec="seconds"))
    log("源目录: %s" % BASE)
    log("输出:   %s" % OUTDIR)
    log("-" * 78)
    log("[A] panel_source_data")
    fig1_params()
    fig2_all()
    fig3_all()
    fig4_all()
    fig5_all()
    fig_s1_s2()
    log("-" * 78)
    log("[B] supplementary_table S9-S15")
    tables_s09()
    tables_s10_s15()

    # 清理已废除的旧口径文件（ΣPIP 最近基因排序，P0 硬伤-4 已废除）
    stale = os.path.join(OUTDIR, "Fig4d_sumpip_nearest_gene.csv")
    if os.path.exists(stale):
        os.remove(stale)
        log("[cleanup] 已删除废除口径文件 Fig4d_sumpip_nearest_gene.csv")

    # manifest（★必须排除自身，否则自指污染）
    SELF = "_manifest.csv"
    hdr = ["file", "family", "label", "content", "source", "bytes", "sha256"]
    rows = [m for m in META if m[0] != SELF]
    rows.sort(key=lambda r: r[0])
    man = [[r[0], r[1], r[2], r[3], r[4], r[5], r[6]] for r in rows]
    write_csv(SELF, hdr, man, "manifest", "-", "Index of all supplementary data files", "-")

    # zip（UTF-8 文件名标志位）
    zip_path = os.path.join(ROOT, "Supplementary_Data.zip")
    allf = sorted(os.listdir(OUTDIR))
    with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as z:
        for f in allf:
            zi = zipfile.ZipInfo("Supplementary_Data/" + f,
                                 date_time=(2026, 9, 28, 12, 0, 0))
            zi.flag_bits |= 0x800
            zi.compress_type = zipfile.ZIP_DEFLATED
            with open(os.path.join(OUTDIR, f), "rb") as fh:
                z.writestr(zi, fh.read())
    log("-" * 78)
    log("CSV 数: %d（panel %d + table %d + manifest 1）"
        % (len(man) + 1,
           sum(1 for r in man if r[1] == "panel_source_data"),
           sum(1 for r in man if r[1] == "supplementary_table")))
    log("zip: %s (%d B)  sha256=%s" % (zip_path, os.path.getsize(zip_path), sha256(zip_path)))

    with open(os.path.join(ROOT, "logs", "_s38_run_p3.log"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(LOG) + "\n")
    print("DONE")


if __name__ == "__main__":
    main()
