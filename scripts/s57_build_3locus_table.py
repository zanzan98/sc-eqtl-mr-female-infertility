# -*- coding: utf-8 -*-
"""s57_build_3locus_table.py -- 任务A 步骤4：三位点（chr1 / chr10 / chr2）比较表

输入（全部为既有产物，只读）：
  tables/31_final_MR_with_method.csv            MR 主结果（Wald ratio 行）
  tables/34_coloc_summary.csv                   coloc.abf / coloc.susie 汇总
  tables/40_coloc_sensitivity_matrix.csv        chr1(+L3) 45 组敏感性
  tables/51_coloc_sensitivity_L2.csv            L2 敏感性
  tables/51b_coloc_sensitivity_L3_recheck.csv   L3 敏感性（复算）
  tables/43_cond_coloc_CDC42.csv                 chr1 条件共定位（v2）
  tables/43c_cond_coloc_negctl.csv               chr1 阴性对照
  tables/54c_cond_coloc_summary_L2L3.csv         L2/L3 条件共定位
  tables/36_finemap_pair.csv                     chr1 精细定位
  tables/52b_finemap_pair_L2.csv / _L3.csv       L2/L3 精细定位（prune 0.9）
  tables/52b_finemap_pair_r099_L2.csv / _L3.csv  L2/L3 精细定位（prune 0.99）
  tables/45c_smr_3p3_combined.csv                chr1 OneK1K-SMR
  tables/56_smr_L2L3_main.csv                    L2/L3 OneK1K-SMR
  tables/47_gtex8_wb_smr_results.csv             chr1 GTEx v8 全血 SMR
  tables/_s52b2_celltype_Neff.csv                各细胞类型真实 eQTL 有效 N

输出：tables/57_三位点比较表.csv
"""
import io, os, sys
import numpy as np
import pandas as pd

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
PROJ = r"D:\endometriosis_project\11_sc_eqtl_mr_project"
TAB = os.path.join(PROJ, "tables")
OUT = os.path.join(TAB, "57_三位点比较表.csv")
LOG = r"D:\endometriosis_project\_s57_three_locus.log"
L = []
def w(s=""):
    L.append(str(s))
    try: print(s)
    except UnicodeEncodeError: print(str(s).encode("ascii", "backslashreplace").decode())

def rd(n):
    p = os.path.join(TAB, n)
    return pd.read_csv(p) if os.path.exists(p) else None

# ---------- 0) 细胞类型真实 N ----------
ne = rd("_s52b2_celltype_Neff.csv")
ne["cell"] = ne["file"].str.extract(r"OneK1K_(.+?)\.cis_qtl_pairs")
NEFF = ne.groupby("cell")["n_eff_af_med"].median().round().astype(int).to_dict()

LOCUS_DEF = {
    "L1": "chr1:22.0-22.9 Mb (GRCh37) / 21.55-22.70 Mb (GRCh38)",
    "L2": "chr10:26.9-27.9 Mb (GRCh37) / 26.50-27.75 Mb (GRCh38)",
    "L3": "chr2:69.3-70.2 Mb (GRCh37) / 68.95-70.10 Mb (GRCh38)",
}

# ---------- 1) MR 主结果 ----------
mr = rd("31_final_MR_with_method.csv")
mr = mr[mr["method"] == "Wald ratio"].copy()

# ---------- 2) coloc.abf / susie ----------
cs = rd("34_coloc_summary.csv")

# ---------- 3) 敏感性矩阵（chr1 + L3）与 L2 ----------
s40 = rd("40_coloc_sensitivity_matrix.csv")
s51 = rd("51_coloc_sensitivity_L2.csv")
# s40 的 locus 列形态 = "<gene>_<cell>"；归一化出 gene/cell
s_all = pd.concat([s40, s51], ignore_index=True, sort=False)

# ---------- 4) 条件共定位 ----------
c43 = rd("43_cond_coloc_CDC42.csv")
c43n = rd("43c_cond_coloc_negctl.csv")
c54 = rd("54c_cond_coloc_summary_L2L3.csv")

# ---------- 5) 精细定位 ----------
f36 = rd("36_finemap_pair.csv")
F2 = {}
for tag, fn in [("L2|0.9", "52b_finemap_pair_L2.csv"), ("L2|0.99", "52b_finemap_pair_r099_L2.csv"),
                ("L3|0.9", "52b_finemap_pair_L3.csv"), ("L3|0.99", "52b_finemap_pair_r099_L3.csv")]:
    d = rd(fn)
    if d is not None:
        F2[tag] = d

# ---------- 6) SMR ----------
s45 = rd("45c_smr_3p3_combined.csv")
s56 = rd("56_smr_L2L3_main.csv")
s47 = rd("47_gtex8_wb_smr_results.csv")
s59 = rd("59_gtex8_L2L3_smr_results.csv")      # 裁定1：GTEx v8 全血 L2/L3 区域取数后 SMR+HEIDI
s45m = s45[s45["analysis"] == "main(5e-8)"].copy() if s45 is not None else None

rows = []
for _, r in cs.iterrows():
    L_id, gene, cell = r["locus_id"], r["gene"], r["cell_type"]
    o = dict(locus=L_id, locus_def=LOCUS_DEF.get(L_id, ""), gene=gene, cell_type=cell,
             celltype_eqtl_N_true=NEFF.get(cell, np.nan))

    # --- MR ---
    m = mr[(mr["locus_id"] == L_id) & (mr["gene"] == gene) & (mr["cell_type"] == cell)]
    if len(m):
        m = m.iloc[0]
        o.update(mr_method=m["method"], mr_nsnp=int(m["nsnp"]), mr_b=float(m["b"]),
                 mr_se=float(m["se"]), mr_p=float(m["p"]), mr_F=float(m["F_stat"]),
                 mr_instrument_grch38=str(m["instrument"]))
    else:
        o.update(mr_method=None, mr_nsnp=None, mr_b=None, mr_se=None, mr_p=None,
                 mr_F=None, mr_instrument_grch38=None)

    # --- coloc.abf / susie ---
    o.update(abf_H3=float(r["abf_H3"]), abf_H4=float(r["abf_H4"]),
             abf_verdict=r["abf_verdict"], abf_nsnp=int(r["nsnp_used_abf"]),
             susie_maxH4=float(r["susie_maxH4"]) if pd.notna(r["susie_maxH4"]) else None,
             susie_n_pairs=r["susie_n_pairs"], susie_verdict=r["susie_verdict"])

    # --- 敏感性（3 窗口 × 3 p12 = 9 组）---
    key = "%s_%s" % (gene, cell)
    ss = s_all[(s_all["gene"] == gene) & (s_all["cell_type"] == cell)]
    if len(ss):
        rob = ss["robust"].astype(str)
        o.update(sens_n=len(ss), sens_H4_min=float(ss["PP.H4"].min()),
                 sens_H4_max=float(ss["PP.H4"].max()),
                 sens_n_robust=int((rob == "稳健").sum()),
                 sens_n_notrobust=int((rob == "敏感性不稳健").sum()),
                 sens_layer=("稳健 9/9" if (rob == "稳健").sum() == 9 else
                             "稳健 %d/9" % (rob == "稳健").sum()),
                 sens_p12_min_H4=float(ss.loc[ss["p12"].astype(str) == "1e-6", "PP.H4"].min())
                 if (ss["p12"].astype(str) == "1e-6").any() else None)
    else:
        o.update(sens_n=None, sens_H4_min=None, sens_H4_max=None, sens_n_robust=None,
                 sens_n_notrobust=None, sens_layer="未做（仅 5 个 strong 对做了 45 组）",
                 sens_p12_min_H4=None)

    # --- 条件共定位 ---
    if L_id == "L1":
        cc = c43[(c43["scope"] == cell)] if c43 is not None else None
        if cc is not None and len(cc):
            lead = cc[cc["cond_snp"] == "rs56318008"]
            if len(lead):
                u = lead[lead["tag"] == "uncond"]
                a = lead[lead["tag"] == "cond_eqtl_only"]
                b = lead[lead["tag"] == "cond_both"]
                nn = c43n[(c43n["scope"] == cell) & (c43n["keep_rule"] == 0.01) &
                          (c43n["tag"] == "cond_both")]
                o.update(cond_done="是", cond_snp="rs56318008",
                         cond_H4_uncond=float(u["PP.H4"].iloc[0]) if len(u) else None,
                         cond_H4_eqtl_only=float(a["PP.H4"].iloc[0]) if len(a) else None,
                         cond_H4_both=float(b["PP.H4"].iloc[0]) if len(b) else None,
                         cond_H3_both=float(b["PP.H3"].iloc[0]) if len(b) else None,
                         cond_negctl_snp="rs12048511",
                         cond_negctl_H4_both=float(nn["PP.H4"].iloc[0]) if len(nn) else None)
            else:
                o.update(cond_done="是(仅 rs2473247 口径)", cond_snp="rs2473247",
                         cond_H4_uncond=None, cond_H4_eqtl_only=None, cond_H4_both=None,
                         cond_H3_both=None, cond_negctl_snp=None, cond_negctl_H4_both=None)
        else:
            o.update(cond_done="未做", cond_snp=None, cond_H4_uncond=None,
                     cond_H4_eqtl_only=None, cond_H4_both=None, cond_H3_both=None,
                     cond_negctl_snp=None, cond_negctl_H4_both=None)
    else:
        cc = c54[(c54["locus"] == L_id)] if c54 is not None else None
        if cc is not None and len(cc):
            k = cc[cc["keep_rule"] == 0.01].iloc[0]
            o.update(cond_done="是", cond_snp=k["cond_snp"],
                     cond_H4_uncond=float(k["H4_uncond"]),
                     cond_H4_eqtl_only=float(k["H4_cond_eqtl_only"]),
                     cond_H4_both=float(k["H4_cond_both"]),
                     cond_H3_both=float(k["H3_cond_both"]),
                     cond_negctl_snp=k["negctl_cond_snp"],
                     cond_negctl_H4_both=float(k["negctl_H4_cond_both"]))
        else:
            o.update(cond_done="未做", cond_snp=None, cond_H4_uncond=None,
                     cond_H4_eqtl_only=None, cond_H4_both=None, cond_H3_both=None,
                     cond_negctl_snp=None, cond_negctl_H4_both=None)

    # --- 精细定位 ---
    if L_id == "L1":
        sc = "%s_%s" % (gene, cell)
        d = f36[f36["scope"] == sc]
        if len(d):
            d = d.iloc[0]
            o.update(fm_prune="0.9（单一口径）",
                     fm_eqtl_converge=(d["converge_eqtl"] if pd.notna(d["converge_eqtl"]) else "见日志"),
                     fm_eqtl_n_cs=d["n_cs_eqtl"],
                     fm_eqtl_cs_lead="%s (%s, %s)" % (d["cs_lead_eqtl_rsid"], d["cs_lead_eqtl_pos37"],
                                                      d["cs_lead_eqtl_gene"]),
                     fm_eqtl_cs_lead_gene=d["cs_lead_eqtl_gene"],
                     fm_out_converge=(d["converge_outcome"] if pd.notna(d["converge_outcome"]) else ""),
                     fm_out_n_cs=d["n_cs_outcome"],
                     fm_out_cs_lead="%s (%s, %s)" % (d["cs_lead_outcome_rsid"], d["cs_lead_outcome_pos37"],
                                                     d["cs_lead_outcome_gene"]),
                     fm_out_cs_lead_gene=d["cs_lead_outcome_gene"],
                     fm_instr_in_set=d["mr_instrument_in_finemap_set"],
                     fm_r2_instr_cs_eqtl=d["r2_instr_vs_cs_lead_eqtl"])
        else:
            o.update(fm_prune="0.9（单一口径）", fm_eqtl_converge="未做", fm_eqtl_n_cs=None,
                     fm_eqtl_cs_lead=None, fm_eqtl_cs_lead_gene=None, fm_out_converge=None,
                     fm_out_n_cs=None, fm_out_cs_lead=None, fm_out_cs_lead_gene=None,
                     fm_instr_in_set=None, fm_r2_instr_cs_eqtl=None)
    else:
        parts = []
        for t in ["%s|0.9" % L_id, "%s|0.99" % L_id]:
            d = F2.get(t)
            if d is None:
                continue
            po = d[d["scope"] == "outcome_only"].iloc[0]
            pe = d[d["scope"] != "outcome_only"]
            if len(pe):
                pe = pe.iloc[0]
                parts.append(dict(prune=t.split("|")[1],
                                  e_cs=pe["cs_lead_eqtl_rsid"], e_pos=pe["cs_lead_eqtl_pos37"],
                                  e_gene=pe["cs_lead_eqtl_gene"], e_pip=pe["max_pip_eqtl"],
                                  o_cs=pe["cs_lead_outcome_rsid"], o_pos=pe["cs_lead_outcome_pos37"],
                                  o_gene=pe["cs_lead_outcome_gene"],
                                  instr=pe["mr_instrument_in_finemap_set"],
                                  r2=pe["r2_instr_vs_cs_lead_eqtl"],
                                  o_conv=po["converge"], o_ncs=po["n_cs"],
                                  o_arg=po["argmax_rsid"], o_argpos=po["argmax_pos37"],
                                  o_arggene=po["argmax_gene"], o_maxpip=po["max_pip"]))
        if parts:
            o.update(fm_prune="0.9 与 0.99 双口径",
                     fm_eqtl_converge="; ".join("prune%s: eqtlCS=%s(PIP %.3f, %s)｜outCS=%s(%s)" %
                                                (p["prune"], p["e_cs"], p["e_pip"], p["e_gene"],
                                                 p["o_cs"], p["o_gene"]) for p in parts),
                     fm_eqtl_n_cs="; ".join("prune%s=%s" % (p["prune"], p["o_ncs"]) for p in parts),
                     fm_eqtl_cs_lead="; ".join("prune%s: %s (%s, %s)" %
                                               (p["prune"], p["e_cs"], p["e_pos"], p["e_gene"]) for p in parts),
                     fm_eqtl_cs_lead_gene="; ".join("%s:%s" % (p["prune"], p["e_gene"]) for p in parts),
                     fm_out_converge="; ".join("prune%s: %s" % (p["prune"], p["o_conv"]) for p in parts),
                     fm_out_n_cs="; ".join("prune%s=%s" % (p["prune"], p["o_ncs"]) for p in parts),
                     fm_out_cs_lead="; ".join("prune%s: %s (%s, %s)" %
                                              (p["prune"], p["o_cs"], p["o_pos"], p["o_gene"]) for p in parts),
                     fm_out_cs_lead_gene="; ".join("%s:%s" % (p["prune"], p["o_gene"]) for p in parts),
                     fm_instr_in_set="; ".join("prune%s=%s" % (p["prune"], p["instr"]) for p in parts),
                     fm_r2_instr_cs_eqtl="; ".join("prune%s=%s" % (p["prune"], p["r2"]) for p in parts),
                     fm_stability=("稳定" if len(set(p["e_gene"] for p in parts)) == 1
                                   and len(set(p["o_gene"] for p in parts)) == 1
                                   else "对 prune 阈值敏感（信号归属不稳定）"))
        else:
            o.update(fm_prune=None, fm_eqtl_converge="未做", fm_eqtl_n_cs=None,
                     fm_eqtl_cs_lead=None, fm_eqtl_cs_lead_gene=None, fm_out_converge=None,
                     fm_out_n_cs=None, fm_out_cs_lead=None, fm_out_cs_lead_gene=None,
                     fm_instr_in_set=None, fm_r2_instr_cs_eqtl=None, fm_stability=None)

    # --- SMR（OneK1K 侧）---
    if L_id == "L1":
        d = s45m[(s45m["cell_type"] == cell) & (s45m["Gene"] == gene)] if s45m is not None else None
        if d is not None and len(d):
            d = d.iloc[0]
            o.update(smr_topSNP=d["topSNP"], smr_b=float(d["b_SMR"]), smr_p=float(d["p_SMR"]),
                     smr_heidi_p=float(d["p_HEIDI"]), smr_heidi_m=d["nsnp_HEIDI"],
                     smr_heidi_verdict=d["heidi_verdict"], smr_N_used=980)
        else:
            o.update(smr_topSNP=None, smr_b=None, smr_p=None, smr_heidi_p=None,
                     smr_heidi_m=None, smr_heidi_verdict="未分析（该细胞该基因 eQTL 未达 5e-8）",
                     smr_N_used=None)
    else:
        d = s56[s56["dataset"].str.contains(L_id)] if s56 is not None else None
        if d is not None and len(d):
            d = d.iloc[0]
            hp = float(d["p_HEIDI"]); m = d["nsnp_HEIDI"]
            o.update(smr_topSNP=d["topSNP"], smr_b=float(d["b_SMR"]), smr_p=float(d["p_SMR"]),
                     smr_heidi_p=hp, smr_heidi_m=m,
                     smr_heidi_verdict=("consistent" if hp >= 0.05 else "heterogeneity"),
                     smr_N_used=int(NEFF.get(cell, np.nan)))
        else:
            o.update(smr_topSNP=None, smr_b=None, smr_p=None, smr_heidi_p=None,
                     smr_heidi_m=None, smr_heidi_verdict="未分析", smr_N_used=None)

    # --- GTEx v8 全血 SMR（bulk；L1 本地已有数据，L2/L3 由裁定1 区域取数补做）---
    o.update(gtex8wb_topSNP=None, gtex8wb_b_SMR=None, gtex8wb_SMR_p=None,
             gtex8wb_HEIDI_p=None, gtex8wb_m=None, gtex8wb_verdict=None,
             gtex8wb_N=None, gtex8wb_note=None)
    if L_id == "L1" and s47 is not None:
        d = s47[(s47["Gene"] == gene) & (s47["analysis"] == "A(5e-8)")]
        note = "GTEx v8 全血 bulk n=670；主分析 --peqtl-smr 5e-8"
        if len(d) == 0:
            d = s47[(s47["Gene"] == gene) & (s47["analysis"] == "B(1e-4)")]
            note = "GTEx v8 全血 bulk n=670；放宽 --peqtl-smr 1e-4（该基因无 SNP 过 5e-8）"
        if len(d):
            d = d.iloc[0]
            hp = float(d["p_HEIDI"])
            o.update(gtex8wb_topSNP=d["topSNP"], gtex8wb_b_SMR=float(d["b_SMR"]),
                     gtex8wb_SMR_p=float(d["p_SMR"]), gtex8wb_HEIDI_p=hp,
                     gtex8wb_m=d["nsnp_HEIDI"], gtex8wb_N=670,
                     gtex8wb_verdict=("consistent" if hp >= 0.05 else "heterogeneity"),
                     gtex8wb_note=note)
        else:
            o.update(gtex8wb_note="GTEx v8 全血未检出该基因 eQTL（或未测）")
    elif s59 is not None:
        d = s59[(s59["Gene"] == gene) & (s59["analysis"] == "A(5e-8)")]
        if len(d):
            d = d.iloc[0]
            hp = float(d["p_HEIDI"])
            o.update(gtex8wb_topSNP=d["topSNP"], gtex8wb_b_SMR=float(d["b_SMR"]),
                     gtex8wb_SMR_p=float(d["p_SMR"]), gtex8wb_HEIDI_p=hp,
                     gtex8wb_m=d["nsnp_HEIDI"], gtex8wb_N=670,
                     gtex8wb_verdict=("consistent" if hp >= 0.05 else "heterogeneity"),
                     gtex8wb_note="GTEx v8 全血 bulk n=670；裁定1 区域定向取数（tabix-over-HTTP，"
                                  "取回 6.11MB=全文件 0.20%）后 cis-SMR+HEIDI；"
                                  "LD=同区域 OneK1K 980 供者子集（与 GTEx 供者不重叠，已知局限）")
        else:
            o.update(gtex8wb_note="GTEx v8 全血区域取数后未生成该基因 SMR 行")
    else:
        o.update(gtex8wb_note="未做：本地 GTEx v8 全血仅覆盖 chr1 三基因；"
                              "YME1L1/ANXA4 按裁定1 区域取数后补做")

    # --- N 口径留痕 ---
    o["N_used_in_sensitivity"] = 980
    o["N_true_celltype"] = NEFF.get(cell, np.nan)
    o["N_mismatch"] = (NEFF.get(cell, 980) != 980)
    rows.append(o)

DF = pd.DataFrame(rows)
DF.to_csv(OUT, index=False, encoding="utf-8-sig")
w("写出 %s（%d 行 × %d 列）" % (OUT, len(DF), DF.shape[1]))
w("")
w("列清单：")
for i, c in enumerate(DF.columns, 1):
    w("  %2d %s" % (i, c))
w("")
w("=== 核心字段速览 ===")
show = ["locus", "gene", "cell_type", "mr_b", "mr_p", "abf_H4", "abf_verdict",
        "susie_maxH4", "sens_layer", "cond_H4_uncond", "cond_H4_both", "cond_negctl_H4_both",
        "sen" if False else "smr_b", "smr_p", "smr_heidi_p", "smr_heidi_verdict",
        "gtex8wb_topSNP", "gtex8wb_b_SMR", "gtex8wb_SMR_p", "gtex8wb_HEIDI_p",
        "gtex8wb_m", "gtex8wb_verdict", "celltype_eqtl_N_true"]
show = [c for c in show if c in DF.columns]
with pd.option_context("display.width", 260, "display.max_columns", 60):
    w(DF[show].to_string(index=False))
io.open(LOG, "w", encoding="utf-8", newline="\n").write("\n".join(L) + "\n")
print("WROTE " + LOG)
