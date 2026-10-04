# -*- coding: utf-8 -*-
"""s32_steiger_egger.py —— 任务二/四：Steiger 方向性检验 + MR-Egger

口径（用户裁定）：
- 对所有 n_index=1 的单 SNP 对：如实声明 MR-Egger 无法运行（工具数 <3，截距不可识别）；
  Steiger 方向性检验单工具在数学上可计算，但方向性结论稳健性有限，一并如实标注。
- 对 LINC00339/CD4_NC（2 SNP 多工具对）：Steiger 可算；MR-Egger 因仅 2 工具仍无法运行
  （截距+斜率两参数 vs 两数据点，完全饱和，无残差自由度）。

Steiger 公式（Hemani et al. 2017 近似，单 SNP）：
  R2_exposure = 2 * beta_exp^2 * eaf_exp * (1 - eaf_exp)
  R2_outcome  = 2 * beta_out^2  * eaf_out * (1 - eaf_out)
  方向性正确 iff R2_exposure > R2_outcome。
  （严格版需按 n 与表型方差缩放，此处用未缩放 R2 相对比较 + 标注近似性）

MR-Egger：工具数 >= 3 才可估计截距；本任务所有对均 <3 工具 → 全部如实声明无法运行。
"""
import os, sys, gzip
import numpy as np
import pandas as pd
from scipy import stats

T = r"D:\endometriosis_project\11_sc_eqtl_mr_project\tables"
RAW = r"D:\endometriosis_project\11_sc_eqtl_mr_project\00_data_raw"
LOG = r"D:\endometriosis_project\_s32_steiger.log"
L = []
def p(s=""):
    L.append(str(s)); print(s)

# ---------- 暴露侧：29b 逐变异 Wald（含 beta_exp/se_exp/af? 无 af，需从 01 表补） ----------
# 29b 表列：gene,cell_type,variant_id,effect_allele,other_allele,beta_exp,se_exp,F,beta_out,se_out,p_out,...
# 无 af 列 → 暴露 af 从 01b_discovery_instruments_grch38.csv 补（按 gene+cell_type+variant_id 匹配）
wd = pd.read_csv(os.path.join(T, "29b_variantlevel_wald.csv"), encoding="utf-8-sig")
wd = wd[wd["status"].isin(["ok", "flipped"])].copy()

inst = pd.read_csv(os.path.join(T, "01b_discovery_instruments_grch38.csv"), encoding="utf-8-sig")
# 01b 列含 af；variant_id 是 GRCh38（1:22113027 形态）。匹配键：gene, cell_type, variant_id
# 注意 01b 的 variant_id 可能带 .0 后缀，统一处理
inst["variant_id_clean"] = inst["variant_id"].astype(str).str.replace(r"\.0$", "", regex=True)
wd["variant_id_clean"] = wd["variant_id"].astype(str).str.replace(r"\.0$", "", regex=True)

afmap = inst.drop_duplicates(["cell_type", "gene", "variant_id_clean"]).set_index(
    ["cell_type", "gene", "variant_id_clean"])["af"].to_dict()

# 补充：敏感阈值工具集（29b_ivw_instruments_sens001.csv）含 af，覆盖多 SNP 对的第二个工具
sens_inst = pd.read_csv(os.path.join(T, "29b_ivw_instruments_sens001.csv"), encoding="utf-8-sig")
sens_inst["variant_id_clean"] = sens_inst["variant_id"].astype(str).str.replace(r"\.0$", "", regex=True)
for _, rr in sens_inst.iterrows():
    key = (rr["cell_type"], rr["gene"], rr["variant_id_clean"])
    if key not in afmap and pd.notna(rr.get("af")):
        afmap[key] = rr["af"]

def get_af_exp(r):
    return afmap.get((r["cell_type"], r["gene"], r["variant_id_clean"]), np.nan)

wd["af_exp"] = wd.apply(get_af_exp, axis=1)

p("=== 暴露侧载入 ===")
p(f"29b 有效对（ok/flipped）: {len(wd)} 行")
p(f"af_exp 可匹配: {wd['af_exp'].notna().sum()}/{len(wd)}")
p("")

# ---------- 结局侧：从缓存提取 eaf ----------
# GCST90483463.h.tsv.gz 列：chromosome, base_pair_location(GRCh38), effect_allele, other_allele,
#   beta, standard_error, effect_allele_frequency, p_value, rsid, ... variant_id(1_22000041_C_T 形态)
# 需要按 染色体:位置 匹配，并核对 effect_allele 方向
outfile = os.path.join(RAW, "gwas", "GCST90483463.h.tsv.gz")

# 结局 n（cases+controls）
N_OUT = 39724 + 650824  # 690548

# 建立 结局 索引：按 variant_id_clean（chr:pos）取 beta/se/eaf/effect_allele
# 29b 的 variant_id 是 GRCh38 "1:22027493"。结局缓存 base_pair_location 也是 GRCh38。
# 用 variant_id_clean 拆 chr:pos 匹配
def parse_chrpos(v):
    v = str(v)
    if ":" in v:
        c, pp = v.split(":")[:2]
        return int(c.replace("chr","")), int(pp)
    return None, None

outcome_cache = {}
need_pos = set()
for v in wd["variant_id_clean"]:
    c, pp = parse_chrpos(v)
    if c: need_pos.add((c, pp))

p(f"需从结局缓存提取的位点数: {len(need_pos)}")
found_out = 0
with gzip.open(outfile, "rt") as f:
    hdr = f.readline().rstrip("\n").split("\t")
    ix = {c: i for i, c in enumerate(hdr)}
    c_chr = ix["chromosome"]; c_pos = ix["base_pair_location"]; c_ea = ix["effect_allele"]
    c_oa = ix["other_allele"]; c_b = ix["beta"]; c_se = ix["standard_error"]
    c_eaf = ix["effect_allele_frequency"]; c_p = ix["p_value"]; c_rs = ix["rsid"]
    for line in f:
        cols = line.rstrip("\n").split("\t")
        ch = int(cols[c_chr]); pos = int(cols[c_pos])
        if (ch, pos) in need_pos:
            outcome_cache[(ch, pos)] = dict(
                ea=cols[c_ea], oa=cols[c_oa], beta=cols[c_b], se=cols[c_se],
                eaf=cols[c_eaf], p=cols[c_p], rs=cols[c_rs])
            found_out += 1
            if found_out >= len(need_pos):
                break
p(f"结局缓存命中: {found_out}/{len(need_pos)}")
p("")

# ---------- 等位基因对齐 + Steiger ----------
# 29b 表已有 effect_allele/other_allele（暴露 effect allele），且 beta_exp/beta_out 已按
# 暴露 effect allele 方向调和（flipped 标记表示结局 EA 与暴露 EA 相反，已翻 beta_out 符号）。
# 但结局 eaf 仍需对齐到暴露 effect allele 方向：
#   若结局缓存的 effect_allele == 暴露 effect_allele → eaf_out 直接用；
#   若结局 effect_allele == 暴露 other_allele（即方向相反）→ eaf_out = 1 - eaf。
rows = []
for _, r in wd.iterrows():
    c, pp = parse_chrpos(r["variant_id_clean"])
    oc = outcome_cache.get((c, pp))
    if oc is None:
        rows.append(dict(**r.to_dict(), eaf_out=np.nan, out_rs="", note="结局位点缺失"))
        continue
    ea_exp = str(r["effect_allele"]); oa_exp = str(r["other_allele"])
    ea_out = oc["ea"]; oa_out = oc["oa"]
    # 对齐 eaf
    eaf_out_raw = float(oc["eaf"])
    if ea_out == ea_exp:
        eaf_out = eaf_out_raw
        align = "same"
    elif ea_out == oa_exp:
        eaf_out = 1.0 - eaf_out_raw
        align = "flip"
    else:
        # 罕见：等位不匹配，标记
        eaf_out = np.nan
        align = "unmatched"
    # Steiger R2（未缩放，近似）
    af_e = r["af_exp"]
    be = r["beta_exp"]; bo = r["beta_out"]
    if (np.isfinite(af_e) and np.isfinite(eaf_out)
            and np.isfinite(be) and np.isfinite(bo)):
        r2e = 2 * be**2 * af_e * (1 - af_e)
        r2o = 2 * bo**2 * eaf_out * (1 - eaf_out)
    else:
        r2e = np.nan; r2o = np.nan
    rows.append(dict(**r.to_dict(), eaf_out=eaf_out, align=align,
                    r2_exposure=r2e, r2_outcome=r2o, out_rs=oc["rs"],
                    note=f"align={align}"))

res = pd.DataFrame(rows)

# 写逐 SNP Steiger 明细
res.to_csv(os.path.join(T, "37_steiger_variantlevel.csv"), index=False, encoding="utf-8-sig")

p("=== Steiger 逐变异明细（前 30 行关键列）===")
disp = res[["gene","cell_type","variant_id","af_exp","eaf_out","align",
            "beta_exp","beta_out","r2_exposure","r2_outcome"]]
p(disp.to_string(index=False))
p("")

# ---------- 汇总：方向性判定 ----------
p("=== Steiger 方向性汇总 ===")
summary = []
for (g, ct), sub in res.groupby(["gene","cell_type"]):
    n = len(sub)
    r2e = sub["r2_exposure"]; r2o = sub["r2_outcome"]
    if n == 1:
        # 单 SNP：方向性可算但稳健性有限
        if np.isfinite(r2e.iloc[0]) and np.isfinite(r2o.iloc[0]):
            direction_ok = "暴露>结局(方向正确)" if r2e.iloc[0] > r2o.iloc[0] else "暴露<结局(方向存疑)"
        else:
            direction_ok = "数据缺失，无法判定"
        summary.append(dict(gene=g, cell_type=ct, nsnp=n,
                            r2_exp=float(r2e.iloc[0]) if n else np.nan,
                            r2_out=float(r2o.iloc[0]) if n else np.nan,
                            steiger=direction_ok, egger="无法运行(工具数<3)"))
    else:
        # 多 SNP：Steiger 用平均 R2 比较
        r2e_mean = r2e.mean(); r2o_mean = r2o.mean()
        direction_ok = "暴露>结局(方向正确)" if r2e_mean > r2o_mean else "暴露<结局(方向存疑)"
        summary.append(dict(gene=g, cell_type=ct, nsnp=n,
                            r2_exp=float(r2e_mean), r2_out=float(r2o_mean),
                            steiger=direction_ok,
                            egger=f"无法运行(仅{n}工具，需>=3)"))

sdf = pd.DataFrame(summary)
sdf.to_csv(os.path.join(T, "38_steiger_egger_summary.csv"), index=False, encoding="utf-8-sig")
p(sdf.to_string(index=False))
p("")

# ---------- MR-Egger 如实声明 ----------
p("=== MR-Egger 结论 ===")
p("所有对工具数均 < 3：")
p("  - 15 个单 SNP 对：n_index=1，MR-Egger 截距不可识别，无法运行。")
p("  - LINC00339/CD4_NC：n_index=2，MR-Egger 回归两参数(截距+斜率) vs 两数据点完全饱和，无残差自由度，无法估计 SE，无法运行。")
p("  → 全部如实声明为「方法学限制：MR-Egger 无法运行」，不伪造任何 Egger 截距/斜率。")
p("")

with open(LOG, "w", encoding="utf-8") as f:
    f.write("\n".join(L))
print("日志已写:", LOG)
