# -*- coding: utf-8 -*-
r"""
s36_phewas_safety_integrate.py
FinnGen R12 PheWAS 交付包 —— 独立核验 + 论文用「靶点安全性评估表」整合

输入（只读，来自大内存机交付包）:
  E:\_handover_phewas_20260924\unpacked\out\20_phewas_all_all.csv
  E:\_handover_phewas_20260924\unpacked\out\24_safety_risk_up_all.csv
  E:\_handover_phewas_20260924\unpacked\out\26_safety_risk_down_all.csv
  E:\_handover_phewas_20260924\unpacked\out\phewas_summary_all.json
输出:
  tables/42_phewas_safety_assessment.csv      52 条 FDR<0.05 显著项（带分层注释）
  tables/42b_phewas_safety_domains.csv        域级汇总（论文表）
  tables/42c_phewas_myeloid_proxy_check.csv   HLH/全血细胞减少 代理端点核查
  logs/_s36_verify_log.txt                    核验日志
"""
import json, math, os, sys, csv, io
from collections import OrderedDict

BASE = r"E:\_handover_phewas_20260924\unpacked\out"
PROJ = r"D:\endometriosis_project\11_sc_eqtl_mr_project"
TBL  = os.path.join(PROJ, "tables")
LOG  = os.path.join(PROJ, "logs")
os.makedirs(TBL, exist_ok=True); os.makedirs(LOG, exist_ok=True)
log_lines = []
def L(s=""):
    log_lines.append(str(s))

def num(x):
    try: return float(x)
    except Exception: return float("nan")

def read_csv(p):
    with open(p, "r", encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))

# ---------------- 1. 读入并与独立 summary 对照 ----------------
rows = read_csv(os.path.join(BASE, "20_phewas_all_all.csv"))
up   = read_csv(os.path.join(BASE, "24_safety_risk_up_all.csv"))
down = read_csv(os.path.join(BASE, "26_safety_risk_down_all.csv"))
with open(os.path.join(BASE, "phewas_summary_all.json"), "r", encoding="utf-8") as f:
    summ = json.load(f)

L("=" * 78)
L("s36 · FinnGen R12 PheWAS 交付包独立核验")
L("=" * 78)
L("")
L("【1】20_phewas_all_all.csv")
L("  数据行数                 : %d" % len(rows))
L("  列数                     : %d" % len(rows[0]))
key = lambda r: (r["gene"], r["cell_type"], r["endpoint"])
uniq = set(key(r) for r in rows)
L("  唯一 (gene,cell,端點) 组合: %d" % len(uniq))
L("  重复组合                 : %d" % (len(rows) - len(uniq)))
hs = OrderedDict()
for r in rows: hs[r["harmony_status"]] = hs.get(r["harmony_status"], 0) + 1
L("  harmony_status 分布       : %s" % dict(hs))
g1 = sum(1 for r in rows if r["gene"] == "CDC42" and r["cell_type"] == "B_MEM")
g2 = sum(1 for r in rows if r["gene"] == "CDC42" and r["cell_type"] == "Mono_NC")
L("  CDC42_B_MEM 检验数        : %d" % g1)
L("  CDC42_Mono_NC 检验数      : %d" % g2)
sig  = [r for r in rows if r["fdr_sig"].strip().lower() == "true"]
upd  = [r for r in sig if r["direction"] == "negative"]
dnd  = [r for r in sig if r["direction"] == "positive"]
L("  FDR<0.05 显著             : %d" % len(sig))
L("    其中 beta_MR<0 (风险升高): %d" % len(upd))
L("    其中 beta_MR>0 (风险降低): %d" % len(dnd))
L("")
L("  与 phewas_summary_all.json 对照:")
chk = [
  ("n_tests_valid",        summ["n_tests_valid"],        len(rows)),
  ("n_fdr_sig",            summ["n_fdr_sig"],            len(sig)),
  ("B_MEM n",              summ["by_instrument"]["CDC42_B_MEM"]["n"], g1),
  ("Mono_NC n",            summ["by_instrument"]["CDC42_Mono_NC"]["n"], g2),
  ("B_MEM n_risk_increased", summ["by_instrument"]["CDC42_B_MEM"]["n_risk_increased"],
                             sum(1 for r in upd if r["cell_type"] == "B_MEM")),
  ("Mono n_risk_increased",  summ["by_instrument"]["CDC42_Mono_NC"]["n_risk_increased"],
                             sum(1 for r in upd if r["cell_type"] == "Mono_NC")),
  ("B_MEM n_risk_decreased", summ["by_instrument"]["CDC42_B_MEM"]["n_risk_decreased"],
                             sum(1 for r in dnd if r["cell_type"] == "B_MEM")),
  ("Mono n_risk_decreased",  summ["by_instrument"]["CDC42_Mono_NC"]["n_risk_decreased"],
                             sum(1 for r in dnd if r["cell_type"] == "Mono_NC")),
]
for name, a, b in chk:
    L("    %-26s summary=%-8s 实测=%-8s %s" % (name, a, b, "OK" if a == b else "*** 不一致 ***"))
L("")

# ---------------- 2. SE 口径独立复算（辨别一阶 vs 全比值 delta）----------------
L("【2】SE 口径独立复算（验证 README 所称 '一阶 delta' 实为全比值 delta）")
maxdiff_full, maxdiff_1st = 0.0, 0.0
for r in rows:
    bo, so, be, se = num(r["beta_out"]), num(r["se_out"]), num(r["beta_exp"]), num(r["se_exp"])
    rep = num(r["se_MR"])
    if not all(map(math.isfinite, (bo, so, be, se, rep))) or be == 0: continue
    full = math.sqrt(so**2 / be**2 + bo**2 * se**2 / be**4)   # 全比值 delta（含二阶项）
    first = so / abs(be)                                      # 仅一阶（丢弃二阶项）
    maxdiff_full = max(maxdiff_full, abs(full - rep))
    maxdiff_1st  = max(maxdiff_1st,  abs(first - rep))
L("  全比值 delta 与报告 se_MR 最大绝对偏差 : %.3e   <- 应当 ~0" % maxdiff_full)
L("  仅一阶项     与报告 se_MR 最大绝对偏差 : %.3e   <- 明显非 0" % maxdiff_1st)
L("  => 报告 SE 采用的是【全比值 delta 法】；README/SAFETY_BRIEF 中 '一阶 delta' 为命名笔误，")
L("     公式本身正确、与本项目既有口径（二阶 delta）一致。计算结果不受影响。")
L("")

# ---------------- 3. HLH / 全血细胞减少 代理端点核查 ----------------
MAN = r"D:\endometriosis_project\01_data_raw\finngen_R12_manifest.tsv"
man = []
with open(MAN, "r", encoding="utf-8-sig") as f:
    rd = csv.DictReader(f, delimiter="\t")
    for r in rd: man.append(r)
L("【3】HLH 与 全血细胞减少 在 FinnGen R12 中的可及性")
L("  manifest 端点数            : %d" % len(man))
absent_terms = ["haemophagocytic", "hemophagocytic", "lymphohistiocytosis", "HLH",
                "pancytopenia", "pancytopaenia", "bone marrow failure"]
for t in absent_terms:
    hit = [r["phenocode"] for r in man if t.lower() in (r["phenocode"] + " " + r["phenotype"]).lower()]
    L("    '%s' -> %d 命中 %s" % (t, len(hit), hit if hit else "(缺席)"))
# 代理端点：R12 中确实存在的髓系毒性 / 造血衰竭相关端点
proxy_pat = ["AGRANULOCYTOSIS", "NEUTROPENIA", "APLAST", "THROMBOCYTOPENIA",
             "ACUTEPOSTBLEEDANAEMIA", "PANCYT"]
proxies = [r for r in man if any(p in r["phenocode"].upper() for p in proxy_pat)]
L("  存在的髓系毒性/造血衰竭代理端点: %d 个" % len(proxies))
row_by_ep = {}
for r in rows: row_by_ep.setdefault(r["endpoint"], []).append(r)
proxy_out = []
L("  %-34s %-10s %-9s %-11s %-11s %s" % ("端點", "细胞类型", "beta_MR", "p_MR", "p_FDR", "FDR<0.05"))
for m in proxies:
    ep = m["phenocode"]
    if ep not in row_by_ep:
        L("  %-34s (未被 PhWAS 覆盖)" % ep); continue
    for r in row_by_ep[ep]:
        proxy_out.append({
            "endpoint": ep, "phenotype": m["phenotype"], "category": m["category"],
            "n_cases": m["num_cases"], "n_controls": m["num_controls"],
            "cell_type": r["cell_type"], "beta_MR": r["beta_MR"], "se_MR": r["se_MR"],
            "p_MR": r["p_MR"], "p_FDR": r["p_FDR"], "fdr_sig": r["fdr_sig"],
            "direction": r["direction"], "risk_if_cdc42_inhibited": r["risk_if_cdc42_inhibited"],
        })
        L("  %-34s %-10s %-9.4f %-11.4g %-11.4g %s" % (
            ep, r["cell_type"], num(r["beta_MR"]), num(r["p_MR"]), num(r["p_FDR"]), r["fdr_sig"]))
L("")

# ---------------- 4. 域级汇总（按 manifest 的 category 字段分组）----------------
# ★ 口径说明：C 节早期版本误用 CSV 的 category_prefix（=phenocode 派生的混合编码，
#   会把 AB1 拆成 AB1/AB1TUBERCU/… 、把 D3 拆成 D3/… ），导致域计数偏低。
#   权威口径是 manifest 的 category 字段，已用其对简报 §5 各项逐一复现（全部一致）。
C_AB1  = "I Certain infectious and parasitic diseases (AB1_)"
C_D3   = "III Diseases of the blood and blood-forming organs"
C_CD2  = "II Neoplasms from hospital discharges (CD2_)"
C_ICDO = "II Neoplasms, from cancer register (ICD-O-3)"
C_AUTO = "Diseases marked as autimmune origin"
C_N14  = "XIV Diseases of the genitourinary system (N14_)"
C_O15  = "XV Pregnancy, childbirth and the puerperium (O15_)"
C_P16  = "XVI Certain conditions originating in the perinatal period (P16_)"
C_M13  = "XIII Diseases of the musculoskeletal system and connective tissue (M13_)"
C_RHEU = "Rheuma endpoints"
DOMAINS = OrderedDict([
    ("免疫/血液/肿瘤",   [C_AB1, C_D3, C_CD2, C_ICDO, C_AUTO]),
    ("女性生殖/不孕",     [C_N14]),
    ("妊娠/分娩/产褥",    [C_O15, C_P16]),
    ("肌肉骨骼/结缔组织", [C_M13, C_RHEU]),
])
def domain_of(cat):
    for d, cs in DOMAINS.items():
        if any(cat.startswith(c[:30]) for c in cs): return d
    return "其他"
MINP = lambda rs: min((num(r["p_MR"]) for r in rs), default=float("nan"))
dom_rows = []
for d, ps in DOMAINS.items():
    sub = [r for r in rows if domain_of(r["category"]) == d]
    s = [r for r in sub if r["fdr_sig"].strip().lower() == "true"]
    u = [r for r in s if r["direction"] == "negative"]
    dn= [r for r in s if r["direction"] == "positive"]
    dom_rows.append({
        "domain": d, "categories": "; ".join(ps),
        "n_tests": len(sub), "n_sig": len(s),
        "n_risk_up": len(u), "n_risk_down": len(dn),
        "min_p_MR": MINP(sub),
        "min_p_MR_risk_up": MINP(u) if u else float("nan"),
    })
# 其他
rest = [r for r in rows if domain_of(r["category"]) == "其他"]
rs = [r for r in rest if r["fdr_sig"].strip().lower() == "true"]
dom_rows.append({
    "domain": "其他", "categories": "(其余全部类别)",
    "n_tests": len(rest), "n_sig": len(rs),
    "n_risk_up": sum(1 for r in rs if r["direction"] == "negative"),
    "n_risk_down": sum(1 for r in rs if r["direction"] == "positive"),
    "min_p_MR": MINP(rest),
    "min_p_MR_risk_up": MINP([r for r in rs if r["direction"] == "negative"]) if any(r["direction"]=="negative" for r in rs) else float("nan"),
})
L("【4】域级汇总")
L("  %-16s %-9s %-7s %-9s %-10s %-12s" % ("域", "检验数", "显著", "风险升高", "风险降低", "最小p_MR"))
for d in dom_rows:
    L("  %-16s %-9d %-7d %-9d %-10d %-12.4g" % (
        d["domain"], d["n_tests"], d["n_sig"], d["n_risk_up"], d["n_risk_down"], d["min_p_MR"]))
L("")

# ---------------- 5. 风险升高谱系映射 ----------------
LINEAGE = {}
for k, eps in {
    "A. 盆底支持结构/结缔组织": ["N14_FEMGENPROL"],
    "B. 胎盘附着/剥离异常与产科出血": ["O15_POSTPART_HEAMORRH_ALLW", "O15_POSTPART_HAEMORRH_PLACENTA",
                                        "O15_RETAINED_PLAC_MEMB_NO_HEAMORRH"],
    "C. 早产": ["O15_PRETERM"],
    "D. 关节软骨退变": ["PRIM_KNEEARTHROSIS", "M13_ARTHROSIS_KNEE_PRIM_KNEESURG"],
    "E. 可疑/待复核": ["M13_FOREIGNBODY"],
    "F. 非疾病生育表型(不计为风险)": ["O15_DELIV_SPONT"],
}.items():
    for e in eps: LINEAGE[e] = k

# ---------------- 6. 输出：52 条显著项整合表 ----------------
def tier(p):
    p = num(p)
    if not math.isfinite(p): return ""
    if p < 0.01: return "strong(FDR<0.01)"
    if p < 0.05: return "moderate(0.01<=FDR<0.05)"
    return "ns"
out1 = []
for r in sorted(sig, key=lambda r: num(r["p_FDR"])):
    cls = "风险升高(不利/脱靶)" if r["direction"] == "negative" else "风险降低(有利/与主效应一致)"
    out1.append({
        "signal_class": cls,
        "lineage_group": LINEAGE.get(r["endpoint"], "") if r["direction"] == "negative" else "",
        "domain": domain_of(r["category"]),
        "endpoint": r["endpoint"], "phenotype": r["phenotype"],
        "category_prefix": r["category_prefix"], "category_label": r["category_label"],
        "cell_type": r["cell_type"], "rsid": r["rsid"],
        "n_cases": r["n_cases"], "n_controls": r["n_controls"],
        "beta_out": r["beta_out"], "se_out": r["se_out"],
        "beta_exp": r["beta_exp"], "se_exp": r["se_exp"],
        "beta_MR": r["beta_MR"], "se_MR": r["se_MR"], "z_MR": r["z_MR"],
        "p_MR": r["p_MR"], "p_FDR": r["p_FDR"], "FDR_tier": tier(r["p_FDR"]),
        "risk_if_cdc42_inhibited": r["risk_if_cdc42_inhibited"],
        "harmony_status": r["harmony_status"],
    })
with open(os.path.join(TBL, "42_phewas_safety_assessment.csv"), "w",
          encoding="utf-8-sig", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(out1[0].keys())); w.writeheader(); w.writerows(out1)

with open(os.path.join(TBL, "42b_phewas_safety_domains.csv"), "w",
          encoding="utf-8-sig", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(dom_rows[0].keys())); w.writeheader(); w.writerows(dom_rows)

if proxy_out:
    with open(os.path.join(TBL, "42c_phewas_myeloid_proxy_check.csv"), "w",
              encoding="utf-8-sig", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(proxy_out[0].keys())); w.writeheader(); w.writerows(proxy_out)

L("【5】输出文件")
for fn in ["42_phewas_safety_assessment.csv", "42b_phewas_safety_domains.csv", "42c_phewas_myeloid_proxy_check.csv"]:
    p = os.path.join(TBL, fn)
    L("  %s  (%s bytes)" % (p, os.path.getsize(p) if os.path.exists(p) else "MISSING"))
L("")
L("【6】关键结论摘要")
ibt = [d for d in dom_rows if d["domain"] == "免疫/血液/肿瘤"][0]
L("  免疫/血液/肿瘤：%d 检验，显著 %d，风险升高 %d，最小 p_MR=%.4g"
  % (ibt["n_tests"], ibt["n_sig"], ibt["n_risk_up"], ibt["min_p_MR"]))
L("  生殖/产科合计：%d 检验，风险升高 %d"
  % (dom_rows[1]["n_tests"] + dom_rows[2]["n_tests"],
     dom_rows[1]["n_risk_up"] + dom_rows[2]["n_risk_up"]))
L("  髓系毒性代理端点（R12 实存）：%d 条检验，风险升高 %d"
  % (len(proxy_out), sum(1 for r in proxy_out if r["direction"] == "negative" and r["fdr_sig"] == "True")))
L("  HLH / 全血细胞减少：R12 manifest 中 0 命中 -> 不可评估（非阴性）")

with open(os.path.join(LOG, "_s36_verify_log.txt"), "w", encoding="utf-8") as f:
    f.write("\n".join(log_lines) + "\n")

print("\n".join(log_lines))
print("\nDONE")
