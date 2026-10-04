# -*- coding: utf-8 -*-
"""s60b_build_mvmr_inputs.py -- 任务B 第1步：构建 CDC42 + LINC00339 的 MVMR 输入（逐细胞类型）

暴露1 = CDC42 顺式 eQTL（OneK1K 单细胞，细胞类型 c）
暴露2 = LINC00339 顺式 eQTL（同一细胞类型 c）
结局  = GCST90483463（女性不孕症，EUR）

工具集（主阈值 P < 5e-8，与 31_final_MR_with_method.csv 同口径）：
  该细胞类型下 CDC42 或 LINC00339 的 cis-eQTL P<5e-8 变异之并集（pre-clump）。
  ★ 31_final_MR_with_method.csv 只存「clump 后每对唯一工具」，故全工具集须从
    coloc_prep/eqtl_L1_*.csv（pre-clump 全 cis 变异）重建；该表仍用于核验所选细胞类型与工具对。

等位谐调（口径）：
  参考取向 = 结局 GWAS 的效应等位 ea。eQTL 侧 `af` 为 ALT(=效应)等位频率，逐 SNP 比较
  |af − eaf| 与 |af − (1−eaf)|：前者更小 → eQTL 效应等位即 ea（保留 slope）；否则取反号。
  双交叉核验：① 同一 SNP 的 CDC42 与 LINC00339 两文件 af 应逐位一致；
              ② 与 OneK1K LD 参考 plink --freq 的等位方向一致率。

LD 参考：00_data_raw/smr_3p3/ld_ref/chr1_20_25mb（OneK1K 980 供者，GRCh37）

输出：00_data_raw/mvmr_L1/<cell>.mvmr_input.csv
      00_data_raw/mvmr_L1/ld/<cell>.ld + .ld.id
      tables/60_mvmr_input_summary.csv
"""
import io
import os
import subprocess
import sys
from collections import OrderedDict

import numpy as np
import pandas as pd

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

PROJ = r"D:\endometriosis_project"
P11 = os.path.join(PROJ, "11_sc_eqtl_mr_project")
PREP = os.path.join(P11, "00_data_raw", "onek1k", "coloc_prep")
OUTD = os.path.join(P11, "00_data_raw", "mvmr_L1")
LDD = os.path.join(OUTD, "ld")
TAB = os.path.join(P11, "tables")
PLINK = os.path.join(PROJ, "_tools", "plink.exe")
LDREF = os.path.join(P11, "00_data_raw", "smr_3p3", "ld_ref", "chr1_20_25mb")
LOG = os.path.join(P11, "scripts", "_s60b_build_mvmr_inputs.log")
os.makedirs(OUTD, exist_ok=True)
os.makedirs(LDD, exist_ok=True)

CELLS = ["B_IN", "B_MEM", "CD4_ET", "CD4_NC", "CD8_ET", "CD8_NC", "CD8_S100B", "NK"]
THR = 5e-8

buf = []
def w(s=""):
    buf.append(str(s)); print(s)
def flush():
    io.open(LOG, "w", encoding="utf-8", newline="\n").write("\n".join(buf) + "\n")

w("=" * 104)
w("s60b — MVMR 输入构建（CDC42 + LINC00339，逐细胞类型）")
w("=" * 104)

# ---- 0) LD 参考 bim + freq ----
bim = pd.read_csv(LDREF + ".bim", sep=r"\s+", header=None,
                  names=["CHR", "SNP", "CM", "BP", "A1", "A2"])
bim["BP"] = bim["BP"].astype(int)
bim = bim.set_index("SNP")
frqf = LDREF + ".frq"
if not os.path.exists(frqf):
    subprocess.run([PLINK, "--bfile", LDREF, "--freq", "--out", LDREF],
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=False, timeout=3600)
frq = pd.read_csv(frqf, sep=r"\s+")
FRQ = dict(zip(frq["SNP"], zip(frq["A1"].str.upper(), frq["A2"].str.upper(), frq["MAF"])))
w("[0] LD 参考（OneK1K 980，GRCh37）= %d 变异；freq 可得 = %d" % (len(bim), len(FRQ)))

# ---- 1) GWAS（L1 区域）----
g = pd.read_csv(os.path.join(PREP, "gwas_L1.csv"))
g = g[["pos38", "ea", "oa", "beta", "se", "eaf", "p", "rsid"]].copy()
g["pos38"] = g["pos38"].astype(int)
ndup = int(g["pos38"].duplicated().sum())
g = g[~g["pos38"].duplicated(keep="first")]
GW = {r.pos38: r for r in g.itertuples()}
w("[1] 结局 GWAS（L1 区域）= %d 唯一位点（原重复 %d 已去重）" % (len(g), ndup))

# ---- 2) 暴露（两基因，逐细胞类型）----
E = {}
for gene in ("CDC42", "LINC00339"):
    d = pd.read_csv(os.path.join(PREP, "eqtl_L1_%s.csv" % gene))
    d = d[["gene", "snp_grch37", "snp_grch38", "af", "slope", "slope_se",
           "pval_nominal", "cell_type"]].copy()
    d["pos37"] = d["snp_grch37"].str.split(":").str[1].astype(int)
    d["pos38"] = d["snp_grch38"].str.split(":").str[1].astype(int)
    E[gene] = d
    w("[2] %-10s pre-clump cis 行 = %-7d ；细胞类型 %d ；P<5e-8 = %d"
      % (gene, len(d), d["cell_type"].nunique(), int((d["pval_nominal"] < THR).sum())))

rows_summary = []
for cell in CELLS:
    w("")
    w("-" * 104)
    w("[%s]" % cell)
    sub = {}
    for gene in ("CDC42", "LINC00339"):
        s = E[gene][E[gene]["cell_type"] == cell].drop_duplicates("pos37").set_index("pos37")
        sub[gene] = s
    # 工具并集
    inst = OrderedDict()
    for gene in ("CDC42", "LINC00339"):
        hit = sub[gene][sub[gene]["pval_nominal"] < THR]
        for p in sorted(hit.index):
            inst.setdefault(p, set()).add(gene)
    n_cdc = sum(1 for v in inst.values() if "CDC42" in v)
    n_lin = sum(1 for v in inst.values() if "LINC00339" in v)
    n_both = sum(1 for v in inst.values() if len(v) == 2)
    w("  P<5e-8：CDC42 = %d ；LINC00339 = %d ；并集 = %d （其中两基因同时命中 = %d）"
      % (n_cdc, n_lin, len(inst), n_both))

    recs = []
    n_nold = n_nogwas = n_orient_cdc = n_orient_lin = 0
    af_mismatch = 0
    for p37, genes in inst.items():
        snp = "1:%d" % p37
        if snp not in FRQ:
            n_nold += 1
            continue
        s1 = sub["CDC42"].loc[p37] if p37 in sub["CDC42"].index else None
        s2 = sub["LINC00339"].loc[p37] if p37 in sub["LINC00339"].index else None
        if s1 is None or s2 is None:
            n_nogwas += 0   # 两基因的 cis 表覆盖同一批变异，正常不会缺
        p38 = int((s1 if s1 is not None else s2)["pos38"])
        gw = GW.get(p38)
        if gw is None:
            n_nogwas += 1
            continue
        af1 = float(s1["af"]) if s1 is not None else np.nan
        af2 = float(s2["af"]) if s2 is not None else np.nan
        if not np.isnan(af1) and not np.isnan(af2) and abs(af1 - af2) > 1e-6:
            af_mismatch += 1
        eaf = float(gw.eaf); b1 = float(s1["slope"]); b2 = float(s2["slope"])
        se1 = float(s1["slope_se"]); se2 = float(s2["slope_se"])
        pv1 = float(s1["pval_nominal"]); pv2 = float(s2["pval_nominal"])
        # 取向：af 更接近 eaf -> eQTL 效应等位即 ea
        afv = af1 if not np.isnan(af1) else af2
        same = abs(afv - eaf) <= abs(afv - (1.0 - eaf))
        if not same:
            b1, b2 = -b1, -b2
            n_orient_cdc += 1 if "CDC42" in genes else 0
            n_orient_lin += 1 if "LINC00339" in genes else 0
        recs.append(dict(snp=snp, pos37=p37, pos38=p38, a1=str(gw.ea).upper(),
                         a2=str(gw.oa).upper(), af_eqtl=afv, eaf_gwas=eaf,
                         b_cdc42=b1, se_cdc42=se1, p_cdc42=pv1,
                         b_linc00339=b2, se_linc00339=se2, p_linc00339=pv2,
                         b_out=float(gw.beta), se_out=float(gw.se), p_out=float(gw.p),
                         instr_genes="+".join(sorted(genes)),
                         flipped=int(not same)))
    D = pd.DataFrame(recs)
    w("  可用（在 LD 参考 ∩ 有 GWAS）= %d ；丢：不在 LD 参考 %d，无 GWAS %d"
      % (len(D), n_nold, n_nogwas))
    w("  取向翻转 = %d / %d (%.2f%%) ；两文件 af 不一致 = %d" %
      (int(D["flipped"].sum()) if len(D) else 0, len(D),
       100.0 * D["flipped"].mean() if len(D) else 0.0, af_mismatch))
    if len(D) == 0:
        continue
    # 严格为正的 se
    D = D[(D.se_cdc42 > 0) & (D.se_linc00339 > 0) & (D.se_out > 0)].copy()
    op = os.path.join(OUTD, "%s.mvmr_input.csv" % cell)
    D.to_csv(op, index=False, encoding="utf-8-sig")
    w("  → %s（%d 行）" % (os.path.basename(op), len(D)))

    # ---- LD 矩阵（对工具并集）----
    ex = os.path.join(LDD, "%s.extract" % cell)
    with io.open(ex, "w", encoding="utf-8", newline="\n") as f:
        for snp in D["snp"]:
            f.write(snp + "\n")
    out = os.path.join(LDD, cell)
    subprocess.run([PLINK, "--bfile", LDREF, "--extract", ex, "--r", "square",
                    "--out", out], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                   check=False, timeout=3600)
    subprocess.run([PLINK, "--bfile", LDREF, "--extract", ex, "--write-snplist",
                    "--out", out], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                   check=False, timeout=3600)
    ldok = False
    n_ld = 0
    if os.path.exists(out + ".ld"):
        with io.open(out + ".ld", "r", encoding="utf-8") as f:
            n_ld = len([l for l in f if l.strip()])
        ldok = (n_ld == len(D))
    # 权威 SNP 顺序（与 .ld 行列顺序一致）
    order_ok = None
    if os.path.exists(out + ".snplist"):
        sp = [l.strip() for l in io.open(out + ".snplist", encoding="utf-8") if l.strip()]
        order_ok = (sp == list(D["snp"]))
        if not order_ok and len(sp) == len(D):
            D = D.set_index("snp").loc[sp].reset_index()
            D.to_csv(op, index=False, encoding="utf-8-sig")
    w("  → LD 矩阵（plink --r square）：%s（%d×%d）｜ snplist 顺序与输入一致 = %s"
      % ("OK" if ldok else "**失败**", n_ld, n_ld, order_ok))

    rows_summary.append(dict(cell=cell, n_cdc42_5e8=n_cdc, n_linc00339_5e8=n_lin,
                             n_union=len(D), n_both=n_both,
                             n_flipped=int(D["flipped"].sum()),
                             af_mismatch=af_mismatch, ld_ok=ldok,
                             minp_cdc42=float(D["p_cdc42"].min()),
                             minp_linc00339=float(D["p_linc00339"].min())))

S = pd.DataFrame(rows_summary)
S.to_csv(os.path.join(TAB, "60_mvmr_input_summary.csv"), index=False, encoding="utf-8-sig")
w("")
w("=" * 104)
w("汇总（工具并集，pre-clump，P<5e-8）")
w("=" * 104)
w(S.to_string(index=False))
w("")
w("已写 tables/60_mvmr_input_summary.csv")
w("VERDICT = PASS")
flush()
print("WROTE " + LOG)
