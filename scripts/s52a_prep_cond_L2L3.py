# -*- coding: utf-8 -*-
"""s52a_prep_cond_L2L3.py -- 任务A：L2/L3 条件共定位输入准备（两阶段）

★ 条件变异的选择规则（由 L1 定标确定）
    chr1 链使用 rs56318008 / 1:22470407 作为条件变异。
    定标检验表明它 **不等于** 「窗口内结局 GWAS p 最小者」（那是 rs61768001 / 1:22465820），
    而等于 **outcome-only susie_rss 精细定位的 argmax**（36_finemap_pair.csv: argmax_snp=rs56318008）。
    => 规则 R2：条件变异 = 本座 outcome-only 精细定位 argmax（由 s53_finemap_L2L3.R 产出）。

阶段 A（无需条件变异）：复刻 s40b harmonise，产出 m_<LOCUS>_<CT>.csv；打印窗口内其他有 cis-eQTL 的基因。
阶段 B（需条件变异）：plink --recode A 取 A1 剂量 -> 带符号 r / r2 vs 条件变异；判据 A/B 等位对齐校验。

用法：
    python s52a_prep_cond_L2L3.py                       # 只跑阶段 A
    python s52a_prep_cond_L2L3.py --cond L2=10:27443623,L3=2:69956721
结果一律落盘；只读；不下载。
"""
import io
import json
import os
import subprocess
import sys
import time

import numpy as np
import pandas as pd

PROJ = r"D:\endometriosis_project\11_sc_eqtl_mr_project"
PREP = os.path.join(PROJ, "00_data_raw", "onek1k", "coloc_prep")
PQ = os.path.join(PROJ, "00_data_raw", "onek1k", "parquet_subset")
COND = os.path.join(PROJ, "00_data_raw", "onek1k", "cond_coloc_L2L3")
PLINK = r"D:\endometriosis_project\_tools\plink.exe"
BFILE = r"D:\endometriosis_project\_onek1k_plink\plink_merged_980_donors"
LOG = r"D:\endometriosis_project\_s52a_prep_cond_L2L3.log"
JSON = os.path.join(COND, "_cond_snps.json")
os.makedirs(COND, exist_ok=True)

L = []


def p(s=""):
    s = str(s)
    L.append(s)
    try:
        print(s)
    except UnicodeEncodeError:
        print(s.encode("ascii", "backslashreplace").decode("ascii"))


T0 = time.time()

# ---- 条件变异来源 ----
COND_SPEC = {}
for i, a in enumerate(sys.argv):
    if a == "--cond" and i + 1 < len(sys.argv):
        for kv in sys.argv[i + 1].split(","):
            k, v = kv.split("=")
            COND_SPEC[k.strip()] = v.strip()
if not COND_SPEC and os.path.exists(JSON):
    COND_SPEC = json.load(io.open(JSON, encoding="utf-8"))
    p("[i] 从 %s 读取条件变异: %s" % (JSON, COND_SPEC))

p("=== s52a L2/L3 条件共定位输入准备（两阶段）===")
p("启动 %s" % time.strftime("%Y-%m-%d %H:%M:%S"))
p("条件变异（阶段B用）: %s" % (COND_SPEC if COND_SPEC else "（未指定，本次仅做阶段 A）"))

# ---- 每个细胞类型的 eQTL 有效样本量（2026-09-28 由 ma_count/af/2 实测，全部落在精确整数网格）----
NEFF_F = os.path.join(PROJ, "tables", "_s52b2_celltype_Neff.csv")
NEFF = {}
if os.path.exists(NEFF_F):
    _ne = pd.read_csv(NEFF_F)
    _ne["cell"] = _ne["file"].str.extract(r"OneK1K_(.+?)\.cis_qtl_pairs")
    NEFF = _ne.groupby("cell")["n_eff_af_med"].median().round().astype(int).to_dict()
N_PANEL = sum(1 for _ in io.open(BFILE + ".fam", encoding="utf-8"))  # LD 参考面板供者数
p("LD 参考面板 N = %d ；eQTL 细胞类型有效 N = %s"
  % (N_PANEL, {k: NEFF.get(k) for k in ("CD4_NC", "Mono_NC", "B_MEM")}))

# GRCh37 蛋白编码基因体（Ensembl GRCh37 REST /overlap/region，2026-09-28 实查）
GENES = {
    "L1": [("LINC00339", 22351681, 22357716), ("CDC42", 22379120, 22419437),
           ("WNT4", 22443798, 22470462)],
    "L2": [("PDSS1", 26986588, 27035727), ("ABI1", 27035522, 27150016),
           ("ANKRD26", 27280843, 27389421), ("YME1L1", 27399383, 27444195),
           ("MASTL", 27443753, 27475853), ("ACBD5", 27484146, 27531059),
           ("PTCHD3", 27687116, 27703297), ("RAB18", 27793197, 27831143)],
    "L3": [("ANTXR1", 69240310, 69476459), ("GFPT1", 69546905, 69614382),
           ("NFU1", 69622882, 69664760), ("AAK1", 69688532, 69901481),
           ("ANXA4", 69871557, 70053596), ("GMCL1", 70056774, 70108528),
           ("SNRNP27", 70120692, 70132707), ("MXD1", 70124820, 70170077),
           ("ASPRV1", 70187226, 70189397), ("PCBP1-AS1", 70189395, 70315978)],
}


def annot(locus, pos37):
    for g, s, e in GENES[locus]:
        if s <= pos37 <= e:
            return g
    best, bd, bside = None, None, None
    for g, s, e in GENES[locus]:
        d, side = (s - pos37, "flank5") if pos37 < s else (pos37 - e, "flank3")
        if bd is None or d < bd:
            bd, best, bside = d, g, side
    return "%s_%s_%dkb" % (bside, best, round(bd / 1000))


# ============================================================ 0) 定标：L1
p("")
p("=" * 78)
p("[0] 定标 —— 用 L1 验证条件变异规则")
g1 = pd.read_csv(os.path.join(PREP, "gwas_L1.csv"))
b1 = pd.read_csv(os.path.join(PREP, "bim_universe_L1.csv"))
e1f = pd.read_csv(os.path.join(PREP, "eqtl_L1_CDC42.csv"))
map1 = e1f[["snp_grch37", "snp_grch38"]].drop_duplicates()
map1["pos37"] = map1["snp_grch37"].str.split(":").str[1].astype(int)
map1["pos38b"] = map1["snp_grch38"].str.split(":").str[1].astype(int)
jb = b1.merge(map1, left_on="snp", right_on="snp_grch37", how="inner")
jg = g1.merge(jb[["snp", "pos37", "pos38b", "a1", "a2"]], left_on="pos38", right_on="pos38b",
              how="inner", suffixes=("", "_b"))
p("  L1 窗口内可对齐变异 = %d" % len(jg))
imin = int(jg["p"].idxmin())
r1row = jg.loc[imin]
p("  R1（窗口内 min GWAS p）= %s  GRCh37=%s  GRCh38=%s  p=%.4e"
  % (r1row["rsid"], r1row["snp"], r1row["pos38"], r1row["p"]))
fm = pd.read_csv(os.path.join(PROJ, "tables", "36_finemap_pair.csv"))
fm0 = fm[fm["scope"] == "outcome_only"]
frs, fpos = str(fm0.iloc[0]["argmax_snp"]), int(fm0.iloc[0]["argmax_pos37"])
p("  R2（outcome-only 精细定位 argmax）= %s  GRCh37=%d  gene=%s"
  % (frs, fpos, fm0.iloc[0]["argmax_gene"]))
p("  实际使用的 chr1 条件变异 = rs56318008 / 1:22470407")
p("  ==> R1 复现 = %s ; R2 复现 = %s"
  % ("是" if str(r1row["snp"]) == "1:22470407" else "否",
     "是" if fpos == 22470407 else "否"))
p("  ==> 结论：条件变异规则采用 R2（outcome-only 精细定位 argmax）")

# ============================================================ 1) 阶段 A
LOCI = [
    dict(locus="L2", gene="YME1L1", ct="CD4_NC", chrom="10",
         lo37=26900000, hi37=27900000, lead37=27443623, lead38=27154694),
    dict(locus="L3", gene="ANXA4", ct="Mono_NC", chrom="2",
         lo37=69300000, hi37=70200000, lead37=69956721, lead38=69729589),
]


def harmonise(locus, gene, ct):
    e_all = pd.read_csv(os.path.join(PREP, "eqtl_%s_%s.csv" % (locus, gene)))
    gw = pd.read_csv(os.path.join(PREP, "gwas_%s.csv" % locus))
    bimU = pd.read_csv(os.path.join(PREP, "bim_universe_%s.csv" % locus))
    gw = gw.copy()
    gw["ea"] = gw["ea"].astype(str).str.upper()
    gw["oa"] = gw["oa"].astype(str).str.upper()
    gw = gw[gw["pos38"].notna() & gw["beta"].notna() & gw["se"].notna() & (gw["se"] > 0)].copy()
    gw["pos38"] = gw["pos38"].astype(int)
    gw["allele_key"] = (gw["pos38"].astype(str) + "_" + np.minimum(gw["ea"], gw["oa"]) + "_"
                        + np.maximum(gw["ea"], gw["oa"]))
    gw = gw[~gw["allele_key"].duplicated()]
    e = e_all[e_all["cell_type"] == ct].copy()
    e = e.merge(bimU, left_on="snp_grch37", right_on="snp", how="left")
    e = e[e["a1"].notna() & e["a2"].notna()].copy()
    e["pos38"] = e["snp_grch38"].astype(str).str.split(":").str[1]
    e = e[e["pos38"].notna()].copy()
    e["pos38"] = e["pos38"].astype(int)
    e["a1"] = e["a1"].astype(str).str.upper()
    e["a2"] = e["a2"].astype(str).str.upper()
    e["allele_key"] = (e["pos38"].astype(str) + "_" + np.minimum(e["a1"], e["a2"]) + "_"
                       + np.maximum(e["a1"], e["a2"]))
    m = e.merge(gw, on="allele_key", suffixes=("_e", "_o"))
    if "pos38_e" in m.columns:
        m = m.rename(columns={"pos38_e": "pos38"})
        m = m.drop(columns=[c for c in ["pos38_o"] if c in m.columns])
    m = m[((m["ea"] == m["a1"]) & (m["oa"] == m["a2"])) |
          ((m["ea"] == m["a2"]) & (m["oa"] == m["a1"]))].copy()
    m["flip"] = (m["ea"] == m["a2"]) & (m["oa"] == m["a1"])
    m["palindromic"] = (np.minimum(m["a1"], m["a2"]) + np.maximum(m["a1"], m["a2"])).isin(["AT", "CG"])
    m["af_mism"] = (m["af"] - m["eaf"]).abs()
    m["af_flip"] = (m["af"] - (1 - m["eaf"])).abs()
    bad = (m["palindromic"]) & (np.minimum(m["af_mism"], m["af_flip"]) > 0.05) & \
          ((m["af_mism"] - m["af_flip"]).abs() < 0.05)
    n_bad = int(bad.sum())
    m = m[~bad].copy()
    fix = (m["palindromic"]) & (m["af_flip"] < m["af_mism"])
    m.loc[fix, "flip"] = ~m.loc[fix, "flip"]
    m["b2"] = np.where(m["flip"], -m["beta"], m["beta"])
    m["maf_e"] = np.minimum(m["af"], 1 - m["af"])
    m["maf_o"] = np.minimum(m["eaf"], 1 - m["eaf"])
    m = m[np.isfinite(m["maf_e"]) & np.isfinite(m["maf_o"]) & (m["maf_e"] > 0) & (m["maf_o"] > 0) &
          np.isfinite(m["slope"]) & np.isfinite(m["slope_se"]) & (m["slope_se"] > 0)].copy()
    m = m[m["rsid"].notna() & (m["rsid"].astype(str) != "")]
    m = m.sort_values("pos38").reset_index(drop=True)
    m["N_out"] = m["n_cases"] + m["n_controls"]
    m["s_out"] = m["n_cases"] / m["N_out"]
    m["snp_grch37"] = m["snp_grch37"].astype(str)
    return m, n_bad


p("")
p("=" * 78)
p("[1] 阶段 A：harmonise + 落盘 m_<LOCUS>_<CT>.csv")
m_tabs = {}
for t in LOCI:
    locus, gene, ct, chrom = t["locus"], t["gene"], t["ct"], t["chrom"]
    p("")
    p("--- [%s] %s / %s / chr%s  (窗口 GRCh37 chr%s:%d-%d)"
      % (locus, gene, ct, chrom, chrom, t["lo37"], t["hi37"]))
    m, n_bad = harmonise(locus, gene, ct)
    p("    harmonise: 交集=%d  flip=%d  palindromic=%d  dropped_af_ambiguous=%d"
      % (len(m), int(m["flip"].sum()), int(m["palindromic"].sum()), n_bad))
    leadid = "%s:%d" % (int(chrom), t["lead37"])
    p("    MR 工具 %s 在交集内 = %s" % (leadid, "是" if (m["snp_grch37"] == leadid).any() else "否"))
    f = os.path.join(COND, "m_%s_%s.csv" % (locus, ct))
    m.to_csv(f, index=False, encoding="utf-8-sig")
    p("    写入 %s rows=%d" % (os.path.basename(f), len(m)))
    m_tabs[locus] = m

# ============================================================ 2) 阶段 B
p("")
p("=" * 78)
if not COND_SPEC:
    p("[2] 阶段 B 跳过（未提供条件变异）。请先运行 s53_finemap_L2L3.R 取 outcome-only argmax，")
    p("    再执行：python s52a_prep_cond_L2L3.py --cond L2=<chr:pos37>,L3=<chr:pos37>")
else:
    p("[2] 阶段 B：plink --recode A + 带符号 r")
    summary_rows = []
    for t in LOCI:
        locus, chrom = t["locus"], t["chrom"]
        raw = os.path.join(COND, "ld_geno_%s.raw" % locus)
        if not os.path.exists(raw):
            cmd = [PLINK, "--bfile", BFILE, "--allow-no-sex", "--chr", chrom,
                   "--from-bp", str(t["lo37"]), "--to-bp", str(t["hi37"]),
                   "--recode", "A", "--out", os.path.join(COND, "ld_geno_%s" % locus),
                   "--memory", "1200"]
            p("    CMD: %s" % " ".join(cmd))
            subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                           check=False, timeout=3600)
        if not os.path.exists(raw):
            p("    !! %s .raw 缺失，跳过" % locus)
            continue

        cond37 = COND_SPEC.get(locus)
        m = m_tabs[locus]
        if cond37 is None:
            p("    !! 未给 %s 的条件变异，跳过" % locus)
            continue
        q = m[m["snp_grch37"] == cond37]
        if len(q) != 1:
            p("    !! 条件变异 %s 在 %s 交集内不唯一/不存在（n=%d），跳过" % (cond37, locus, len(q)))
            continue
        cond_rs, cond_rows = str(q.iloc[0]["rsid"]), q.iloc[0]
        p("")
        p("    --- [%s] 条件变异 = %s (%s) 注释=%s" % (locus, cond_rs, cond37, annot(locus, int(cond37.split(":")[1]))))
        p("        eQTL z_e=%+.3f (slope=%+.5f se=%.5f p=%.3e) | 结局 z_o=%+.3f (b2=%+.5f se=%.5f p=%.3e)"
          % (cond_rows["slope"] / cond_rows["slope_se"], cond_rows["slope"], cond_rows["slope_se"],
             cond_rows["pval_nominal"], cond_rows["b2"] / cond_rows["se"], cond_rows["b2"],
             cond_rows["se"], cond_rows["p"]))

        rw = pd.read_csv(raw, sep=r"\s+")
        snp_cols = [c for c in rw.columns
                    if c not in ("FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE")]
        cnt = {}
        for c in snp_cols:
            sid, al = c.rsplit("_", 1)
            cnt[sid] = al
        bim_ref = pd.read_csv(BFILE + ".bim", sep="\t", header=None, dtype={0: str},
                              names=["chr", "snp", "cm", "pos", "a1", "a2"])
        bim_ref = bim_ref[(bim_ref["chr"] == chrom) & (bim_ref["pos"] >= t["lo37"]) &
                          (bim_ref["pos"] <= t["hi37"])].copy()
        bim_ref["a1"] = bim_ref["a1"].str.upper(); bim_ref["a2"] = bim_ref["a2"].str.upper()
        map_a1 = dict(zip(bim_ref["snp"], bim_ref["a1"]))
        map_a2 = dict(zip(bim_ref["snp"], bim_ref["a2"]))
        pos_of = dict(zip(bim_ref["snp"], bim_ref["pos"]))
        ids = [c.rsplit("_", 1)[0] for c in snp_cols]
        flip = np.array([0 if cnt[s] == map_a1.get(s) else (1 if cnt[s] == map_a2.get(s) else -1)
                         for s in ids], dtype=np.int8)
        p("        recode 列 = %d （==bim A1: %d ; ==A2: %d ; 不可判定: %d）"
          % (len(ids), int((flip == 0).sum()), int((flip == 1).sum()), int((flip < 0).sum())))
        if (flip < 0).any():
            p("        !! 存在不可对齐列，跳过")
            continue
        X = rw[snp_cols].astype(np.float64).values.copy()
        X[:, np.where(flip == 1)[0]] = 2.0 - X[:, np.where(flip == 1)[0]]

        idx_of = {s: j for j, s in enumerate(ids)}
        rows = []
        for _, qq in m.iterrows():
            j = idx_of.get(str(qq["snp_grch37"]))
            if j is None:
                continue
            rows.append((np.nanmean(X[:, j]) / 2.0, float(qq["af"]), float(qq["eaf"])))
        chk = np.array(rows, dtype=float)
        eA = np.abs(chk[:, 0] - chk[:, 1])
        eB_align = np.abs(chk[:, 0] - chk[:, 2])
        eB_flip = np.abs(chk[:, 0] - (1.0 - chk[:, 2]))
        p("        [判据A·原口径] A1剂量频率 vs eQTL af: n=%d max=%.4f median=%.4f >0.05=%d"
          % (len(chk), eA.max(), np.median(eA), int((eA > 0.05).sum())))
        p("        [判据B·诊断] vs 结局 eaf: 对齐 mean|dev|=%.4f 翻转 mean|dev|=%.4f 比值=%.4f"
          % (eB_align.mean(), eB_flip.mean(), eB_align.mean() / eB_flip.mean()))

        # ---- 判据A′（统计口径，2026-09-28 起）----
        # eQTL 与 LD 面板的 N 可能不同（每细胞类型 N_eff 不同）。若 eQTL 样本 ⊂ 面板，
        #   Var(af_eQTL - af_panel) = p(1-p) * [1/(2*N_e) - 1/(2*N_panel)]
        # 该式是嵌套假设下的方差；若两者只是部分重叠，真实方差更大 => 本式给出 SD 下界
        # => |z| 上界 => 判据保守。
        ne = NEFF.get(t["ct"], N_PANEL)
        pe = np.clip(chk[:, 1], 1e-6, 1 - 1e-6)
        var_nest = pe * (1 - pe) * (1.0 / (2.0 * ne) - 1.0 / (2.0 * N_PANEL))
        p("        [判据A′·统计] N_eQTL=%d  N_panel=%d" % (ne, N_PANEL))
        if ne >= N_PANEL or var_nest.max() <= 0:
            p("            N_eQTL >= N_panel => 期望差为 0，要求逐变异精确一致")
            passA = bool(eA.max() <= 1e-4)
            p("            max|dev|=%.6f  （须 <= 1e-4）  => %s"
              % (eA.max(), "通过" if passA else "不通过"))
            zmax = np.nan; zmean = np.nan; zmed = np.nan
        else:
            with np.errstate(invalid="ignore", divide="ignore"):
                zabs = eA / np.sqrt(var_nest)
            zabs = zabs[np.isfinite(zabs)]
            zmax, zmean, zmed = float(zabs.max()), float(zabs.mean()), float(np.median(zabs))
            # 1300~2400 个（强相关）变异，Bonferroni 阈 ~3.9；取 6 作为宽松上限，
            # 并加均值/中位数齐性检查（若存在真实的样本集不重叠或系统性等位错误，
            # 偏差会整体抬高，而非仅个别变异超出）。
            passA = bool(zmax <= 6.0 and zmean <= 1.5 and zmed <= 1.0)
            p("            max|z|=%.3f  mean|z|=%.3f  median|z|=%.3f  "
              "（阈值 max<=6, mean<=1.5, median<=1.0） => %s"
              % (zmax, zmean, zmed, "通过" if passA else "不通过"))
        if (not passA) or eB_align.mean() >= eB_flip.mean():
            p("        !! 等位对齐/同源性判定失败，跳过")
            continue
        p("        -> 等位对齐通过（判据A′：eQTL 与 LD 面板频率差在 N 差可解释的抽样误差内）")

        sd = np.nanstd(X, axis=0, ddof=1)
        af_ld = np.nanmean(X, axis=0) / 2.0
        if cond37 not in ids:
            p("        !! 条件变异不在剂量矩阵中，跳过")
            continue
        jj = ids.index(cond37)
        xk = X[:, jj] - np.nanmean(X[:, jj])
        num = np.nansum((X - np.nanmean(X, axis=0)) * xk[:, None], axis=0)
        den = (len(xk) - 1) * sd * sd[jj]
        with np.errstate(invalid="ignore", divide="ignore"):
            r = num / den
        r[jj] = 1.0
        out = pd.DataFrame(dict(snp_grch37=ids, pos37=[pos_of.get(s) for s in ids],
                                a1=[map_a1.get(s) for s in ids], a2=[map_a2.get(s) for s in ids],
                                af_ld=af_ld, sd_a1=sd, r_to_cond=r, r2_to_cond=r ** 2))
        out["in_intersect"] = out["snp_grch37"].isin(set(m["snp_grch37"]))
        out["n_eqtl_neff"] = int(ne)
        out["n_ld_panel"] = int(N_PANEL)
        outf = os.path.join(COND, "r_to_%s_%s.csv" % (locus, cond_rs))
        out.to_csv(outf, index=False, encoding="utf-8-sig")
        p("        写入 %s rows=%d" % (os.path.basename(outf), len(out)))
        p("        |r|（不含条件变异）：max=%.4f  >0.99=%d  >0.9=%d  >0.8=%d"
          % (np.nanmax(np.abs(r[np.arange(len(r)) != jj])),
             int((np.abs(r) > 0.99).sum()) - 1, int((np.abs(r) > 0.9).sum()),
             int((np.abs(r) > 0.8).sum())))
        summary_rows.append(dict(locus=locus, gene=t["gene"], cell_type=t["ct"],
                                 cond_snp37=cond37, cond_rsid=cond_rs,
                                 cond_annot=annot(locus, int(cond37.split(":")[1])),
                                 n_eqtl_neff=int(ne), n_ld_panel=int(N_PANEL),
                                 max_dev_af=float(eA.max()),
                                 max_abs_z_neff=float(zmax) if np.isfinite(zmax) else None,
                                 judgeA2_pass=bool(passA),
                                 n_intersect=len(m), n_ld=len(ids),
                                 max_abs_r_excl=float(np.nanmax(np.abs(r[np.arange(len(r)) != jj])))))
    if summary_rows:
        pd.DataFrame(summary_rows).to_csv(os.path.join(COND, "_cond_input_summary.csv"),
                                          index=False, encoding="utf-8-sig")
        json.dump(COND_SPEC, io.open(JSON, "w", encoding="utf-8"))
        p("")
        p("    条件变异已记录至 %s" % JSON)

# ============================================================ 3) 窗口内其他基因
p("")
p("=" * 78)
p("[3] 窗口内 cis-eQTL 基因一览（parquet，同一细胞类型）")
gene_rows = []
for t in LOCI:
    f = os.path.join(PQ, "OneK1K_%s.cis_qtl_pairs.chr%s.parquet" % (t["ct"], t["chrom"]))
    df = pd.read_parquet(f, columns=["phenotype_id", "variant_id", "pval_nominal"])
    df["pos37"] = df["variant_id"].astype(str).str.split(":").str[1].astype(int)
    df = df[(df["pos37"] >= t["lo37"]) & (df["pos37"] <= t["hi37"])]
    g = df.groupby("phenotype_id").agg(minp=("pval_nominal", "min"), n=("pos37", "size")).reset_index()
    g = g.sort_values("minp")
    p("  [%s] %s，窗口内 cis-eQTL 基因数 = %d（列为 p 升序，前 25）"
      % (t["locus"], t["ct"], len(g)))
    for _, q in g.head(25).iterrows():
        sub = df[(df["phenotype_id"] == q["phenotype_id"]) & (df["pval_nominal"] == q["minp"])]
        lp = int(sub.iloc[0]["pos37"])
        p("     %-16s minp=%.3e n_var=%-5d lead37=%d (%s)%s"
          % (q["phenotype_id"], q["minp"], q["n"], lp, annot(t["locus"], lp),
             "  <== MR 基因" if q["phenotype_id"] == t["gene"] else ""))
    for _, q in g.iterrows():
        sub = df[(df["phenotype_id"] == q["phenotype_id"]) & (df["pval_nominal"] == q["minp"])]
        lp = int(sub.iloc[0]["pos37"])
        gene_rows.append(dict(locus=t["locus"], cell_type=t["ct"], gene=q["phenotype_id"],
                              min_p=q["minp"], n_var=q["n"], lead37=lp,
                              annot=annot(t["locus"], lp), is_target=(q["phenotype_id"] == t["gene"]),
                              p_lt_5e8=bool(q["minp"] < 5e-8)))
    del df
pd.DataFrame(gene_rows).to_csv(os.path.join(COND, "_window_genes.csv"),
                               index=False, encoding="utf-8-sig")

p("")
p("总耗时 %.1f s" % (time.time() - T0))
io.open(LOG, "w", encoding="utf-8", newline="\n").write("\n".join(L) + "\n")
print("WROTE " + LOG)
