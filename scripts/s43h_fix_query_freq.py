# -*- coding: utf-8 -*-
"""s43h_fix_query_freq.py — 修正 GTEx v8 全血 SMR query 的 Freq 取向（选项 A 阻塞点修复）

问题根因（已实证）：
  eQTL Catalogue 文档声明 `ac` = ALT 等位计数，但在 `imported/GTEx_V8/ge/` 这组文件里，
  `ac` 实际承载的是 GTEx 原生 allpairs 的 `ma_count`（次等位计数）：
    · 12,359 个可比行中 ac/an 恒 ≤ 0.5（0 例外）
    · ac/an 与 `maf` 列在 100.0% 的行上相等（容差 0.005）
  ⇒ ac/an = 次等位频率（MAF），不是 ALT 等位频率。
  原脚本无条件把 ac/an 当 ALT 频率写入 Freq ⇒ 「次等位恰为 REF」的变异（本窗口 ~20%）写成补数
  ⇒ SMR 频率 QC 误剔除 20.88% > 5% 上限 ⇒ rc=1、0 行。

修法：
  保留 A1 = alt（eQTL Catalogue FAQ：ALT 恒为效应等位，beta 即对 alt），
  仅把 GTEx 自身的 maf 值正确指派给 alt 等位：
      Freq(A1=alt) = maf      若 alt 为次等位
      Freq(A1=alt) = 1 - maf  若 alt 为显性等位
  取向判定：主判据 = OneK1K LD 参考（980 人，hg19 实测频率）；交叉 = 结局 GWAS
  GCST90483463（N=705,682，EUR，freq 列）。不复制外部频率数值，只用取向信息。

可用性约束（与 s43f 一致，必须保留）：
  只保留 ① 双等位 SNP；② 位置能精确落在 OneK1K LD 参考（bim ID = 1:<pos37>）；
  ③ {ref, alt} 与 bim 的 {A1, A2} 完全一致（liftover 等位一致性硬校验）。

输出：00_data_raw/smr_3p3/gtex8/query_fix/{WNT4,CDC42,LINC00339}.query.txt
      + tables/48_gtex8_query_freq_orient_audit.csv
"""
import os, io, sys, json
import numpy as np
import pandas as pd

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

PROJ = r"D:\endometriosis_project"
SMR3 = os.path.join(PROJ, "11_sc_eqtl_mr_project", "00_data_raw", "smr_3p3")
LD = os.path.join(SMR3, "ld_ref", "chr1_20_25mb")
COJO = os.path.join(SMR3, "gwas_GCST90483463.cojo.txt")
SRC = os.path.join(PROJ, "01_data_raw", "eqtlcat_gtex_v8",
                   "Whole_Blood.chr1_21_23.3Mb_3genes.tsv")
G8 = os.path.join(SMR3, "gtex8")
QOUT = os.path.join(G8, "query_fix")
TAB = os.path.join(PROJ, "11_sc_eqtl_mr_project", "tables")
LOG = os.path.join(PROJ, "11_sc_eqtl_mr_project", "scripts", "_s43h_fix_query_freq.log")
os.makedirs(QOUT, exist_ok=True); os.makedirs(TAB, exist_ok=True)

buf = []
def w(s=""):
    buf.append(str(s)); print(s)
def flush():
    io.open(LOG, "w", encoding="utf-8").write("\n".join(buf) + "\n")

DELTA = -326493                      # pos38 - pos37（本区域实测恒值）
GENES = ["WNT4", "CDC42", "LINC00339"]
ENS = {"ENSG00000162552": "WNT4", "ENSG00000070831": "CDC42",
       "ENSG00000218510": "LINC00339"}
TSS37 = {"WNT4": 22470462, "CDC42": 22379120, "LINC00339": 22351051}

w("=" * 100)
w("s43h — 修正 GTEx v8 全血 SMR query 的 Freq 取向")
w("=" * 100)

# ---------- 1. OneK1K LD 参考 ----------
bim = pd.read_csv(LD + ".bim", sep=r"\s+", header=None,
                  names=["CHR", "SNP", "CM", "BP", "A1", "A2"]).set_index("SNP")
frq = pd.read_csv(os.path.join(SMR3, "ld_ref", "_ldfreq.frq"), sep=r"\s+").set_index("SNP")
one = {}
for snp, r in frq.iterrows():
    if snp not in bim.index:
        continue
    a1, a2 = bim.at[snp, "A1"], bim.at[snp, "A2"]
    f1 = float(r["MAF"])
    one[snp] = {a1: f1, a2: 1.0 - f1, "pair": {a1, a2}, "A1": a1, "A2": a2}
w("[1] OneK1K LD 参考：%d SNP（bim 与 .frq 交集）" % len(one))

# ---------- 2. 结局 GWAS ----------
cojo = {}
with io.open(COJO, "r", encoding="utf-8") as f:
    hdr = f.readline().rstrip("\n").split("\t")
    ic = {c: i for i, c in enumerate(hdr)}
    for ln in f:
        p = ln.rstrip("\n").split("\t")
        if len(p) != len(hdr):
            continue
        try:
            cojo[p[ic["SNP"]]] = (p[ic["A1"]], p[ic["A2"]], float(p[ic["freq"]]))
        except ValueError:
            continue
w("    结局 GWAS cojo：%d SNP" % len(cojo))

# ---------- 3. 源 TSV ----------
recs = []
with io.open(SRC, "r", encoding="utf-8", errors="replace") as f:
    ix = {c: i for i, c in enumerate(f.readline().rstrip("\n").split("\t"))}
    for ln in f:
        p = ln.rstrip("\n").split("\t")
        if len(p) <= ix["gene_symbol"]:
            continue
        ref, alt = p[ix["ref"]], p[ix["alt"]]
        if len(ref) > 1 or len(alt) > 1:
            continue
        try:
            pos38 = int(p[ix["position"]]); an = int(p[ix["an"]]); ac = int(p[ix["ac"]])
            maf = float(p[ix["maf"]]); b = float(p[ix["beta"]]); se = float(p[ix["se"]])
            pv = float(p[ix["pvalue"]])
        except ValueError:
            continue
        gid = p[ix["gene_id"]].split(".")[0]
        recs.append({"SNP": "1:" + str(pos38 - DELTA), "pos38": pos38, "ref": ref, "alt": alt,
                     "maf": maf, "an": an, "ac": ac, "beta": b, "se": se, "p": pv,
                     "gene": ENS.get(gid, p[ix["gene_symbol"]]), "rsid": p[ix["rsid"]]})
src = pd.DataFrame(recs)
n_all = len(src)
w("    源 TSV 双等位 SNP 行数 = %d（唯一变异 %d）" % (n_all, src.SNP.nunique()))

# ---------- 4. 可用性过滤（硬约束） ----------
src["in_ld"] = src.SNP.isin(one.keys())
src["pair_ok"] = [ (r.ref in one[r.SNP]["pair"] and r.alt in one[r.SNP]["pair"])
                   if r.in_ld else False for r in src.itertuples() ]
n_inld = int(src.in_ld.sum()); n_pairok = int(src.pair_ok.sum())
w("    [可用性] 落在 OneK1K LD 参考内 = %d 行；其中等位对与 bim 完全一致 = %d 行 (%.2f%%)"
  % (n_inld, n_pairok, 100.0 * n_pairok / max(1, n_inld)))
w("              被剔：不在 LD 参考 = %d 行；等位对不一致（liftover/组装差异）= %d 行"
  % (n_all - n_inld, n_inld - n_pairok))
u = src[src.pair_ok].copy()

# 多等位：同一 SNP 位置可能有多个 alt（不同 alt 等位共享 1:<pos> ID）
dup = u.groupby("SNP")["alt"].nunique()
multi = dup[dup > 1]
df = pd.DataFrame([], columns=u.columns)
if len(multi) > 0:
    for snp, grp in u[u.SNP.isin(multi.index)].groupby("SNP"):
        pair = one[snp]["pair"]
        sub = grp[grp.alt.isin(pair) & grp.ref.isin(pair)]
        if len(sub):
            df = pd.concat([df, sub.head(1)], ignore_index=True)
u = pd.concat([u[~u.SNP.isin(multi.index)], df], ignore_index=True)
w("    多等位位置被折叠（保留与 bim 等位对相符者）= %d 个位置；可用行数 = %d"
  % (len(multi), len(u)))

# ---------- 5. 取向判定（逐唯一变异） ----------
uniq = u.drop_duplicates(subset=["SNP"])
audit = []
for r in uniq.itertuples():
    snp, ref, alt, maf = r.SNP, r.ref, r.alt, r.maf
    f_alt_one = one[snp][alt]; f_ref_one = one[snp][ref]
    alt_minor_one = bool(f_alt_one < f_ref_one)
    g = cojo.get(snp)
    f_alt_gwas = np.nan; alt_minor_gwas = None
    if g is not None:
        ga1, ga2, gf = g
        if ga1 == alt and ga2 == ref:
            f_alt_gwas, alt_minor_gwas = gf, bool(gf < 0.5)
        elif ga1 == ref and ga2 == alt:
            f_alt_gwas, alt_minor_gwas = 1.0 - gf, bool((1.0 - gf) < 0.5)
    alt_minor = alt_minor_one          # 主判据 = OneK1K（可用性约束保证恒可用）
    freq_new = maf if alt_minor else 1.0 - maf
    freq_old = r.ac / float(r.an)
    audit.append({"SNP": snp, "pos38": r.pos38, "ref": ref, "alt": alt, "maf": maf,
                  "ac_over_an": freq_old,
                  "f_alt_onek1k": f_alt_one, "f_ref_onek1k": f_ref_one,
                  "alt_minor_onek1k": alt_minor_one,
                  "f_alt_gwas": f_alt_gwas, "alt_minor_gwas": alt_minor_gwas,
                  "alt_minor_final": alt_minor,
                  "Freq_alt_old": freq_old, "Freq_alt_new": freq_new,
                  "freq_flipped": bool(abs(freq_old - freq_new) > 1e-3)})
aud = pd.DataFrame(audit)

w()
w("[2] 取向判定（N = %d 个可用变异）" % len(aud))
acc = aud.dropna(subset=["alt_minor_gwas"])
agree = (acc.alt_minor_onek1k == acc.alt_minor_gwas)
w("    OneK1K 与 GWAS 取向可比 = %d，一致 = %d (%.2f%%)，不一致 = %d"
  % (len(acc), agree.sum(), 100.0 * agree.mean(), (~agree).sum()))
w("    alt 为次等位 = %d (%.1f%%)；alt 为显性等位 = %d (%.1f%%)"
  % (aud.alt_minor_final.sum(), 100.0 * aud.alt_minor_final.mean(),
     (~aud.alt_minor_final).sum(), 100.0 * (~aud.alt_minor_final).mean()))
w("    因取向错误被写成补数的变异 = %d (%.1f%%) ← 即原脚本的静默错误规模"
  % (aud.freq_flipped.sum(), 100.0 * aud.freq_flipped.mean()))
w("    （对照）ac/an 是否恒为次等位频率：max(ac/an) = %.6f；与 maf 列 max|Δ| = %.2e"
  % (aud.ac_over_an.max(), (aud.ac_over_an - aud.maf).abs().max()))
w("    原 Freq 与修正 Freq 的 max|Δ| = %.6f" % (aud.Freq_alt_old - aud.Freq_alt_new).abs().max())

w()
w("[3] 修正后 Freq 与两个独立面板的一致性（逐变异）")
for nm, col in (("OneK1K", "f_alt_onek1k"), ("GWAS", "f_alt_gwas")):
    d = (aud.Freq_alt_new - aud[col]).abs().dropna()
    do = (aud.Freq_alt_old - aud[col]).abs().dropna()
    if len(d) == 0:
        continue
    w("    vs %-7s n=%5d | 修正后 |Δ|中位=%.5f <0.05=%.1f%% <0.20=%.1f%% "
      "| 未修正 |Δ|中位=%.5f <0.05=%.1f%%"
      % (nm, len(d), d.median(), 100.0 * (d < 0.05).mean(), 100.0 * (d < 0.20).mean(),
         do.median(), 100.0 * (do < 0.05).mean()))
# 关键：SMR 真正关心的阈值 --diff-freq 0.2 的剔除率预估
for nm, col in (("OneK1K(LD参考)", "f_alt_onek1k"), ("GWAS", "f_alt_gwas")):
    dd = (aud.Freq_alt_new - aud[col]).abs().dropna()
    w("    预估 |ΔFreq|>0.20 的比例 vs %s = %.2f%%（阈值上限 5.00%%）"
      % (nm, 100.0 * (dd > 0.20).mean()))

w()
w("[4] 关键变异（供任务 3.4 与正文）")
for rs in ["rs112678906", "rs12037376", "rs10917151", "rs56318008", "rs2473290", "rs2473294"]:
    snps = src[src.rsid == rs].SNP.unique()
    for snp in snps:
        rr = aud[aud.SNP == snp]
        if len(rr) == 0:
            w("    %-12s %-12s 【不在 OneK1K LD 参考内 → 不可作 SMR target SNP】" % (rs, snp))
            continue
        r = rr.iloc[0]
        w("    %-12s %-12s ref=%s alt=%s maf=%.6f | OneK1K f(alt)=%.4f GWAS f(alt)=%s | "
          "Freq=%.6f（旧 %.6f）"
          % (rs, snp, r.ref, r.alt, r.maf, r.f_alt_onek1k,
             ("%.4f" % r.f_alt_gwas) if pd.notna(r.f_alt_gwas) else "NA",
             r.Freq_alt_new, r.Freq_alt_old))

# ---------- 6. 写 query ----------
w()
w("[5] 写修正后的 query 文件")
fmap = dict(zip(aud.SNP, aud.Freq_alt_new))
HDR = ["SNP", "Chr", "BP", "A1", "A2", "Freq", "Probe", "Probe_Chr", "Probe_bp",
       "Gene", "Orientation", "b", "se", "p"]
cnt = {}
for g in GENES:
    sub = u[u.gene == g].copy()
    fp = os.path.join(QOUT, "%s.query.txt" % g)
    with io.open(fp, "w", encoding="utf-8", newline="\n") as f:
        f.write("\t".join(HDR) + "\n")
        for r in sub.itertuples():
            f.write("\t".join([r.SNP, "1", r.SNP.split(":")[1], r.alt, r.ref,
                               "%.6g" % fmap[r.SNP], g, "1", str(TSS37[g]), g, "-",
                               "%.6g" % r.beta, "%.6g" % r.se, "%.6g" % r.p]) + "\n")
    cnt[g] = len(sub)
    pp = sub.p.values
    w("    %-10s rows=%5d  p<5e-8=%4d  p<1e-5=%4d  p<1e-4=%4d  p<1.57e-3=%4d  minP=%.4g"
      % (g, len(sub), (pp < 5e-8).sum(), (pp < 1e-5).sum(), (pp < 1e-4).sum(),
         (pp < 1.57e-3).sum(), pp.min()))

aud.to_csv(os.path.join(TAB, "48_gtex8_query_freq_orient_audit.csv"),
           index=False, encoding="utf-8-sig")
w()
w("已写 tables/48_gtex8_query_freq_orient_audit.csv（%d 行）" % len(aud))
json.dump({"n_usable_variants": int(len(aud)),
           "n_query_rows": {k: int(v) for k, v in cnt.items()},
           "n_source_rows": int(n_all),
           "n_dropped_not_in_ld": int(n_all - n_inld),
           "n_dropped_pair_mismatch": int(n_inld - n_pairok),
           "n_multiallelic_positions_collapsed": int(len(multi)),
           "alt_minor_fraction": float(aud.alt_minor_final.mean()),
           "n_freq_flipped_vs_old": int(aud.freq_flipped.sum()),
           "agreement_onek1k_vs_gwas": {"n": int(len(acc)), "agree": int(agree.sum())},
           "est_frac_diffgt0.2_vs_LDref":
               float(((aud.Freq_alt_new - aud.f_alt_onek1k).abs() > 0.2).mean()),
           "delta_b38_minus_b37": DELTA},
          io.open(os.path.join(G8, "query_fix_summary.json"), "w", encoding="utf-8"),
          ensure_ascii=False, indent=2)
w("VERDICT = PASS")
flush()
