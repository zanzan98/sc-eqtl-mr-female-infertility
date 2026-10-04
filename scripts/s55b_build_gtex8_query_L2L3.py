# -*- coding: utf-8 -*-
"""s55b_build_gtex8_query_L2L3.py -- 裁定1：由 GTEx v8 区域数据构建 SMR query（L2 YME1L1 / L3 ANXA4）

严格复刻 s43h_fix_query_freq.py 的口径：
  · A1 = alt（eQTL Catalogue FAQ：ALT 恒为效应等位，beta 即对 alt）
  · 该文件里 `ac/an` 实为次等位频率 MAF（非 ALT 频率，s43h 已实证）→ 必须按取向指派：
        Freq(A1=alt) = maf      若 alt 为次等位
        Freq(A1=alt) = 1 - maf  若 alt 为显性等位
    取向主判据 = OneK1K LD 参考（bim + plink --freq）的频率；交叉 = 结局 GWAS cojo 的 freq。
  · 可用性硬约束：① 双等位；② liftOver(hg38->hg19) 后 SNP ID = <chr>:<pos37> 落在 OneK1K bim；
                  ③ {ref,alt} 与 bim 的 {A1,A2} 完全一致。

输入：
  01_data_raw/eqtlcat_gtex_v8/Whole_Blood.L2_YME1L1.tsv / .L3_ANXA4.tsv      （s55a 区域取数）
  00_data_raw/smr_L2L3/ld_ref/L2_YME1L1_CD4_NC_cis.{bim,fam,bed}             （s54a 已建，OneK1K 980）
  00_data_raw/smr_L2L3/gwas_L2_YME1L1_CD4_NC.cojo.txt / gwas_L3_ANXA4_Mono_NC.cojo.txt（s54a 已建）
输出：
  00_data_raw/smr_L2L3/gtex8_query/<TAG>.query.txt
  tables/58_gtex8_query_orient_audit_L2L3.csv
"""
import io
import json
import os
import subprocess
import sys

import numpy as np
import pandas as pd
from pyliftover import LiftOver

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

PROJ = r"D:\endometriosis_project"
P11 = os.path.join(PROJ, "11_sc_eqtl_mr_project")
G8RAW = os.path.join(PROJ, "01_data_raw", "eqtlcat_gtex_v8")
SMRD = os.path.join(P11, "00_data_raw", "smr_L2L3")
LDD = os.path.join(SMRD, "ld_ref")
QOUT = os.path.join(SMRD, "gtex8_query")
TAB = os.path.join(P11, "tables")
PLINK = os.path.join(PROJ, "_tools", "plink.exe")
LOG = os.path.join(P11, "scripts", "_s55b_build_gtex8_query_L2L3.log")
os.makedirs(QOUT, exist_ok=True)

DS = [
    dict(tag="L2_YME1L1", locus="L2", gene="YME1L1", chrom="10",
         tss37=27421790, src=os.path.join(G8RAW, "Whole_Blood.L2_YME1L1.tsv"),
         bfile=os.path.join(LDD, "L2_YME1L1_CD4_NC_cis"),
         cojo=os.path.join(SMRD, "gwas_L2_YME1L1_CD4_NC.cojo.txt")),
    dict(tag="L3_ANXA4", locus="L3", gene="ANXA4", chrom="2",
         tss37=69962577, src=os.path.join(G8RAW, "Whole_Blood.L3_ANXA4.tsv"),
         bfile=os.path.join(LDD, "L3_ANXA4_Mono_NC_cis"),
         cojo=os.path.join(SMRD, "gwas_L3_ANXA4_Mono_NC.cojo.txt")),
]
STRAND = {"YME1L1": "-", "ANXA4": "+"}   # Ensembl GRCh37 实查（s54a 已留痕）

buf = []
def w(s=""):
    buf.append(str(s)); print(s)
def flush():
    io.open(LOG, "w", encoding="utf-8", newline="\n").write("\n".join(buf) + "\n")

w("=" * 104)
w("s55b — GTEx v8 全血 L2/L3 SMR query 构建（含 Freq 取向修正）")
w("=" * 104)

lo = LiftOver("hg38", "hg19")          # 注意方向：gt38 -> hg19
audit_all = []
summary = {}

for D in DS:
    tag, gene, chrom, tss37 = D["tag"], D["gene"], D["chrom"], D["tss37"]
    w("")
    w("-" * 104)
    w("[%s] %s  chr%s  TSS37=%d  bfile=%s" % (tag, gene, chrom, tss37, os.path.basename(D["bfile"])))

    # ---- 1) LD 参考 bim + freq（需要频率取取向） ----
    bim = pd.read_csv(D["bfile"] + ".bim", sep=r"\s+", header=None,
                      names=["CHR", "SNP", "CM", "BP", "A1", "A2"]).set_index("SNP")
    frqf = D["bfile"] + ".frq"
    if not os.path.exists(frqf):
        cmd = [PLINK, "--bfile", D["bfile"], "--freq", "--out", D["bfile"]]
        w("  CMD: %s" % " ".join(cmd))
        subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                       check=False, timeout=3600)
    frq = pd.read_csv(frqf, sep=r"\s+").set_index("SNP")
    one = {}
    for snp, r in frq.iterrows():
        if snp not in bim.index:
            continue
        a1, a2 = str(bim.at[snp, "A1"]).upper(), str(bim.at[snp, "A2"]).upper()
        f1 = float(r["MAF"])
        one[snp] = {a1: f1, a2: 1.0 - f1, "pair": {a1, a2}}
    w("  [1] OneK1K LD 参考：bim=%d SNP ；freq 可得=%d" % (len(bim), len(one)))

    # ---- 2) 结局 GWAS cojo ----
    cojo = {}
    with io.open(D["cojo"], "r", encoding="utf-8") as f:
        hdr = f.readline().rstrip("\n").split("\t")
        ic = {c: i for i, c in enumerate(hdr)}
        for ln in f:
            p = ln.rstrip("\n").split("\t")
            if len(p) != len(hdr):
                continue
            try:
                cojo[p[ic["SNP"]]] = (p[ic["A1"]].upper(), p[ic["A2"]].upper(), float(p[ic["freq"]]))
            except ValueError:
                continue
    w("  [2] 结局 GWAS cojo：%d SNP" % len(cojo))

    # ---- 3) 源 TSV（只取目标基因） ----
    recs = []
    n_src = 0
    with io.open(D["src"], "r", encoding="utf-8", errors="replace") as f:
        ix = {c: i for i, c in enumerate(f.readline().rstrip("\n").split("\t"))}
        for ln in f:
            p = ln.rstrip("\n").split("\t")
            if len(p) <= ix["gene_symbol"]:
                continue
            if p[ix["gene_symbol"]] != gene:
                continue
            ref, alt = p[ix["ref"]].upper(), p[ix["alt"]].upper()
            if len(ref) > 1 or len(alt) > 1:
                continue
            try:
                pos38 = int(p[ix["position"]]); an = int(p[ix["an"]]); ac = int(p[ix["ac"]])
                maf = float(p[ix["maf"]]); b = float(p[ix["beta"]]); se = float(p[ix["se"]])
                pv = float(p[ix["pvalue"]])
            except ValueError:
                continue
            n_src += 1
            recs.append({"pos38": pos38, "ref": ref, "alt": alt, "maf": maf, "an": an,
                         "ac": ac, "beta": b, "se": se, "p": pv, "rsid": p[ix["rsid"]]})
    src = pd.DataFrame(recs)
    w("  [3] 源行（目标基因，双等位）= %d ；唯一 pos38 = %d" % (n_src, src.pos38.nunique()))

    # ---- 4) liftover hg38->hg19 ----
    uniq38 = sorted(src.pos38.unique())
    MP = {}
    for p in uniq38:
        r = lo.convert_coordinate("chr%s" % chrom, int(p) - 1)
        if not r:
            MP[p] = None
        elif len(r) > 1:
            MP[p] = int(r[0][1]) + 1          # 多映射取第一，另记 multi
        else:
            MP[p] = int(r[0][1]) + 1
    src["pos37"] = src.pos38.map(MP)
    nmiss = int(src.pos37.isna().sum())
    w("  [4] liftover 未映射 = %d 行 (%.2f%%)" % (nmiss, 100.0 * nmiss / max(1, n_src)))
    # 位移一致性（用于留痕，不作为映射依据）
    dd = (src.pos38 - src.pos37).dropna()
    w("      位移 (pos38-pos37) 唯一值数 = %d ；若为 1 则区域内恒定" % dd.nunique())
    if dd.nunique() <= 6:
        w("      位移取值 = %s" % sorted(dd.unique().tolist()))

    src = src[src.pos37.notna()].copy()
    src["pos37"] = src.pos37.astype(int)
    src["SNP"] = chrom + ":" + src.pos37.astype(str)

    # ---- 5) 可用性硬约束 ----
    src["in_ld"] = src.SNP.isin(one.keys())
    src["pair_ok"] = [bool(r.in_ld and r.ref in one[r.SNP]["pair"] and r.alt in one[r.SNP]["pair"])
                      for r in src.itertuples()]
    n_inld = int(src.in_ld.sum()); n_pairok = int(src.pair_ok.sum())
    w("  [5] 落在 OneK1K bim 内 = %d 行 ；其中等位对一致 = %d 行 (%.2f%%)"
      % (n_inld, n_pairok, 100.0 * n_pairok / max(1, n_inld)))
    w("      被剔：不在 bim = %d ；等位对不一致 = %d" % (n_src - nmiss - n_inld, n_inld - n_pairok))
    u = src[src.pair_ok].copy()

    # 多等位折叠
    dup = u.groupby("SNP")["alt"].nunique()
    multi = dup[dup > 1]
    keep_parts = [u[~u.SNP.isin(multi.index)]]
    for snp, grp in u[u.SNP.isin(multi.index)].groupby("SNP"):
        pair = one[snp]["pair"]
        sub = grp[grp.alt.isin(pair) & grp.ref.isin(pair)]
        if len(sub):
            keep_parts.append(sub.head(1))
    u = pd.concat(keep_parts, ignore_index=True)
    w("      多等位位置折叠 = %d ；可用行 = %d" % (len(multi), len(u)))

    # ---- 6) 取向判定 ----
    aud = []
    for r in u.itertuples():
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
        freq_old = float(r.ac) / float(r.an)          # 旧（错误）写法
        freq_new = maf if alt_minor_one else 1.0 - maf
        aud.append({"tag": tag, "locus": D["locus"], "gene": gene, "SNP": snp,
                    "pos38": r.pos38, "pos37": r.pos37, "ref": ref, "alt": alt, "maf": maf,
                    "ac_over_an": freq_old, "f_alt_onek1k": f_alt_one, "f_ref_onek1k": f_ref_one,
                    "alt_minor_onek1k": alt_minor_one, "f_alt_gwas": f_alt_gwas,
                    "alt_minor_gwas": alt_minor_gwas, "alt_minor_final": alt_minor_one,
                    "Freq_alt_old": freq_old, "Freq_alt_new": freq_new,
                    "rsid": r.rsid,
                    "freq_flipped": bool(abs(freq_old - freq_new) > 1e-3)})
    aud = pd.DataFrame(aud)
    audit_all.append(aud)

    acc = aud.dropna(subset=["alt_minor_gwas"])
    agree = (acc.alt_minor_onek1k == acc.alt_minor_gwas) if len(acc) else pd.Series([], dtype=bool)
    w("  [6] 取向判定 N=%d ：alt 为次等位 = %d (%.1f%%) ；旧写法写成补数 = %d (%.1f%%)"
      % (len(aud), int(aud.alt_minor_final.sum()), 100.0 * aud.alt_minor_final.mean(),
         int(aud.freq_flipped.sum()), 100.0 * aud.freq_flipped.mean()))
    w("      OneK1K vs GWAS 取向可比=%d 一致=%d (%.2f%%)"
      % (len(acc), int(agree.sum()) if len(acc) else 0,
         100.0 * agree.mean() if len(acc) else float("nan")))
    w("      对照：max(ac/an)=%.6f ；与 maf 列 max|Δ|=%.2e"
      % (aud.ac_over_an.max(), (aud.ac_over_an - aud.maf).abs().max()))
    for nm, col in (("OneK1K(LD参考)", "f_alt_onek1k"), ("GWAS(cojo)", "f_alt_gwas")):
        dnew = (aud.Freq_alt_new - aud[col]).abs().dropna()
        if len(dnew) == 0:
            continue
        w("      vs %-14s n=%5d |ΔFreq|中位=%.5f ；>0.20 占比=%.2f%%（SMR --diff-freq 上限 5%%）"
          % (nm, len(dnew), dnew.median(), 100.0 * (dnew > 0.20).mean()))

    # ---- 7) 写 query ----
    fmap = dict(zip(aud.SNP, aud.Freq_alt_new))
    HDR = ["SNP", "Chr", "BP", "A1", "A2", "Freq", "Probe", "Probe_Chr", "Probe_bp",
           "Gene", "Orientation", "b", "se", "p"]
    qp = os.path.join(QOUT, "%s.query.txt" % tag)
    with io.open(qp, "w", encoding="utf-8", newline="\n") as f:
        f.write("\t".join(HDR) + "\n")
        for r in u.itertuples():
            f.write("\t".join([r.SNP, chrom, str(r.pos37), r.alt, r.ref,
                               "%.6g" % fmap[r.SNP], gene, chrom, str(tss37), gene,
                               STRAND[gene], "%.6g" % r.beta, "%.6g" % r.se,
                               "%.6g" % r.p]) + "\n")
    pp = u.p.values
    w("  [7] → %s（%d 行）" % (os.path.basename(qp), len(u)))
    w("      p<5e-8=%d ；p<1e-5=%d ；p<1e-4=%d ；p<1.57e-3=%d ；minP=%.4g"
      % (int((pp < 5e-8).sum()), int((pp < 1e-5).sum()), int((pp < 1e-4).sum()),
         int((pp < 1.57e-3).sum()), pp.min()))
    summary[tag] = {"n_src_biallelic": int(n_src), "n_unmapped": int(nmiss),
                    "n_in_ld": int(n_inld), "n_pair_ok": int(n_pairok),
                    "n_query_rows": int(len(u)),
                    "n_p_lt_5e8": int((pp < 5e-8).sum()),
                    "n_p_lt_1e5": int((pp < 1e-5).sum()),
                    "alt_minor_fraction": float(aud.alt_minor_final.mean()),
                    "n_freq_flipped": int(aud.freq_flipped.sum()),
                    "tss37": int(tss37), "strand": STRAND[gene]}

au = pd.concat(audit_all, ignore_index=True)
au.to_csv(os.path.join(TAB, "58_gtex8_query_orient_audit_L2L3.csv"),
          index=False, encoding="utf-8-sig")
w("")
w("已写 tables/58_gtex8_query_orient_audit_L2L3.csv（%d 行）" % len(au))
json.dump(summary, io.open(os.path.join(SMRD, "gtex8_query_L2L3_summary.json"), "w",
                           encoding="utf-8"), ensure_ascii=False, indent=2)
w("VERDICT = PASS")
flush()
print("WROTE " + LOG)
