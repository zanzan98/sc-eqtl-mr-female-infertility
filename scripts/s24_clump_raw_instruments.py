# -*- coding: utf-8 -*-
"""s24_clump_raw_instruments.py —— 对三座 raw 工具变量做 PLINK clumping（per gene×cell_type）

裁定 (a)：clumping r²<0.001 / 10 Mb（主）。
★ 忠实使用 PLINK 1.90 `--clump`，参考面板 = OneK1K 980 供者真实基因型（GRCh37）。
★ 每个 gene×cell_type 单独 clump → 工具集按细胞类型分别确定。

输出：tables/28_raw_instruments_clumped.csv（含 clump_index 标记、nsnp）
用法：
  python s24_clump_raw_instruments.py [--r2 0.001] [--kb 10000] [--tag main]
"""
import argparse, os, subprocess, time
import pandas as pd

PROJ = r"D:\endometriosis_project\11_sc_eqtl_mr_project"
PRE = r"D:\endometriosis_project\_onek1k_plink\plink_merged_980_donors"
PLINK = r"D:\endometriosis_project\_tools\plink.exe"
IV = os.path.join(PROJ, "tables", "27_raw_instruments_preclump.csv")
WORK = r"D:\endometriosis_project\_clump_raw"
LOG = r"D:\endometriosis_project\_s24_clump.log"

L = []
def p(s=""):
    L.append(str(s)); print(s)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--r2", type=float, default=0.001)
    ap.add_argument("--kb", type=int, default=10000)
    ap.add_argument("--tag", default="main")
    ap.add_argument("--p1", type=float, default=5e-8)
    a = ap.parse_args()

    os.makedirs(WORK, exist_ok=True)
    T0 = time.time()
    p("=== s24 PLINK clumping（r2<%g ; window %d kb ; tag=%s） ===" % (a.r2, a.kb, a.tag))
    p("时间 %s" % time.strftime("%H:%M:%S"))

    iv = pd.read_csv(IV, encoding="utf-8-sig")
    p("候选工具表 = %d 行 ; gene×cell_type 组合 = %d" % (len(iv), iv.groupby(["gene", "cell_type"]).ngroups))

    keep = []
    stats = []
    for (gene, ct), g in iv.groupby(["gene", "cell_type"], sort=True):
        g = g[g["variant_id_grch38"].notna()].copy()
        if g.empty:
            stats.append(dict(gene=gene, cell_type=ct, n_input=0, n_index=0, note="no_grch38"))
            continue
        chrm = str(g["chrom"].iloc[0])
        # 若同一位点在不同细胞类型重复出现，PLINK 的 clump 输入要求 SNP 唯一 → 去重取最小 p
        g = g.sort_values("pval_nominal").drop_duplicates(subset=["variant_id_grch37"])
        stem = "%s_%s_%s" % (a.tag, gene, ct)
        inf = os.path.join(WORK, stem + ".clumpin")
        with open(inf, "w", encoding="utf-8") as fh:
            fh.write("SNP\tP\n")
            for r in g.itertuples(index=False):
                fh.write("%s\t%.6g\n" % (r.variant_id_grch37, r.pval_nominal))
        outf = os.path.join(WORK, stem)
        cmd = [PLINK, "--bfile", PRE, "--chr", chrm,
               "--clump", inf, "--clump-snp-field", "SNP", "--clump-field", "P",
               "--clump-p1", "%.3g" % a.p1, "--clump-r2", "%g" % a.r2,
               "--clump-kb", str(a.kb), "--allow-no-sex",
               "--out", outf]
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
        clumped = outf + ".clumped"
        if not os.path.exists(clumped):
            stats.append(dict(gene=gene, cell_type=ct, n_input=len(g), n_index=0,
                              note="no_clumped_out rc=%d" % r.returncode))
            p("   !! %s/%s 无 .clumped（rc=%d）" % (gene, ct, r.returncode))
            continue
        try:
            cl = pd.read_csv(clumped, sep=r"\s+", comment=None)
        except Exception as e:
            stats.append(dict(gene=gene, cell_type=ct, n_input=len(g), n_index=0, note="read_err %r" % e))
            continue
        if cl.empty or "SNP" not in cl.columns:
            stats.append(dict(gene=gene, cell_type=ct, n_input=len(g), n_index=0, note="empty_clumped"))
            p("   !! %s/%s .clumped 为空" % (gene, ct))
            continue
        idx = set(cl["SNP"].astype(str))
        sub = g[g["variant_id_grch37"].isin(idx)].copy()
        sub["clump_index"] = 1
        sub["n_input_p5e8"] = len(g)
        sub["n_index"] = len(idx)
        keep.append(sub)
        dup = len(idx) - len(sub)
        stats.append(dict(gene=gene, cell_type=ct, n_input=len(g), n_index=len(idx),
                          note=("ok" if dup == 0 else "index_not_in_input=%d" % dup)))
        p("   %-10s %-12s 输入 %4d -> 独立 %3d   (%.1f s)"
          % (gene, ct, len(g), len(idx), time.time() - T0))

    res = pd.concat(keep, ignore_index=True) if keep else pd.DataFrame()
    outp = os.path.join(PROJ, "tables", "28_raw_instruments_clumped_%s.csv" % a.tag)
    if len(res):
        res = res.sort_values(["chrom", "gene", "cell_type", "pval_nominal"]).reset_index(drop=True)
        res.to_csv(outp, index=False, encoding="utf-8-sig")
        p("")
        p("已写出 %s（%d 行；%d 个 gene×cell_type 组合）"
          % (outp, len(res), res.groupby(["gene", "cell_type"]).ngroups))
    st = pd.DataFrame(stats)
    st.to_csv(os.path.join(PROJ, "tables", "28b_clump_stats_%s.csv" % a.tag),
              index=False, encoding="utf-8-sig")
    p("")
    p("=== clump 汇总 ===")
    p(st.to_string(index=False))
    p("")
    p("独立工具数分布：")
    p(st.groupby("n_index").size().rename("gene×cell_type 数").to_string())
    p("")
    p("★ 多 SNP（nsnp>=2，可 IVW）= %d 个组合 ; 单 SNP = %d 个"
      % (int((st["n_index"] >= 2).sum()), int((st["n_index"] == 1).sum())))
    p("总耗时 %.1f min" % ((time.time() - T0) / 60))
    open(LOG, "w", encoding="utf-8").write("\n".join(L) + "\n")


if __name__ == "__main__":
    main()
