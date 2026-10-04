# -*- coding: utf-8 -*-
"""s43f_build_gtex8_wnt4.py — 由 GTEx v8 全血区域数据构建 SMR 输入（BESD）

输入：01_data_raw/eqtlcat_gtex_v8/Whole_Blood.chr1_21_23.3Mb_3genes.tsv（s43e 流式产出）
      LD 参考 bim：00_data_raw/smr_3p3/ld_ref/chr1_20_25mb.bim（hg19 位置，ID=`1:<pos37>`）

关键：GTEx v8（eQTL Catalogue harmonised）坐标为 **GRCh38**，而 LD 参考与 GWAS 为 **hg19/GRCh37**。
      → 必须把 GTEx variant 的 b38 位置回标为 hg19，SNP ID 才能与 bim 匹配。
      做法：**不做常数偏移假设**，而是把 LD 参考 bim 内全部 chr1:20–25 Mb 变异（9,623 个）
      正向 liftover（hg19→hg38，链文件已缓存），建立精确的 b38→hg19 反向映射；
      并统计偏移 δ 的分布以验证其是否恒定。

输出：00_data_raw/smr_3p3/gtex8/query/<GENE>.query.txt  （SMR query 14 列）
"""
import os, io, sys, csv, json
from collections import Counter
import numpy as np

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

ROOT = r"D:\endometriosis_project"
PROJ = os.path.join(ROOT, "11_sc_eqtl_mr_project")
SRC = os.path.join(ROOT, "01_data_raw", "eqtlcat_gtex_v8",
                   "Whole_Blood.chr1_21_23.3Mb_3genes.tsv")
BIM = os.path.join(PROJ, "00_data_raw", "smr_3p3", "ld_ref", "chr1_20_25mb.bim")
OUTD = os.path.join(PROJ, "00_data_raw", "smr_3p3", "gtex8")
QDIR = os.path.join(OUTD, "query")
LOG = os.path.join(PROJ, "scripts", "_s43f_build_gtex8_wnt4.log")

# 三基因 GRCh38 坐标（Ensembl REST 实查，见 00_data_raw/smr_3p3/strand.json）
GENES = {
    "WNT4":      {"ens": "ENSG00000162552", "start38": 22117313, "end38": 22143969, "strand": "-"},
    "CDC42":     {"ens": "ENSG00000070831", "start38": 22052627, "end38": 22101360, "strand": "+"},
    "LINC00339": {"ens": "ENSG00000218510", "start38": 22024558, "end38": 22031223, "strand": "+"},
}
PRIMARY = "WNT4"

buf = []
def w(s=""):
    buf.append(str(s)); print(s)
def flush():
    os.makedirs(os.path.dirname(LOG), exist_ok=True)
    with io.open(LOG, "w", encoding="utf-8") as f:
        f.write("\n".join(buf) + "\n")

os.makedirs(QDIR, exist_ok=True)
w("=" * 94)
w("s43f_build_gtex8_wnt4.py — GTEx v8 全血 → SMR 输入（BESD query）")
w("=" * 94)

# ---- 1. 读 LD 参考 bim（chr1 20-25 Mb，hg19）
w()
w("[1] LD 参考 bim 读取与 hg19→hg38 正向 liftover（建立 b38→hg19 精确反向映射）")
bim = []
with io.open(BIM, "r", encoding="utf-8", errors="replace") as f:
    for ln in f:
        p = ln.rstrip("\n").split("\t")
        if len(p) < 6:
            continue
        bim.append((p[0], p[1], int(p[3]), p[4], p[5]))   # chr, id, pos, A1, A2
w("    bim 变异数 = %d  bp %d..%d" % (len(bim), min(b[2] for b in bim), max(b[2] for b in bim)))

from pyliftover import LiftOver
lo = LiftOver("hg19", "hg38")
MAP38_37 = {}     # pos38 -> pos37
DELTA = []
for _, _, pos37, _, _ in bim:
    r = lo.convert_coordinate("chr1", pos37 - 1)
    if not r or len(r) > 1:
        continue
    pos38 = int(r[0][1]) + 1
    MAP38_37.setdefault(pos38, pos37)
    DELTA.append(pos38 - pos37)
DELTA = np.array(DELTA)
w("    可正向映射 = %d / %d" % (len(MAP38_37), len(bim)))
w("    ★ 偏移 δ = pos38 − pos37 : min=%d  max=%d  median=%d  唯一值数=%d  众数=%d"
  % (DELTA.min(), DELTA.max(), int(np.median(DELTA)), len(set(DELTA.tolist())),
     Counter(DELTA.tolist()).most_common(1)[0][0]))
# 区域内（21.4–23.5 Mb hg19）δ 是否恒定
sub = DELTA[(np.array([b[2] for b in bim]) >= 21400000) & (np.array([b[2] for b in bim]) <= 23500000)]
w("    区域内(hg19 21.4–23.5 Mb) δ: min=%d max=%d 唯一值=%s"
  % (sub.min(), sub.max(), sorted(set(sub.tolist()))[:6]))
DELTA_MODE = Counter(DELTA.tolist()).most_common(1)[0][0]

# ---- 2. 读 GTEx v8 区域数据
w()
w("[2] GTEx v8 全血区域数据读取")
rows = {g: [] for g in GENES}
hdr = None
with io.open(SRC, "r", encoding="utf-8", errors="replace") as f:
    hdr = f.readline().rstrip("\n").split("\t")
    ix = {c: i for i, c in enumerate(hdr)}
    for ln in f:
        p = ln.rstrip("\n").split("\t")
        if len(p) != len(hdr):
            continue
        sym = p[ix["gene_symbol"]]
        if sym not in GENES:
            continue
        rows[sym].append(p)
for g in GENES:
    w("    %-10s 原始行 = %d" % (g, len(rows[g])))

# ---- 3. 构建 query
QCOLS = ["SNP", "Chr", "BP", "A1", "A2", "Freq", "Probe", "Probe_Chr",
         "Probe_bp", "Gene", "Orientation", "b", "se", "p"]
summary = {}
for g, meta in GENES.items():
    tss38 = meta["end38"] if meta["strand"] == "-" else meta["start38"]
    # ★ δ = pos38 − pos37（s43f 实测定值）；故 pos37 = pos38 − δ
    tss37 = tss38 - DELTA_MODE
    out = []
    n_nomap = 0
    n_bad = 0
    n_dup = 0
    seen = set()
    for p in rows[g]:
        pos38 = int(p[ix["position"]])
        pos37 = MAP38_37.get(pos38)
        if pos37 is None:
            n_nomap += 1
            continue
        ref = p[ix["ref"]].upper()
        alt = p[ix["alt"]].upper()
        if not ref or not alt or len(ref) > 1 or len(alt) > 1:   # 仅保留 SNP（与 LD 参考同类）
            n_bad += 1
            continue
        snp = "1:%d" % pos37
        if snp in seen:
            n_dup += 1
            continue
        seen.add(snp)
        an = float(p[ix["an"]]); ac = float(p[ix["ac"]])
        freq = (ac / an) if an else float("nan")
        out.append({
            "SNP": snp, "Chr": 1, "BP": pos37, "A1": alt, "A2": ref, "Freq": "%.6g" % freq,
            "Probe": g, "Probe_Chr": 1, "Probe_bp": tss37, "Gene": g,
            "Orientation": meta["strand"], "b": p[ix["beta"]], "se": p[ix["se"]],
            "p": p[ix["pvalue"]],
        })
    out.sort(key=lambda d: int(d["BP"]))
    qf = os.path.join(QDIR, "%s.query.txt" % g)
    with io.open(qf, "w", encoding="utf-8", newline="\n") as f:
        f.write("\t".join(QCOLS) + "\n")
        for d in out:
            f.write("\t".join(str(d[c]) for c in QCOLS) + "\n")
    pv = np.array([float(d["p"]) for d in out]) if out else np.array([])
    n_1e3 = int((pv < 1.57e-3).sum()) if len(pv) else 0
    n_5e8 = int((pv < 5e-8).sum()) if len(pv) else 0
    n_1e5 = int((pv < 1e-5).sum()) if len(pv) else 0
    n_1e4 = int((pv < 1e-4).sum()) if len(pv) else 0
    argmin = out[int(np.argmin(pv))] if len(pv) else None
    summary[g] = dict(n_raw=len(rows[g]), n_query=len(out), n_nomap=n_nomap, n_nonSNP=n_bad,
                      n_dup=n_dup,
                      tss37=tss37, minP=float(pv.min()) if len(pv) else float("nan"),
                      argmin=argmin, n_5e8=n_5e8, n_1e5=n_1e5, n_1e4=n_1e4, n_1e3=n_1e3)

w()
w("[3] 逐基因 query 构建结果（仅保留与 LD 参考同类、且 b38→hg19 可映射的双等位 SNP）")
w("%-10s %7s %8s %8s %7s %7s %10s %8s %8s %8s %8s" %
  ("gene", "raw", "query", "未映射", "非SNP", "同ID去重", "minP", "p<5e-8", "p<1e-5", "p<1e-4", "p<1.57e-3"))
for g in ("WNT4", "CDC42", "LINC00339"):
    s = summary[g]
    w("%-10s %7d %8d %8d %7d %7d %10.3e %8d %8d %8d %8d"
      % (g, s["n_raw"], s["n_query"], s["n_nomap"], s["n_nonSNP"], s["n_dup"], s["minP"],
         s["n_5e8"], s["n_1e5"], s["n_1e4"], s["n_1e3"]))
w()
w("★ WNT4 TSS: GRCh38 %d → hg19 %d (δ=%d)" % (GENES["WNT4"]["end38"], summary["WNT4"]["tss37"], DELTA_MODE))
w("   WNT4 最显著位点: %s  b=%s se=%s p=%s"
  % (summary["WNT4"]["argmin"]["SNP"] if summary["WNT4"]["argmin"] else "NA",
     summary["WNT4"]["argmin"]["b"] if summary["WNT4"]["argmin"] else "",
     summary["WNT4"]["argmin"]["se"] if summary["WNT4"]["argmin"] else "",
     "%.4e" % summary["WNT4"]["minP"]))
w()
w("★★ HEIDI 可行性预判（SMR 需 eQTL p<1.57e-3 的 SNP ≥3 个，且剪枝后仍 ≥3）")
for g in ("WNT4", "CDC42", "LINC00339"):
    s = summary[g]
    w("   %-10s p<1.57e-3 的 SNP 数 = %d  → %s"
      % (g, s["n_1e3"], "可供 HEIDI" if s["n_1e3"] >= 3 else "**HEIDI 可能因 m<3 无法运行**"))
w()
w("★ SMR target SNP 阈值预判（默认 --peqtl-smr 5e-8）")
for g in ("WNT4", "CDC42", "LINC00339"):
    s = summary[g]
    need = "5e-8 可用" if s["n_5e8"] else ("需放宽至 1e-5" if s["n_1e5"] else
                                        ("需放宽至 1e-4" if s["n_1e4"] else "**无 SNP 过 1e-4，须用 --target-snp 强制**"))
    w("   %-10s minP=%.3e  → %s" % (g, s["minP"], need))

json.dump({g: {k: (v if k != "argmin" else (v["SNP"] if v else None))
               for k, v in summary[g].items()} for g in summary},
          io.open(os.path.join(OUTD, "gtex8_build_summary.json"), "w", encoding="utf-8"),
          ensure_ascii=False, indent=1, default=str)
w()
w("VERDICT = %s" % ("PASS" if summary["WNT4"]["n_query"] > 0 else "FAIL"))
flush()
