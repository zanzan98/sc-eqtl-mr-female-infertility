# -*- coding: utf-8 -*-
"""s17 —— Replication MR（1M-scBloodNL 作为暴露，结局沿用 GCST90483463）

约定与红线：
  * 复现 = 用同一结局重跑 Wald ratio，暴露换成 scBloodNL 的 cis-eQTL。
  * Discovery 与 Replication 分开呈现；**不做 meta 合并、不做跨数据集 IVW**（用户指令 4）。
  * 复现判据（用户指令 2，需在报告中显式声明）：
      强复现 = 方向一致 且 P < 0.05
      弱复现 = 方向一致 且 0.05 <= P < 0.10
      未复现 = 方向相反 或 P >= 0.10
  * 细胞类型为**多对一映射**（OneK1K 14 精细 -> scBloodNL lowerres 10 种）。

数据侧事实（本轮实测）：
  * scBloodNL 归档按「条件目录 / 细胞类型文件」组织：./<cond>/<celltype>_expression_eQTLsFDR-ProbeLevel.txt.gz
  * 22 列 MatrixEQTL 风格；SNPName=rsID, SNPChrPos=GRCh37 位置, HGNCName=基因符号,
    AlleleAssessed=被评估等位, Beta (SE)=两套数据(v2;v3)各自的 beta/SE, OverallZScore/PValue/FDR。
  * 坐标系 GRCh37，与 OneK1K 一致 -> 直接按位置/等位匹配，无需 liftover。

★ 匹配键口径（Yan 2026-09-22 裁定，硬约束）：
  「eQTL 工具变量与基因型矩阵的匹配一律用 rsID 键，绝对不用 chr:pos 键；
    若需统一 build，用 chain-file liftover（PLINK --update-map），不要手算坐标偏移。」
  实测约束冲突与处置：
   * **OneK1K 侧无任何 rsID**：`plink_merged_980_donors.bim` 第 2 列是 `chr:pos`（GRCh37，5,326,925 变异）；
     `01b_discovery_instruments_grch38.csv` 亦无 rsID 列（仅 variant_id/chr37/pos37/chr38/pos38）。
     → 暴露↔复现两侧**同在 GRCh37**，chr:pos 键**不构成错配**（这是 build 一致的情形，不是跨 build）。
   * 为避免"看似可用实则错配"，本脚本对每个工具变量额外做**跨源 rsID 交叉验证**：
     工具 GRCh38 位置 → 结局 harmonised 文件的 `rsid` 列（非循环，独立来源）→ 得 `expected_rsid`；
     与 scBloodNL 的 `SNPName` 比对，写入 `rsid_consistent` 列。**不一致即报错级警示。**
   * **共定位阶段**（coloc）所有 eQTL/LD/GWAS 三方合并**一律以 rsID 为键**（见 08 文档）。
"""
import argparse, gzip, io, os, sys, tarfile
import numpy as np
import pandas as pd
from scipy.stats import norm

ROOT = r"D:\endometriosis_project\11_sc_eqtl_mr_project"
T = os.path.join(ROOT, "tables")
COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}


def log(msgs, s):
    msgs.append(str(s))


def load_expected_rsid(outcome_gz, targets_by_chrom, msgs=None):
    """从结局 harmonised 文件取 rsID：键 = (chromosome, base_pair_location=GRCh38)。

    这是**非循环**的 rsID 来源——不依赖 scBloodNL 自身，故可独立验证其 SNPName。
    harmonised 表头实测：chromosome, base_pair_location, effect_allele, other_allele,
    beta, standard_error, effect_allele_frequency, p_value, **rsid**, rs_id, ...
    """
    want = {c: set(v) for c, v in targets_by_chrom.items()}
    remaining = sum(len(v) for v in want.values())
    got = {}
    with gzip.open(outcome_gz, "rt", encoding="utf-8", errors="replace") as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        ic = hdr.index("chromosome")
        ip = hdr.index("base_pair_location")
        ir = hdr.index("rsid")
        for line in fh:
            fr = line.rstrip("\n").split("\t")
            if len(fr) <= ir:
                continue
            c = fr[ic].strip()
            if c not in want:
                continue
            try:
                p = int(fr[ip])
            except ValueError:
                continue
            if p in want[c]:
                got[(c, p)] = fr[ir].strip()
                want[c].discard(p)
                remaining -= 1
                if remaining <= 0:
                    break
    if msgs is not None:
        log(msgs, "结局 harmonised rsID 回填 = %d / %d 个工具 GRCh38 位置"
            % (len(got), sum(len(v) for v in targets_by_chrom.values())))
    return got


class LimitedReader(io.RawIOBase):
    def __init__(self, fh, limit):
        self.fh = fh
        self.n = 0
        self.limit = limit

    def read(self, n=-1):
        if n is None or n < 0:
            n = 1 << 20
        n = min(n, self.limit - self.n)
        if n <= 0:
            return b''
        b = self.fh.read(n)
        self.n += len(b)
        return b

    def readable(self):
        return True


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--archive", required=True, help="scBloodNL eqtls_genome_wide tar.gz 路径")
    ap.add_argument("--condition", default="UT", help="刺激条件目录名，如 UT / 24hCA")
    ap.add_argument("--limit-bytes", type=float, default=0.0,
                    help="只读前 N 字节（0 = 完整文件）；用于前缀就地验证")
    ap.add_argument("--plan", default=os.path.join(T, "22_replication_test_plan.csv"))
    ap.add_argument("--out", required=True)
    ap.add_argument("--log", required=True)
    ap.add_argument("--reuse-hits", default="",
                    help="复用已缓存的原始命中表（跳过归档扫描，用于快速迭代分析口径）")
    ap.add_argument("--outcome", default=os.path.join(
        ROOT, "00_data_raw", "gwas", "GCST90483463.h.tsv.gz"),
        help="结局 harmonised 文件（GRCh38），用于非循环地回填/校验 rsID")
    ap.add_argument("--no-rsid-check", action="store_true",
                    help="跳过 rsID 跨源交叉验证（默认执行）")
    ap.add_argument("--tier1-out", default="",
                    help="额外单独输出 Tier1 明细表（用户指令 5：Tier1 必须单独列出）")
    a = ap.parse_args()

    msgs = []
    plan = pd.read_csv(a.plan, dtype={"variant_id_grch37": str, "variant_id_grch38": str})
    pool = pd.read_csv(os.path.join(T, "18_replication_candidate_pool.csv"),
                       dtype={"variant_id_grch37": str, "variant_id_grch38": str})
    plan = plan.merge(pool[["gene", "cell_type", "effect_allele", "other_allele", "af"]],
                      on=["gene", "cell_type"], how="left")
    assert plan["effect_allele"].notna().all(), "等位列合并失败"
    log(msgs, "=== s17 Replication MR :: condition=%s ===" % a.condition)
    log(msgs, "计划行数 = %d ; 基因 = %s ; 粗分群 = %s"
        % (len(plan), sorted(plan.gene.unique()), sorted(plan.cell_type_scBloodNL.unique())))

    # 目标集合：(GRCh37 位置, 基因) -> 多行索引
    targets = {}
    for _, r in plan.iterrows():
        chrom, pos = r["variant_id_grch37"].split(":")
        targets.setdefault((chrom, int(pos), r["gene"]), []).append(r)
    log(msgs, "待查 (chr,pos,gene) 组合 = %d" % len(targets))

    # ---- rsID 跨源交叉验证表（非循环来源 = 结局 harmonised 的 rsid 列）----
    targets_by_chrom = {}
    for v in sorted(set(plan["variant_id_grch38"].dropna().astype(str))):
        c38, p38 = v.split(":")
        targets_by_chrom.setdefault(c38, []).append(int(p38))
    expected_rsid = {}
    if not a.no_rsid_check:
        if os.path.exists(a.outcome):
            log(msgs, "rsID 交叉验证来源 = %s" % os.path.basename(a.outcome))
            try:
                expected_rsid = load_expected_rsid(a.outcome, targets_by_chrom, msgs)
            except Exception as ex:
                log(msgs, "!! rsID 回填失败（%s: %s）→ rsid_consistent 全部置 NA" % (type(ex).__name__, ex))
        else:
            log(msgs, "!! 结局文件不存在，跳过 rsID 交叉验证：%s" % a.outcome)
    else:
        log(msgs, "--no-rsid-check：跳过 rsID 交叉验证")

    need_cells = sorted(plan.cell_type_scBloodNL.unique())
    log(msgs, "需扫描的 scBloodNL 细胞类型文件 = %s" % need_cells)

    fsize = os.path.getsize(a.archive)
    limit = fsize if a.limit_bytes <= 0 else min(fsize, a.limit_bytes)
    log(msgs, "归档 = %s (%d B) ; 读取上限 = %d B (%.2f GB)"
        % (os.path.basename(a.archive), fsize, limit, limit / 1e9))

    COL = {"pval": "PValue", "snp": "SNPName", "chr": "SNPChr", "pos": "SNPChrPos",
           "probe": "ProbeName", "cistrans": "CisTrans", "assessed": "AlleleAssessed",
           "z": "OverallZScore", "nsamp": "DatasetsNrSamples", "hgnc": "HGNCName",
           "betase": "Beta (SE)", "fdr": "FDR", "meta": "Meta-Beta (SE)"}

    hits = []
    scanned = {}
    if a.reuse_hits:
        h0 = pd.read_csv(a.reuse_hits, dtype={"chrom": str})
        hits = h0.to_dict("records")
        log(msgs, "复用缓存命中表 %s（%d 行），跳过归档扫描" % (a.reuse_hits, len(h0)))
    else:
        try:
            with open(a.archive, "rb") as fh:
                rd = LimitedReader(fh, limit)
                gz = gzip.GzipFile(fileobj=rd, mode="rb")
                tf = tarfile.open(fileobj=gz, mode="r|")
                for m in tf:
                    base = os.path.basename(m.name)
                    parts = m.name.strip("./").split("/")
                    if len(parts) < 2 or parts[0] != a.condition:
                        continue
                    ctype = parts[1].split("_")[0]
                    if base.split("_")[0] not in need_cells:
                        continue
                    log(msgs, "")
                    log(msgs, "--- 扫描 %s ---" % m.name)
                    f = tf.extractfile(m)
                    g = gzip.GzipFile(fileobj=f, mode="rb")
                    hdr = None
                    ix = {}
                    nrow = 0
                    nhit = 0
                    for raw in io.BufferedReader(g, buffer_size=1 << 22):
                        line = raw.decode("utf-8", errors="replace").rstrip("\n")
                        if hdr is None:
                            hdr = line.split("\t")
                            for k, c in COL.items():
                                if c in hdr:
                                    ix[k] = hdr.index(c)
                            continue
                        nrow += 1
                        fr = line.split("\t")
                        if len(fr) < len(hdr):
                            continue
                        try:
                            pos = int(float(fr[ix["pos"]]))
                        except Exception:
                            continue
                        hgnc = fr[ix["hgnc"]].strip()
                        key = (fr[ix["chr"]].strip(), pos, hgnc)
                        if key not in targets:
                            continue
                        nhit += 1
                        hits.append(dict(
                            condition=a.condition, file_celltype=ctype, gene=hgnc,
                            snp=fr[ix["snp"]], chrom=key[0], pos_grch37=pos,
                            assessed=fr[ix["assessed"]].strip(),
                            cistrans=fr[ix["cistrans"]].strip(),
                            probe=fr[ix["probe"]], hgnc=hgnc,
                            z=fr[ix["z"]], nsamples=fr[ix["nsamp"]],
                            beta_se=fr[ix["betase"]], meta_beta_se=fr[ix["meta"]],
                            pval=fr[ix["pval"]], fdr=fr[ix["fdr"]]))
                    scanned[m.name] = (nrow, nhit)
                    log(msgs, "    行数=%d  命中目标行=%d" % (nrow, nhit))
                tf.close()
        except Exception as ex:
            log(msgs, "流式扫描结束（%s: %s）" % (type(ex).__name__, ex))

    log(msgs, "")
    log(msgs, "各成员扫描统计:")
    for k, v in scanned.items():
        log(msgs, "   %-64s rows=%9d hits=%d" % (k, v[0], v[1]))

    if not hits:
        log(msgs, "")
        log(msgs, "!! 条件 '%s' 下未找到任何目标行（可能该条件目录不在读取上限内）" % a.condition)
        pd.DataFrame().to_csv(a.out, index=False, encoding="utf-8-sig")
        with open(a.log, "w", encoding="utf-8") as fh:
            fh.write("\n".join(msgs))
        return

    h = pd.DataFrame(hits)
    h["pval"] = pd.to_numeric(h["pval"], errors="coerce")
    h["z"] = pd.to_numeric(h["z"], errors="coerce")
    h["fdr"] = pd.to_numeric(h["fdr"], errors="coerce")
    h.to_csv(a.out.replace(".csv", "_rawhits.csv"), index=False, encoding="utf-8-sig")

    log(msgs, "")
    log(msgs, "原始命中行 = %d ; 唯一 (gene,pos) = %d"
        % (len(h), h.groupby(["gene", "pos_grch37"]).ngroups))

    # ---- 合并到计划 ----
    rows = []
    for _, r in plan.sort_values(["tier", "gene", "cell_type"]).iterrows():
        chrom, pos = r["variant_id_grch37"].split(":")
        pos = int(pos)
        sub = h[(h["gene"] == r["gene"]) & (h["pos_grch37"] == pos)
                & (h["file_celltype"] == r["cell_type_scBloodNL"])
                & (h["cistrans"].str.lower() == "cis")]
        rec = dict(
            tier=r["tier"], locus_id=r["locus_id"], gene=r["gene"],
            cell_type_discovery=r["cell_type"],
            cell_type_replication=r["cell_type_scBloodNL"],
            variant_id_grch37=r["variant_id_grch37"], variant_id_grch38=r["variant_id_grch38"],
            rsid=sub["snp"].iloc[0] if len(sub) else "",
            discovery_b=r["discovery_b"], discovery_se=r["discovery_se"],
            discovery_p=r["discovery_p"], discovery_qval=r["discovery_qval"],
            discovery_F=r["discovery_F"],
            beta_outcome=r["beta_outcome"], se_outcome=r["se_outcome"],
            p_outcome=r["p_outcome"], outcome_GWS=bool(r["outcome_GWS"]),
        )
        if len(sub) == 0:
            rec.update(replication_b=np.nan, replication_se=np.nan, replication_p=np.nan,
                       direction_consistent=False, replication_status="no_record_in_replication",
                       interpretation_flag="NOT_TESTABLE_no_record", instrument_qc="absent",
                       n_probe=0, scb_beta=np.nan, scb_se=np.nan, scb_z=np.nan,
                       scb_F=np.nan, scb_pval=np.nan, scb_fdr=np.nan,
                       allele_check="NA", n_datasets=0)
            rows.append(rec)
            continue

        b = sub["beta_se"].str.strip()
        nd = 0
        bs = []
        ses = []
        for item in str(b.iloc[0]).split(";"):
            item = item.strip()
            if "(" not in item:
                continue
            try:
                bb = float(item.split("(")[0])
                ss = float(item.split("(")[1].rstrip(")"))
            except Exception:
                continue
            if np.isfinite(bb) and np.isfinite(ss) and ss > 0:
                bs.append(bb)
                ses.append(ss)
        nd = len(bs)
        if nd == 0:
            rec.update(replication_b=np.nan, replication_se=np.nan, replication_p=np.nan,
                       direction_consistent=False, replication_status="no_usable_beta_se",
                       interpretation_flag="NOT_TESTABLE_no_beta_se", instrument_qc="unknown",
                       n_probe=len(sub), scb_beta=np.nan, scb_se=np.nan,
                       scb_z=float(sub["z"].iloc[0]), scb_F=(float(sub["z"].iloc[0]) ** 2),
                       scb_pval=float(sub["pval"].iloc[0]),
                       scb_fdr=float(sub["fdr"].iloc[0]), allele_check="NA", n_datasets=0)
            rows.append(rec)
            continue

        wts = [1.0 / (s * s) for s in ses]
        b_meta = sum(w * x for w, x in zip(wts, bs)) / sum(wts)
        se_meta = float(np.sqrt(1.0 / sum(wts)))

        # ---- 等位对齐 ----
        ea = str(r["effect_allele"]).upper()
        oa = str(r["other_allele"]).upper()
        ass = str(sub["assessed"].iloc[0]).upper()
        if ass == ea:
            flip, chk = False, "same"
        elif ass == COMP.get(ea, ""):
            flip, chk = False, "complement_of_effect"
        elif ass == oa:
            flip, chk = True, "swapped"
        elif ass == COMP.get(oa, ""):
            flip, chk = True, "complement_of_other"
        else:
            flip, chk = None, "unresolved(%s_%s_vs_%s)" % (ea, oa, ass)
        if flip is None:
            rec.update(replication_b=np.nan, replication_se=np.nan, replication_p=np.nan,
                       direction_consistent=False, replication_status="allele_unresolved",
                       interpretation_flag="NOT_TESTABLE_allele_unresolved", instrument_qc="unknown",
                       n_probe=len(sub), scb_beta=b_meta, scb_se=se_meta,
                       scb_z=float(sub["z"].iloc[0]), scb_F=(float(sub["z"].iloc[0]) ** 2),
                       scb_pval=float(sub["pval"].iloc[0]),
                       scb_fdr=float(sub["fdr"].iloc[0]), allele_check=chk, n_datasets=nd)
            rows.append(rec)
            continue
        if flip:
            b_meta = -b_meta

        if not (b_meta == b_meta) or abs(b_meta) < 1e-300:
            rec.update(replication_b=np.nan, replication_se=np.nan, replication_p=np.nan,
                       direction_consistent=False, replication_status="zero_exposure_effect",
                       interpretation_flag="NOT_TESTABLE_zero_exposure_effect", instrument_qc="unknown",
                       n_probe=len(sub), scb_beta=b_meta, scb_se=se_meta,
                       scb_z=float(sub["z"].iloc[0]), scb_F=(float(sub["z"].iloc[0]) ** 2),
                       scb_pval=float(sub["pval"].iloc[0]),
                       scb_fdr=float(sub["fdr"].iloc[0]), allele_check=chk, n_datasets=nd)
            rows.append(rec)
            continue

        rep_b = float(r["beta_outcome"]) / b_meta
        bo = float(r["beta_outcome"])
        seo = float(r["se_outcome"])
        # ★ 与 s06 完全同口径：二阶 delta 法（含暴露端不确定度），保证 Discovery/Replication 可比
        rep_se = float(np.sqrt(seo ** 2 / b_meta ** 2 + (bo ** 2) * (se_meta ** 2) / (b_meta ** 4)))
        rep_p = float(2 * norm.sf(abs(rep_b / rep_se))) if rep_se > 0 else np.nan
        # 复现队列中该工具变量的 eQTL 是否成立
        scb_z = float(sub["z"].iloc[0]) if np.isfinite(sub["z"].iloc[0]) else np.nan
        scb_p = float(sub["pval"].iloc[0])
        scb_q = float(sub["fdr"].iloc[0])
        instr_qc = "valid" if (scb_q < 0.05) else ("weak_nominal" if scb_p < 0.05 else "not_detected")
        dc = (np.sign(rep_b) == np.sign(r["discovery_b"]))
        if dc and rep_p < 0.05:
            cls = "strong_replication"
        elif dc and rep_p < 0.10:
            cls = "weak_replication"
        else:
            cls = "not_replicated"
        flag = ("informative" if instr_qc == "valid"
                else ("PARTIALLY_informative_instrument_nominal_only" if instr_qc == "weak_nominal"
                      else "UNINFORMATIVE_instrument_not_detected_in_replication"))
        rec.update(replication_b=rep_b, replication_se=rep_se, replication_p=rep_p,
                   direction_consistent=bool(dc), replication_status=cls,
                   interpretation_flag=flag, instrument_qc=instr_qc,
                   n_probe=len(sub), scb_beta=b_meta, scb_se=se_meta,
                   scb_z=scb_z, scb_F=(scb_z ** 2 if np.isfinite(scb_z) else np.nan),
                   scb_pval=scb_p, scb_fdr=scb_q, allele_check=chk, n_datasets=nd)
        rows.append(rec)

    res = pd.DataFrame(rows)

    # ---- 用户指令 3：提供固定 schema 的列名别名 ----
    res["cell_type"] = res["cell_type_discovery"]
    res["replication_tier"] = res["tier"]

    # ---- rsID 跨源交叉验证（硬约束：以 rsID 为准，chr:pos 仅作定位）----
    def _exp_rsid(v):
        if not isinstance(v, str) or ":" not in v:
            return ""
        c, p = v.split(":")
        return expected_rsid.get((c, int(p)), "")
    res["expected_rsid_outcome"] = [_exp_rsid(v) for v in res["variant_id_grch38"]]
    def _rsid_match(r):
        e = str(r["expected_rsid_outcome"]).strip()
        g = str(r["rsid"]).strip()
        if not e or not g or e.lower() in ("nan", "none"):
            return "NA"
        return "match" if e == g else "MISMATCH"
    res["rsid_consistent"] = res.apply(_rsid_match, axis=1)
    res["build_exposure"] = "GRCh37"      # OneK1K
    res["build_replication"] = "GRCh37"   # 1M-scBloodNL（SNPChrPos 实测 hg19）
    res["build_outcome"] = "GRCh38"       # GWAS Catalog harmonised base_pair_location
    res["match_key_used"] = "chr:pos(GRCh37=GRCh37) + rsID 交叉验证"
    res["method"] = "Wald ratio"

    n_mm = int((res["rsid_consistent"] == "MISMATCH").sum())
    n_ok = int((res["rsid_consistent"] == "match").sum())
    n_na = int((res["rsid_consistent"] == "NA").sum())
    log(msgs, "")
    log(msgs, "=== rsID 跨源交叉验证（scBloodNL SNPName vs 结局 harmonised rsid）===")
    log(msgs, "   match=%d  MISMATCH=%d  NA=%d" % (n_ok, n_mm, n_na))
    if n_mm:
        log(msgs, "   ★★ 警示：存在 rsID 不一致行，须先解决匹配键问题再看结果：")
        for _, r in res[res["rsid_consistent"] == "MISMATCH"].iterrows():
            log(msgs, "      %s %s  预期=%s  复现=%s"
                % (r["gene"], r["variant_id_grch38"], r["expected_rsid_outcome"], r["rsid"]))

    res.to_csv(a.out, index=False, encoding="utf-8-sig")

    # ---- 用户指令 5：Tier1 必须单独列出 ----
    t1 = res[res["tier"] == "Tier1"].copy()
    if a.tier1_out:
        t1.to_csv(a.tier1_out, index=False, encoding="utf-8-sig")
    log(msgs, "")
    log(msgs, "=== ★ Tier1（必须单独列出，逐行，不并入汇总）===")
    if len(t1) == 0:
        log(msgs, "   （计划中无 Tier1 行）")
    else:
        for _, r in t1.iterrows():
            log(msgs, "   %-8s / %-10s -> %-10s  rsid=%s(%s)"
                % (r["gene"], r["cell_type_discovery"], r["cell_type_replication"],
                   r["rsid"], r["rsid_consistent"]))
            log(msgs, "       discovery : b=%+.4f  p=%.4g  q=%.4g  F=%.1f  结局位点 P=%.4g (%s)"
                % (r["discovery_b"], r["discovery_p"], r["discovery_qval"], r["discovery_F"],
                   r["p_outcome"],
                   "结局 GWS" if bool(r["outcome_GWS"]) else "结局仅提示性，非 GWS"))
            log(msgs, "       replication: b=%s  se=%s  p=%s  方向一致=%s"
                % (("%+.4f" % r["replication_b"]) if np.isfinite(r["replication_b"]) else "NA",
                   ("%.4f" % r["replication_se"]) if np.isfinite(r["replication_se"]) else "NA",
                   ("%.4g" % r["replication_p"]) if np.isfinite(r["replication_p"]) else "NA",
                   r["direction_consistent"]))
            log(msgs, "       判定=%s ; 工具QC=%s ; 解释标记=%s ; 等位=%s"
                % (r["replication_status"], r["instrument_qc"],
                   r["interpretation_flag"], r["allele_check"]))
            # ---- 用户裁定的强制限定语：按「判定 × 工具QC」自动选词 ----
            # ★ 措辞来自用户 2026-09-22 指令，逐字使用，不得改写、不得省略、不得反向归因。
            if r["replication_status"] == "strong_replication" and r["instrument_qc"] == "weak_nominal":
                log(msgs, "       ★★ 强制限定语：Discovery 信号在 1M-scBloodNL 中获得名义复现，"
                          "但工具质量偏弱，需共定位进一步支持。")
            elif (r["replication_status"] == "not_replicated"
                  and r["instrument_qc"] == "not_detected"):
                log(msgs, "       ★★ 强制限定语：发现未获支持（**严禁写作「因果被否证」**）。"
                          "理由：复现数据集中该工具变异未被检测到或工具强度不足，"
                          "而非真实生物学效应不存在。")
            elif (r["replication_status"] == "not_replicated" and bool(r["direction_consistent"])
                  and r["instrument_qc"] == "weak_nominal"):
                log(msgs, "       ★★ 强制限定语：功效不足，非方向冲突。")

    log(msgs, "")
    log(msgs, "=== 复现结果（Discovery vs Replication 分列呈现，未做 meta 合并）===")
    show = res[["tier", "gene", "cell_type_discovery", "cell_type_replication", "rsid",
                "discovery_b", "discovery_p", "replication_b", "replication_p",
                "direction_consistent", "replication_status", "allele_check", "n_datasets"]]
    log(msgs, show.to_string(index=False))
    log(msgs, "")
    log(msgs, "=== 暴露端（eQTL）在复现队列中的状态 ===")
    show2 = res[["gene", "cell_type_replication", "rsid",
                 "scb_beta", "scb_se", "scb_z", "scb_F", "scb_pval", "scb_fdr",
                 "instrument_qc", "interpretation_flag"]]
    log(msgs, show2.to_string(index=False))
    log(msgs, "")
    log(msgs, "分类计数:")
    for k, v in res["replication_status"].value_counts().items():
        log(msgs, "   %-26s %d" % (k, v))
    log(msgs, "")
    log(msgs, "工具变量 QC 计数:")
    for k, v in res["instrument_qc"].value_counts().items():
        log(msgs, "   %-26s %d" % (k, v))
    log(msgs, "")
    log(msgs, "解释标记计数:")
    for k, v in res["interpretation_flag"].value_counts().items():
        log(msgs, "   %-52s %d" % (k, v))
    log(msgs, "")
    log(msgs, "按 Tier 分层:")
    for tier, g in res.groupby("tier"):
        log(msgs, "   %-12s n=%d : %s" % (tier, len(g), dict(g["replication_status"].value_counts())))
    log(msgs, "")
    log(msgs, "★ 判据声明：强复现=方向一致且P<0.05；弱复现=方向一致且0.05<=P<0.10；"
              "未复现=方向相反或P>=0.10。scBloodNL 样本量约 119（v2=79;v3=40），"
              "远低于 OneK1K（982），功效按 Z ∝ sqrt(N) 约为其 0.35 倍。")
    log(msgs, "★ 细胞类型为多对一映射（OneK1K 精细亚型 -> scBloodNL lowerres），"
              "故本步骤不能确认精细亚型特异性。")
    log(msgs, "★ 方法口径：本分析为单 SNP 工具变量，**主方法为 Wald ratio**；"
              "IVW（多 SNP）将在 OneK1K raw tensorQTL 全量数据到位后作为**敏感性分析**补跑。"
              "本表 nsnp 全为 1，属设计使然，非缺陷。")
    log(msgs, "★ 匹配键口径：暴露（OneK1K, GRCh37）与复现（scBloodNL, GRCh37）**同在 GRCh37**，"
              "chr:pos 键不构成跨 build 错配；每个工具另做 rsID 跨源交叉验证"
              "（来源 = 结局 harmonised 的 rsid 列，非循环）。共定位阶段三方合并一律按 rsID 键。")
    log(msgs, "★ 无 meta 合并、无跨数据集 IVW：Discovery 与 Replication 仅分列呈现（用户指令 4）。")

    with open(a.log, "w", encoding="utf-8") as fh:
        fh.write("\n".join(msgs))
    print("done")


if __name__ == "__main__":
    main()
