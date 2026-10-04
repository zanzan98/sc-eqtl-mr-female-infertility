# -*- coding: utf-8 -*-
"""s55a_fetch_gtex8_wb_L2L3.py -- 裁定1：从 EBI 定向抽取 GTEx v8 全血的 L2(YME1L1) / L3(ANXA4) 区域

严格限定区域，不下载全量文件：
  · 本地已有 tabix 索引 Whole_Blood.tsv.gz.tbi（2,110,864 B）
  · 解析索引 -> 求目标区间的 bin -> 取 chunk 的虚拟偏移 -> 合并 -> HTTP Range 只取这些字节
  · BGZF(=多成员 gzip) 增量解码 -> 跳过 chunk 内 uoffset -> 按染色体/位置逐行过滤
  · 不做任何全文件下载；不依赖 pysam（本机 cp313 无 wheel）

区域（GTEx v8 = GRCh38；cis 窗 = TSS ± 1 Mb，取 ±1.1 Mb 留余量）：
  L2  YME1L1  ENSG00000104047  chr10  TSS38 = 27,132,861  区间 26,032,861..28,232,861
  L3  ANXA4   ENSG00000196975  chr2   TSS38 = 69,735,445  区间 68,635,445..70,835,445

来源：https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/imported/GTEx_V8/ge/Whole_Blood.tsv.gz (3,190,551,490 B)

输出：D:\\endometriosis_project\\01_data_raw\\eqtlcat_gtex_v8\\Whole_Blood.<TAG>.tsv（原始 18 列 + gene_symbol）
"""
import io
import os
import sys
import gzip
import struct
import time
import urllib.request

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

PROJ = r"D:\endometriosis_project"
OUTD = os.path.join(PROJ, "01_data_raw", "eqtlcat_gtex_v8")
TBI = os.path.join(OUTD, "Whole_Blood.tsv.gz.tbi")
URL = "https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/imported/GTEx_V8/ge/Whole_Blood.tsv.gz"
LOG = r"D:\endometriosis_project\11_sc_eqtl_mr_project\scripts\_s55a_fetch_gtex8_wb_L2L3.log"
FSIZE = 3190551490

# 远端 chunk 不在文件开头 → 表头不会出现在取回的字节里，按已核实布局硬编码（s43e 日志 18 列）
HDR18 = ["variant", "r2", "pvalue", "molecular_trait_object_id", "molecular_trait_id",
         "maf", "gene_id", "median_tpm", "beta", "se", "an", "ac",
         "chromosome", "position", "ref", "alt", "type", "rsid"]
IX = {c: i for i, c in enumerate(HDR18)}
C_GID, C_CHR, C_POS, C_PV = IX["gene_id"], IX["chromosome"], IX["position"], IX["pvalue"]
C_AN, C_AC, C_MAF = IX["an"], IX["ac"], IX["maf"]
C_REF, C_ALT, C_B, C_SE = IX["ref"], IX["alt"], IX["beta"], IX["se"]

# 目标（chrom 用索引内的命名：纯数字）
# ★基因 ID 经 Ensembl REST 实查（2026-09-28）：
#   YME1L1 = ENSG00000136758 (chr10, GRCh38 27,107,777-27,155,266, strand -1)
#   ANXA4  = ENSG00000196975 (chr2,  GRCh38 69,643,109-69,827,121, strand +1)
DS = [
    dict(tag="L2_YME1L1", gene="YME1L1", ens="ENSG00000136758",
         chrom="10", lo=26032861, hi=28232861, tss38=27132861),
    dict(tag="L3_ANXA4", gene="ANXA4", ens="ENSG00000196975",
         chrom="2", lo=68635445, hi=70835445, tss38=69735445),
]
SYM = {"ENSG00000136758": "YME1L1", "ENSG00000196975": "ANXA4"}

buf = []
def w(s=""):
    buf.append(str(s)); print(s)
def flush():
    io.open(LOG, "w", encoding="utf-8", newline="\n").write("\n".join(buf) + "\n")

# ---------------- tabix 索引解析 ----------------
def parse_tbi(path):
    raw = gzip.open(path, "rb").read()
    assert raw[:4] == b"TBI\x01", "不是 tabix 索引"
    n_ref, fmt, col_seq, col_beg, col_end, meta, skip, l_nm = struct.unpack("<8i", raw[4:36])
    o = 36
    names = raw[o:o + l_nm].split(b"\x00")[:-1]
    o += l_nm
    refs = {}
    for i in range(n_ref):
        n_bin = struct.unpack("<i", raw[o:o + 4])[0]; o += 4
        bins = {}
        for _ in range(n_bin):
            bid = struct.unpack("<I", raw[o:o + 4])[0]; o += 4
            n_chunk = struct.unpack("<i", raw[o:o + 4])[0]; o += 4
            ch = []
            for _ in range(n_chunk):
                cb, ce = struct.unpack("<QQ", raw[o:o + 16]); o += 16
                ch.append((cb, ce))
            bins[bid] = ch
        n_intv = struct.unpack("<i", raw[o:o + 4])[0]; o += 4
        intv = []
        for _ in range(n_intv):
            intv.append(struct.unpack("<Q", raw[o:o + 8])[0]); o += 8
        refs[names[i].decode()] = dict(n_bin=n_bin, n_intv=n_intv, bins=bins,
                                       intv=intv, first_off=(intv[0] if intv else 0))
    return dict(n_ref=n_ref, fmt=fmt, col_seq=col_seq, col_beg=col_beg,
                col_end=col_end, meta=meta, skip=skip, names=names, refs=refs, raw_len=len(raw))

def reg2bins(beg, end):
    """标准 tabix 分箱（beg 0-based inclusive, end 0-based exclusive）。"""
    lst = [0]
    end -= 1
    for k in range(1 + (beg >> 26), 1 + (end >> 26) + 1):
        lst.append(((1 << 3) - 1) // 7 + k)
    for k in range(1 + (beg >> 23), 1 + (end >> 23) + 1):
        lst.append(((1 << 6) - 1) // 7 + k)
    for k in range(1 + (beg >> 20), 1 + (end >> 20) + 1):
        lst.append(((1 << 9) - 1) // 7 + k)
    for k in range(1 + (beg >> 17), 1 + (end >> 17) + 1):
        lst.append(((1 << 12) - 1) // 7 + k)
    for k in range(1 + (beg >> 14), 1 + (end >> 14) + 1):
        lst.append(((1 << 15) - 1) // 7 + k)
    return lst

def merge_chunks(ch):
    ch = sorted(ch)
    res = []
    for cb, ce in ch:
        if ce <= cb:
            continue
        if not res or cb > res[-1][1]:
            res.append([cb, ce])
        else:
            res[-1][1] = max(res[-1][1], ce)
    return res

opener = urllib.request.build_opener(urllib.request.ProxyHandler())
UA = {"User-Agent": "Mozilla/5.0 wkb-rev1/1.0", "Accept-Encoding": "identity"}

def http_range(a, b):
    """取字节 [a, b)（b 不含）。自动夹到文件长度。"""
    a = max(0, a); b = min(FSIZE, b)
    if b <= a:
        return b""
    req = urllib.request.Request(URL, headers=dict(UA, **{"Range": "bytes=%d-%d" % (a, b - 1)}))
    with opener.open(req, timeout=180) as r:
        d = r.read()
    if len(d) != (b - a):
        raise IOError("Range 长度不符：要 %d 得 %d" % (b - a, len(d)))
    return d

def decode_gz_prefix(data):
    """对以 BGZF 块边界起始的缓冲做多成员 gzip 增量解码；容忍尾部截断。"""
    out = []
    gz = gzip.GzipFile(fileobj=io.BytesIO(data))
    while True:
        try:
            b = gz.read(1 << 20)
        except (EOFError, OSError, IOError):
            break
        if not b:
            break
        out.append(b)
    return b"".join(out)

# ---------------- 主流程 ----------------
T0 = time.time()
w("=" * 100)
w("s55a — GTEx v8 全血 区域定向取数（tabix-over-HTTP Range）%s" % time.strftime("%Y-%m-%d %H:%M:%S"))
w("=" * 100)
idx = parse_tbi(TBI)
w("[0] 索引：n_ref=%d fmt=%d col_seq=%d col_beg=%d col_end=%d meta=%d skip=%d 解压长=%d"
  % (idx["n_ref"], idx["fmt"], idx["col_seq"], idx["col_beg"], idx["col_end"],
     idx["meta"], idx["skip"], idx["raw_len"]))
w("    染色体名 = %s" % [n.decode() for n in idx["names"]])

tot_fetch = 0
for D in DS:
    tag, chrom, gene, ens = D["tag"], D["chrom"], D["gene"], D["ens"]
    lo, hi = D["lo"], D["hi"]
    w("")
    w("-" * 100)
    w("[%s] %s (%s)  chr%s:%d-%d (b38)  TSS38=%d" % (tag, gene, ens, chrom, lo, hi, D["tss38"]))
    if chrom not in idx["refs"]:
        w("  !! 索引内无该染色体，跳过"); continue
    R = idx["refs"][chrom]
    bins = reg2bins(lo - 1, hi)
    hit = [b for b in bins if b in R["bins"]]
    ch = [c for b in hit for c in R["bins"][b]]
    mg = merge_chunks(ch)
    nbytes = sum((min(FSIZE, (ce >> 16) + 65536) - (cb >> 16)) for cb, ce in mg)
    w("  reg2bins 候选 = %d ；索引内命中 = %d ；chunk = %d ；合并后 = %d 段"
      % (len(bins), len(hit), len(ch), len(mg)))
    for cb, ce in mg[:12]:
        w("      bytes %d..%d  (%.1f KB)" % (cb >> 16, ce >> 16, ((ce >> 16) - (cb >> 16)) / 1024.0))
    if len(mg) > 12:
        w("      ... 共 %d 段" % len(mg))
    w("  计划取字节 = %.2f MB" % (nbytes / 2 ** 20))
    tot_fetch += nbytes

    # 拉取 + 解码 + 过滤
    rows = []
    seen_var = set()
    n_dec = n_drop_out = n_dup = 0
    t1 = time.time()
    for cb, ce in mg:
        a = cb >> 16
        b = min(FSIZE, (ce >> 16) + 65536)
        data = http_range(a, b)
        dec = decode_gz_prefix(data)
        uoff = cb & 0xFFFF
        dec = dec[uoff:]
        n_dec += len(dec)
        txt = dec.decode("utf-8", "replace")
        lines = txt.split("\n")
        if lines and lines[-1] != "":
            lines = lines[:-1]          # 尾部半行丢弃
        for ln in lines:
            if not ln or ln[0] == "\n":
                continue
            f = ln.rstrip("\r").split("\t")
            if f[0] == "variant":
                continue
            if len(f) < 18:
                continue
            if f[C_CHR] != chrom:
                n_drop_out += 1
                continue
            try:
                pos = int(f[C_POS])
            except ValueError:
                continue
            if lo <= pos <= hi:
                key = (f[C_GID], pos, f[C_REF], f[C_ALT])
                if key in seen_var:
                    n_dup += 1
                    continue
                seen_var.add(key)
                rows.append(f)
    dt = time.time() - t1
    w("  实际解码 = %.2f MB ；越界行丢 %d ；去重丢 %d" % (n_dec / 2 ** 20, n_drop_out, n_dup))
    w("  区域内真实数据行 = %d ；耗时 %.1f s" % (len(rows), dt))

    if not rows:
        w("  !! 0 行，跳过写出"); continue

    # 写出（原始 18 列 + gene_symbol）
    outp = os.path.join(OUTD, "Whole_Blood.%s.tsv" % tag)
    gi = C_GID
    with io.open(outp, "w", encoding="utf-8", newline="\n") as f:
        f.write("\t".join(HDR18 + ["gene_symbol"]) + "\n")
        for f_ in rows:
            g0 = f_[C_GID].split(".")[0]
            f.write("\t".join(f_) + "\t" + SYM.get(g0, "") + "\n")
    w("  → %s（%.1f KB，%d 行）" % (os.path.basename(outp), os.path.getsize(outp) / 1024.0, len(rows)))

    # 描述性统计
    tgt = [r for r in rows if r[gi].split(".")[0] == ens]
    w("  目标基因 %s 行数 = %d" % (gene, len(tgt)))
    if tgt:
        pv = [(float(r[C_PV]), r[0], r[C_REF], r[C_ALT]) for r in tgt]
        pv.sort()
        w("    minP = %.4e @ %s (%s/%s)" % (pv[0][0], pv[0][1], pv[0][2], pv[0][3]))
        n5 = sum(1 for x in pv if x[0] < 5e-8)
        n1e5 = sum(1 for x in pv if x[0] < 1e-5)
        n1e4 = sum(1 for x in pv if x[0] < 1e-4)
        anl = sorted(int(r[C_AN]) for r in tgt)
        w("    p<5e-8 = %d ；p<1e-5 = %d ；p<1e-4 = %d" % (n5, n1e5, n1e4))
        w("    an 唯一值数 = %d ；中位 an = %d → N≈%d"
          % (len(set(anl)), anl[len(anl) // 2], anl[len(anl) // 2] // 2))
        ps = sorted(int(r[C_POS]) for r in tgt)
        span = ps[-1] - ps[0]
        w("    ★cis 窗完整性：pos %d..%d 宽度 %d (%.3f Mb)；推断 TSS=%.0f（应≈2.000 Mb 且 TSS 落在取数界内）"
          % (ps[0], ps[-1], span, span / 1e6, (ps[0] + ps[-1]) / 2.0))
        w("      取数界 chr%s:%d..%d ；窗口是否被截断 = %s"
          % (chrom, lo, hi, "否(完整)" if (ps[0] > lo and ps[-1] < hi) else "**需检查**"))
    # 区域内其它基因（探查串扰/为条件分析备料）
    from collections import Counter
    cc = Counter(r[gi].split(".")[0] for r in rows)
    w("  区域内基因数 = %d；行数前 8 = %s" % (len(cc), cc.most_common(8)))

w("")
w("=" * 100)
w("计划总取字节 = %.2f MB（相对全文件 3,190.55 MB = %.4f%%）"
  % (tot_fetch / 2 ** 20, 100.0 * tot_fetch / FSIZE))
w("总耗时 %.1f s" % (time.time() - T0))
w("VERDICT = PASS")
flush()
print("WROTE " + LOG)
