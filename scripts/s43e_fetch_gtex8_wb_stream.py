# -*- coding: utf-8 -*-
"""s43e_fetch_gtex8_wb_stream.py — 选项 A：流式抓取 eQTL Catalogue GTEx v8 全血（区域过滤，提前停）

用户约束（严格遵守）：
  · 不加载整个文件；不做 BGZF 二进制手工解析；
  · 用 Python gzip 模块逐行读取；仅保留目标区域行；
  · 内存受限（可用 ~3 GB），必须流式。

关键前提（本轮实测）：该文件带 tabix 索引（.tbi n_ref=23, col_seq=13, col_beg=14），
  故按 (chromosome, position) 升序排列 → 流式读到 chr1 越过区域上界即可停止，
  无需下完 2.97 GB（实测速率 0.77 MB/s，全量需 ~66 min）。

抓取范围：chromosome == "1" 且 21,000,000 <= position <= 23,300,000（GRCh38）
          且 molecular_trait_id ∈ {WNT4, CDC42, LINC00339}
  范围理由：GTEx 的 cis 窗 = TSS ±1 Mb。三基因 TSS(b38)：
    LINC00339 22,024,558 / CDC42 22,052,627 / WNT4 22,143,969
    → 并集 = chr1:21,024,558–23,143,969。取 21.0–23.3 Mb 覆盖全部（含 SMR 所需完整 cis 窗）。
  ※ 用户原述 22.03–22.14 Mb 仅覆盖三基因基因体，会截断 cis 窗与 WNT4 的 7/42 显著对，故加宽，
    已在报告中显式说明。

输出：00_data_raw/eqtlcat_gtex_v8/Whole_Blood.chr1_21_23.3Mb_3genes.tsv（UTF-8, LF, 制表符）
"""
import os, io, sys, gzip, time, json, urllib.request, urllib.error, ctypes

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

PROJ = r"D:\endometriosis_project\11_sc_eqtl_mr_project"
OUTD = r"D:\endometriosis_project\01_data_raw\eqtlcat_gtex_v8"
OUT = os.path.join(OUTD, "Whole_Blood.chr1_21_23.3Mb_3genes.tsv")
LOG = os.path.join(PROJ, "scripts", "_s43e_fetch_gtex8_wb_stream.log")

URL = "https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/imported/GTEx_V8/ge/Whole_Blood.tsv.gz"
CHR = "1"
LO, HI = 21000000, 23300000
GENES = {
    "ENSG00000162552": "WNT4",
    "ENSG00000070831": "CDC42",
    "ENSG00000218510": "LINC00339",
}
WATCH_RSID = {"rs12037376", "rs10917151", "rs112678906",
              "rs56318008", "rs2473290", "rs2501299", "rs2473294"}

buf = []
def w(s=""):
    buf.append(str(s)); print(s)

def flush():
    os.makedirs(os.path.dirname(LOG), exist_ok=True)
    with io.open(LOG, "w", encoding="utf-8") as f:
        f.write("\n".join(buf) + "\n")

class MS(ctypes.Structure):
    _fields_ = [("dwLength", ctypes.c_ulong), ("dwMemoryLoad", ctypes.c_ulong),
                ("ullTotalPhys", ctypes.c_ulonglong), ("ullAvailPhys", ctypes.c_ulonglong),
                ("ullTotalPageFile", ctypes.c_ulonglong), ("ullAvailPageFile", ctypes.c_ulonglong),
                ("ullTotalVirtual", ctypes.c_ulonglong), ("ullAvailVirtual", ctypes.c_ulonglong),
                ("ullAvailExtendedVirtual", ctypes.c_ulonglong)]
def mem():
    st = MS(); st.dwLength = ctypes.sizeof(MS)
    ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(st))
    return st.ullAvailPhys / 1048576.0, st.dwMemoryLoad

class CountingReader:
    """包住 HTTP 响应，统计已消费字节数（供进度/提前停判断）。"""
    def __init__(self, fp):
        self.fp = fp; self.n = 0
    def read(self, size=-1):
        b = self.fp.read(size)
        self.n += len(b); return b
    def readable(self): return True
    def seekable(self): return False
    def tell(self): return self.n
    def close(self):
        try: self.fp.close()
        except Exception: pass

os.makedirs(OUTD, exist_ok=True)
OP = urllib.request.build_opener(urllib.request.ProxyHandler())
UA = {"User-Agent": "Mozilla/5.0 wkb-3p3/1.0", "Accept-Encoding": "identity"}

w("=" * 92)
w("s43e_fetch_gtex8_wb_stream.py — GTEx v8 全血 流式区域抓取（提前停）")
w("URL  = %s" % URL)
w("范围 = chromosome=='%s' 且 %d <= position <= %d ; 基因 = %s"
  % (CHR, LO, HI, ",".join(GENES.values())))
w("=" * 92)
av, load = mem()
w("起始内存: 可用 %.0f MB (占 %d%%)" % (av, load))

idx = None
n_lines = 0
n_kept = 0
n_region_rows = 0          # 区域内全部基因的行数（用于判定排序）
last_pos = -1
per_gene = {}
watch_hits = {}
stop_reason = "stream_end"
t0 = time.time()

resp = None
reader = None
gz = None
try:
    resp = OP.open(urllib.request.Request(URL, headers=UA), timeout=120)
    w("HTTP %s  Content-Length=%s" % (resp.status, resp.headers.get("Content-Length")))
    reader = CountingReader(resp)
    gz = gzip.GzipFile(fileobj=reader, mode="rb")
    fout = io.open(OUT, "w", encoding="utf-8", newline="\n")
    header_out = False
    while True:
        line = gz.readline()
        if not line:
            stop_reason = "gzip_stream_end"
            break
        s = line.decode("utf-8", "replace")
        if s.endswith("\n"):
            s = s[:-1]
        if s.endswith("\r"):
            s = s[:-1]
        if not s:
            continue

        if idx is None:
            # 第一行是表头（无 '#' 前缀）
            cols = s.split("\t")
            idx = {c: i for i, c in enumerate(cols)}
            need = ("chromosome", "position", "molecular_trait_id", "variant",
                    "ref", "alt", "beta", "se", "pvalue", "maf", "median_tpm", "rsid")
            missing = [c for c in need if c not in idx]
            if missing:
                w("!! 表头缺列: %s ；实际列 = %s" % (missing, cols))
                break
            w("表头列 (%d): %s" % (len(cols), " | ".join(cols)))
            # 输出表头（在末尾追加 gene_symbol 列）
            fout.write("\t".join(cols) + "\tgene_symbol\n")
            header_out = True
            continue

        n_lines += 1
        f = s.split("\t")
        chrom = f[idx["chromosome"]]
        pos = int(f[idx["position"]])
        gid = f[idx["molecular_trait_id"]]

        # 提前停：chr1 内越过上界
        if chrom == CHR and pos > HI:
            stop_reason = "passed_region_upper_bound"
            break
        # 提前停：离开 chr1（文件按 chrom,pos 升序）
        if chrom != CHR and last_pos > 0:
            stop_reason = "left_chromosome_1"
            break
        last_pos = pos

        if chrom == CHR and LO <= pos <= HI:
            n_region_rows += 1
            sym = GENES.get(gid)
            if sym:
                n_kept += 1
                fout.write(s + "\t" + sym + "\n")
                d = per_gene.setdefault(sym, {"n": 0, "pmin": 1.0, "pmin_var": "",
                                              "tpm": f[idx["median_tpm"]],
                                              "pos_min": pos, "pos_max": pos})
                d["n"] += 1
                d["pos_max"] = pos
                pv = float(f[idx["pvalue"]])
                if pv < d["pmin"]:
                    d["pmin"] = pv
                    d["pmin_var"] = f[idx["variant"]]
        rs = f[idx["rsid"]]
        if rs in WATCH_RSID:
            watch_hits.setdefault(rs, {"n": 0, "genes": set(), "rows": []})
            watch_hits[rs]["n"] += 1
            if GENES.get(gid):
                watch_hits[rs]["genes"].add(GENES[gid])
            if len(watch_hits[rs]["rows"]) < 4:
                watch_hits[rs]["rows"].append({
                    "gene": GENES.get(gid, gid), "variant": f[idx["variant"]],
                    "pos": pos, "ref": f[idx["ref"]], "alt": f[idx["alt"]],
                    "beta": f[idx["beta"]], "se": f[idx["se"]], "pvalue": f[idx["pvalue"]],
                    "maf": f[idx["maf"]],
                })

        if n_lines % 2000000 == 0:
            av2, l2 = mem()
            w("  ... 已读 %d 行 / 消费 %.1f MB / 区域内 %d 行 / 保留 %d 行 / 内存可用 %.0f MB (占%d%%) / %.1f min"
              % (n_lines, reader.n / 2**20, n_region_rows, n_kept, av2, l2, (time.time() - t0) / 60.0))
    fout.close()
except Exception as e:
    stop_reason = "EXC:%s:%s" % (type(e).__name__, str(e)[:160])
finally:
    try:
        if gz: gz.close()
    except Exception:
        pass

dt = time.time() - t0
w()
w("-" * 92)
w("停止原因      = %s" % stop_reason)
w("网络消费字节  = %d B (%.1f MB)  用时 %.1f s → 平均 %.2f MB/s"
  % (reader.n if reader else 0, (reader.n if reader else 0) / 2**20, dt,
     (reader.n if reader else 0) / 2**20 / max(dt, 1e-6)))
w("已读数据行    = %d" % n_lines)
w("chr1 区域内行 = %d" % n_region_rows)
w("保留(3 基因)  = %d" % n_kept)
w()
w("逐基因（GTEx v8 全血，区域 %s:%d-%d b38）:" % (CHR, LO, HI))
for sym in ("WNT4", "CDC42", "LINC00339"):
    d = per_gene.get(sym)
    if not d:
        w("  %-10s 0 行  ← 该基因在此区域无任何被检验的 cis 变异" % sym)
        continue
    w("  %-10s n=%-7d median_tpm=%-10s minP=%.3e @ %s ; pos %d..%d"
      % (sym, d["n"], d["tpm"], d["pmin"], d["pmin_var"], d["pos_min"], d["pos_max"]))
w()
w("关注 rsID 命中情况:")
for rs in sorted(watch_hits):
    h = watch_hits[rs]
    w("  %-12s 行数=%-5d 命中基因=%s" % (rs, h["n"], sorted(h["genes"]) if h["genes"] else "(非三基因)"))
    for r in h["rows"]:
        w("        gene=%-10s %s pos=%d %s/%s beta=%s se=%s p=%s maf=%s"
          % (r["gene"], r["variant"], r["pos"], r["ref"], r["alt"], r["beta"], r["se"],
             r["pvalue"], r["maf"]))

av3, load3 = mem()
w()
w("结束内存: 可用 %.0f MB (占 %d%%)" % (av3, load3))
w("输出: %s  (%.1f KB)" % (OUT, (os.path.getsize(OUT) / 1024.0) if os.path.exists(OUT) else -1))
w("VERDICT = %s" % ("PASS" if n_kept > 0 and stop_reason in
                    ("passed_region_upper_bound", "left_chromosome_1", "gzip_stream_end") else "CHECK"))
flush()
