# -*- coding: utf-8 -*-
"""s42a_tabix_probe.py -- HTTP Range + tabix 单点取数（可行性验证，含逐块 BGZF 解压）

关键点：FinnGen R12 summary_stats 为 BGZF（每块是完整 gzip member，块长在头 16..18 字节）。
       末端块常被 Range 截断 → 必须逐块解压并容忍截断，不能用 gzip.decompress 整段。
"""
import gzip
import io
import os
import struct
import ssl
import sys
import time
import zlib
import urllib.request

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

LOG = r"D:\endometriosis_project\_s42a_tabix_probe.log"
BASE = "https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release/"
CTX = ssl.create_default_context()
L = []
_fh = None


def p(s=""):
    L.append(str(s))
    print(s)
    if _fh is not None:
        _fh.write(str(s) + "\n")
        _fh.flush()


def ranged(url, start, end, timeout=300):
    req = urllib.request.Request(url, headers={"Range": "bytes=%d-%d" % (start, end),
                                              "User-Agent": "tabix-probe/1.0"})
    with urllib.request.urlopen(req, timeout=timeout, context=CTX) as r:
        return r.read()


def bgzf_read(buf):
    """逐 BGZF 块解压，容忍垃圾/截断。返回解压后的 bytes。"""
    out = []
    i, n = 0, len(buf)
    nblk = 0
    while i + 18 <= n:
        if buf[i] != 0x1f or buf[i + 1] != 0x8b:
            i += 1
            continue
        bsize = struct.unpack("<H", buf[i + 16:i + 18])[0] + 1
        if bsize < 18 or i + bsize > n:
            break                      # 末块被截断
        blk = buf[i:i + bsize]
        i += bsize
        nblk += 1
        try:
            out.append(gzip.decompress(blk))
        except Exception:
            pass
    return b"".join(out), nblk


def reg2bins(beg, end, min_shift=14, depth=5):
    bins, l = [], 0
    t, s = 0, min_shift + depth * 3
    for lvl in range(depth + 1):
        bins += list(range(t + (beg >> s), t + (end >> s) + 1))
        s -= 3
        t += 1 << ((lvl + 1) * 3)
    return bins


EP = "finngen_R12_AB1_ACTINOMYCOSIS"
TARGETS = [(22135618, "G", "A"), (22096228, "G", "A"), (22143914, "C", "T")]

_fh = io.open(LOG, "w", encoding="utf-8", newline="\n")
p("=== s42a tabix 可行性验证 %s ===" % time.strftime("%Y-%m-%d %H:%M:%S"))
p("python=%s" % sys.version.split()[0])

# ---------- 1) .tbi ----------
tbi_url = BASE + EP + ".gz.tbi"
t0 = time.time()
tbi = ranged(tbi_url, 0, 20 * 10 ** 6)
t1 = time.time()
p("")
p("[1] .tbi：%d bytes  %.1f s  (%.0f KB/s)" % (len(tbi), t1 - t0, len(tbi) / 1024 / max(t1 - t0, 1e-9)))
raw = gzip.decompress(tbi)
p("    解压后 %d bytes" % len(raw))
magic = raw[:4]
n_ref, fmt, col_seq, col_beg, col_end, meta, skip, l_nm = struct.unpack("<8i", raw[4:36])
names = raw[36:36 + l_nm].decode().split("\x00")[:-1]
p("    magic=%r n_ref=%d fmt=%d col_seq=%d col_beg=%d col_end=%d meta=%s l_nm=%d"
  % (magic, n_ref, fmt, col_seq, col_beg, col_end, meta, l_nm))
p("    names 前 3 = %s ; '1' 索引 = %s" % (names[:3], names.index("1") if "1" in names else -1))

pos = 36 + l_nm
ref_intv = {}
ref_end_bytes = []
for i in range(n_ref):
    n_bin = struct.unpack("<i", raw[pos:pos + 4])[0]; pos += 4
    chunks_of = {}
    for _ in range(n_bin):
        binid, n_chunk = struct.unpack("<ii", raw[pos:pos + 8]); pos += 8
        chunks_of[binid] = [struct.unpack("<QQ", raw[pos + 16 * k: pos + 16 * (k + 1)])
                            for k in range(n_chunk)]
        pos += 16 * n_chunk
    n_intv = struct.unpack("<i", raw[pos:pos + 4])[0]; pos += 4
    intv = [struct.unpack("<Q", raw[pos + 8 * k: pos + 8 * (k + 1)])[0] for k in range(n_intv)]
    pos += 8 * n_intv
    ref_intv[i] = (chunks_of, intv)
    ref_end_bytes.append(pos)
p("    解析完成 pos=%d / %d" % (pos, len(raw)))
i1 = names.index("1")
p("    chr1: bins=%d n_intv=%d  解压后结束字节=%d（.tbi 压缩比 %.1fx）"
  % (len(ref_intv[i1][0]), len(ref_intv[i1][1]), ref_end_bytes[i1], len(raw) / max(len(tbi), 1)))

# ---------- 2) 目标窗口 ----------
chunks_of, intv = ref_intv[i1]
p("")
p("[2] 目标窗口与 virtual offset")
want = []
for pos38, ref, alt in TARGETS:
    beg0 = pos38 - 1
    win = beg0 >> 14
    bins = reg2bins(beg0, beg0)
    cks = [c for b in bins if b in chunks_of for c in chunks_of[b]]
    v0 = intv[win] if win < len(intv) else 0
    v1 = intv[win + 1] if win + 1 < len(intv) else 0
    p("    %d win=%d 命中chunks=%d intv[win]=%d intv[win+1]=%d  (blk=%d/%d)"
      % (pos38, win, len(cks), v0, v1, v0 >> 16, v1 >> 16))
    want.append(dict(pos38=pos38, ref=ref, alt=alt, win=win, chunks=cks, v0=v0, v1=v1))

# ---------- 3) 取 .gz Range ----------
vals = []
for w in want:
    vals += [w["v0"], w["v1"]] + [x for c in w["chunks"] for x in c]
vals = [v for v in vals if v > 0]
cstart, cend = (min(vals) >> 16), (max(vals) >> 16)
p("")
p("[3] .gz 字节范围 [%d, %d] = %.0f KB" % (cstart, cend, (cend - cstart) / 1024))
gz_url = BASE + EP + ".gz"
t0 = time.time()
buf = ranged(gz_url, cstart, cend + 131071)     # 多取 2 个块作保险
t1 = time.time()
p("    实取 %d bytes  %.1f s  (%.0f KB/s)" % (len(buf), t1 - t0, len(buf) / 1024 / max(t1 - t0, 1e-9)))
p("    首 4 字节 = %s（1f8b0804 = BGZF 块起点）" % buf[:4].hex())
txt, nblk = bgzf_read(buf)
p("    成功解压块数 = %d  解压后 %d bytes" % (nblk, len(txt)))
lines = txt.decode("utf-8", "replace").split("\n")
p("    行数 = %d" % len(lines))
p("    表头 = %s" % lines[0][:250])
hdr = lines[0].lstrip("#").split("\t")
p("    列名 = %s" % hdr)

# ---------- 4) 目标行 ----------
p("")
p("[4] 目标位置命中")
ip = hdr.index("pos") if "pos" in hdr else 1
hit = 0
for w in want:
    found = [ln for ln in lines[1:] if ln and ln.split("\t")[ip] == str(w["pos38"])]
    p("    pos=%s 命中 %d 行" % (w["pos38"], len(found)))
    for ln in found[:1]:
        p("      %s" % ln[:220])
    hit += len(found)
p("    合计 = %d / %d" % (hit, len(want)))
p("")
p("总耗时 %.1f s" % (time.time() - t0))
p("VERDICT = %s" % ("PASS" if hit == len(want) else "FAIL(未全命中)"))
_fh.close()
io.open(LOG + ".done", "w", encoding="utf-8", newline="\n").write("done\n")
print("WROTE " + LOG)
