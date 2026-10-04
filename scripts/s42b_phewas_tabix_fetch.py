# -*- coding: utf-8 -*-
"""s42b_phewas_tabix_fetch.py -- 任务 3.2 数据层：用 HTTP Range + tabix 取**全端点**面板

为什么需要它
-----------
本机没有 FinnGen R12 本地汇总统计，也没有原始 4,938 行全量表（后者在大内存机/已失联的 E 盘）。
R12 PheWeb 的 variant API 只对 ~256/2,470 端点返回结果 → 无法重算全面板 FDR。
R12 公开桶 (storage.googleapis.com/finngen-public-data-r12) 的 summary_stats 是 **BGZF** 且带 **.tbi**，
支持 HTTP Range → 可按 tabix 只取目标位置的极小字节范围。

做法
----
对每个端点：
  ① 递增 Range 取 .tbi 前缀，解析出 chr1 的 bins 与 **linear index (intv)**；
  ② 由 intv[win] / intv[win+1] 算出含目标位点的 BGZF 字节范围（~0.1-0.2 MB）；
  ③ Range 取该范围，**逐 BGZF 块** gzip 解压（块独立，容忍末块截断），按 pos 抽 3 行；
  ④ 增量追加 CSV（崩溃可续跑）。

目标 SNP（GRCh38，Ensembl 已核验）：1:22135618:G:A、1:22096228:G:A、1:22143914:C:T

产出：00_data_raw/finngen_r12_tabix/phewas_3snp_all_endpoints.csv
      （endpoint, source_file, pos38, ref, alt, rsids, nearest_genes, pval, mlogp, beta, sebeta, af_alt, ...）
      以及 phewas_3snp_fetch_status.csv（逐端点状态，失败必须留痕）

用法：python s42b_phewas_tabix_fetch.py [--limit N] [--workers 1]
"""
import argparse
import csv
import gzip
import http.client
import io
import json
import os
import socket
import struct
import ssl
import sys
import threading
import time
import urllib.error
import urllib.request

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

PROJ = r"D:\endometriosis_project\11_sc_eqtl_mr_project"
OUTD = os.path.join(PROJ, "00_data_raw", "finngen_r12_tabix")
MANIFEST = r"D:\endometriosis_project\01_data_raw\finngen_R12_manifest.tsv"
BASE = "https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release/"
LOG = os.path.join(OUTD, "_fetch_log.txt")
OUT = os.path.join(OUTD, "phewas_3snp_all_endpoints.csv")
STAT = os.path.join(OUTD, "phewas_3snp_fetch_status.csv")

TARGETS = [(22135618, "G", "A", "rs12037376"),
           (22096228, "G", "A", "rs10917151"),
           (22143914, "C", "T", "rs56318008")]
WINS = sorted((p - 1) >> 14 for p, _, _, _ in TARGETS)

CTX = ssl.create_default_context()
COLS = ["endpoint", "source_file", "pos38", "ref", "alt", "rsids", "nearest_genes",
        "pval", "mlogp", "beta", "sebeta", "af_alt", "af_alt_cases", "af_alt_controls"]
L = []
_lk = threading.Lock()


def p(s=""):
    with _lk:
        L.append(str(s))
        print(s)


HOST = "storage.googleapis.com"
_TL = threading.local()


def _proxy_addr():
    """取本机 HTTPS 代理（Windows 下 urllib 会读 IE 注册表；实测为 127.0.0.1:10808）。

    ★ 为什么必须走代理：本机对 googleapis.com 的**系统 DNS 解析被阻断**
      （`Resolve-DnsName storage.googleapis.com` / `socket.getaddrinfo` 均 gaierror 11001），
      只有经代理隧道才能访问 GCS 桶。直连会直接死在 DNS 上。
    """
    raw = (os.environ.get("HTTPS_PROXY") or os.environ.get("https_proxy")
           or urllib.request.getproxies().get("https") or "")
    if not raw:
        return None
    raw = raw.split("://")[-1].rstrip("/")
    if "@" in raw:
        raw = raw.rsplit("@", 1)[1]
    if ":" in raw:
        h, p = raw.rsplit(":", 1)
        return h, int(p)
    return raw, 80


PROXY = _proxy_addr()


def _tunnel(timeout=180):
    """建立（并复用）经代理 CONNECT 到 GCS 的 TLS 连接。

    标定实测：单请求 196 KB 用 urllib 需 2.2-2.8 s（每请求都要重新 CONNECT + TLS 握手）；
    每端点 2 次请求 = 2 次握手 → 复用后显著省时。
    ★ 直连回退：若未配置代理则直接 socket 连接（本地 DNS 可用时）。
    """
    c = getattr(_TL, "c", None)
    if c is not None:
        return c
    if PROXY:
        s = socket.create_connection(PROXY, timeout=timeout)
        s.sendall(("CONNECT %s:443 HTTP/1.1\r\nHost: %s:443\r\n"
                   "Proxy-Connection: keep-alive\r\n\r\n" % (HOST, HOST)).encode())
        buf = b""
        while b"\r\n\r\n" not in buf:
            d = s.recv(4096)
            if not d:
                raise IOError("proxy_closed_during_connect")
            buf += d
        if b" 200" not in buf.split(b"\r\n")[0]:
            s.close()
            raise IOError("proxy_connect_failed:%r" % buf[:60])
        sock = CTX.wrap_socket(s, server_hostname=HOST)
    else:
        sock = CTX.wrap_socket(socket.create_connection((HOST, 443), timeout=timeout),
                               server_hostname=HOST)
    c = http.client.HTTPSConnection(HOST, timeout=timeout, context=CTX)
    c.sock = sock                      # ★ 直接注入已建好的 TLS socket → request() 不再 connect()
    _TL.c = c
    return c


def ranged(url, start, end, timeout=180, tries=4):
    """HTTP Range 取 [start, end]；经代理复用连接，失败则重建连接重试。"""
    path = url.split(HOST, 1)[1]
    hdr = {"Range": "bytes=%d-%d" % (start, end),
           "User-Agent": "wkb-phewas/1.0",
           "Accept-Encoding": "identity"}
    last = None
    for k in range(tries):
        try:
            c = _tunnel(timeout)
            c.request("GET", path, headers=hdr)
            r = c.getresponse()
            body = r.read()
            if r.status not in (200, 206):
                raise IOError("HTTP %d (len=%d)" % (r.status, len(body)))
            return body
        except Exception as e:
            last = e
            try:
                _TL.c.close()
            except Exception:
                pass
            _TL.c = None
            time.sleep(1.0 + 1.5 * k)
    raise last


def bgzf_read(buf):
    out, i, n = [], 0, len(buf)
    while i + 18 <= n:
        if buf[i] != 0x1F or buf[i + 1] != 0x8B:
            i += 1
            continue
        bsize = struct.unpack("<H", buf[i + 16:i + 18])[0] + 1
        if bsize < 18 or i + bsize > n:
            break
        try:
            out.append(gzip.decompress(buf[i:i + bsize]))
        except Exception:
            pass
        i += bsize
    return b"".join(out)


def reg2bins(beg, end, min_shift=14, depth=5):
    bins, t, s = [], 0, min_shift + depth * 3
    for lvl in range(depth + 1):
        bins += list(range(t + (beg >> s), t + (end >> s) + 1))
        s -= 3
        t += 1 << ((lvl + 1) * 3)
    return bins


def parse_tbi(raw, max_refs=None):
    """解析 tabix；只解析前 max_refs 条 ref 即返回（chr1 是 ref 0 → 无需下完整份 .tbi）。
    不完整则抛 ValueError。返回 {'names','ref','end','n_ref'}"""
    if raw[:4] != b"TBI\x01":
        raise ValueError("not_tbi")
    n_ref, fmt, col_seq, col_beg, col_end, meta, skip, l_nm = struct.unpack("<8i", raw[4:36])
    names = raw[36:36 + l_nm].decode("utf-8", "replace").split("\x00")[:-1]
    limit = n_ref if max_refs is None else min(n_ref, max_refs)
    pos = 36 + l_nm
    ref = {}
    for i in range(limit):
        if pos + 4 > len(raw):
            raise ValueError("truncated_at_bins")
        n_bin = struct.unpack("<i", raw[pos:pos + 4])[0]
        pos += 4
        chunks_of = {}
        for _ in range(n_bin):
            if pos + 8 > len(raw):
                raise ValueError("truncated_in_bins")
            binid, n_chunk = struct.unpack("<ii", raw[pos:pos + 8])
            pos += 8
            if pos + 16 * n_chunk > len(raw):
                raise ValueError("truncated_in_chunks")
            chunks_of[binid] = [struct.unpack("<QQ", raw[pos + 16 * k:pos + 16 * (k + 1)])
                                for k in range(n_chunk)]
            pos += 16 * n_chunk
        if pos + 4 > len(raw):
            raise ValueError("truncated_at_intv")
        n_intv = struct.unpack("<i", raw[pos:pos + 4])[0]
        pos += 4
        if pos + 8 * n_intv > len(raw):
            raise ValueError("truncated_in_intv")
        intv = [struct.unpack("<Q", raw[pos + 8 * k:pos + 8 * (k + 1)])[0] for k in range(n_intv)]
        pos += 8 * n_intv
        ref[i] = (chunks_of, intv)
    return dict(names=names, ref=ref, end=pos, n_ref=n_ref)


def load_tbi_slices(gz_url, slices=(163840, 327680, 786432, 1572864, 3145728, 6291456)):
    """递增前缀取 .tbi，**只解析 chr1（ref 0）成功后即返回**；返回 (parsed, used_bytes)

    ★ gz_url 形如 .../finngen_R12_X.gz → .tbi 为 gz_url + '.tbi'
    ★ 只解析 ref 0：chr1 段解压后约 0.46 MB（压缩约 0.15 MB）→ 前缀 ~0.2 MB 即可；
      早期版本解析全部 23 条 ref，被迫下完整份 1.77 MB 的 .tbi。
    ★ 前缀下界已实测标定（_calib_s42b.py）：131,072 B 解析失败、**163,840 B 成功** → 取 163,840。
    """
    tbi_url = gz_url + ".tbi"
    cum = bytearray()
    for s in slices:
        if len(cum) >= s:
            continue
        need = s - len(cum)
        chunk = ranged(tbi_url, len(cum), len(cum) + need - 1)
        if not chunk:
            break
        cum += chunk
        try:
            # ★ 必须用逐块解压：gzip.decompress 对**截断的**多成员流直接抛 EOFError，
            #   会让 196 KB 前缀永远解析失败，从而被迫下完整份 1.77 MB 的 .tbi。
            raw = bgzf_read(bytes(cum))          # ★ bgzf_read 返回单个 bytes（不是元组）
            return parse_tbi(raw, max_refs=1), len(cum)
        except ValueError:
            continue
        except Exception:
            continue
    raise ValueError("tbi_unparsed")


def data_range(intv, win_lo, win_hi):
    """★ 必须从**最小窗口**取起点、从**最大窗口之后**取终点，
    否则会漏掉较早窗口（如 win 1348）的记录。"""
    n = len(intv)
    start_i = None
    for i in range(min(win_lo, n - 1), -1, -1):
        if intv[i] > 0:
            start_i = i
            break
    if start_i is None:
        return None
    end_i = None
    for i in range(min(win_hi, n - 1) + 1, n):
        if intv[i] > 0:
            end_i = i
            break
    start = intv[start_i] >> 16
    end = (intv[end_i] >> 16) + 2 * 65536 if end_i is not None else start + 6 * 65536
    return start, end


def fetch_endpoint(ep):
    """返回 (rows, status)；rows 为目标位点命中的行"""
    url = BASE + "finngen_R12_" + ep + ".gz"
    st = dict(endpoint=ep, status="ok", tbi_bytes=0, data_bytes=0, n_hit=0, note="")
    try:
        tb, used = load_tbi_slices(url)
        st["tbi_bytes"] = used
    except Exception as e:
        st["status"] = "tbi_fail"
        st["note"] = "%s:%s" % (type(e).__name__, str(e)[:60])
        return [], st
    i1 = tb["names"].index("1") if "1" in tb["names"] else None
    if i1 is None:
        st["status"] = "no_chr1"
        return [], st
    chunks_of, intv = tb["ref"][i1]
    rng = data_range(intv, min(WINS), max(WINS))
    if rng is None or rng[0] <= 0:
        st["status"] = "no_index"
        return [], st
    try:
        buf = ranged(url, rng[0], rng[1])
    except Exception as e:
        st["status"] = "data_fail"
        st["note"] = "%s:%s" % (type(e).__name__, str(e)[:60])
        return [], st
    st["data_bytes"] = len(buf)
    txt = bgzf_read(buf).decode("utf-8", "replace")
    want = {str(pos): (ref, alt, rs) for pos, ref, alt, rs in TARGETS}
    rows = []
    ncol_seen = 0
    for ln in txt.split("\n"):
        if not ln:
            continue
        f = ln.split("\t")
        # ★ 列数容错：case/control 端点为 13 列（…af_alt, af_alt_cases, af_alt_controls）；
        #   **定量端点（category = "Quantitative endpoints"，如 BMI_IRN/HEIGHT_IRN/WEIGHT_IRN，
        #   num_controls = 0）只有 11 列，没有 af_alt_cases / af_alt_controls**。
        #   早期版本硬写 `len(f) < 13: continue` → 这 3 个端点的 3 个目标位点被**静默丢弃**
        #   （status 记 partial_hit / hit=0），这也是"静默失败"型缺陷。
        if len(f) < 11:
            continue
        if f[1] in want:
            ncol_seen = max(ncol_seen, len(f))
            vals = (f[1:14] + [""] * 13)[:13]
            r = dict(zip(COLS, [ep, os.path.basename(url)] + vals))
            rows.append(r)
    st["n_hit"] = len(rows)
    if ncol_seen and ncol_seen < 13:
        st["note"] = "cols=%d(quantitative)" % ncol_seen
    if len(rows) < len(TARGETS):
        st["status"] = "partial_hit"
        st["note"] = "hit=%d" % len(rows)
    return rows, st


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--workers", type=int, default=1)
    args = ap.parse_args()

    os.makedirs(OUTD, exist_ok=True)
    fh = io.open(LOG, "w", encoding="utf-8", newline="\n")

    def pl(s=""):
        with _lk:
            L.append(str(s))
            print(s)
            fh.write(str(s) + "\n")
            fh.flush()

    man = []
    with io.open(MANIFEST, encoding="utf-8") as f:
        rd = csv.DictReader(f, delimiter="\t")
        for r in rd:
            man.append(r["phenocode"])
    done = set()
    if os.path.exists(STAT):
        with io.open(STAT, encoding="utf-8") as f:
            for r in csv.DictReader(f):
                if r.get("status") in ("ok", "no_chr1", "no_index"):
                    done.add(r["endpoint"])
    todo = [e for e in man if e not in done]
    # ★ --limit 必须在**过滤已完成之后**生效：否则重复 `--limit N` 第二轮 todo=0，
    #   分块续跑会直接空转（早期版本即犯此错）。
    if args.limit:
        todo = todo[:args.limit]
    pl("=== s42b 全端点取数 %s ===" % time.strftime("%Y-%m-%d %H:%M:%S"))
    pl("端点总数 = %d ；已完成（可续跑）= %d ；本轮计划 = %d" % (len(man), len(done), len(todo)))

    new_file = not os.path.exists(OUT)
    fout = io.open(OUT, "a", encoding="utf-8", newline="\n")
    w = csv.DictWriter(fout, fieldnames=COLS)
    if new_file:
        w.writeheader()
    new_stat = not os.path.exists(STAT)
    fs = io.open(STAT, "a", encoding="utf-8", newline="\n")
    ws = csv.DictWriter(fs, fieldnames=["endpoint", "status", "tbi_bytes", "data_bytes",
                                        "n_hit", "note"])
    if new_stat:
        ws.writeheader()

    tb_sum = db_sum = 0
    n_ok = n_bad = 0
    t0 = time.time()
    wlock = threading.Lock()

    def work(ep):
        try:
            return fetch_endpoint(ep)
        except Exception as e:
            return [], dict(endpoint=ep, status="exc", tbi_bytes=0, data_bytes=0,
                            n_hit=0, note="%s" % type(e).__name__)

    from concurrent.futures import ThreadPoolExecutor, as_completed
    nw = max(1, int(args.workers))
    pl("并发 worker = %d" % nw)
    with ThreadPoolExecutor(max_workers=nw) as ex:
        futs = {ex.submit(work, ep): ep for ep in todo}
        for k, fu in enumerate(as_completed(futs), 1):
            ep = futs[fu]
            rows, st = fu.result()
            with wlock:
                for r in rows:
                    w.writerow(r)
                ws.writerow(st)
                fout.flush(); fs.flush()
                tb_sum += st.get("tbi_bytes", 0) or 0
                db_sum += st.get("data_bytes", 0) or 0
                if st["status"] == "ok":
                    n_ok += 1
                else:
                    n_bad += 1
            if k % 25 == 0 or k == len(todo):
                el = time.time() - t0
                pl("  %5d/%d  用时 %.0fs  吞吐 %.0f KB/s  ETA %.0f min  ok=%d 其他=%d"
                   % (k, len(todo), el, (tb_sum + db_sum) / 1024 / max(el, 1e-9),
                      (len(todo) - k) * el / k / 60.0, n_ok, n_bad))
    pl("")
    pl("完成：ok=%d  非 ok=%d" % (n_ok, n_bad))
    pl("累计流量：.tbi %.1f MB ；.gz %.1f MB" % (tb_sum / 1e6, db_sum / 1e6))
    pl("总耗时 %.1f s" % (time.time() - t0))
    fout.close(); fs.close(); fh.close()
    print("WROTE " + LOG)


if __name__ == "__main__":
    main()
