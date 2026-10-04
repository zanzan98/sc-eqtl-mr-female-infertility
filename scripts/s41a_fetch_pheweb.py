# -*- coding: utf-8 -*-
"""s41a_fetch_pheweb.py -- 任务 3.2 步骤1：抓取 FinnGen R12 三个 SNP 的全端点关联结果

三个 SNP（GRCh38，已由 Ensembl REST 核验）：
  rs12037376  1:22135618  ref=G  （CDC42 / B_MEM 工具；多等位 G/A/C，R12 中为 G-A）
  rs10917151  1:22096228  ref=G  （CDC42 / Mono_NC 工具；多等位 G/A/C/T）
  rs56318008  1:22143914  ref=C  （WNT4 可信集 lead / 条件变异；C/T）

抓取路径：GET https://r12.finngen.fi/api/variant/{chr}-{pos}-{ref}-{alt}
  说明：Windows curl(schannel) 因吊销列表脱机报 CRYPT_E_REVOCATION_OFFLINE；
        Python urllib 走自带 CA，可正常访问（已实测 200）。

产出（目录 00_data_raw/finngen_r12_pheweb/）：
  var_<rsid>.json            原始 JSON（缓存，避免重复请求）
  var_<rsid>.tsv             扁平化后的 results 表
  _fetch_log.txt

纪律：只读远程 + 只新增本地产物；不使用 Bash 管道写中文；一切结果落盘再读。
"""
import io
import json
import os
import ssl
import sys
import time
import urllib.error
import urllib.request

import pandas as pd

PROJ = r"D:\endometriosis_project\11_sc_eqtl_mr_project"
COND = os.path.join(PROJ, "00_data_raw", "onek1k", "cond_coloc")
OUT = os.path.join(PROJ, "00_data_raw", "finngen_r12_pheweb")
LOG = os.path.join(OUT, "_fetch_log.txt")
CTX = ssl.create_default_context()

# rsid -> (chr, pos38, ref_alleles_to_try)
TARGETS = {
    "rs12037376": ("1", 22135618, ["G-A", "G-C", "A-G", "C-G"]),
    "rs10917151": ("1", 22096228, ["G-A", "G-C", "G-T", "A-C", "C-T", "A-T"]),
    "rs56318008": ("1", 22143914, ["C-T", "T-C"]),
}
L = []


def p(s=""):
    L.append(str(s))
    print(s)


def fetch(url, timeout=180):
    req = urllib.request.Request(url, headers={"User-Agent": "phewas-3.2/1.0"})
    with urllib.request.urlopen(req, timeout=timeout, context=CTX) as r:
        return r.status, r.read()


def free_gb():
    import subprocess
    try:
        out = subprocess.check_output(
            ["powershell", "-NoProfile", "-Command",
             "(Get-CimInstance Win32_OperatingSystem).FreePhysicalMemory"],
            text=True, stderr=subprocess.DEVNULL, timeout=60).strip()
        return int(out) / 1024.0 / 1024.0
    except Exception:
        return float("nan")


T0 = time.time()
os.makedirs(OUT, exist_ok=True)
p("=== s41a 抓取 FinnGen R12 PheWeb（任务 3.2 步骤1）===")
p("time %s   可用内存 %.2f GB（闸门：>14 GB 系统已用即停）" % (time.strftime("%Y-%m-%d %H:%M:%S"), free_gb()))

# ---- 0) 从既有 harmonise 表取三个 SNP 的 bim 等位与 af（本地核验，不请求网络） ----
p("")
p("[0] 本机等位核验（m_*.csv，bim A1/A2）")
mB = pd.read_csv(os.path.join(COND, "m_B_MEM.csv"))
mM = pd.read_csv(os.path.join(COND, "m_Mono_NC.csv"))
POS37 = {  # rsid -> GRCh37 "1:pos"
    "rs12037376": "1:22462111",
    "rs10917151": "1:22422721",
    "rs56318008": "1:22470407",
}
for rs, sid in POS37.items():
    qb = mB[mB["snp_grch37"] == sid]
    qm = mM[mM["snp_grch37"] == sid]
    if len(qb):
        r0 = qb.iloc[0]
        p("  %-12s bim a1=%s a2=%s | eqtl af=%.4f slope=%.6f se=%.6f | pos38=%s"
          % (rs, r0["a1"], r0["a2"], r0["af"], r0["slope"], r0["slope_se"], r0["pos38"]))
    else:
        p("  %-12s !! 不在 B_MEM 交集内" % rs)
    if len(qm) and not len(qb):
        r1 = qm.iloc[0]
        p("     （Mono_NC: a1=%s a2=%s slope=%.6f se=%.6f）" % (r1["a1"], r1["a2"], r1["slope"], r1["slope_se"]))

# ---- 1) 抓取 ----
p("")
p("[1] 抓取 action=variant（每次 ~0.9 MB，共 3 个 SNP）")
meta = {}
for rs, (ch, pos, tries) in TARGETS.items():
    jf = os.path.join(OUT, "var_%s.json" % rs)
    if os.path.exists(jf) and os.path.getsize(jf) > 10000:
        p("  %-12s 已缓存，跳过（%d bytes）" % (rs, os.path.getsize(jf)))
        meta[rs] = json.load(io.open(jf, encoding="utf-8"))["_meta_used"]
        continue
    got = None
    for al in tries:
        url = "https://r12.finngen.fi/api/variant/%s-%d-%s" % (ch, pos, al)
        try:
            st, body = fetch(url)
            p("  %-12s %s -> status=%s bytes=%d" % (rs, al, st, len(body)))
            if st == 200 and len(body) > 10000:
                got = (al, body)
                break
        except urllib.error.HTTPError as e:
            p("  %-12s %s -> HTTP %s" % (rs, al, e.code))
        except Exception as e:
            p("  %-12s %s -> %s: %s" % (rs, al, type(e).__name__, e))
        time.sleep(0.5)
    if got is None:
        p("  !! %s 全部等位组合均失败" % rs)
        continue
    al, body = got
    js = json.loads(body)
    js["_meta_used"] = al
    io.open(jf, "w", encoding="utf-8", newline="\n").write(json.dumps(js))
    meta[rs] = al
    p("     -> 落地 %s（keys=%s）" % (jf, sorted(js.keys())))
    time.sleep(0.5)

# ---- 2) 扁平化 ----
p("")
p("[2] 扁平化 results")
summary = []
for rs in TARGETS:
    jf = os.path.join(OUT, "var_%s.json" % rs)
    if not os.path.exists(jf):
        p("  %-12s 无缓存，跳过" % rs)
        continue
    js = json.load(io.open(jf, encoding="utf-8"))
    res = js.get("results", [])
    p("  %-12s results 条数 = %d   used_alleles=%s" % (rs, len(res), js.get("_meta_used")))
    if res:
        p("     字段 = %s" % sorted(res[0].keys()))
        p("     首条 = %s" % json.dumps(res[0], ensure_ascii=False)[:400])
    rows = []
    for it in res:
        d = dict(it)
        for k in ("pheno", "phenocode", "phenostring", "category"):
            if k in d and isinstance(d[k], dict):
                d[k] = json.dumps(d[k], ensure_ascii=False)
        rows.append(d)
    df = pd.DataFrame(rows)
    tf = os.path.join(OUT, "var_%s.tsv" % rs)
    df.to_csv(tf, sep="\t", index=False, encoding="utf-8")
    summary.append((rs, len(df), tf))
    p("     -> %s  rows=%d cols=%d" % (tf, len(df), df.shape[1]))

# ---- 3) 与既有结果交叉核对（回归测试锚点） ----
p("")
p("[3] 与既有 42_phewas_safety_assessment.csv 交叉核对（回归锚点）")
old = pd.read_csv(os.path.join(PROJ, "tables", "42_phewas_safety_assessment.csv"))
p("  既有显著信号行 = %d" % len(old))
for rs in ["rs12037376", "rs10917151"]:
    tf = os.path.join(OUT, "var_%s.tsv" % rs)
    if not os.path.exists(tf):
        continue
    df = pd.read_csv(tf, sep="\t")
    sub = old[old["rsid"] == rs]
    if not len(sub):
        continue
    key = "phenocode" if "phenocode" in df.columns else ("pheno" if "pheno" in df.columns else None)
    if key is None:
        p("  %s: 无 phenocode 列，跳过核对" % rs)
        continue
    mm = sub.merge(df, left_on="endpoint", right_on=key, how="left")
    hit = mm["beta"].notna().sum()
    dev_b = (mm["beta"] - mm["beta_out"]).abs()
    p("  %-12s 既有 %d 行，在 R12 API 命中 %d 行；|Δbeta_out| max = %.3e（>1e-6 的个数 %d）"
      % (rs, len(sub), int(hit), dev_b.max(), int((dev_b > 1e-6).sum())))

p("")
p("产物目录 = %s" % OUT)
p("总耗时 %.1f s   末次可用内存 %.2f GB" % (time.time() - T0, free_gb()))
io.open(LOG, "w", encoding="utf-8", newline="\n").write("\n".join(L) + "\n")
print("WROTE " + LOG)
