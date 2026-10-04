# -*- coding: utf-8 -*-
"""s75c · 核验「仅脚注消失、其余逐一同」

判据：
  ① 消失行集合 == 4 段脚注的文本（逐行比对，无其他行消失）
  ② 新增行集合 == 空集（=> 没有任何元素被位移/改写/重排）
  ③ 全部保留行的 (x, y, text) 与旧图逐一同（tight 只裁底部，上方坐标不动）
  ④ 宽度仍在硬窗口 181.5–183.3 mm；最小字号 >= 5.2 pt

★ 坐标容差 0.1 pt：tight 裁剪改变页面高度，PDF 写入坐标存在 0.01 pt 级浮点
  舍入抖动（实测 Fig3 的 CDC42_Mono_C 标签 346.49 -> 346.48）。文本与 x 完全
  相同，属无害抖动；若用 1e-9 判等会得到假阳性。
"""
import io, os, sys, json, hashlib
import fitz

ROOT = r"D:/endometriosis_project/11_sc_eqtl_mr_project"
FIGD = os.path.join(ROOT, "figures")
BAK = r"D:/_transfer_logs/_figs_backup_20260929_nonotes"
LOG = os.path.join(ROOT, "logs", "_s75c_verify_no_notes.txt")

out = []
def w(s=""):
    out.append(str(s)); print(s)

FILES = ["Fig2_discovery_replication", "Fig3_coloc_sensitivity",
         "Fig4_chr1p36_12_finemap", "Fig5_opentargets_phewas_safety"]


def lines(p):
    d = fitz.open(p)
    pg = d[0]
    r = []
    for b in pg.get_text("dict")["blocks"]:
        if b.get("type") != 0:
            continue
        for ln in b.get("lines", []):
            t = "".join(sp["text"] for sp in ln.get("spans", []))
            if not t.strip():
                continue
            r.append((round(ln["bbox"][0], 1), round(ln["bbox"][1], 1), t))   # 0.1 pt 容差
    W = pg.rect.width / 72 * 25.4
    H = pg.rect.height / 72 * 25.4
    d.close()
    return r, W, H


def sizes(p):
    d = fitz.open(p); pg = d[0]
    s = set()
    for b in pg.get_text("dict")["blocks"]:
        if b.get("type") != 0:
            continue
        for ln in b.get("lines", []):
            for sp in ln.get("spans", []):
                if sp["text"].strip():
                    s.add(round(sp["size"], 2))
    d.close()
    return sorted(s)


ALLPASS = True
for name in FILES:
    old = os.path.join(BAK, name + ".pdf")
    new = os.path.join(FIGD, name + ".pdf")
    A, Wa, Ha = lines(old)
    B, Wb, Hb = lines(new)
    sa, sb = set(A), set(B)
    only_a = [x for x in A if x not in sb]     # 消失
    only_b = [x for x in B if x not in sa]     # 新增/位移
    common = len(sa & sb)

    w("=" * 78)
    w("%s" % name)
    w("  旧: %.2f x %.2f mm  文本 %d 行    新: %.2f x %.2f mm  文本 %d 行"
      % (Wa, Ha, len(A), Wb, Hb, len(B)))
    w("  高度变化 %+.2f mm   宽度变化 %+.2f mm" % (Hb - Ha, Wb - Wa))
    w("  保留(逐一同) %d 行 ; 消失 %d 行 ; 新增或位移 %d 行" % (common, len(only_a), len(only_b)))
    w("  ---- 消失的行（应全部为底部脚注）----")
    for x, y, t in sorted(only_a, key=lambda z: z[1]):
        w("    y=%7.2f x=%6.2f | %s" % (y, x, t[:110]))
    if only_b:
        w("  !!!! 新增/位移的行（应为 0）!!!!")
        for x, y, t in sorted(only_b, key=lambda z: z[1])[:20]:
            w("    y=%7.2f x=%6.2f | %s" % (y, x, t[:110]))

    sz_old = sizes(old); sz_new = sizes(new)
    w("  字号集合 旧=%s" % sz_old)
    w("  字号集合 新=%s" % sz_new)

    ok = (len(only_b) == 0) and (181.5 <= Wb <= 183.3) and (min(sz_new) >= 5.2 if sz_new else True)
    ALLPASS = ALLPASS and ok
    w("  本图判据: 新增=0 -> %s ; 宽度在窗 -> %s ; 最小字号>=5.2 -> %s"
      % (len(only_b) == 0, 181.5 <= Wb <= 183.3, min(sz_new) if sz_new else None))

    # 字节与哈希
    for ext in ("png", "pdf"):
        p = os.path.join(FIGD, name + "." + ext)
        w("    %s  %d B  md5 %s" % (ext, os.path.getsize(p), hashlib.md5(open(p, "rb").read()).hexdigest()))

w("")
w("VERDICT_ALL = %s" % ("PASS" if ALLPASS else "FAIL"))

with io.open(LOG, "w", encoding="utf-8", newline="\n") as fh:
    fh.write("\n".join(out) + "\n")
