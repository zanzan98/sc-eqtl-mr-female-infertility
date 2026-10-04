# -*- coding: utf-8 -*-
"""s117_fig2_bgcolor_audit.py —— Fig2(a) 背景浅色系的客观审计（只读）。

问题：背景点要改成「浅色系彩色轮换」，但阳性点已占用 绿 / 蓝 两个色相族：
    CDC42     #519D78  中绿
    LINC00339 #8BCF8B  浅绿
    YME1L1    #C4E9CA  极浅绿
    ANXA4     #6CBAD8  浅蓝
    WNT4      #367DB0  深蓝
    region/nr #3D9F3C  深绿（显著性阈值线 + 红叉标记）
若背景再用浅绿/浅蓝，会与阳性点混。

本脚本对一批候选背景色，逐一计算与上述 6 个「前景语义色」的
  · CIE76 dE（原色）
  · CIE76 dE（protan / deutan / tritan 模拟）
  · |dL*|（明度差）
取对 6 色的最小值，判定是否安全。

口径（对齐项目既有纪律）：
  安全 = min dE >= 15（GB12 共现集口径） 且 min{protan,deutan} dE >= 12
  附加 = 背景色 L* 必须高于全部前景色（背景--前景的明度分层）
"""
import io
import os
import sys

sys.path.insert(0, r"C:/Users/28144/.workbuddy/skills/cvd-safe-palette-design/scripts")
import cvd_sim as CV            # noqa: E402

LOG = r"D:/_transfer_logs/_p0_scratch/_s117_bgcolor_audit.txt"
L = []


def w(s=""):
    L.append(str(s))


FG = [
    ("CDC42", "#519D78"),
    ("LINC00339", "#8BCF8B"),
    ("YME1L1", "#C4E9CA"),
    ("ANXA4", "#6CBAD8"),
    ("WNT4", "#367DB0"),
    ("region/nr", "#3D9F3C"),
    ("ns grey", "#B8B8B8"),
]

# 候选背景色池
CAND = [
    ("user-1  A6CEE3 浅蓝", "#A6CEE3"),
    ("user-2  FDBF6F 浅橙", "#FDBF6F"),
    ("user-3  B2DF8A 浅绿", "#B2DF8A"),
    ("user-4  CAB2D6 浅紫", "#CAB2D6"),
    ("user-5  FFD9A0 浅杏", "#FFD9A0"),
    ("user-6  FBB4AE 浅粉", "#FBB4AE"),
    ("alt-1   F2E29B 浅黄", "#F2E29B"),
    ("alt-2   E0C9A6 浅卡其", "#E0C9A6"),
    ("alt-3   E5D5E8 淡薰衣草", "#E5D5E8"),
    ("alt-4   F5D6C6 浅桃", "#F5D6C6"),
    ("alt-5   D9E2C0 浅黄绿", "#D9E2C0"),
    ("alt-6   CFD8E8 浅灰蓝", "#CFD8E8"),
    ("alt-7   EAD9A0 浅麦", "#EAD9A0"),
    ("alt-8   DCC7DE 浅藕荷", "#DCC7DE"),
    ("alt-9   F0C8B4 浅陶", "#F0C8B4"),
    ("grey    B0B0B0 现状灰A", "#B0B0B0"),
    ("grey    D5D5D5 现状灰B", "#D5D5D5"),
]

w("=" * 108)
w("Fig2(a) 背景浅色系 vs 前景语义色 —— 可辨性审计")
w("口径：安全 = min dE >= 15(原色) 且 min dE >= 12 (protan/deutan)；且背景 L* 高于全部前景")
w("=" * 108)
w("")
w("前景语义色:")
for n, c in FG:
    w("   %-11s %s   L* = %6.2f" % (n, c, CV.lab(c)[0]))
w("")
w("%-24s %-9s %6s %7s %7s %7s %6s  %-11s %s"
  % ("候选背景色", "hex", "L*", "min dE", "protan", "deutan", "tritan", "最近前景", "判定"))
w("-" * 108)

rows = []
for name, c in CAND:
    d_orig = []
    d_pr = []
    d_de = []
    d_tr = []
    for n, f in FG:
        d_orig.append((CV.dE76(c, f), n))
        d_pr.append(CV.dE76(c, f, "protan"))
        d_de.append(CV.dE76(c, f, "deutan"))
        d_tr.append(CV.dE76(c, f, "tritan"))
    mn, near = min(d_orig)
    Lc = CV.lab(c)[0]
    lmin = min(CV.lab(f)[0] for _, f in FG)
    darker_ok = Lc > lmin
    safe = (mn >= 15 and min(d_pr) >= 12 and min(d_de) >= 12 and darker_ok)
    rows.append((safe, name, c, Lc, mn, min(d_pr), min(d_de), min(d_tr), near, darker_ok))
    w("%-24s %-9s %6.2f %7.2f %7.2f %7.2f %6.2f  %-11s %s"
      % (name, c, Lc, mn, min(d_pr), min(d_de), min(d_tr), near,
         ("OK" if safe else "!! 风险") + ("" if darker_ok else " (更暗)")))
w("")
w("=" * 108)
w("汇总")
w("=" * 108)
ok = [r for r in rows if r[0]]
bad = [r for r in rows if not r[0]]
w("安全候选 %d / %d :" % (len(ok), len(rows)))
for r in ok:
    w("   %-24s %s   min dE %5.2f (最近 %s)" % (r[1], r[2], r[4], r[8]))
w("")
w("**风险/不合格候选 %d :**" % len(bad))
for r in bad:
    w("   %-24s %s   min dE %5.2f (最近 %s)%s"
      % (r[1], r[2], r[4], r[8], "" if r[9] else "  [明度也未高于前景]"))
w("")

# 候选之间（背景互相）的可辨性
w("=" * 108)
w("参考：候选池内部两两最小 dE（仅列安全候选之间）")
w("=" * 108)
okc = [(r[1], r[2]) for r in ok]
for i in range(len(okc)):
    for j in range(i + 1, len(okc)):
        w("   %-24s vs %-24s  dE = %5.2f" % (okc[i][0], okc[j][0],
                                             CV.dE76(okc[i][1], okc[j][1])))

with io.open(LOG, "w", encoding="utf-8", newline="\n") as fh:
    fh.write("\n".join(L) + "\n")
print("done -> %s" % LOG)
