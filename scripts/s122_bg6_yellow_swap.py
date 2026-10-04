# -*- coding: utf-8 -*-
"""s122_bg6_yellow_swap.py —— 给 Fig2(a) 的「黄色」换色（只算不改）。

用户反馈：BG6[0] = `#FFFF99` 看不清。
原因（量）：`#FFFF99` 的 CIE L* = 98.07，几乎贴白底；背景点只有 ms=1.05 pt，
             抗锯齿会把这种高亮低反差的小点"冲淡"成一片白黄。
原则：换成**明度更低、更实**的浅色，同时必须仍满足
  C1 对 6 个前景语义色 min dE >= 15
  C2 与其余 5 个背景色 min dE >= 15（不能和它们撞）
  C3 对白底 dE >= 25（在白纸上看得见）
"""
import io
import os
import sys

sys.path.insert(0, r"C:/Users/28144/.workbuddy/skills/cvd-safe-palette-design/scripts")
import cvd_sim as CV            # noqa: E402

LOG = r"D:/_transfer_logs/_p0_scratch/_s122_yellow_swap.txt"
L = []


def w(s=""):
    L.append(str(s))


FG = [("CDC42", "#519D78"), ("LINC00339", "#8BCF8B"), ("YME1L1", "#C4E9CA"),
      ("ANXA4", "#6CBAD8"), ("WNT4", "#367DB0"), ("region/nr", "#3D9F3C")]

# BG6 里除黄色之外的 5 色
KEEP = ["#CAB2D6", "#FDBF6F", "#FDDAEC", "#FFD9A0", "#FBB4AE"]

POOL = [
    ("原黄（要换掉）  #FFFF99", "#FFFF99"),
    ("Pastel1 #FFFFCC", "#FFFFCC"),
    ("user    #F2E29B", "#F2E29B"),
    ("Pastel1 #FED9A6", "#FED9A6"),
    ("Pastel1 #E5D8BD", "#E5D8BD"),
    ("user    #E0C9A6", "#E0C9A6"),
    ("Pastel1 #DECBE4", "#DECBE4"),
    ("user    #E5D5E8", "#E5D5E8"),
]
# 自生成：暖黄区（hue 40-72），L* 78-88，S 0.40-0.70
seen = set()
for hue in range(40, 73, 4):
    for Lt in (78, 81, 84, 87):
        for S in (0.40, 0.50, 0.60, 0.70):
            h = CV.hsl_hex_L(float(hue), float(Lt), float(S))
            if h not in seen:
                seen.add(h)
                POOL.append(("gen hue=%d L*=%d S=%.2f" % (hue, Lt, S), h))

w("=" * 112)
w("候选评估（KEEP 5 色 = %s）" % ", ".join(KEEP))
w("=" * 112)
w("%-28s %-9s %6s %8s %8s %8s  %s"
  % ("候选", "hex", "L*", "对前景", "对白底", "与其余5色", "判定"))
w("-" * 112)

rows = []
for name, h in POOL:
    Lc = CV.lab(h)[0]
    d_fg = min(CV.dE76(h, f) for _, f in FG)
    d_wh = CV.dE76(h, "#FFFFFF")
    d_kp = min(CV.dE76(h, k) for k in KEEP)
    fails = []
    if d_fg < 15:
        fails.append("前景 %.2f<15" % d_fg)
    if d_kp < 15:
        fails.append("与其余5色 %.2f<15" % d_kp)
    if d_wh < 25:
        fails.append("对白底 %.2f<25" % d_wh)
    tag = "通过" if not fails else "剔除: " + "; ".join(fails)
    rows.append((len(fails) == 0, h, Lc, d_fg, d_wh, d_kp, name))
    w("%-28s %-9s %6.2f %8.2f %8.2f %8.2f  %s" % (name, h, Lc, d_fg, d_wh, d_kp, tag))

ok = [r for r in rows if r[0]]
w("")
w("=" * 112)
w("通过的 %d 个，按「明度从低到高」排（越低越实、越看得清）：" % len(ok))
w("=" * 112)
ok.sort(key=lambda r: r[2])
for r in ok:
    w("   %-9s  L*=%6.2f  对前景 %6.2f  对白底 %6.2f  与其余5色 %6.2f   [%s]"
      % (r[1], r[2], r[3], r[4], r[5], r[6]))

w("")
w("★ 推荐（取通过者中 L* 最低、同时贴近原作浅色观感的 3 个）：")
for r in ok[:3]:
    w("   %s   L*=%5.2f   对前景 %5.2f   对白底 %5.2f" % (r[1], r[2], r[3], r[4]))

with io.open(LOG, "w", encoding="utf-8", newline="\n") as fh:
    fh.write("\n".join(L) + "\n")
print("done -> %s" % LOG)
