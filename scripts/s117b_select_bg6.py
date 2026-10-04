# -*- coding: utf-8 -*-
"""s117b_select_bg6.py —— 为 Fig2(a) 选 6 色「背景浅色系」（只看不改）。

设计原则（不交给优化器自由生成，只在**成熟设计家族**内筛选）：
  候选池 = ColorBrewer Pastel1 / Paired 的浅色 + 用户点名的 2 个自定义浅色。
  这些色同属「柔和浅色」家族，观感统一；自动生成的灰浊色一律不用。

硬约束：
  C1  与 6 个前景语义色（CDC42/LINC00339/YME1L1/ANXA4/WNT4/region-nr）的 min CIE76 dE >= 15
  C3  与白底 dE >= 12（否则浅到在白纸上看不见）
  C4  入选 6 色彼此 min dE >= 15（相邻染色体能看出变化）
报告项（不作硬门槛，因背景点直径仅 ~1 pt、前景点 ~4 pt，尺寸本身提供分层）：
  R1  CIE L*（背景通常应比前景亮，但用户点名的浅紫 L*=75.6 / 浅粉 79.7 略低，不因此剔除）
  R2  protan / deutan 下与前景的 min dE（红绿色盲下的残余风险）

用户点名的 5 色（#FDBF6F #CAB2D6 #FBB4AE #FFD9A0 + 被 C1 判冲突的 #A6CEE3/#B2DF8A）优先入选，
只有在违反 C1/C3/C4 时才让位。
"""
import io
import os
import sys

sys.path.insert(0, r"C:/Users/28144/.workbuddy/skills/cvd-safe-palette-design/scripts")
import cvd_sim as CV            # noqa: E402

LOG = r"D:/_transfer_logs/_p0_scratch/_s117b_bg6_select.txt"
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
]

PRIORITY = ["#FDBF6F", "#FBB4AE", "#CAB2D6", "#FFD9A0"]

POOL = [
    ("Paired   #A6CEE3", "#A6CEE3", "浅蓝  <- 用户点名"),
    ("Paired   #B2DF8A", "#B2DF8A", "浅绿  <- 用户点名"),
    ("Paired   #FB9A99", "#FB9A99", "浅红"),
    ("Paired   #FDBF6F", "#FDBF6F", "浅橙  <- 用户点名"),
    ("Paired   #CAB2D6", "#CAB2D6", "浅紫  <- 用户点名"),
    ("Paired   #FFFF99", "#FFFF99", "浅黄"),
    ("Pastel1  #FBB4AE", "#FBB4AE", "浅粉  <- 用户点名"),
    ("Pastel1  #DECBE4", "#DECBE4", "浅紫"),
    ("Pastel1  #FED9A6", "#FED9A6", "浅橙"),
    ("Pastel1  #FFFFCC", "#FFFFCC", "浅黄"),
    ("Pastel1  #E5D8BD", "#E5D8BD", "浅米"),
    ("Pastel1  #FDDAEC", "#FDDAEC", "淡粉"),
    ("Pastel1  #F2F2F2", "#F2F2F2", "浅灰"),
    ("user     #FFD9A0", "#FFD9A0", "浅杏  <- 用户点名"),
    ("user     #F2E29B", "#F2E29B", "浅黄  (替代浅绿)"),
    ("user     #E5D5E8", "#E5D5E8", "淡薰衣草 (替代浅蓝)"),
    ("user     #E0C9A6", "#E0C9A6", "浅卡其"),
]

w("=" * 118)
w("候选池筛选（只在 ColorBrewer Pastel1 / Paired + 用户自定义色 这一设计家族内挑）")
w("=" * 118)
w("%-18s %-9s %6s %8s %8s %8s %8s  %s"
  % ("来源", "hex", "L*", "对前景", "对白底", "protan", "deutan", "判定/备注"))
w("-" * 118)

ok = []
for src, h, note in POOL:
    Lc = CV.lab(h)[0]
    d_fg = min(CV.dE76(h, f) for _, f in FG)
    d_wh = CV.dE76(h, "#FFFFFF")
    d_pr = min(CV.dE76(h, f, "protan") for _, f in FG)
    d_de = min(CV.dE76(h, f, "deutan") for _, f in FG)
    fails = []
    if d_fg < 15:
        fails.append("C1 对前景 %.2f<15" % d_fg)
    if d_wh < 12:
        fails.append("C3 对白底 %.2f<12" % d_wh)
    tag = "通过" if not fails else "剔除: " + "; ".join(fails)
    w("%-18s %-9s %6.2f %8.2f %8.2f %8.2f %8.2f  %s   %s"
      % (src, h, Lc, d_fg, d_wh, d_pr, d_de, tag, note))
    if not fails:
        ok.append((src, h, note, Lc, d_fg, d_wh, d_pr, d_de))

w("")
w("通过 C1/C3 共 %d 色" % len(ok))
w("")

# --- 选择：优先放入用户点名色，再用池内余量最大的补足，全部受 C4 约束 ---
sel = []
byh = {o[1]: o for o in ok}


def try_add(o, why):
    if not sel:
        sel.append(o)
        w("   + %s  %s（%s）" % (o[1], o[2], why))
        return True
    d = min(CV.dE76(o[1], s[1]) for s in sel)
    if d >= 15:
        sel.append(o)
        w("   + %s  %s（%s）  与已选 min dE = %.2f" % (o[1], o[2], why, d))
        return True
    w("   - %s  %s 跳过（%s）：与已选 min dE 仅 %.2f < 15" % (o[1], o[2], why, d))
    return False


for h in PRIORITY:
    if h in byh:
        try_add(byh[h], "用户点名")
rest = sorted((o for o in ok if o not in sel), key=lambda x: -x[4])
i = 0
while len(sel) < 6 and i < len(rest):
    try_add(rest[i], "池内余量第%d" % (i + 1))
    i += 1

w("")
w("=" * 118)
w("选出 %d 色：" % len(sel))
w("=" * 118)
for i, (src, h, note, Lc, d_fg, d_wh, d_pr, d_de) in enumerate(sel, 1):
    w("  %d. %-9s  L*=%5.2f  对前景 %5.2f  对白底 %5.2f  protan %5.2f  deutan %5.2f   [%s]"
      % (i, h, Lc, d_fg, d_wh, d_pr, d_de, src.strip()))
w("")
w("两两 min dE 矩阵：")
w("        " + "".join("%9s" % s[1] for s in sel))
for a in sel:
    w("%-8s%s" % (a[1], "".join("%9.1f" % (0 if a is b else CV.dE76(a[1], b[1])) for b in sel)))
mn = min(CV.dE76(sel[i][1], sel[j][1]) for i in range(len(sel)) for j in range(i + 1, len(sel)))
w("")
w("★ 彼此 min dE = %.2f  (要求 >= 15)  -> %s" % (mn, "OK" if mn >= 15 else "不达"))
w("★ 对前景 min dE = %.2f ；对白底 min dE = %.2f ；protan min = %.2f ；deutan min = %.2f"
  % (min(s[4] for s in sel), min(s[5] for s in sel),
     min(s[6] for s in sel), min(s[7] for s in sel)))
w("   (protan/deutan 列为**报告项**：红绿色盲下背景浅色会与前景浅绿/浅蓝趋同，"
  "实际区分依赖尺寸 0.85 vs 3.6 pt 与位置分层)")

# 相邻排列：greedy 链最大化相邻色差
rest = list(sel)
chain = [max(rest, key=lambda x: x[3])]     # 起点：最亮
rest.remove(chain[0])
while rest:
    nxt = max(rest, key=lambda x: CV.dE76(chain[-1][1], x[1]))
    chain.append(nxt)
    rest.remove(nxt)
adj = [CV.dE76(chain[i][1], chain[i + 1][1]) for i in range(len(chain) - 1)]
w("")
w("★ 相邻排列（相邻 dE 越大越能看出染色体分块）：")
for c in chain:
    w("     %s   %s" % (c[1], c[2]))
w("     相邻 dE = %s   min = %.2f" % (", ".join("%.1f" % a for a in adj), min(adj)))
w("")
w("★ 交付用 py 字面量（按相邻排列顺序）：")
w("   BG6 = [%s]" % ", ".join('"%s"' % c[1] for c in chain))

with io.open(LOG, "w", encoding="utf-8", newline="\n") as fh:
    fh.write("\n".join(L) + "\n")
print("done -> %s" % LOG)
