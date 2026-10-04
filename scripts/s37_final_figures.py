#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""s37_final_figures.py -- 论文最终图表 Fig1-Fig5 组装

规格：Nature 双栏 183 mm（实际标定到 182.4 mm，落点恒在硬窗口 181.5-183.3 内）；
      Arial；**最小字号 5.2 pt**（Nature 图内文字下限 5 pt；本脚本无 <5.2 pt 的文字）；
      图内文字全部 ASCII（_fig_style.save 做字形核验；并已设 axes.unicode_minus=False
      以避免刻度负号写成 U+2212）。
数据：tables/ 下已封表；不改动任何源表。
"""
import os, sys, csv, math, re
from collections import defaultdict, OrderedDict

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import _fig_style as FS                      # noqa: E402
import numpy as np                            # noqa: E402
import matplotlib                             # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt               # noqa: E402
import matplotlib.gridspec as gridspec        # noqa: E402
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle

# ★图内脚注（panel notes）开关：FIG_NOTES=0 时不画底部说明文字（2026-09-29 用户裁定）
SHOW_NOTES = os.environ.get("FIG_NOTES", "1") != "0"  # noqa: E402

FAM = FS.setup()
BASE = os.path.join(FS.ROOT, "tables")
OUT = []
TARGET_MM = 182.4          # 标定目标（居中，留出硬窗口两侧余量；Nature 双栏上限 183 mm）
LO_MM, HI_MM = 181.5, 183.3    # P0 验收窗口（用户要求 181.5-183.3 mm）
CAL = {}                   # 自动标定得到的 scale，写入 _fig_scale.json 复现
CENTER_RESID = {}          # 2026-10-01：各行面板组居中后的残差(mm)，供核验脚本读取
ROW_GAP_RESID = {}         # 2026-10-01：各相邻行带的实测可见空白(mm)
ROWSPAN_RESID = {}         # 2026-10-01：Fig4 专用 —— 三排面板组的「轴区」左右缘与列间距
ALIGN_RESID = {}           # 2026-10-01（第九轮）：Fig4 专用 —— A/C/E 刻度文字左缘对齐残差(mm)
ROW0_HR_RESID = {}         # 2026-10-01（第九轮）：Fig4 专用 —— 动态 height_ratios 与实际行高(mm)
                           #   （供 _s_int_rowspan 类探针读取，核验「三排等宽 + A/C/E 齐平」）
INK_DX_RESID = {}          # 2026-10-02（第九轮追补）：Fig4 专用 —— 「墨迹左缘」二次微调实测量
ACENTER_RESID = {}         # 2026-10-02（第十一轮）：Fig4 专用 —— (a) 绘图区中心对齐 (c) 的位移实测量
LETTER_ANCHOR_OVR = {}     # 2026-10-02（第十一轮）：Fig4 专用 —— B/D/F 字母改用 corner 锚的实测量
ESHIFT_RESID = {}          # 2026-10-02（第十二轮）：Fig4 专用 —— (e) 绘图区左移实测量（字母 E 不动）
FYLIM_RESID = {}           # 2026-10-02（第十二轮）：Fig4 专用 —— (f) y 轴范围内缩 + 图例下移实测量
# 2026-10-02（Fig5 第一轮细改）：四个留痕字典
A5LEG_RESID = {}           # (a) 图例：右缘是否越出面板右缘 + 下移量
D5LAB_RESID = {}           # (d) y 刻度标签：相邻净空（重叠判据）
LETTER5_ALIGN = {}         # 字母 A/F 对齐到 D 的位移实测量
F5LEG_RESID = {}           # (f) 图例下移实测量
SHIFT5_RESID = {}          # 字母 E 左上位移的实测量（第二轮）
D5BAR_RESID = {}           # (d) 五根柱的配色 + 阈值线色（第二轮）
AXALIGN5_RESID = {}        # (f) 绘图区右移到与 (d) 左缘对齐的位移量（第三轮）
# ★★ 2026-10-01 修正（重要，Fig2 细改第 2 轮发现）：
#   历史版本**只写不读** `_fig_scale.json` ⇒ 任何在新进程里用 `S.CAL.get(name, 1.0)`
#   的下游脚本（典型：`s120_fig2_preview.py` 生成预览裁剪矩形）都会**静默退化成
#   scale=1.0**，于是拿**错误的 axes 位置**去定位/裁剪 PNG（实测 s120 的 crop 矩形
#   就是这样算错的，`_preview_log.txt` 里 scale=1.00000 却对着 scale=0.99507 的 PNG）。
#   现在在 import 期把已落盘的标定读回 CAL；读不到就保持空字典（与历史行为一致）。
#   ★ **对 fit_save 零影响**：它从 s=1.0 起迭代收敛，从不以 CAL 作种子 ⇒ 不改变任何
#     已交付图件的像素。
if True:
    try:
        import json as _json
        _cs = os.path.join(FS.FIGDIR, "_fig_scale.json")
        if os.path.exists(_cs):
            with open(_cs, "r", encoding="utf-8") as _fh:
                CAL.update(_json.load(_fh).get("scales", {}))
    except Exception as _e:                     # 任何异常都不应阻断出图
        print("[warn] _fig_scale.json 读取失败: %r" % (_e,))

# ============================================================================
# ★★ P0（2026-10-01）：幅高压缩参数化。
#    **默认值 = 原硬编码值**，故用默认值时出图与历史逐字节一致。
#    改这些键即可压缩幅高，**不触任何数据**。
#    实测机制：tight 内容高 = (top - bottom) * H_nom + 文字外溢(常数)
#      ⇒ ① hspace 只重分配「行高 / 行间距」，**不改变总高**；
#         ② 真正压缩总高的是 top / bottom / H_nom 三个键。
#      ⇒ 压缩时通常要「同时」降 H_nom 与降 hspace，以保住每个面板的行高。
# ============================================================================
LAY = {
    # ★ 2026-10-01 本轮：h 7.287 -> 6.70（**只压幅高，不动数据/排版比例**）。
    #   原因：字母锚定把左侧墨迹收窄 ⇒ tight 宽 191.9 -> 177.3 mm ⇒ fit_save 把 scale
    #        从 0.95053 顶到 1.02881 ⇒ 幅高被等比放大到 179.8 mm，**突破 P0 的 ≤170 mm**。
    #   标定：`_s161_fig3_height_probe.py` 实测 h 与预测幅高近似线性
    #         （7.287->179.8、7.00->174.4、6.70->168.6 mm）⇒ 取 6.70（余量 1.4 mm）。
    #   ★ hr / hspace 不变 ⇒ 行高按同一比例缩小（面板略扁），宽高比恢复 >1.07。
    # ★ 2026-10-01（逐图细改 · Fig3 · 「D 要跟 AB 对齐」）：
    #   h 6.700 -> **6.550**。原因：字母列对齐后三个字母都落到 A 的标准锚点，
    #   而 A 的锚点在画布坐标 x≈+3.2 mm，**B/D 原本的锚点是 −4.5 / −1.0 mm（在画布外）**
    #   ⇒ 对齐后最左墨迹由"字母"变成 B 行标签 `ANXA4_Mono_NC  (suggestive)`（≈3.47 mm），
    #   内容变窄 ⇒ fit_save 把 scale 从 1.03449 顶到 **1.05452** ⇒ 幅高 168.19 -> **170.56 mm**，
    #   越过 P0 的 ≤170 mm。标定（`_s180_fig3_h_probe.py`，三点实测）：
    #       h 6.700 -> 170.56（OVER） / 6.550 -> **167.72（OK）** / 6.400 -> 164.89
    #   取 6.550（面板高度只降 2.2%，数据与排版比例不变）。
    # ★★ 2026-10-01 21:2x（D 再上移 1.50 -> 3.05 mm）：D 是第 3 行最上墨迹，行带工具把
    #   第 3 行整行下压 1.566 mm 以维持 4.74 mm 行带 ⇒ 幅高 168.66 -> **170.22 mm**，又越限。
    #   ⇒ h 6.550 -> **6.500**（经验斜率 18.9 mm / 单位 h，预计 -0.95 mm）。
    # ★★ 21:3x 定稿：D 再上移 0.625 mm（3.05 -> 3.675，为拿到 ≈1.2 mm 的净空）
    #   ⇒ h 6.500 -> **6.475**（抵消该上移带来的幅高增量，守住 ≤170 mm）。
    # ★★ 2026-10-01 22:0x（D 的 Y 轴对齐 A）引起**幅高来源变更**，h 必须重标定：
    #   `align_axes_left` 把 axd 搬到与 axa 同一列后，D 面板右侧的 `PP.H4 = 0.8` 注记
    #   （`axd.text(1.015, ..., transform=axd.transAxes)`）也随之右移到 **182.25 mm**，
    #   成了全图最右墨迹 ⇒ `fit_save` 在 **scale=1.0** 就命中宽度窗口（182.80 mm），
    #   不再像以前那样放大到 1.05452 ⇒ 整个图变小、幅高掉到 163.45 mm。
    #   ⇒ 改为直接调 h 恢复版面高度。实测（scale=1.0，`_s203_fig3_scale_probe.py`）：
    #       h 6.475 -> 163.45 / 6.700 -> 167.30 / 6.800 -> 168.99 （≈16.9 mm / 单位 h）
    #   ⇒ 取 **6.810**（预计 ≈169.2 mm，留 ~0.8 mm 余量）。
    # ★★ 2026-10-01 22:3x（用户第六轮逐图 · Fig3-D 虚线缩短）引起**第二次幅高来源变更**：
    #   D 的注记原是全图最右墨迹（182.21 mm）；缩短虚线后注记左移到 ~165.8 mm
    #   ⇒ 右界缩到 ~175.2 mm ⇒ tight 宽变窄 ⇒ `fit_save` 把 scale 从 **1.0 顶到 1.04514**
    #   （实测 182.12 × **172.21 mm**，宽 OK 但**幅高超 170**）。
    #   ⇒ 依经验斜率 ≈16.9 mm/单位 h（scale=1.0 口径），本状态下等效 ≈17.7 mm/单位
    #      ⇒ 需降 ~3.4 mm ⇒ h 6.810 -> **6.620**（预计 ≈168.8 mm）。
    "fig3": dict(h=6.620, hr=[1.0, 0.90, 0.52], hspace=0.35,
                 top=0.980, bottom=0.075),
    # ★ 2026-10-01 本轮：h 6.138 -> 5.730（**只压幅高**）。
    #   原因：本图改前两处行带的可见空白只有 ~1.2 / ~0 mm（行1 与行2 的墨迹几乎相接）
    #        ⇒ 撑到全局标准的 4.74 mm 会把两行往下推 ⇒ 幅高 168.27 -> 178.3 mm，
    #        **突破 P0 的 ≤170 mm**。
    #   标定：每减 0.1 的 h 约减 2.52 mm 幅高（= 25.4·(top−bottom)·Δh·s，s=1.08971）
    #        ⇒ 需 −10.3 mm ⇒ Δh ≈ −0.41。
    #   ★★ 实测两轮微调：h=5.730 → 169.97 mm（余量仅 0.03 mm，太险）；
    #      取 **5.640** ⇒ 再降 ≈2.3 mm ⇒ ≈167.7 mm（余量 ≈2.3 mm）。
    #      （行带改用 `get_tightbbox` 口径后比旧的 raster 口径多要 ~0.26 mm/边界，
    #        故 5.730 的旧标定值不够。）
    # ★★★ 2026-10-01（逐图**第九轮** · Fig4 —— **退回「原版」基准重做**）
    #   用户原话：「这一版我不满意，**还是在原版的基础上修改**，图AX轴斜45度摆放，
    #   缩小三排图之间的间距，**稍微放大图A**，**图BX轴太长请缩短**，**F图X轴同样请稍微缩短**，
    #   使**图A、C、E的Y轴图标最左端对齐**，并且**三个字母也要最左端对齐**」
    #   ⇒ 第八轮「(a) 长到与 (e) 等宽 66.90 × (b) 缩到与 (f) 等长 66.90」的方案**整体否决**
    #     （那一版里 (a) 已到 66.90 无法"再稍微放大"、(b) 已到 66.90 谈不上"太长"，
    #      只有**原版**「(a) 42.38 / (b) 86.32」才与本条的两句话自洽）。
    #   经 AskUserQuestion 三项裁定：A 的 X 轴 → **约 46 mm**（+9%，"稍微放大"）；
    #     「Y 轴图标最左端」= **Y 轴刻度文字左缘**（非轴线、非"两者都要"）；
    #     B → **66 mm**、F → **62 mm**（两者**不等长**，用户手工指定）。
    #
    #   ---- 原版实测基线（`_s234`/`_s236`/`_s237` 从 `Fig4_pre_r8.png` 像素反解，mm）----
    #     (a) [25.062, 67.438] 42.376 正方 | colorbar [69.385, 71.671] 2.286
    #     (b) [95.124, 181.443] **86.319**（X 轴线实测长）
    #     (c) [25.019, 114.259] 89.240 | (d) [136.781, 181.443] 44.662
    #     (e) [25.062,  91.991] 66.930 | (f) ≈[114.5, 181.4] 66.9
    #     行 0/1/2 = 42.2 / 40.2 / 37.2 mm；行带 **4.74 mm × 2**（= Fig2 全局标准）
    #     ★ 三排 **y 刻度文字左缘**：A 13.928 / C 20.617 / E 3.768 ⇒ **极差 16.85 mm**
    #       （这就是用户看到的错位）；而三排**轴线**左缘 A/C/E **同为 25.06**
    #       ⇒ 原版错的是**文字列**，轴线本已齐平。
    #
    #   ---- 本轮几何（全部按**绝对 mm** 指定 ⇒ 与 fit_save 的 scale 解耦，见 fig4 内 `_mm2x/_mm2y`）----
    #     · 行 0 行高 **= 46.0 mm = (a) 的正方形边长**（`hr` 由 fig4 内按 scale 动态反解）；
    #     · (a) 46×46、(b) X 轴 **66**、(f) X 轴 **62**；三排**右缘统一**到列 1 右缘 R；
    #     · (a)、(c) 水平左移，使 **A/C/E 的 y 刻度文字左缘共线**（基准 = 三者中最靠左者，
    #       即原版的 (e)）⇒ 三排**轴线**将按各自文字宽度错开（这是用户明确选择的代价）；
    #     · 三个字母 A/C/E 的**左缘共线**（新工具 `FS.align_letter_left_edges`）；
    #     · 行带 4.74 → **3.00 mm**（用户「缩小三排图之间的间距」）。
    #     · (b) 因行 0 由 70.9 → 46 mm 而变矮 ⇒ ylim 与三处注记分层需重排（见 fig4 内注释）。
    #
    #   ★★★ `ink_dx_mm` —— 「墨迹左缘」二次微调（第九轮追补，2026-10-02）：
    #     `FS.align_ytick_label_left()` 对齐的是**文字逻辑包围盒**左缘（三排实测均 8.442 mm，
    #     极差 0.000）。但用户要的是「最左端」= 人眼看到的**墨迹**左缘，二者相差
    #     **各排首字符的左侧边距（side bearing）**：A 排首字符 `r`（bearing 最小）、
    #     C 排 `0`/`1`、E 排 `(`（bearing 最大）。
    #     实测（成品 PNG，`scripts/_s240_fig4_r9_leftedge_verify.py`，logo 口径 = 紧贴轴线的那段文字）：
    #         A = 8.805 mm | C = 9.059 mm | E = 9.821 mm  ⇒ 极差 **1.016 mm**（≈0.77 个字高，肉眼可辨）
    #     ⇒ 令 C、E 两组再左移 (bearing_i − bearing_A)：
    #         c: 9.059 − 8.805 = **−0.254 mm**   e: 9.821 − 8.805 = **−1.016 mm**
    #     ★ 只平移 (c)(e) 两组绘图区（(a) 不动）；(b)(d)(f) 不在组内 ⇒ 三排右缘仍统一。
    #     ★★ 维护：若各排 y 刻度**文字内容/字号**变更 ⇒ 必须重跑 `_s240` 重测这两个值，
    #        `_s240` 会把「三排极差」直接报出来（>0.10 mm 即视为回归）。
    #   ★★★ 2026-10-02（第十轮）右列三图 (b)(d)(f)：
    #     用户原话「图D Y轴与B和F对齐，图D X轴长度与F长度一致」。
    #     ⇒ **列 1 左缘 `_L1 = _Rx − _bw` 成为三图共用的 Y 轴位置**；(b) 保持 66 mm，
    #        (d)、(f) 各 62 mm（`d_mm` / `f_mm`）。
    #     实测：三图 Y 轴 = 115.113 / 115.113 / 115.113 mm（极差 0.000）；
    #            (d) X 轴 62.000 = (f)；右端 181.113 / 177.113 / 177.113（差 4.0 mm）。
    #   ★★★ 2026-10-02（第十一轮）用户原话：
    #        「图D、F的最右端与图B对齐。D、F都在图的左上角。B、D、F对齐。
    #          图A往右移动一点，与图C居中」
    #     ⇒ ① **(d)、(f) 加宽到与 (b) 等宽 66 mm**（`d_mm`/`f_mm` 62 → **66**）
    #          ⇒ 三图**左右缘全部齐平** = `[_L1, _Rx]` = `[115.113, 181.113]`。
    #          ★ 这与第十轮的「右缘差 4 mm」是**用户改主意**（第十轮 AskUserQuestion 曾裁定
    #            「B 保持 66 mm」并接受右缘差 4 mm；第十一轮要求右缘也对齐 ⇒ 只能把 (d)(f)
    #            加宽到 66，而不是把 (b) 缩到 62，也不是把 (d)(f) 右移 —— 右移会破坏
    #            第十轮的「三图 Y 轴共线」）。
    #       ② **B/D/F 三个面板字母左缘共线**：做法 = 三者统一改用
    #          `FS.letter_offset_for(..., anchor="corner")`（锚 = **面板绘图区左上角**，
    #          Δx 仍 = 全局标准 6.41 mm）。因三图左缘已同为 115.113 ⇒ 字母左缘天然全等。
    #          ★ 这是对全局标准「有数值刻度就用 `anchor="word"`」的**显式例外**（用户点名要求），
    #            与 `ink_dx_mm` / `align_ytick_label_left` 同为 Fig4 的登记例外。
    #          ★ 纵向不动（沿用同行左侧字母的 dy）⇒ 保住既有硬要求 A/B、C/D、E/F 同行共线。
    #       ③ **(a) 绘图区中心 ≡ (c) 绘图区中心**（+13.099 mm，见 fig4() 内 `_d_center_mm`）；
    #          色条随行。★ 字母 A/C/E 由 `x_px=` 钉回原公共左缘 ⇒ **字母列不动**。
    #
    #   ★★★ 2026-10-02（第十二轮）用户原话：
    #        「图E的图再往左移动，字母E不动，图往左移动5mm。
    #          F图Y轴减少高度，删掉刻度17.5，最高刻度15即可，
    #          并且右上角图标往下移，要在删除后的Y轴内。」
    #     ⇒ ① **(e) 绘图区整体左移 5.0 mm**（`e_shift_mm = −5.0`），**面板字母 E 不动**：
    #          复用第十一轮「钉住对齐基准」手法 —— `_chosen_ace_px` 在平移前算出，
    #          之后 `align_letter_left_edges(..., x_px=...)` 把 A/C/E 三字母钉回该公共左缘。
    #          ★ 副作用（与第十一轮 A 排同口径）：(e) 的 y 刻度文字随绘图区左移 5 mm
    #            （9.059 → 4.059 mm），故 **E 排刻度文字不再与 C 排共线**。
    #       ② **(f) 的 y 轴上界 17.5 → 15**（`f_ymax = 15.0`；数据最大柱 = 12 ⇒ 安全），
    #          刻度改为 [0, 2.5, …, 15]。
    #       ③ **(f) 右上角图例下移**，使其完全落在新的 0–15 轴内（`f_leg_anchor_y`）：
    #          原 `loc="upper right"`（轴分数顶 0.9762）⇒ 下移到 `f_leg_anchor_y`。
    #
    #   ★★★ 2026-10-02（第十二轮 · 用户澄清「**我想缩减 Y 轴长度**」，修正 ②）：
    #        用户原话：「F图我的意思是，保持原来Y轴到15的距离，直接删掉15-17.5的长度，
    #                  我想缩减Y轴长度。」
    #        ⇒ **不是**只把 ylim 改成 0–15（那会让柱变高、与「保持原来距离」相悖），
    #          而是**把 (f) 面板物理变矮**：0→15 段的**物理长度保持不变**（比例尺不变），
    #          直接**删掉 15–17.5 这一段的高度**。
    #        实现：`axf` 的 y0 不变（**底缘仍与 (e) 对齐**），高度 × `f_ymax / f_ymax_prev`
    #              = 15.0 / 17.5 = **0.857143** ⇒ 38.536 → **33.031 mm**，顶由 46.212 → **40.708 mm**。
    #          ★ 附带的必然结果（已知、须向用户披露）：F 字母锚在本面板**绘图区左上角**，
    #            面板顶下移 ⇒ **F 字母随之下降 5.5 mm，E/F 不再同水平线**（各字母仍在自己
    #            面板的左上角，属标准面板标注）；(f) 的 y 标题（旋转 3 行、按面板居中）
    #            也随之下移约 2.75 mm。
    #          ★ 本行带顶界由 (e) 决定（(e) 未变）⇒ `set_row_ink_gaps` 实测不变、幅面不变。
    "fig4": dict(h=6.18, hr=[1.0, 0.90, 0.85], hspace=0.30,
                 top=0.985, bottom=0.075, row_gap=3.00,
                 row0_mm=46.0, a_mm=46.0, b_mm=66.0, d_mm=66.0, f_mm=66.0,
                 cb_pad_mm=2.00, cb_w_mm=2.40,
                 e_shift_mm=-5.0, f_ymax=15.0, f_ymax_prev=17.5, f_leg_anchor_y=0.93,
                 ink_dx_mm={"c": -0.254, "e": -1.016}),
    # 2026-10-02（Fig5 第一轮细改，用户原话）：
    #   「图A的图标往下和左移动，图标最右端不要超过图的最右端。」
    #     ⇒ (a) 图例改 **右对齐到绘图区右缘**（`loc="upper right"`, `bbox_to_anchor=(1.0, a_leg_anchor_y)`）
    #       ⇒ 右缘恒 = 绘图区右缘（约束由构造保证），整体左移；并缩字号 5.7→5.2 + 收紧间距，
    #         使其不再横越整个面板宽度。`a_leg_anchor_y` 控制下移量。
    #   「图D的Y轴标签重叠了，可适当拉宽Y轴长度，不要使文字重叠。」
    #     ⇒ (d) 的 y 刻度标签由 **3 行折行改 2 行**（把 "(N tests)" 并到域名末行，文字一字不改），
    #       并收紧 `linespacing` 到 1.15；配合本轮 scale 上抬后的轴高增长消除重叠。
    #   「以字母D为标准，字母A、F和D纵向对齐。」
    #     ⇒ 面板字母锚定**之后**再加一步**事后水平平移**：把 A、F 的左缘搬到 D 的左缘。
    #   「图F右上角的图标适当往下移动。」
    #     ⇒ (f) 图例 `bbox_to_anchor=(1.0, f_leg_anchor_y)`。
    # ★★★ 2026-10-02（Fig5 第二轮）用户：「图A最上端的图标往下移动3mm。」
    #   换算（`_s280` 实测，scale 0.9842）：(a) 轴高 = 62.046 mm、`borderaxespad` 0.5×5.2pt = 0.917 mm
    #   ⇒ 图例顶距轴顶 = 0.105×62.046 − 0.917 = 5.598 mm（与实测一致）。
    #   下移 3 mm ⇒ Δanchor = −3.0 / 62.046 = −0.0483476 ⇒ 1.105 → 1.05665。
    #   ★ 可行性已实测：图例底现距轴顶 +0.758 mm、最上柱顶距轴顶 2.820 mm
    #     ⇒ 下移 3 mm 后图例底（轴顶 −2.242 mm）仍比最上柱顶（−2.820 mm）高 **0.578 mm**，不压柱。
    "fig5": dict(h=7.297, outer_hr=[1.0, 1.10], outer_hspace=0.45,
                 gsr_hr=[2.0, 1.0], gsr_hspace=0.70, bot_hspace=0.45,
                 top=0.960, bottom=0.085, left=0.170, right=0.995,
                 a_leg_anchor_y=1.05665, f_leg_anchor_y=0.900,
                 d_lab_linespacing=1.15, letter_align_to="d",
                 letter_align_members=("a", "f"),
                 # ★★★ 2026-10-02（Fig5 第二轮）用户：「字母E往左上角移动5mm。」
                 #   口径（AskUserQuestion 已确认）= 水平、垂直**各** 5 mm（负 dx = 左，正 dy = 上）。
                 letter_shift_mm={"e": (-5.0, 5.0)}),
}


def measure_mm(fig, dpi=600):
    """用内存 PNG 实测 bbox_inches='tight' 后的物理尺寸。

    ★ 必须与 FS.save 的 dpi 一致（600）：150 dpi 下测得的宽度与 600 dpi 落盘宽度
      曾出现最大 ~2 mm 偏差（Fig5 实测 181.0 vs 标定值），导致验收窗口误判。
    """
    import io
    from PIL import Image
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=dpi, bbox_inches="tight", pad_inches=0.02,
                facecolor="white")
    buf.seek(0)
    im = Image.open(buf)
    return im.size[0] / dpi * 25.4, im.size[1] / dpi * 25.4


def fit_save(build, name, target=TARGET_MM, tol=0.45, iters=10):
    """把整幅图的 tight 宽度自动标定到 target mm 附近，再交给 FS.save 做字形核验。

    P0 修正：
      ① 内部停止判据 = |w - target| <= tol（target=182.4, tol=0.45 -> 181.95-182.85），
         并要求同时落在硬验收窗口 [LO_MM, HI_MM] 内 —— 只判「在窗口内」会让落点偶然
         贴在 183.3 上界（Fig5 实测），再导出一次就可能越界；居中落点更稳。
      ② 落盘后再按实测宽度闭环修正，最多 3 轮。
    """
    s = 1.0
    w_mm = float("nan")

    def _good(x):
        return LO_MM <= x <= HI_MM and abs(x - target) <= tol

    for _ in range(iters):
        fig = build(s)
        w_mm, _h = measure_mm(fig, dpi=600)
        plt.close(fig)
        if _good(w_mm):
            break
        s *= target / w_mm
    paths, wmm, hmm = FS.save(build(s), name)
    tries = 0
    while not _good(wmm) and tries < 3:
        tries += 1
        s *= target / wmm
        paths, wmm, hmm = FS.save(build(s), name)
        print("   [re-fit %d] %s -> %.2f mm" % (tries, name, wmm))
    CAL[name] = round(s, 5)
    # ★ 2026-10-01：**每次标定即落盘**。否则 _fig_scale.json 会停留在上一次整轮运行的
    #   旧值，下游脚本（预览/核验）读回的 scale 与磁盘上的 PNG 不匹配（见上方 CAL 处注释）。
    #   因为 import 期已把旧 JSON 读进 CAL，这里写回不会丢掉其它图的记录。
    try:
        import json as _json
        with open(os.path.join(FS.FIGDIR, "_fig_scale.json"), "w", encoding="utf-8") as _fh:
            _json.dump({"target_mm": TARGET_MM, "scales": CAL,
                        "p0_height_cap_mm": 170.0, "p1_panel_label_pt": 8},
                       _fh, indent=1, ensure_ascii=False)
    except Exception as _e:
        print("   [warn] _fig_scale.json 写入失败: %r" % (_e,))
    if not _good(wmm):
        print("   [!! WIDTH OUT OF RANGE] %s = %.2f mm (want %.1f-%.1f, aim %.1f)"
              % (name, wmm, LO_MM, HI_MM, target))
    return paths, wmm, hmm


def load(f):
    with open(os.path.join(BASE, f), encoding="utf-8-sig", newline="") as fh:
        return list(csv.DictReader(fh))


def num(x, d=float("nan")):
    try:
        return float(x)
    except Exception:
        return d


def txt_on(cmap, norm, v, dark="#222222", light="white", cut=0.45):
    """按单元格真实亮度选字色（深底白字 / 浅底深字），避免固定阈值下的低对比。"""
    r, g, b, _a = cmap(norm(v))
    lum = 0.2126 * r + 0.7152 * g + 0.0722 * b
    return light if lum < cut else dark


GENE_C = {"CDC42": FS.C["locus1"], "LINC00339": FS.C["locus2"],
          "WNT4": FS.C["locus3"], "YME1L1": FS.C["gene_ym"], "ANXA4": FS.C["gene_an"]}
# ★ 2026-10-01（逐图第七轮 · Fig4 要求 2）：Fig4(c) 的「落在三个基因体之外」的点色。
#   取 GB12 中**全项目未被占用**的一档浅钢蓝 #96C2D4（GB12 #10），
#   既是真彩色（不再是中性灰），又足够内敛、不与 WNT4 的 #367DB0 抢读。
#   ★ 不用 #6CBAD8（本图 (f) 已用作 rs10917151 的杆色，语义会撞车）。
INTERGENIC_C = "#96C2D4"

# ---------------------------------------------------------------- Fig2(a) 背景点色盘
# 2026-10-01（图件细改 · Fig2 · 问题2）：Fig2(a) 的**非显著变异背景点**由原来的
#   两档灰 (#B0B0B0 / #D5D5D5) 改为「浅色系彩色」，按染色体序号 i % 6 轮换。
#
# ★ 保留用户点名的 4 色：#FDBF6F 浅橙 / #CAB2D6 浅紫 / #FFD9A0 浅杏 / #FBB4AE 浅粉。
# ★ 剔除用户原列表中的 #A6CEE3(浅蓝) 与 #B2DF8A(浅绿)：实测与**前景阳性点同族**
#     —— #A6CEE3 vs ANXA4 #6CBAD8  dE = 14.01
#     —— #B2DF8A vs LINC00339 #8BCF8B  dE = 13.19
#     两者都 < 15，会把背景点与阳性点混在一起（这正是本次要解决的问题，不能用同样近的色）。
# ★ 改用 #FFFF99(浅黄) 与 #FDDAEC(淡粉) 补足 6 色，仍取自 ColorBrewer Pastel1 / Paired
#   的柔和浅色家族，保证观感统一、不出现荧光感。
#
# ★★ 2026-10-01 追加修订（用户反馈「黄色看不清」）：
#    原 #FFFF99 的 CIE L* = 98.07，几乎贴白底；背景点只有 ~1 pt，抗锯齿会把这种
#    高亮低反差的小点冲淡成"一片白黄"（ΔE 高 ≠ 视觉可辨）。
#    —— 换成 **#E4D467 浅金黄 L* = 84.26**：明度降 13.8，在浅色系里仍属"浅"，
#       但对比明显变实。审计：对前景 min dE 39.33 / 对白底 57.85 / 与其余 5 色 min dE 22.81，全部通过。
#       候选依据 scripts/s122_bg6_yellow_swap.py（暖黄区 40-72 度 x L*78-87 x S0.4-0.7 共 160 个候选）。
#    —— 备选（若仍不满意，改这一个字面量即可）：更淡 #E5DE6E / 偏橙 #E7D27A / 更亮 #D9D926。
#    —— 不可选：浅蓝/浅绿/青系（与 ANXA4、LINC00339、WNT4 冲突，dE 全部 < 15）。
#
# 客观审计（scripts/s117_fig2_bgcolor_audit.py、scripts/s117b_select_bg6.py，
#           色度学复用 cvd-safe-palette-design/scripts/cvd_sim.py）：
#   · 对 6 个前景语义色的 min CIE76 dE >= 22.81（门槛 15）
#   · 对白底 min dE >= 18.38（保证在白纸上看得见，不是"白上白"）
#   · 6 色彼此 min dE = 18.20（染色体分块可辨）
#   ★ protan 下 min dE = 4.85 是**已知残余风险**（红绿色盲下浅紫与浅绿趋同），
#     实际区分依赖尺寸（背景约 1.0 pt vs 前景 4.6 pt）与"背景在下、阳性在上"的位置分层。
BG6 = ["#E4D467", "#CAB2D6", "#FDBF6F", "#FDDAEC", "#FFD9A0", "#FBB4AE"]

# ★★★ 2026-10-01（用户）：「图 ABC 中**所有点放大一倍**，并且点加上与之**同色系颜色的外轮廓**，
#     外轮廓颜色要比点的颜色**加深 2 个度**」。
#   三项口径（用户逐条选定）：
#     · **范围** = 全部点，**含 (a) 的 8612 个背景点**；
#     · **加深** = 1 个度 = CIE L* 降 10 ⇒ 共 **−20**（`FS.darken()`：保色相 H、保 HSL 饱和度 S）；
#     · **放大** = **面积 ×2 ⇒ 直径 ×√2 = 1.41421**（即 ms 直接乘 1.41421）。
#   ★ 描边宽度：原本没有 mew 的（背景点、显著点、(b) 圆点/菱形、(c) 发现圆点）按
#     「直径的 10%，clamp 到 [0.30, 0.70] pt」统一给（落在 Nature 0.25–1 pt 窗口内）；
#     原本**已经是"纯描边型"**的标记 —— (c) 的白底复核方块（mew 0.7）与 [NR] 叉（mew 1.0）——
#     **保留原 mew**，因为那已经是它们唯一的线重，再动会改变它们的可辨性。
#   ★ 图例：`(a)` 的 `markerscale` 由 1.6 改为 **1.6/1.41421 = 1.1314**，使**图例标记的绝对
#     尺寸保持不变**（用户要放大的是图上的点，不是图例）。
#   ★ 实测代价（`scripts/_s139_fig2_point_outline_plan.py`）：背景点 ms 1.05→1.485 后，
#     8612 个点在 (a) 绘图区（145.06×47.88 mm）上的**面积覆盖率 13.4% → 26.9%**
#     （若取"直径 ×2"会是 53.4%）。这是用户已知并明确选择的结果。
P_UP = 1.41421                     # ms 放大系数（面积 ×2）
P_UP_BG2 = 1.41421                 # 背景 6 色点**第二轮**再放大（用户：面积「再」×2 ⇒ 直径再 ×√2）
P_DL = 20.0                        # 轮廓加深总量（L* 降 20 = 2 个度）
#   ★ 背景点合计 = 1.05 × 1.41421 × 1.41421 = **2.1002 pt**（相对最初 = 面积 ×4）。
#     实测代价（`_s139`）：8612 个点在 (a) 绘图区上的面积覆盖率 **13.4% → 26.9% → 53.4%**。
#     用户已知此数并明确要求"再放大 2 倍"，故照办。回退：`P_UP_BG2 → 1.0` 即回到上一版。


def _mew_for(ms):
    """描边宽度 = 数值直径的 10%，clamp 到 [0.30, 0.70] pt。"""
    return min(0.70, max(0.30, 0.10 * ms))


GENE_CD = {k: FS.darken(v, P_DL) for k, v in GENE_C.items()}    # 基因色 -> 深 2 度
BG6_D = [FS.darken(c, P_DL) for c in BG6]                       # 背景 6 色 -> 深 2 度
GENE_YM_D = FS.darken(FS.C["gene_ym"], P_DL)                    # (b) IVW 菱形
NR_D = FS.darken(FS.C["nr"], P_DL)                              # (c) [NR] 叉
NEUTRAL_D = FS.darken("#333333", P_DL)                          # (c) 图例代理

# Fig2 的 (b)/(c) 列间距（GridSpec wspace = 相对轴宽的比例，**不是 mm**）。
# 2026-10-01 细改前为 0.54（绝对空白带 30.1 mm < (c) 最宽 y 标签 37.56 mm ⇒ 标签压进 (b) 8.66 mm）；
# 现取 0.90（绝对空白带约 43.8 mm ⇒ 可辨间隙约 5 mm）。提为模块常量便于后续微调与复现。
FIG2_WSPACE = 0.90

# Fig2 (a) 面板的**整体水平位移**（figure 分数的横坐标偏移，负值 = 左移）。
# 2026-10-01 用户要求「(a) 整体往左移、跟下面两组图居中」。
#   量因（scripts/s118_fig2_revise.py 实测）：(a) 的绘图区 [22.75, 163.78] mm 与
#   (b) 左缘 / (c) 右缘本就**等宽且左右对齐**；它看着偏右，是因为
#   左侧 36.0 mm 全被 (b) 的长 y 标签占掉（最长 'LINC00339 | CD4_NC (IVW, r2<0.01)' = 34.3 mm），
#   右侧只剩 5.5 mm ⇒ 整块内容在纸面上「左重右轻」。
#   解法：把 (a) 左移 17.75 mm（17.75 / (7.2*0.99507*25.4) = 0.0975）,
#   使 (a) 在成图里左右留白各约 18 mm（真正居中）。
#   ★ 代价：(a) 的左缘不再与 (b) 左缘对齐，上下错开 17.75 mm —— 这是「纸面居中」与
#     「上下严格对齐」之间的取舍，用户明确选了前者。★ 改回 0.0 即恢复原来的严格对齐。
#   ★ 必须在 subplots_adjust 之后设置，否则会被 gridspec 重新布局覆盖。
FIG2_A_SHIFT = 0.0975

# ★★★ 2026-10-01（用户）：「调整图与图之间的距离，**下移图A**，使图A与下面两张图的距离是目前的 1/3」。
#   量因（`scripts/_s136_fig2_row_gap.py`，PNG 像素口径逐行扫"最长全白横带"）：
#     改前 A 行最低墨迹 y = 58.80 mm（= A 的 x 轴标题 'Chromosome' 底），
#          B/C 行最高墨迹 y = 73.07 mm（= C 字母上缘）  ⇒ **可见空白带 = 14.22 mm**。
#     目标 = 14.22 / 3 = **4.74 mm**。
#   ★★ 实现方式（**关键取舍，务必先读**）：
#     可见空白带 = 轴间距 − 常数。常数 = A 轴下溢出 6.50 + B/C 轴上溢出 2.99 = **9.49 mm**（与文字无关地固定）。
#     ⇒ 要让可见空白带 = 4.74，轴间距须 = 4.74 + 9.49 = **14.24 mm**（改前 23.71）。
#     ★ **不能只改 hspace**：图高 H 固定时，hspace 只把"行高 ↔ 行间距"重新分配 ——
#       降 hspace 会让**两个面板都被拉高**，面板内布局（刻度间距）与字母的 dy 全部失准。
#     ⇒ 本方案**同时降 hspace 与降图高**，令**两个行高(mm)一字不动**、只压掉行间距：
#         u = 0.72·H / (Σratio + hspace·avg_ratio)，Σratio = 0.92+1.25 = 2.17，avg_ratio = 1.085
#         改前 H = 7.3·scale = 7.3×1.02351 in = 189.78 mm ⇒ u = 52.04
#           ⇒ 行高 0.92u = 47.88 / 1.25u = 65.05、行间距 0.42×1.085u = 23.71 ✓（与实测逐项吻合，模型可信）
#         要 u 不变 ⇒ 0.72·H = 47.88 + 14.24 + 65.05 = 127.17 mm ⇒ H = 176.62 mm
#           ⇒ 新图高比 = 7.3 × 176.62 / 189.78 = **6.7938**
#         hspace：gap = hspace × 1.085 × u = 14.24 ⇒ hspace = 14.24 / (1.085×52.04) = **0.2521**
#   ★ 收益：**面板高度与横向布局完全不动** ⇒ 三个字母的 dx/dy（第六轮定案值）**继续精确有效**，
#     即"字母 ↔ 图 的距离"这条新标准**不会被本次改动扰动**。
#     幅高由 149.18 → 约 139.7 mm；幅宽不变 182.24 mm（fit_save 只标定幅宽 ⇒ scale 仍 1.02351）。
#     视觉上：A 整块在 figure 坐标里**上下缘同步下移 11.58 mm**，即用户所说的「下移图A」。
#   ★ 回退：`FIG2_H_RATIO → 7.3`、`FIG2_HSPACE → 0.42` 即完全恢复。
FIG2_H_RATIO = 6.7938
FIG2_HSPACE = 0.2521

# Fig2 三个面板字母的**左上偏移**（`panel_label` 的 dx / dy，单位 = axes 分数）。
# ★★★ 2026-10-01 第六轮定案（用户：「三个字母跟图的距离再调近一点，比现在的距离再调近 1/3」）。
#   **口径沿用第五轮，只把目标位移整体乘 2/3**：Δx 9.62 → **6.4133 mm**；Δy 3.63 → **2.42 mm**。
#   Δx = 字母右缘→该面板 Y 轴最上面那个词的左缘；Δy = 字母下缘→该词的上缘。**三块面板统一（含 C）**。
#
#   ★★ 两个必须知道的坑（第五轮实测发现，仍适用）：
#   ① **A 的 y 刻度里有一个 `'14'`（还有 `'-2'`）落在 ylim 之外 ⇒ 不被绘制、也不进 tight bbox**
#      ⇒ 求"最上面那个词"**必须过滤掉落在绘图区之外的刻度标签**，否则会选到看不见的 '14'。
#   ② (b)/(c) 的 y 标签**右对齐**，三块面板的"最上面那个词"宽度不同 ⇒ **dx 必然各不相同**。
#
#   ★ 基准（第五轮定案后的实测值，`scripts/_s135_letter_word_distance.py`，画布 mm）：
#       面板 A：词 '12'                           左缘  1.20  上缘 165.86 ；字母右缘 −8.42  下缘 169.49
#       面板 B：词 'CDC42 | B_IN'                 左缘  9.12  上缘  92.58 ；字母右缘 −0.50  下缘  96.21
#       面板 C：词 'CDC42 | B_IN -> B (m)  [NR]'  左缘 91.30  上缘  93.22 ；字母右缘 81.68  下缘  96.85
#       三者 Δx / Δy = 9.62~9.63 / 3.62~3.63 mm（已统一）⇒ 本轮目标 = 各 × 2/3。
#   ★ 本轮定值（字母宽 2.03 mm 不变；Δx_t = 6.4133，Δy_t = 2.42）：
#       A：目标左缘 =  1.20 − 6.4133 − 2.03 =  −7.2433 → dx = (−7.2433 −   5.15)/145.06 = −0.08544
#          目标下缘 = 165.86 + 2.42         = 168.28    → dy = 1 + (168.28 − 167.01)/ 47.88 = 1.02652
#       B：目标左缘 =  9.12 − 6.4133 − 2.03 =   0.6767 → dx = ( 0.6767 −  23.40)/ 50.02 = −0.45428
#          目标下缘 =  92.58 + 2.42         =  95.00    → dy = 1 + ( 95.00 −  95.41)/ 65.05 = 0.99370
#       C：目标左缘 = 91.30 − 6.4133 − 2.03 =  82.8567 → dx = (82.8567 − 118.44)/ 50.02 = −0.71138
#          目标下缘 =  93.22 + 2.42         =  95.64    → dy = 1 + ( 95.64 −  95.41)/ 65.05 = 1.00354
#   ★ 因为 dx/dy 是**轴分数**，`fit_save` 之后整体缩放时三个面板的间隙**同步缩放 ⇒ 恒保持相等**。
#   ★ **幅宽不变**：A 字母左缘 −7.24 mm 仍落在 (b) 标签左缘 −12.10 mm **之内**，
#     左右界都由别的东西决定 ⇒ `fit_save` 的 scale 仍 1.02351、axw/axh 不变
#     ⇒ 上面按 mm 反算出的 dx/dy 会**精准落地**（无需迭代）。
#   ★ **幅高会下降**：A 字母下移 1.21 mm，而它仍是全图最高墨迹 ⇒ tight bbox 上界下移，
#     实测 150.41 → ≈149.2 mm（正常，非 bug）。
#   ★ 口径史（全部可回退）：距绘图区 8 mm（−0.1595）→ 齐绘图区左缘（0.0）→
#     齐顶行左缘（−0.2853/−0.5426）→ 齐面板 y 文字块最左缘（−0.7096/−0.7755）→
#     以 C 为基准统一 (Δx,Δy)（−0.10754/−0.51839/−0.7755）→ **本轮：整体移近 1/3**。
FIG2_A_LABEL_DX = -0.08544
FIG2_A_LABEL_DY = 1.02652
FIG2_B_LABEL_DX = -0.45428
FIG2_B_LABEL_DY = 0.99370
FIG2_C_LABEL_DX = -0.71138
FIG2_C_LABEL_DY = 1.00354

# P0 审稿意见修正：coloc.abf strong_shared 二级命名（与正文 §3.4 / §D2 完全一致）
ROBUST = ["CDC42_B_MEM", "CDC42_Mono_NC"]        # 45 组敏感性 9/9 均 PP.H4 > 0.8
ROBUST_SET = set(ROBUST)
TIER_C = {"robust_strong": FS.C["robust_strong"],
          "suggest_strong": FS.C["suggest_strong"],
          "moderate_shared": FS.C["moderate_sh"],
          "no_shared": FS.C["no_shared"]}


def coloc_tier(r):
    """coloc.abf 判定 -> 二级命名：robust_strong / suggest_strong / moderate_shared / no_shared。"""
    v = r["abf_verdict"]
    if v == "strong_shared":
        return "robust_strong" if ("%s_%s" % (r["gene"], r["cell_type"])) in ROBUST_SET \
            else "suggest_strong"
    return v


# ============================================================ Fig1
def fig1(scale=1.0):
    fig = plt.figure(figsize=(7.2 * scale, 4.5 * scale))
    ax = fig.add_axes([0.0, 0.0, 1.0, 1.0])
    ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.axis("off")

    def box(x, y, w, h, title, body, fc="#F4F4F4", ec="#4D4D4D", ts=7.0, bs=6.3):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.5,rounding_size=1.0",
                                    fc=fc, ec=ec, lw=0.6))
        ax.text(x + w / 2, y + h - 1.4, title, ha="center", va="top",
                fontsize=ts, fontweight="bold")
        ax.text(x + w / 2, y + h - 5.6, body, ha="center", va="top",
                fontsize=bs, linespacing=1.45)

    def arr(x1, y1, x2, y2):
        ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>",
                                     mutation_scale=6.5, lw=0.6, color="#4D4D4D",
                                     shrinkA=0, shrinkB=0))

    box(2, 84, 44, 13, "Exposure: OneK1K sc-eQTL",
        "14 immune cell types; 980 donors\ncis-eQTL re-run (TensorQTL), GRCh37", fc=FS.C["tint_cool"])
    box(54, 84, 44, 13, "Replication: 1M-scBloodNL",
        "unstimulated; 9 MERGED cell types; N = 119\n"
        "direction-consistent nominal replication only", fc=FS.C["tint_cool"])
    box(14, 67, 72, 11.5, "Instrument selection (cis-eQTL)",
        "P < 5e-8 ; PLINK clumping r2 < 0.001 within 10 Mb\n"
        "single instrument per gene-cell pair (26/26 pairs)", fc="#F7F7F7")
    box(14, 50, 72, 11.5, "Outcome: female infertility GWAS",
        "GCST90483463 ; 40,024 cases / 665,658 controls (EUR)\n"
        "sensitivity outcome: GCST90483469", fc="#F7F7F7")
    methods = [
        (1.5, "MR", "Wald ratio\n(primary)\nIVW / Egger\nnot runnable"),
        (17.8, "Colocalization", "coloc.abf\n(primary)\nsusie\ndiagnostic"),
        (34.1, "Conditional", "condition on\nthe region lead\ncoloc + PheWAS"),
        (50.4, "Fine-mapping", "susie_rss\n(external LD)\ncredible-set\nlead placement"),
        (66.7, "SMR-HEIDI", "cis SMR + HEIDI\nper gene and\ncell type"),
        (83.0, "Safety", "Open Targets\nPheWAS safety\npower analysis"),
    ]
    for x, t, b in methods:
        box(x, 25, 14.5, 18.5, t, b, fc="#FFFFFF", ts=6.4, bs=6.0)
    box(1.5, 4, 97, 16.5, "Safety scan of the chr1p36.12 region: FinnGen R12 PheWAS",
        "2,469 endpoints x 2 instruments = 4,938 tests ; FDR < 0.05: 52 (12 risk-increasing, 40 risk-decreasing)\n"
        "immune / blood / tumour domain: no risk-increasing signal within available power (968 tests; median min. detectable OR 3.15)\n"
        "conditioning on the region lead removes all 52 signals -> the whole profile is region-driven\n"
        "both instruments in high LD with the WNT4 credible-set lead -> signals attributed to the REGION",
        fc=FS.C["tint_warm"])

    for x in (24, 76):
        arr(x, 84, x, 78.7)
    arr(50, 67, 50, 61.7)
    arr(50, 50, 50, 43.8)
    for x in (8.75, 25.05, 41.35, 57.65, 73.95, 90.25):
        arr(x, 50, x, 43.8)
    arr(49.5, 25, 49.5, 20.8)
    return fig


# ============================================================ Fig2
def fig2(scale=1.0):
    disc = load("10_discovery_MR_main.csv")
    allf = load("31_final_MR_with_method.csv")
    fin = [r for r in allf if r["thr"].startswith("P<5e-8")]
    ivw = [r for r in allf if not r["thr"].startswith("P<5e-8")][:1]
    rep = load("24_replication_MR_UT.csv")

    fig = plt.figure(figsize=(7.2 * scale, FIG2_H_RATIO * scale))
    # 2026-10-01（图件细改 · Fig2 · 问题1）：wspace 0.54 -> **0.90**（提为模块常量 FIG2_WSPACE）。
    #   诊断（scripts/s116_fig2_diagnose.py）实测：(c) 面板的 y 刻度标签**左端已越出列间空白带**，
    #   压进 (b) 面板绘图区 8.66 mm —— 不是"看起来挤"，是真遮挡。
    #   根因不是"wspace 设小了"（0.54 已是轴宽的 54%，用户以为的 0.15-0.20 反而更小），
    #   而是 (c) 的最宽标签达 37.56 mm（'LINC00339 | CD4_NC -> CD4T (m)  [NR]'），
    #   而列间绝对空白带只有 0.54/(2+0.54) x 141.7 = 30.1 mm。
    #   ★ wspace 是「相对轴宽」的比例而不是 mm ⇒ 要把绝对空白带做到 >= 40 mm，
    #     解 141.7 s/(2+s) >= 40 得 s >= 0.80；取 **0.90** 留出约 5 mm 的可辨间隙
    #     （实测：间隙 = 空白带 - 刻度 pad 1.23 - 标签宽 37.56）。
    #   ★ 代价：单个面板宽由约 57 mm 降到约 48.6 mm（-15%），幅宽仍由 fit_save 标定 182.4 mm。
    gs = gridspec.GridSpec(2, 2, height_ratios=[0.92, 1.25], hspace=FIG2_HSPACE,
                           wspace=FIG2_WSPACE)
    axm = fig.add_subplot(gs[0, :]); axf = fig.add_subplot(gs[1, 0]); axr = fig.add_subplot(gs[1, 1])

    # --- (a) Manhattan (gene TSS positions, GRCh38) ---
    byc = defaultdict(list)
    for r in disc:
        byc[int(r["chr"])].append(r)
    order = sorted(byc); off, cur, GAP = {}, 0.0, 6e7
    for c in order:
        off[c] = cur
        cur += max(int(x["tss_pos_grch38"]) for x in byc[c]) + GAP
    tick, tlab = [], []
    for i, c in enumerate(order):
        xs = [off[c] + int(x["tss_pos_grch38"]) for x in byc[c]]
        ys = [-math.log10(max(num(x["p"], 1.0), 1e-15)) for x in byc[c]]
        # 2026-10-01（问题2）：背景点改为**按染色体轮换的浅色系**（原来只有两档灰）。
        #   ms 0.85 -> 1.05 pt：浅色点对比度低于灰点，略放大以保证在 600 dpi 下不至于淡到看不见；
        #   8612 个点的密度下，再大就会糊成色块。
        _ms_bg = 1.05 * P_UP * P_UP_BG2
        axm.plot(xs, ys, ls="none", marker="o", ms=_ms_bg,
                 mew=_mew_for(_ms_bg), mec=BG6_D[i % len(BG6_D)],
                 color=BG6[i % len(BG6)], zorder=1)
        tick.append((min(xs) + max(xs)) / 2); tlab.append(str(c))
    for r in disc:
        if r["sig_P_FDR"] == "True":
            x = off[int(r["chr"])] + int(r["tss_pos_grch38"])
            # 2026-10-01（问题2）：阳性点 ms 3.6 -> 4.6 pt，并加 0.4 pt 白色细描边，
            #   使其在浅色背景点上被"抠"出来（Nature：描边须落在 0.25-1 pt）。
            axm.plot([x], [-math.log10(num(r["p"]))], ls="none", marker="o", ms=4.6 * P_UP,
                     mew=_mew_for(4.6 * P_UP), mec=GENE_CD.get(r["gene"], "#000000"),
                     color=GENE_C.get(r["gene"], "#000000"),
                     zorder=3, label=r["gene"])
    thr = -math.log10(0.05 / len(disc))
    axm.axhline(thr, color=FS.C["region"], ls="--", lw=0.6, zorder=2)
    axm.text(axm.get_xlim()[1], thr + 0.12, "Bonferroni  p = 0.05/%d" % len(disc),
             ha="right", va="bottom", fontsize=6.2, color=FS.C["region"])
    h, l = axm.get_legend_handles_labels()
    seen, hh, ll = set(), [], []
    for a, b in zip(h, l):
        if b not in seen:
            seen.add(b); hh.append(a); ll.append(b)
    # ★ 2026-10-01（用户）：「右上角图标的点重叠了要上下分开一点」。
    #   实测量化（scripts/_s142_fig2_legend_probe.py）：图例标记 7.36 pt = **2.60 mm 高**，
    #   而图例行距只有 **2.81 mm** ⇒ 两个圆形标记的 bbox 间隙仅 **0.21 mm** ——
    #   数值上"不算相交"，但深色描边已彼此贴住，肉眼就是**粘成一团**。
    #   修法：加大 `labelspacing`（单位 = 字号）0.5 -> **0.95**，行距约 2.81 -> 约 3.68 mm
    #   ⇒ 标记间隙约 **1.08 mm**。★ 不动 `markerscale`：用户要的是"分开"，不是"变小"。
    axm.legend(hh, ll, loc="upper right", ncol=2, fontsize=6.2,
               markerscale=1.6 / P_UP,      # 图例标记绝对尺寸保持不变
               labelspacing=0.95)
    axm.set_xticks(tick); axm.set_xticklabels(tlab, fontsize=5.4)
    axm.set_xlabel("Chromosome"); axm.set_ylabel("-log10(p)  (MR, single-SNP Wald)")
    axm.set_xlim(-2e6, cur); FS.finalize(axm)
    FS.panel_label(axm, "a", dx=FIG2_A_LABEL_DX, dy=FIG2_A_LABEL_DY)

    # --- (b) Forest: 15 Discovery pairs ---
    ordg = {"CDC42": 0, "LINC00339": 1, "YME1L1": 2, "ANXA4": 3}
    fin = sorted(fin, key=lambda r: (r["locus_id"], ordg.get(r["gene"], 9), r["cell_type"]))
    ylab, yp = [], []
    yy = 0
    for r in fin:
        b, se = num(r["b"]), num(r["se"])
        axf.errorbar([b], [yy], xerr=[[1.96 * se], [1.96 * se]], fmt="o", ms=3.0 * P_UP,
                     mew=_mew_for(3.0 * P_UP), mec=GENE_CD.get(r["gene"], "#333333"),
                     color=GENE_C.get(r["gene"], "#333333"), elinewidth=0.6, capsize=1.4)
        ylab.append("%s | %s" % (r["gene"], r["cell_type"])); yp.append(yy); yy -= 1
    for r in ivw:
        b, se = num(r["b"]), num(r["se"])
        axf.errorbar([b], [yy], xerr=[[1.96 * se], [1.96 * se]], fmt="D", ms=3.0 * P_UP,
                     mew=_mew_for(3.0 * P_UP), mec=GENE_YM_D,
                     color=FS.C["gene_ym"], elinewidth=0.6, capsize=1.4)
        ylab.append("%s | %s (IVW, r2<0.01)" % (r["gene"], r["cell_type"])); yp.append(yy); yy -= 1
    axf.axvline(0, color="#7F7F7F", lw=0.5, ls="-")
    axf.set_yticks(yp); axf.set_yticklabels(ylab, fontsize=5.9)
    axf.set_xlabel("MR effect (beta, 95% CI)"); axf.set_ylim(yy + 0.6, 0.9)
    for s in ("left",):
        pass
    axf.spines["left"].set_visible(False); axf.tick_params(axis="y", length=0)
    FS.finalize(axf); FS.panel_label(axf, "b", dx=FIG2_B_LABEL_DX, dy=FIG2_B_LABEL_DY)

    # --- (c) Replication (P0: cell-type mapping annotated; not-reportable marked) ---
    MERGED = {"B", "monocyte", "CD4T", "CD8T"}   # 1M-scBloodNL 注释为合并层级
    rep = sorted(rep, key=lambda r: (ordg.get(r["gene"], 9), r["cell_type_discovery"]))
    # P0：x 轴范围只覆盖「可报告」的数据（[-0.168, 0.382]）；NR 行的伪影值(-56.2)不入图，
    # 因此无需再用 -0.9 的极宽范围压缩分辨率，也不需要在面板内塞入长文字。
    axr.set_xlim(-0.28, 0.62)
    ylab2, yp2, yy = [], [], 0
    for r in rep:
        db, ok = num(r["discovery_b"]), (r["effect_reportable"] == "TRUE")
        col = GENE_C.get(r["gene"], "#333333")
        cold = GENE_CD.get(r["gene"], NEUTRAL_D)     # 同色系深 2 度（描边）
        rc = r["cell_type_replication"]
        if ok:
            rb = num(r["replication_b"])
            axr.plot([db, rb], [yy, yy], lw=0.5, color="#9E9E9E", zorder=1)
            axr.plot([db], [yy], marker="o", ms=3.0 * P_UP, mew=_mew_for(3.0 * P_UP),
                     color=col, mec=cold, zorder=3)
            # 复核方块：保持白底（"空心=复核"的语义不变），把它的**边线**换成同色系深 2 度；
            #   mew 保留 0.7 —— 这条边线就是该标记唯一的线重，不再按尺寸比例缩放。
            axr.plot([rb], [yy], marker="s", ms=3.0 * P_UP, mfc="white", mew=0.7, zorder=3,
                     color=cold)
        else:
            # P0：工具不可检出 -> 只画红叉，绝不画连线/方块（不暗示方向一致）
            # P1（2026-10-01）：mew 1.1 -> 1.0 pt（Nature：描边 0.25-1 pt @最终尺寸）。
            # 2026-10-01 加点改版：ms 4.6 -> 6.5（面积 x2）、色 #3D9F3C -> 深 2 度；
            #   ★ mew 保留 1.0：叉是纯线型标记，"外轮廓"对它无意义，线重即其全部视觉权重。
            axr.plot([db], [yy], marker="x", ms=4.6 * P_UP, mew=1.0, color=NR_D, zorder=4)
        ylab2.append("%s | %s -> %s%s%s"
                     % (r["gene"], r["cell_type_discovery"], rc,
                        " (m)" if rc in MERGED else "",
                        "" if ok else "  [NR]"))
        yp2.append(yy); yy -= 1
    axr.axvline(0, color="#7F7F7F", lw=0.5)
    axr.set_yticks(yp2); axr.set_yticklabels(ylab2, fontsize=5.6)
    # NR 行的轴标签整体标红，替代原先塞在面板内的 "not reportable" 文字（消除压叠）
    for lbl, r in zip(axr.get_yticklabels(), rep):
        if r["effect_reportable"] != "TRUE":
            lbl.set_color(FS.C["nr"])
    axr.set_xlabel("MR effect (beta)")
    axr.spines["left"].set_visible(False); axr.tick_params(axis="y", length=0)
    axr.plot([], [], marker="o", ls="none", ms=3.2 * P_UP, color="#333333",
             mec=NEUTRAL_D, mew=_mew_for(3.2 * P_UP), label="Discovery")
    axr.plot([], [], marker="s", ls="none", ms=3.2 * P_UP, mfc="white", mew=0.7,
             color=NEUTRAL_D, label="Replication")
    axr.plot([], [], marker="x", ls="none", ms=4.6 * P_UP, mew=1.0, color=NR_D,
             label="instrument not\ndetectable")
    axr.legend(loc="lower right", fontsize=5.2, handletextpad=0.35, labelspacing=0.35,
               borderpad=0.35)
    FS.finalize(axr); FS.panel_label(axr, "c", dx=FIG2_C_LABEL_DX, dy=FIG2_C_LABEL_DY)
    fig.subplots_adjust(bottom=0.16)
    # 2026-10-01：必须在 subplots_adjust **之后**把 (a) 左移，否则会被重新布局覆盖。
    _pa = axm.get_position()
    axm.set_position([_pa.x0 - FIG2_A_SHIFT, _pa.y0, _pa.width, _pa.height])
    if SHOW_NOTES:
        fig.text(0.5, 0.004,
                 "Cell-type mapping (discovery -> replication). 1M-scBloodNL annotation is MERGED (m):\n"
                 "B_MEM / B_IN -> B ;  Mono_C / Mono_NC -> monocyte ;  CD4_NC / CD4_ET -> CD4T ;\n"
                 "CD8_NC / CD8_ET / CD8_S100B -> CD8T ;  NK -> NK.  A merged-monocyte result is NOT\n"
                 "replication of Mono_NC.  Replication is direction-consistent NOMINAL replication, driven\n"
                 "by instrument detectability at N = 119.  [NR] = the instrument was undetectable in the\n"
                 "replication cohort (CDC42|B_MEM F = 2.01 ; CDC42|B_IN F = 7.2e-05): no effect estimate\n"
                 "from these rows is reportable, and no direction-consistency claim is made for them.\n"
                 "Both instruments act as REGION tags: each regulates CDC42 and LINC00339 in the same cell\n"
                 "types, and neither regulates WNT4 in either layer tested (see Fig. 4f).",
                 ha="center", va="bottom", fontsize=5.2, color="#555555")
    return fig


# ============================================================ Fig3
def fig3(scale=1.0):
    col = load("34_coloc_summary.csv")
    sen = load("40_coloc_sensitivity_matrix.csv")
    cc = load("43_cond_coloc_CDC42.csv")          # task 3.1 conditional colocalization

    _p = LAY["fig3"]
    fig = plt.figure(figsize=(7.2 * scale, _p["h"] * scale))
    gs = gridspec.GridSpec(3, 2, height_ratios=_p["hr"], hspace=_p["hspace"], wspace=0.34)
    axa = fig.add_subplot(gs[0, :]); axb = fig.add_subplot(gs[1, 0])
    axc = fig.add_subplot(gs[1, 1]); axd = fig.add_subplot(gs[2, :])

    # (a) PP.H4 by pair —— P0：strong_shared 拆为「稳健 / 提示性」两档
    ordg = {"CDC42": 0, "LINC00339": 1, "YME1L1": 2, "ANXA4": 3}
    col = sorted(col, key=lambda r: (ordg.get(r["gene"], 9), r["cell_type"]))
    y = np.arange(len(col))[::-1]
    axa.barh(y, [num(r["abf_H4"]) for r in col],
             color=[TIER_C[coloc_tier(r)] for r in col],
             height=0.72, edgecolor="#4D4D4D", lw=0.3)
    axa.axvline(0.8, color=FS.C["region"], ls="--", lw=0.6)
    # ★★ 2026-10-01 22:1x（用户第五轮逐图 · Fig3-A）：
    #    需求原话「图A中绿色的 PP.H4=0.8 往上移动，不要跟绿色的柱状图重叠，看不清字了」。
    #    实测（`_s208_fig3a_label_probe.py` / `_s209_fig3a_top_zoom.py`，**精确色分离**：
    #    区域色 GB12 `#3D9F3C` 在面板 A 里只出现在 x=0.8 的虚线与本注记上，
    #    与柱色 `#367DB0`(robust)/`#519D78`(suggest)/`#AADCA9`(moderate)/`#BFBFBF`(no_shared)
    #    都不同 ⇒ 可把注记从柱里干净地分离出来）：
    #      · 注记墨迹 bbox = x 137.63–148.98、y **7.155–8.721**（高 1.566 mm）
    #      · 最上一根柱（CDC42 | B_IN，suggest_strong `#519D78`，值 0.864）= y **7.428–9.695**
    #        ⇒ 二者**竖向重叠 1.228 mm**（注记下缘 8.721 深入柱内）—— 就是"看不清字"的成因；
    #          后半段（x > 145.588）落在柱外，所以只有前 4 个字被吃掉。
    #      · 面板 A 的 ylim ≈ (−1.096, 15.096)，1 数据单位 = 3.2386 mm
    #        ⇒ 最上柱顶到绘图区上缘只有 **2.384 mm** 净空，而注记高 1.566 mm ⇒ 必须贴着上缘放。
    #    ⇒ 把锚点从「data y = len(col) − 0.55（va=top）」改为**轴分数坐标 y = 0.995（va=top）**，
    #      即贴着绘图区上缘、留在 legend 下方；实测与柱顶净空 ≈0.5 mm、与上缘 ≈0.26 mm。
    #    ★ 用**轴分数**而非 data 坐标：不依赖 ylim 的具体数值，日后重标定不会失效。
    #    ★ 注记仍完全落在 axes **内部** ⇒ tight bbox / 幅面 182.80 × 169.16 mm 不变（`_s168` 复核）。
    axa.text(0.805, 0.995, "PP.H4 = 0.8", transform=axa.get_xaxis_transform(),
             fontsize=6.2, color=FS.C["region"], va="top")
    axa.set_yticks(y); axa.set_yticklabels(["%s | %s" % (r["gene"], r["cell_type"]) for r in col],
                                           fontsize=6.0)
    axa.set_xlabel("coloc.abf  PP.H4"); axa.set_xlim(0, 1.05)
    axa.spines["left"].set_visible(False); axa.tick_params(axis="y", length=0)
    hs = [plt.Rectangle((0, 0), 1, 1, fc=TIER_C[k]) for k in
          ["robust_strong", "suggest_strong", "moderate_shared", "no_shared"]]
    # ★★ 2026-10-01 22:3x（用户第六轮逐图 · Fig3-A）：
    #   需求原话「图A最上面的图标往下移，靠近最上面的绿色横行柱状」。
    #   实测（`_s212_fig3_precise.py`）：legend 的 4 个 handle 在 y **1.312–2.667**、
    #   内容行墨迹到 y **3.006**（盒底 ≈3.77 mm），而绘图区上缘在 **5.00 mm**、
    #   最上面那根柱（CDC42 | B_IN，suggest_strong）顶在 **7.428 mm**
    #   ⇒ 原来 legend 盒底距柱顶 **≈3.7 mm**，中间一大片空白。
    #   ★ 关键：legend 右缘只到 **x≈119.7 mm**，而 A 面板那条 "PP.H4 = 0.8" 注记在
    #     x 137.6–149.0 / y 5.334–6.858 ⇒ **两者 x 不重叠**，所以 legend 可以一路下移
    #     到柱顶之上而不撞注记（这是本次能"贴近柱子"的前提）。
    #   ⇒ 锚点 axes 分数 **1.005 → 0.9468**。标定关系（实测反推）：
    #        盒底_y_mm = anchor_y_mm − 0.95（borderaxespad）；anchor_y_mm = 57.49 − f×52.45
    #        ⇒ f=0.9468 ⇒ anchor 7.83 mm ⇒ 盒底 **6.88 mm**，距柱顶净空 **0.55 mm**。
    #   ★ 代价：legend 会略微越过绘图区上缘（盒高 3.26 mm > 柱顶上方净空 2.43 mm），
    #     但 A 面板无上/右脊柱，视觉上无碍；且最上墨迹由 legend 变为**字母 A**（y 2.667）
    #     ⇒ tight bbox 上边界下移、幅高略降（见日志 §17）。
    axa.legend(hs, ["strong, robust (n=2)", "strong, suggestive (n=3)",
                    "moderate_shared (n=5)", "no_shared (n=5)"],
               loc="lower left", bbox_to_anchor=(0.0, 0.9468), ncol=4, fontsize=5.4,
               columnspacing=0.9, handletextpad=0.35)
    FS.finalize(axa)

    # (b) sensitivity heatmap 5 loci x 9 combos ; robust group first, separated
    loci, wl, pl = [], ["500Kb", "1Mb", "2Mb"], ["1e-4", "1e-5", "1e-6"]
    for r in sen:
        if r["locus"] not in loci:
            loci.append(r["locus"])
    loci = [x for x in ROBUST if x in loci] + [x for x in loci if x not in ROBUST]
    M = np.full((len(loci), len(wl) * len(pl)), np.nan)
    for r in sen:
        i = loci.index(r["locus"]); j = wl.index(r["window"]) * len(pl) + pl.index(r["p12"])
        M[i, j] = num(r["PP.H4"])
    im = axb.imshow(M, cmap=FS.cmap_for("seqA"), vmin=0.3, vmax=1.0, aspect="auto")
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            v = M[i, j]
            axb.text(j, i, "%.2f" % v, ha="center", va="center", fontsize=5.6,
                     color=txt_on(im.cmap, im.norm, v))
            if v < 0.8:
                axb.add_patch(Rectangle((j - .5, i - .5), 1, 1, fill=False,
                                        ec=FS.C["sig"], lw=1.0))
    axb.set_yticks(range(len(loci)))
    axb.set_yticklabels(["%s  (%s)" % (l, "robust" if l in ROBUST_SET else "suggestive")
                         for l in loci], fontsize=6.0)
    axb.set_xticks(range(len(wl) * len(pl)))
    axb.set_xticklabels(["%s\n%s" % (w, p) for w in wl for p in pl], fontsize=5.4)
    axb.set_xlabel("window (upper row) / p12 (lower row)")
    axb.axvline(2.5, color="#4D4D4D", lw=0.6)
    axb.axvline(5.5, color="#4D4D4D", lw=0.6)
    axb.axhline(len(ROBUST) - 0.5, color="#2B2B2B", lw=1.0)
    cb = fig.colorbar(im, ax=axb, fraction=0.045, pad=0.03)
    cb.set_label("PP.H4", fontsize=6.2); cb.ax.tick_params(labelsize=5.8)
    axb.tick_params(axis="both", length=2)

    # (c) PP.H4 vs p12 (mean over windows)
    style = {"CDC42_B_MEM": ("strong", FS.C["robust_strong"]),
             "CDC42_Mono_NC": ("strong", FS.C["risk_up"]),
             "CDC42_B_IN": ("unstable", FS.C["locus1"]), "CDC42_Mono_C": ("unstable", FS.C["nominal"]),
             "ANXA4_Mono_NC": ("unstable", FS.C["gene_ym"])}
    xp = [0, 1, 2]
    for loc in loci:
        ys = []
        for p in pl:
            v = [num(r["PP.H4"]) for r in sen if r["locus"] == loc and r["p12"] == p]
            ys.append(np.mean(v))
        tag, cl = style.get(loc, ("", "#333333"))
        # ★ 2026-10-01（全局「点的大小与风格」标准，见 _fig_style.py）：
        #   面积 ×2（ms ×√2）+ **同色系深 2 度**描边（L* −20）+ mew = 直径的 10%。
        _ms3 = FS.pt_up(3.0)
        axc.plot(xp, ys, marker="o", ms=_ms3, lw=0.8, color=cl,
                 mec=FS.edge_color(cl), mew=FS.mew_for(_ms3),
                 ls=("-" if tag == "strong" else "--") if tag else "-",
                 label="%s (%s)" % (loc, "robust, 9/9" if loc in ROBUST_SET
                                    else "suggestive, < 0.8 at 1e-6"))
    axc.axhline(0.8, color=FS.C["region"], ls=":", lw=0.6)
    axc.text(0.02, 0.815, "PP.H4 = 0.8", fontsize=6.0, color=FS.C["region"])
    axc.set_xticks(xp); axc.set_xticklabels(pl)
    axc.set_xlabel("p12 prior"); axc.set_ylabel("PP.H4 (mean over windows)")
    axc.set_ylim(0.3, 1.03)
    # ★ 图例标记**绝对尺寸不变**（图上点放大了 √2 ⇒ markerscale 同除 √2），见 _fig_style.py
    axc.legend(fontsize=5.6, loc="lower left", ncol=1,
               markerscale=FS.markerscale_keep(1.0))
    FS.finalize(axc)

    # (d) conditional colocalization (task 3.1): conditioning on the WNT4 credible-set lead
    #     abolishes the shared signal; conditioning on an alternative lead does not.
    Q1 = "rs56318008"
    Q2 = {"B_MEM": "rs2473247", "Mono_NC": "rs12048511"}

    def h4(cell, snp, tag):
        for r in cc:
            if (r["scope"] == cell and r["keep_rule"] == "0.01" and
                    r["cond_snp"] == snp and r["tag"] == tag):
                return num(r["PP.H4"])
        return float("nan")

    cells = ["B_MEM", "Mono_NC"]
    xg = np.arange(len(cells))
    w = 0.26
    series = [
        ("unconditional", [h4(c, Q1, "uncond") for c in cells], FS.C["robust_strong"]),
        ("+ WNT4 credible-set lead", [h4(c, Q1, "cond_both") for c in cells], FS.C["suggest_strong"]),
        ("+ alternative lead", [h4(c, Q2[c], "cond_both") for c in cells], FS.C["cyan"]),
    ]
    for k, (lab, vals, cl) in enumerate(series):
        axd.bar(xg + (k - 1) * w, vals, width=w, color=cl, edgecolor="#4D4D4D", lw=0.3,
                label=lab)
        for i, v in enumerate(vals):
            axd.text(xg[i] + (k - 1) * w, v + 0.018, "%.2f" % v, ha="center", va="bottom",
                     fontsize=5.2)
    # ★★ 2026-10-01 22:3x（用户第六轮逐图 · Fig3-D）：
    #   需求原话「图D穿过柱状图的绿色虚线最两端缩短，绿色的 PP.H4=0.8 也随着虚线的缩短往左边移动」。
    #   实测（`_s212_fig3_precise.py`）：虚线短划 x **28.872 → 169.803**（长 140.931 mm，
    #   `axhline` 默认铺满绘图区 28.700–170.610）；而 D 的两组柱只占
    #   x **39.667–92.204** 与 **107.191–159.685** ⇒ 左右各空出 ~10.8 / ~10.9 mm。
    #   ⇒ 两端各内缩到「柱组外缘 + 4.0 mm」：
    #        xmin = (39.667 − 4.0 − 28.700) / 141.910 = **0.04910**
    #        xmax = (159.685 + 4.0 − 28.700) / 141.910 = **0.95122**
    #   ⇒ 线段 x 34.667–163.685 mm（长 129.02 mm，比原来短 11.91 mm、两端各缩 5.95 mm 视觉量）。
    #   ★ 注记随之左移，**保持原来的 2.13 mm 间距**：
    #        原 172.740 = 170.610 + 2.13 ⇒ 新 163.685 + 2.13 = 165.815 mm
    #        ⇒ axes 分数 = (165.815 − 28.700) / 141.910 = **0.9664**（原 1.015）
    #   ★★ 副作用（须记）：注记原是全图**最右墨迹**（182.21 mm），左移后右界缩到 ~175.2 mm
    #      ⇒ tight bbox 变窄 ⇒ `fit_save` 会放大整图去凑 182.4 mm 宽 ⇒ 幅高随之抬高，
    #        需要重标定 `LAY["fig3"]["h"]`（见日志 §17）。
    axd.axhline(0.8, color=FS.C["region"], ls="--", lw=0.6, xmin=0.04910, xmax=0.95122)
    axd.set_xticks(xg); axd.set_xticklabels(["CDC42 | B_MEM", "CDC42 | Mono_NC"], fontsize=6.0)
    axd.set_ylim(0, 1.30); axd.set_xlim(-0.55, 1.55)
    axd.set_ylabel("coloc.abf  PP.H4")
    axd.text(0.9664, 0.615, "PP.H4 = 0.8", transform=axd.transAxes, ha="left", va="center",
             fontsize=5.2, color=FS.C["region"], clip_on=False)
    axd.legend(fontsize=5.2, loc="upper center", bbox_to_anchor=(0.5, 1.005), ncol=3,
               columnspacing=1.4, handletextpad=0.4)
    FS.finalize(axd)

    fig.subplots_adjust(bottom=LAY["fig3"]["bottom"], top=LAY["fig3"]["top"])
    # ★★ 全局标准（2026-10-01 用户定案）：「各行面板组居中于整幅」。
    #   ★ 必须在 subplots_adjust **之后**调用，否则会被 gridspec 重新布局覆盖。
    #   ★ (b) 的 colorbar 是**独立 axes**（cb.ax），必须与 (b) 放进同一组一起平移。
    _G3 = [[axa], [axb, cb.ax, axc], [axd]]
    CENTER_RESID["fig3"] = FS.center_row_groups(fig, _G3)
    # ★★ 2026-10-01 21:5x（用户第三轮逐图 · Fig3）：「图D的Y轴与图A对齐」。
    #   A（gs[0,:]）与 D（gs[2,:]）本是**同一 GridSpec 列**，宽度都是 149.65 mm；
    #   但 `center_row_groups()` 逐**行**独立平移 ⇒ 实测（成品 PNG 底部横轴脊柱）：
    #       A 绘图区 [28.700, 178.350]   D 绘图区 [14.140, 163.790]
    #   ⇒ 左缘相差 **14.560 mm**，用户看到的就是"D 的 Y 轴单独往左凸出去"。
    #   ⇒ 居中之后把 axd 搬回 axa 那一列（宽度相同 ⇒ 右缘同时对齐，两面板完全叠合）。
    #   ★ 本操作是「各行面板组居中」的**显式例外**：第 3 行不再自居中（用户明确要求优先）。
    #   ★ 必须在面板字母锚定**之前**：字母锚定读的是 axes 的最终位置。
    FS.align_axes_left(axa, [axd], fig=fig, verbose=True)
    # ★★ 全局标准（2026-10-01）：「面板字母 ↔ 本面板 Y 轴最上面那个词」的距离 = Fig2 标准
    #    （Δx = 6.41 mm / Δy = 2.42 mm）。必须在**最终**水平布局之后锚定，
    #    因为 letter_offset_for() 要真实渲染才量得到 axes 的最终位置与文字 bbox。
    #   ★★ 2026-10-01（用户逐图细改 · Fig3）：「D 要跟 A、B 对齐」。
    #      A（gs[0,:]）、B（gs[1,0]）、D（gs[2,:]）的**绘图区左缘相同**，是同一列；
    #      但「锚到各自 Y 轴最上面那个词」会让三者 x 各不相同（实测 A 8.30 / B 0.51 /
    #      D 4.10 mm —— v10 里 D 更是跑到 0.72 mm，用户就是看到这个才提的）。
    #      ⇒ 改用 `align_letter_column()`：三者在标准 dy 下共用一条竖线（取各面板标准锚点中
    #        最靠右的那个候选，既贴近标签列又不会被某个宽标签拖到画布边缘）。
    #      C 在右列、没有同列伙伴 ⇒ 仍走逐面板锚定。
    #   ★★ 2026-10-01 20:5x（用户第二轮逐图微调，三处「毫米级」手工位移）：
    #      · **B 上移 2.905 mm** —— 原 B 字母 y0=72.504、C 字母 y0=69.599（同属第 2 行、
    #        但两面板"最上面那个词"的 y 位置不同）⇒ 用户要 B「跟 C 在同一水平线」。
    #      · **D 上移 1.50 mm** —— 原 D 字母下缘 139.396 与其旋转 y 轴标题 `coloc.abf PP.H4`
    #        上缘 139.860 只差 **0.464 mm**（视觉上"连在一起"）⇒ 上移后净空 ≈1.96 mm。
    #        ★ 注意：D 字母是第 3 行的**最上墨迹**，上移后 `set_row_ink_gaps` 会把第 3 行
    #          整体下压同样距离以维持 4.74 mm 行带 ⇒ 幅高会随之升高（实测见日志）。
    #      · **C 右移 1.50 mm** —— 原 C 字母右缘 106.682 距本面板 y 刻度左缘 113.045 有
    #        6.363 mm（≈标准 6.41）⇒ 用户要「往右移一点点」。
    #   ★★ 2026-10-01 21:1x（用户第三轮逐图微调 · Fig3-D 再上移，落到 Y 轴线左上角）：
    #      需求原话「D字母再往上移。在Y轴线的左上角。」
    #      实测（成品 PNG 墨迹口径，`_s189`/`_s191`/`_s192`）：
    #        · D 面板 **Y 轴脊柱**（左竖线）= x 14.139 mm、顶点 y **137.162 mm**（长 27.86 mm）
    #        · 参照 B 面板：字母墨迹下缘 72.052、热图顶 74.042 ⇒ **字母下缘高于绘图区上缘 1.99 mm**
    #        · 改前 D 字母墨迹 = x 5.588–7.324、y **136.188–138.220** ⇒ 下缘压在 Y 轴线顶点
    #          **下方 1.058 mm**（即字母横跨轴线顶端，不是"左上角"）
    #      ⇒ 位移 **1.50 → 3.05 mm**（再上移 1.55 mm）⇒ D 字母下缘落到轴线顶点**上方 ≈0.45 mm**，
    #        即真正落在「Y 轴线的左上角」。
    #      ★ 代价：D 仍是第 3 行最上墨迹 ⇒ 行带工具把第 3 行下压同样距离，幅高约 +0.9 mm
    #        （仍须 <= 170 mm 的 P0 硬约束，见下方 h 标定与日志 §14）。
    #      ★★ 21:3x 复核后定稿：3.05 mm 时实测**间隙只有 0.550 mm**（`_s194`），论文缩印后
    #        几乎贴住轴线 ⇒ 再上移 0.625 mm（总 3.675），实测间隙 ≈1.2 mm ≈ 3.4 pt，
    #        与 B 面板（1.95 mm）同量级但更贴合「就在轴线左上角」。h 同步 6.500 -> 6.475。
    FS.align_letter_column([(axa, "a"), (axb, "b"), (axd, "d")],
                           shifts_mm={"B": (0.0, 2.905), "D": (0.0, 3.675)}, verbose=True)
    FS.panel_label_anchored(axc, "c", shift_mm=(1.50, 0.0))
    # ★★ 全局标准（2026-10-01）：「图与图之间的距离」= Fig2 的可见空白带 4.74 mm。
    #    ★ 必须排在面板字母**之后**：字母是下行最上方的墨迹，先放字母再量才准。
    #    ★ 本图两个边界的文字外溢量不同（B/C 行下方是两行刻度 + 轴标题，D 行上方只有图例）
    #      ⇒ 单一 hspace 数学上无法让两个边界同时 = 4.74 ⇒ 逐边界测量并竖直平移。
    ROW_GAP_RESID["fig3"] = FS.set_row_ink_gaps(fig, _G3, [4.74, 4.74])
    if SHOW_NOTES:
        fig.text(0.5, 0.005,
                 "Panel (d), conditional colocalization: PP.H4 at the two robust CDC42 loci after conditioning on the\n"
                 "WNT4 outcome credible-set lead rs56318008 (both datasets) versus on an alternative lead (B_MEM\n"
                 "rs2473247 / Mono_NC rs12048511). Conditioning on the WNT4 lead abolishes the shared signal\n"
                 "(0.99 -> 0.06-0.11), i.e. the shared component is carried by the region lead, not by CDC42.\n"
                 "Two-tier naming of coloc.abf strong_shared (identical to Results 3.3 / 3.4 and to the summary\n"
                 "table): robust = PP.H4 > 0.8 in all 9 window x p12 settings; suggestive = drops below 0.8 only at\n"
                 "the strictest prior p12 = 1e-6 (red boxes in b). Only the two robust loci (CDC42_B_MEM,\n"
                 "CDC42_Mono_NC) enter the primary conclusions.",
                 ha="center", va="bottom", fontsize=5.2, color="#555555")
    return fig


# ============================================================ Fig4
def fig4(scale=1.0):
    ld = load("23_chr1_LD_matrix.csv")
    V = load("35_finemap_variant.csv")
    o = [r for r in V if r["scope"] == ""]

    # P0：删除「ΣPIP 按最近基因相加」的归属排序，改为 credible-set lead 落点计数
    fp = load("36_finemap_pair.csv")
    pairs = [r for r in fp if r["scope"] != "outcome_only"]
    cnt = {"WNT4": 0, "LINC00339": 0, "CDC42": 0}
    for r in pairs:
        g = r["cs_lead_outcome_nopal_gene"]
        if g == "WNT4":
            cnt["WNT4"] += 1
        elif "LINC00339" in g:
            cnt["LINC00339"] += 1
        elif g == "CDC42":
            cnt["CDC42"] += 1

    _p = LAY["fig4"]
    # ★★★ 2026-10-01（逐图第九轮 · Fig4）：行 0 的行高**必须正好等于 (a) 的正方形边长**
    #   （`row0_mm` = 46 mm）。理由：本轮 (a) 不再由列宽锁定（列宽 66.9 ≫ 46），而是
    #   **行高锁定**；若行 0 的格子高于 46，则 (a) 上/下会留出纯白（还白占幅高预算）；
    #   若低于 46，(a) 就长不到 46。
    #   GridSpec 机制（见 MEMORY「行高只由 h 与 top/bottom 决定」）：
    #       Σ行高 = (top − bottom)·H_nom / (1 + 2·hspace/3)      （hspace=0.30 ⇒ 系数 0.8333）
    #       hr 只做**重分配** ⇒ 令 hr[0] ≡ 1.0，把剩余行高按原版行 1 : 行 2 = 0.95 : 0.88 分。
    #   ⇒ 行 1 / 行 2 会随 h 一起长高，(a) 恒定 46 mm。改 h 只影响第 2/3 排 + 总幅高。
    _FH_nom = _p["h"] * scale * 25.4                      # 名义幅高 mm
    _FW_nom = 7.2 * scale * 25.4                          # 名义幅宽 mm
    _sum_row_mm = ((_p["top"] - _p["bottom"]) * _FH_nom
                   / (1.0 + 2.0 * _p["hspace"] / 3.0))
    _r0_mm = float(_p["row0_mm"])
    _rem_mm = max(_sum_row_mm - _r0_mm, 1.0)
    _hr = [1.0, (_rem_mm * 0.95 / 1.83) / _r0_mm, (_rem_mm * 0.88 / 1.83) / _r0_mm]
    ROW0_HR_RESID["fig4"] = dict(hr=_hr, sum_row_mm=_sum_row_mm, row0_mm=_r0_mm,
                                 rem_mm=_rem_mm)

    def _mm2x(_mm):
        """绝对 mm → figure 分数（x）。★ 用**名义幅宽**换算 ⇒ 落盘后就是该 mm。"""
        return _mm / _FW_nom

    def _mm2y(_mm):
        """绝对 mm → figure 分数（y）。"""
        return _mm / _FH_nom

    def _x2mm(_frac):
        return _frac * _FW_nom

    def _y2mm(_frac):
        return _frac * _FH_nom

    fig = plt.figure(figsize=(7.2 * scale, _p["h"] * scale))
    gs = gridspec.GridSpec(3, 2, height_ratios=_hr, hspace=_p["hspace"], wspace=0.34)
    axa = fig.add_subplot(gs[0, 0]); axb = fig.add_subplot(gs[0, 1])
    axc = fig.add_subplot(gs[1, :]); axd = fig.add_subplot(gs[1, 1])
    axc.remove(); axd.remove()
    gs2 = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs[1, :], wspace=0.34,
                                           width_ratios=[2.0, 1.0])
    axc = fig.add_subplot(gs2[0]); axd = fig.add_subplot(gs2[1])
    axe = fig.add_subplot(gs[2, 0]); axf = fig.add_subplot(gs[2, 1])

    # (a) LD heatmap
    labels = [r[list(r.keys())[0]] for r in ld]
    M = np.array([[num(r[k], 0.0) for k in list(ld[0].keys())[1:]] for r in ld])
    im = axa.imshow(M, cmap=FS.cmap_for("seqB"), vmin=0, vmax=1)
    axa.set_xticks(range(len(labels))); axa.set_yticks(range(len(labels)))
    short = [l.split(" / ")[0] for l in labels]
    # ★ 2026-10-01（逐图第八轮 · Fig4 要求 1）用户原话「图AX轴标签斜45度摆放」：
    #   由 90°（竖直）改 45°，`ha="right"` + `rotation_mode="anchor"` 让文字自左下向右上、
    #   右端锚在刻度上。★ 附带好处：竖直外溢由 8.6 mm 降到约 7.5 mm，**不增加行带高度**。
    axa.set_xticklabels(short, rotation=45, ha="right", rotation_mode="anchor", fontsize=5.4)
    axa.set_yticklabels([l.split(" / ")[0] for l in labels], fontsize=5.4)
    for i in range(len(labels)):
        for j in range(len(labels)):
            axa.text(j, i, "%.2f" % M[i, j], ha="center", va="center", fontsize=5.2,
                     color=txt_on(im.cmap, im.norm, M[i, j], dark="#333333"))
    # ★★ 2026-10-01（逐图第八轮 · Fig4 要求 1）**关键修正**：色条必须画在**独立的 cax**
    #   上，不能用 `fig.colorbar(im, ax=axa)`。原因：后者会让 matplotlib 把父轴宽度压成
    #   `1 − fraction − pad = 0.924` 倍（fraction=0.046 / pad=0.03），实测 (a) 的格子宽由
    #   列宽 60.57 mm 掉到 **55.97 mm** ⇒ 即使行 0 行高足够，(a) 也**永远长不到 (e) 的
    #   60.57 mm**（这正是上一轮「(a) 与 (e) 不等宽」的隐藏原因）。
    #   改用 cax 后：色条**不再侵占 (a)**，(a) 的正方形边长 = 列宽 = (e) 的 X 轴长。
    #   cax 的位置在下方定位块里按「(a) 右缘 + pad」显式给出（允许伸进两列之间的白带）。
    cax = fig.add_axes([0.5, 0.5, 0.01, 0.1])      # 占位，稍后由定位块接管
    cb = fig.colorbar(im, cax=cax)
    cb.set_label("r2 (1000G EUR)", fontsize=6.0); cb.ax.tick_params(labelsize=5.6)
    axa.tick_params(length=1.5)

    # (b) gene / variant track (GRCh37)
    genes = [("CDC42", 22388297, 22397199, GENE_C["CDC42"]),
             ("LINC00339", 22352108, 22355978, GENE_C["LINC00339"]),
             ("WNT4", 22444975, 22470451, GENE_C["WNT4"])]
    for k, (g, s, e, c) in enumerate(genes):
        axb.add_patch(Rectangle((s, k - 0.16), e - s, 0.32, fc=c, ec="#3A3A3A", lw=0.4))
        axb.text((s + e) / 2, k + 0.34, g, ha="center", fontsize=6.2)
    # P0（2026-10-01）：把两处 instrument 注记由 3 行并为 2 行（"rsID" / "细胞型 + r2"）。
    # ★ 2026-10-01（逐图第八轮 · Fig4 要求 2）用户原话「图B的X轴缩短，使其长度和图F一致，
    #   且调整箭头下字的大小，不要跟箭头重叠」——第二句的落点（第一句见下方定位块）：
    #     ① 字号**保持 5.2 pt**（项目下限 / Nature 图内文字下限 5 pt；曾经试到 4.8 pt，
    #        `_fig_overlap_check` 立即报出「低于 5 pt 的字符数 = 133」⇒ 已回退，不得下调）；
    #        改为把每处注记**拆成 3 行**（rsID / instrument 类型 / r2），最长行
    #        "r2 = 0.899 vs CS lead" 由 ~22 mm 收窄到 ~20 mm；
    #     ② 三处 y 分层为 −0.95 / −1.95 / −3.15（相邻层差 1.0 数据单位 ≈ 9.8 mm）；
    #     ③ **按 x 邻域改变对齐**，这是消除「字压箭头」的关键（三根箭头杆分别在
    #        x = 22422721 / 22462111 / 22470407，后两根只差 8296 bp ≈ 3.4 mm）：
    #          · rs10917151（最左）ha="center" —— 文本块 x ∈ [22.398, 22.447] Mb，
    #            左右各 16 mm / 12 mm 内无别的杆；
    #          · rs12037376（居中）ha="right"  —— 文本整体退到标记左侧
    #            （x ∈ [22.413, 22.462] Mb），**避开 rs56318008 那根垂到 y=−2.97 的杆**；
    #          · rs56318008（最右）ha="right"  —— 居中对齐会**越出 x 上限 22.49 Mb** 约 2 mm
    #            （文本宽 20 mm，而标记右缘只剩 7.9 mm），故同样右对齐到标记。
    #     ④ 逐个复核（文本块 y 区间 vs 每根杆的 y 区间）：
    #          杆1 [−0.28,−0.77] 杆2 [−0.28,−1.77] 杆3 [−0.28,−2.97]
    #          文1 [−1.17,−1.87] 文2 [−2.17,−2.87] 文3 [−3.37,−4.07]
    #        ⇒ 与自身杆**只在端点相接**、与其他任何杆**y 区间完全不相交** ⇒ 0 处重叠；
    #        `_fig_overlap_check`（Text 对象口径）复核 = **0 处**。
    # ★★ 2026-10-01（逐图第九轮 · Fig4）**(b) 因行 0 由 70.9 → 46 mm 而变矮**
    #   ⇒ 三处注记必须重排，否则 3 行文本块（5.2 pt × 3 行 × linespacing 1.25 ≈ 6.8 mm）
    #     在更小的"每个数据单位对应 mm"下会互相挤压、并越过 ylim 下界。
    #   算法（不缩字号，只调**结构**）：
    #     · 保持 3 行/块（r8 已用它把最长行压到 ≈19 mm，避免压到相邻箭头杆）；
    #     · 文本块高 ≈ 6.8 mm = **1.145 数据单位**（46 mm / 7.75 单位 = 5.94 mm/单位）
    #       ⇒ 相邻层差必须 ≥ 1.145 + 1.22 mm 净空 = **1.35 单位**；
    #     · 取 −0.95 / −2.30 / −3.65；ylim 下界 −4.35 → **−5.25**：
    #         文 1 [−1.17, −2.315]、文 2 [−2.52, −3.665]、文 3 [−3.87, −5.015]
    #       ⇒ 块间净空各 1.22 mm、文 3 距 ylim 下界 1.40 mm。
    #     · 箭头杆 y 区间 [−0.28, my+0.18] = [−0.77 / −2.12 / −3.47]
    #       ⇒ 与**自身**杆只在端点相接、与**其他**杆的 y 区间完全不相交 ⇒ 0 处重叠。
    marks = [(22422721, ["rs10917151", "Mono_NC instr.", "r2 = 0.847 vs CS lead"],
              -0.95, "center"),
             (22462111, ["rs12037376", "B_MEM instr.", "r2 = 0.899 vs CS lead"],
              -2.30, "right"),
             (22470407, ["rs56318008", "outcome CS lead", "WNT4 gene body"],
              -3.65, "right")]
    for p, _lines, my, ha in marks:
        axb.plot([p, p], [-0.28, my + 0.18], lw=0.4, color="#9E9E9E", zorder=1)
        # ★ 全局「点的大小与风格」标准（_fig_style.py）：面积 ×2 + 同色系深 2 度描边
        _msv = FS.pt_up(4.0)
        axb.plot([p], [my], marker="v", ms=_msv, color="#2B2B2B",
                 mec=FS.edge_color("#2B2B2B"), mew=FS.mew_for(_msv), zorder=2)
        axb.text(p, my - 0.22, "\n".join(_lines), ha=ha, va="top", fontsize=5.2,
                 linespacing=1.25)
    axb.set_xlim(22325000, 22490000); axb.set_ylim(-5.25, 2.5)
    axb.set_yticks([]); axb.set_xlabel("GRCh37 position on chr1 (bp)")
    axb.set_xticks([22350000, 22400000, 22450000])
    axb.set_xticklabels(["22.35 Mb", "22.40 Mb", "22.45 Mb"], fontsize=5.8)
    axb.spines["left"].set_visible(False)
    FS.finalize(axb)

    # (c) PIP along position
    xs = [int(r["pos37"]) for r in o]; ys = [num(r["PIP"], 0.0) for r in o]
    # ★★ 2026-10-01（逐图第七轮 · Fig4 要求 2，用户原话「图C下面的灰色小点面积再放大2倍，
    #   颜色换成彩色的」+ 追问后选「按所属基因区域着色」）：
    #   ① **面积再 ×2**：在既有全局 P_UP=√2 之上再乘一个 √2 ⇒ 直径 1.6·√2·√2 = 3.2
    #      （面积 = 4 倍原始）；
    #   ② 由单一中性灰 #7F7F7F 改为**按所属基因区域着色**：落在三个基因体内的点用该基因色
    #      （与 (a)(b)(c) 基因配色、下方 axvspan 色带、以及 (d) 的条形同源），
    #      基因体之外统一归入 `INTERGENIC_C`。
    #   ★ 区间口径：(b) 里同一组 GRCh37 坐标（硬编码自注释，与所绘完全一致）。
    #   ★ 实测分布（n=541）：LINC00339 体内 7 / CDC42 体内 7 / WNT4 体内 6 / 体外 521。
    #     ——「只有 20 个点在基因体内」是数据事实，不是绘图遗漏。
    _msc = FS.pt_up(1.6) * math.sqrt(2.0)
    _REG = [("LINC00339", 22352108, 22355978), ("CDC42", 22388297, 22397199),
            ("WNT4", 22444975, 22470451)]

    def _bucket_of(_pos):
        for _g, _s, _e in _REG:
            if _s <= _pos <= _e:
                return _g
        return "intergenic"

    def _col_of(_k):
        return INTERGENIC_C if _k == "intergenic" else GENE_C[_k]

    _buckets = OrderedDict((g, ([], [])) for g, _, _ in _REG)
    _buckets["intergenic"] = ([], [])
    for _x, _y in zip(xs, ys):
        _bk = _bucket_of(_x)
        _buckets[_bk][0].append(_x); _buckets[_bk][1].append(_y)
    for _k, (_bx, _by) in _buckets.items():
        if not _bx:
            continue
        _c = _col_of(_k)
        axc.plot(_bx, _by, ls="none", marker="o", ms=_msc, mec=FS.edge_color(_c),
                 mew=FS.mew_for(_msc), color=_c, zorder=2)
    # ★★ 2026-10-01 修正 —— **颜色键冲突（本轮自查发现，必须消除）**：
    #   原实现把 lead 点画成 `FS.C["sig"]`，但 `GB12_C["sig"] == GB12_C["locus1"] == "#519D78"`
    #   **恰好就是 CDC42 的桶色**（`GENE_C["CDC42"] = FS.C["locus1"]`）。
    #   ⇒ 图例写着「#519D78 = in CDC42」，而被箭头标出的 lead 点 rs56318008（PIP 0.9142）
    #     实际落在 **WNT4** 基因体内（pos37 22470407 ≤ 22470451），却也是 #519D78
    #     ⇒ 读者必然读成「in CDC42」，与图例/文字自相矛盾。
    #   改法：lead 点与其注记**统一改用它自己所属桶的颜色**（此处 = WNT4 蓝 #367DB0），
    #   突出手段回到「更大的标记 + 箭头 + 文字标签」，而不是引入第五种颜色 —
    #   这样四色图例仍然**完备且互斥**（与用户裁定的「按所属基因区域着色」完全一致）。
    _top = sorted(o, key=lambda r: -num(r["PIP"]))[0]
    _tc = _col_of(_bucket_of(int(_top["pos37"])))
    _mst = FS.pt_up(4.2)
    axc.plot([int(_top["pos37"])], [num(_top["PIP"])], marker="o", ms=_mst,
             mec=FS.edge_color(_tc), mew=FS.mew_for(_mst),
             color=_tc, zorder=3)
    axc.annotate("%s (PIP %.4f)\n%s" % (_top["rsid"], num(_top["PIP"]), _top["gene_annot"]),
                 xy=(int(_top["pos37"]), num(_top["PIP"])),
                 xytext=(int(_top["pos37"]) - 130000, 0.66), fontsize=5.8, color=_tc,
                 arrowprops=dict(arrowstyle="->", lw=0.5, color=_tc))
    for g, s, e, c in genes:
        axc.axvspan(s, e, color=c, alpha=0.32, lw=0)
    axc.set_ylim(-0.03, 1.06); axc.set_xlim(min(xs), max(xs))
    axc.set_xticks([22000000, 22300000, 22600000, 22900000])
    axc.set_xticklabels(["22.0", "22.3", "22.6", "22.9"], fontsize=5.8)
    axc.set_xlabel("GRCh37 position on chr1 (Mb)")
    # ★★ 2026-10-01（逐图第八轮 · Fig4 高度预算）**换行 = 消除纯白浪费**：
    #   原单行 y 标签 "posterior inclusion prob. (outcome)" 共 36 字符 ≈ 40.5 mm 长，
    #   **比行 1 的轴高（27 mm）还长** ⇒ 旋转 90° 后上下各外溢 8.6 mm，
    #   把 (c) 的 tightbbox 由 27 mm 撑到 **44.2 mm**（`_s231` 实测 `yaxis tightbbox` 即是元凶）
    #   ⇒ 在 600 dpi 交付 PNG 上就是 (c) 上下两段 ~8.6 mm 的**纯白**，白占 ~17 mm 幅高。
    #   **文字一字未改，仅把它折成 2 行**（最长行 "prob. (outcome)" ≈ 18 mm < 轴高）。
    axc.set_ylabel("posterior inclusion\nprob. (outcome)")
    axc.text(0.015, 0.965,
             "gene bodies (Mb): LINC00339 22.352-22.356 | CDC42 22.388-22.397 | WNT4 22.445-22.470",
             transform=axc.transAxes, fontsize=5.4, color="#444444", va="top", ha="left")
    # ★ 2026-10-01（逐图第七轮 · Fig4 要求 2）：颜色现在承载「变异落在哪个基因体内」这一语义
    #   ⇒ 必须给键。放在 (c) 左侧空白区（PIP 0.15-0.9 的中上部**全空**，唯一的点标注在
    #   x≈0.31-0.63 / y≈0.66，故本图例取 x 0.010-0.18、顶端 y 0.885，二者不相交）。
    # ★★ 2026-10-01（逐图第八轮 · Fig4）**图例句柄必须 `set_in_layout(False)`**：
    #   原实现用 `axc.plot([], [], ...)` 造空折线当句柄 ⇒ 这些**空 Line2D** 是 axes 的
    #   child 且 `in_layout=True`，被 `Axes.get_default_bbox_extra_artists()` 收进
    #   tightbbox 的并集；而空数据的 `get_window_extent` 是**退化/异常**值 ⇒ 实测把 (c) 的
    #   tightbbox 由真实的 ~34 mm 撑到 **44.20 mm**（上下各多 8.61 mm 空白），
    #   (d) 40.70 mm。这些空白会原样进入 `bbox_inches='tight'` 的交付 PNG ⇒ **白占 ~10 mm 幅高**。
    #   改用**脱离 axes 的独立 Line2D** 并显式 `set_in_layout(False)`，图例外观不变。
    _hs = []
    for _k in ("LINC00339", "CDC42", "WNT4", "intergenic"):
        _c = INTERGENIC_C if _k == "intergenic" else GENE_C[_k]
        _hd = matplotlib.lines.Line2D([], [], ls="none", marker="o", ms=2.6, color=_c,
                                      mec=FS.edge_color(_c), mew=0.35)
        _hd.set_in_layout(False)
        _hs.append(_hd)
    axc.legend(_hs, ["in LINC00339", "in CDC42", "in WNT4", "outside gene bodies"],
               loc="upper left", bbox_to_anchor=(0.010, 0.885), fontsize=5.2,
               frameon=False, handletextpad=0.35, labelspacing=0.30,
               borderpad=0.0, borderaxespad=0.0, handlelength=0.9)
    FS.finalize(axc)

    # (d) credible-set lead placement across the 13 chr1 pairs (P0: replaces summed-PIP bars)
    order = ["WNT4", "LINC00339", "CDC42"]
    vals = [cnt[g] for g in order]
    axd.barh(np.arange(3)[::-1], vals, color=[GENE_C[g] for g in order],
             height=0.6, edgecolor="#4D4D4D", lw=0.3)
    for i, g in enumerate(order):
        axd.text(vals[i] + 0.28, 2 - i, "%d" % vals[i], va="center", fontsize=8.0,
                 fontweight="bold", color=(FS.C["region"] if vals[i] == 0 else "#222222"))
    axd.set_yticks(np.arange(3)[::-1])
    axd.set_yticklabels(["WNT4\n(gene body)", "LINC00339\n(flank, 9 kb)", "CDC42\n(gene body)"],
                        fontsize=5.8)
    axd.set_xlabel("chr1 pairs whose credible-set lead\nfalls in this interval  (n = 13)")
    axd.set_xlim(0, 15.0)
    axd.set_xticks([0, 5, 10])
    axd.tick_params(axis="x", labelsize=5.8)
    axd.spines["left"].set_visible(False); axd.tick_params(axis="y", length=0)
    axd.text(0.97, 0.06, "CDC42 never carried the\noutcome credible set",
             transform=axd.transAxes, ha="right", va="bottom", fontsize=5.2, color=FS.C["region"])
    FS.finalize(axd)

    # (e) cis-SMR + HEIDI (task 3.3): b_SMR per cell type, one row per gene
    smr = [r for r in load("45c_smr_3p3_combined.csv") if r["analysis"] == "main(5e-8)"]
    gw = [r for r in load("47_gtex8_wb_smr_results.csv") if r["analysis"].startswith("A(")]
    rows = []
    for g, cl in (("CDC42", GENE_C["CDC42"]), ("LINC00339", GENE_C["LINC00339"])):
        pts = [(num(r["b_SMR"]), 1.96 * num(r["se_SMR"]), r["heidi_verdict"])
               for r in smr if r["Gene"] == g]
        rows.append(("%s\n(OneK1K, %d cells)" % (g, len(pts)), pts, cl))
    rows.append(("WNT4\n(GTEx v8 whole blood)",
                 [(num(r["b_SMR"]), 1.96 * num(r["se_SMR"]), "n/a") for r in gw],
                 GENE_C["WNT4"]))
    for j, (lab, pts, cl) in enumerate(rows):
        yv = len(rows) - 1 - j
        # ★ 全局「点的大小与风格」标准（_fig_style.py）：面积 ×2（ms ×√2）。
        #   实心点 → 加同色系深 2 度描边、mew = 直径 10%；
        #   "空心 = HEIDI heterogeneity"是**纯描边型**标记 ⇒ 保留原 mew=0.9，
        #   只把边线色换成同色系深 2 度（与 Fig2 (c) 的白底复核方块同一处理）。
        _mse = FS.pt_up(3.0)
        for b, e, verd in pts:
            if verd == "heterogeneity":
                axe.errorbar([b], [yv], xerr=[[e], [e]], fmt="o", ms=_mse, mfc="white",
                             mew=0.9, color=FS.edge_color(cl), elinewidth=0.5, capsize=1.2)
            else:
                axe.errorbar([b], [yv], xerr=[[e], [e]], fmt="o", ms=_mse,
                             mec=FS.edge_color(cl), mew=FS.mew_for(_mse),
                             color=cl, elinewidth=0.5, capsize=1.2)
    axe.axvline(0, color="#7F7F7F", lw=0.5)
    axe.set_yticks(range(len(rows)))
    axe.set_yticklabels([r[0] for r in rows][::-1], fontsize=5.6)
    axe.set_xlabel("b_SMR (95% CI)")
    axe.set_ylim(-0.75, len(rows) - 0.25)
    axe.set_xlim(-0.40, 0.62)
    axe.spines["left"].set_visible(False); axe.tick_params(axis="y", length=0)
    axe.text(0.985, 0.05,
             "filled = HEIDI consistent\nopen = HEIDI heterogeneity\n"
             "all OneK1K rows P_SMR FDR < 0.05\nWNT4: n.s. (P_SMR = 0.148, m = 11)",
             transform=axe.transAxes, ha="right", va="bottom", fontsize=5.2, color="#555555",
             linespacing=1.35)
    FS.finalize(axe)

    # (f) instrument neighbouring-gene check (task 3.4): both instruments are REGION tags
    rt = load("49b_task34_region_tag_summary.csv")
    genes3 = ["CDC42", "LINC00339", "WNT4"]
    insts = [("rs12037376", FS.C["robust_strong"]), ("rs10917151", FS.C["cyan"])]
    xg2 = np.arange(len(genes3)); w2 = 0.32
    for k, (s, cl) in enumerate(insts):
        vals = []
        for g in genes3:
            v = 0
            for r in rt:
                if r["SNP_rsID"] == s and r["gene"] == g:
                    v = int(r["OneK1K_n_sig"])
            vals.append(v)
        axf.bar(xg2 + (k - 0.5) * w2, vals, width=w2, color=cl, edgecolor="#4D4D4D", lw=0.3,
                label=s)
        for i, v in enumerate(vals):
            axf.text(xg2[i] + (k - 0.5) * w2, v + 0.3, "%d" % v, ha="center", va="bottom",
                     fontsize=5.4)
    axf.set_xticks(xg2); axf.set_xticklabels(genes3, fontsize=6.0)
    # ★★★ 2026-10-02（第十二轮 · 用户「F图Y轴减少高度，删掉刻度17.5，最高刻度15即可」）：
    #   y 轴上界 17.5 → **15.0**（数据最大柱 = 12 ⇒ 不裁剪），并**显式**给出刻度
    #   [0, 2.5, …, 15]，杜绝自动定位器再生成 17.5。
    #   ★ 只动**数据范围**，不动面板的物理行高 ⇒ E/F 底缘仍对齐（既有硬要求）。
    _fymax = float(_p.get("f_ymax", 15.0))
    axf.set_ylim(0, _fymax)
    _ftick_step = 2.5
    _fticks = [round(i * _ftick_step, 1) for i in range(int(_fymax / _ftick_step) + 1)]
    if _fticks[-1] < _fymax - 1e-9:
        _fticks.append(_fymax)
    axf.set_yticks(_fticks)
    # ★ 2026-10-01（逐图第八轮 · Fig4 高度预算）同上：原来 2 行、最长行 29 字符 ≈ 32.8 mm
    #   > 轴高 25 mm ⇒ (f) 顶沿再外溢 3.9 mm。折成 3 行（最长 "significant cis eQTL" 22 字符
    #   ≈ 24.9 mm ≤ 轴高）⇒ (f) 的顶沿回到字母 F 的锚点，行 2 的墨迹顶由 24.98 → 22.79 mm。
    axf.set_ylabel("cell types with a\nsignificant cis eQTL\n(of 14)")
    axf.set_xlabel("candidate gene in the region")
    # ★★★ 2026-10-02（第十二轮 · 用户「右上角图标往下移，要在删除后的Y轴内」）：
    #   y 轴上界降到 15 后，图例必须整体落在 **0–15** 之内。图例是**轴分数**锚定 ⇒ 只改
    #   `bbox_to_anchor` 的 y 即可：原 `loc="upper right"`（顶 = 轴分数 0.9762，数据 17.08
    #   ⇒ 在旧轴上贴顶、按数据口径"越过 15"）；本轮把锚点顶降到 `f_leg_anchor_y`
    #   ⇒ 图例整体下移并**完全落在新轴内**（实测值记录在 `FYLIM_RESID`）。
    #   ★ 位置在 WNT4 组上方（该组柱高 = 0 ⇒ 无遮挡），且不压 "not in the OneK1K panel" 注记。
    axf.legend(fontsize=5.2, loc="upper right",
               bbox_to_anchor=(1.0, float(_p.get("f_leg_anchor_y", 0.93))),
               ncol=1, handletextpad=0.35, labelspacing=0.3)
    axf.text(2.0, 2.8, "not in the\nOneK1K panel", ha="center", va="bottom", fontsize=5.2,
             color=FS.C["region"], linespacing=1.3)
    FS.finalize(axf)

    fig.subplots_adjust(bottom=LAY["fig4"]["bottom"], top=LAY["fig4"]["top"])
    # ★★★ 2026-10-01（逐图**第九轮** · Fig4）：水平重排 —— **全部按绝对 mm 指定**。
    #   目标（用户原话 + AskUserQuestion 三项裁定）：
    #     ① (a) = **46 × 46 mm**（"稍微放大图A"；原版 42.376）
    #     ② (b) X 轴 = **66 mm**（"图BX轴太长请缩短"；原版 86.319）
    #     ③ (f) X 轴 = **62 mm**（"F图X轴同样请稍微缩短"；原版 ≈66.9）
    #     ④ 三排**右缘**统一到列 1 右缘 R（右边界齐平）
    #     ⑤ (a)+(色条) 与 (c) 水平左移，使 **A/C/E 的 y 刻度文字左缘共线**
    #   为什么用绝对 mm 而不是 GridSpec 列宽：
    #     `_mm2x/_mm2y` 以**名义幅宽/幅高**为分母 ⇒ 落盘时该 axes 就是这么多 mm，
    #     与 `fit_save` 的 scale 无关（tight 裁剪不缩放内容，只裁边）。
    #     ⇒ 三个"用户点名的 mm 值"第一次就是精确的，不必与标定耦合计较。
    _pb = axb.get_position()                      # 行 0 未被 aspect 影响 ⇒ 行 0 格子真值
    _pd = axd.get_position()                      # 行 1 的 GridSpec 格（与 (c) 同行、同纵向范围）
    _pe = axe.get_position(); _pf = axf.get_position()
    _row0_y0, _row0_y1 = _y2mm(_pb.y0), _y2mm(_pb.y1)
    _row1_y0, _row1_y1 = _y2mm(_pd.y0), _y2mm(_pd.y1)
    _row2_y0, _row2_y1 = _y2mm(_pe.y0), _y2mm(_pe.y1)
    _c0x0, _c0x1 = _x2mm(_pe.x0), _x2mm(_pe.x1)   # 列 0（= 原版 (e) 那列）
    _c1x0, _Rx = _x2mm(_pf.x0), _x2mm(_pf.x1)     # 列 1 右缘 = 全图右界锚
    _a_mm = float(_p["a_mm"]); _bw = float(_p["b_mm"])
    _dw = float(_p["d_mm"]); _fw = float(_p["f_mm"])
    _cb_pad = float(_p["cb_pad_mm"]); _cb_w = float(_p["cb_w_mm"])
    _row0_h = _row0_y1 - _row0_y0
    # ★★★ 第十轮：**列 1 左缘 = 三图共用的 Y 轴位置**（(b) 的宽度定义它）。
    _L1 = _Rx - _bw
    # (a) 先放在列 0 左缘；x0 的真正取值由下方「刻度文字左缘对齐」统一平移决定。
    _ax0 = _c0x0
    axa.set_position([_mm2x(_ax0), _mm2y(_row0_y0), _mm2x(_a_mm), _mm2y(_row0_h)])
    # (a) 的**独立色条**：紧贴 (a) 右缘。★ 必须用独立 cax（`fig.colorbar(im, ax=axa)` 会把
    #   父轴宽度压成 1−fraction−pad = 0.924 倍 ⇒ (a) 永远长不到目标边长）。
    cax.set_position([_mm2x(_ax0 + _a_mm + _cb_pad), _mm2y(_row0_y0),
                      _mm2x(_cb_w), _mm2y(_row0_h)])
    # (b) 宽 66 mm、**左缘 = 列 1 左缘 `_L1`**、纵向与 (a) 完全同一范围（保证字母 A/B 同水平线）。
    axb.set_position([_mm2x(_L1), _mm2y(_row0_y0), _mm2x(_bw), _mm2y(_row0_h)])
    # (d) 宽 62 mm（= 与 (f) 一致，用户第十轮要求）、**左缘 = `_L1`**（与 (b)(f) 的 Y 轴齐平）。
    axd.set_position([_mm2x(_L1), _mm2y(_row1_y0),
                      _mm2x(_dw), _mm2y(_row1_y1 - _row1_y0)])
    # (f) 宽 66 mm、**左缘 = `_L1`**；纵向**底缘 = 行 2 底**（与 (e) 对齐）。
    # ★★★ 2026-10-02（第十二轮澄清）：**物理缩减 (f) 的 Y 轴长度** ——
    #   0→15 段比例尺不变，直接删掉 15–17.5 的高度 ⇒ 高度 × (f_ymax / f_ymax_prev)。
    #   ★ 保持 y0 不动 ⇒ (f) 的 X 轴仍与 (e) 同高（这条硬要求不变）；
    #     代价是 (f) **顶**低于 (e) 顶（预期、已向用户披露）。
    _f_ratio = (float(_p.get("f_ymax", 15.0))
                / float(_p.get("f_ymax_prev", 17.5)))
    axf.set_position([_mm2x(_L1), _mm2y(_row2_y0),
                      _mm2x(_fw), _mm2y((_row2_y1 - _row2_y0) * _f_ratio)])
    ROWSPAN_RESID["fig4"] = dict(
        fx0=_c0x0, fx1=_Rx, gap=_c1x0 - _c0x1,
        a_mm=_a_mm, b_mm=_bw, d_mm=_dw, f_mm=_fw, col1_left=_L1,
        a0=_ax0, b0=_L1, d0=_L1, f0=_L1, row0_h=_row0_h,
        b_right=_L1 + _bw, d_right=_L1 + _dw, f_right=_L1 + _fw,
        rows=[(axa, axb), (axc, axd), (axe, axf)])
    # ★★ 全局标准（2026-10-01）：「各行面板组居中于整幅」——
    #   ★ Fig4 例外：三排共用同一 GridSpec 列 ⇒ **必须整块一起平移**（单一 group），
    #     否则逐排居中会把三排重新推开。本图 group 内没有其它面板 ⇒ 实际位移 ≈ 0。
    CENTER_RESID["fig4"] = FS.center_row_groups(
        fig, [[axa, cb.ax, axb, axc, axd, axe, axf]])
    # ★★★ 2026-10-01（逐图第九轮 · Fig4）用户原话「使图A、C、E的Y轴图标最左端对齐」，
    #   口径经裁定 = **Y 轴刻度文字左缘**。做法：以三者中最靠左者（原版的 (e)，其
    #   "(GTEx v8 whole blood)" 是最宽的标签）为基准，把 (a)+色条 与 (c) 整体水平左移。
    #   ⚠ 代价：三排**轴线**从此按各自标签宽度错开（x0_i = L + w_i），这是用户明确选择的；
    #     它是「各行组居中」的显式例外（同 `align_axes_left` 的地位，见 `_s168` 的 MARGIN_BASE）。
    ALIGN_RESID["fig4"] = FS.align_ytick_label_left(
        [([axa, cb.ax], axa), ([axc], axc), ([axe], axe)])
    # ★★★ 第九轮追补（2026-10-02）：上一步对齐的是**文字逻辑包围盒**左缘；用户原话要的是
    #   「最左端」= 肉眼可见的**墨迹**左缘。二者相差 = 各排首字符的**左侧边距**：
    #     A 排首字符 `r`（bearing 小） / C 排 `0`、`1` / E 排 `(`（bearing 最大）。
    #   成品 PNG 实测（`scripts/_s240_fig4_r9_leftedge_verify.py`）：
    #     A = 8.805 | C = 9.059 | E = 9.821 mm ⇒ 极差 1.016 mm（≈0.77 个字高，肉眼可辨）。
    #   ⇒ 令 C、E 两组各再左移 (bearing_i − bearing_A)（值取自 `LAY["fig4"]["ink_dx_mm"]`，
    #     是**实测常量** ⇒ 文字变更须重测，见 LAY 内注释）。
    #   ★ 只动 (c)(e) 两组；(b)(d)(f) 不在组内 ⇒ 三排右缘齐平(181.128)不受影响。
    #   ★ 必须在 `align_letter_left_edges()` **之前**做：字母是 transAxes 文字，会随面板一起搬，
    #     所以字母的共线要对齐**最终**的 axes 位置。
    _ink = _p.get("ink_dx_mm") or {}
    for _key, _seq in (("c", (axc,)), ("e", (axe,))):
        _d_mm = float(_ink.get(_key, 0.0))
        if abs(_d_mm) > 1e-9:
            _dfr = _mm2x(_d_mm)
            for _ax in _seq:
                _pp = _ax.get_position()
                _ax.set_position([_pp.x0 + _dfr, _pp.y0, _pp.width, _pp.height])
        INK_DX_RESID.setdefault("fig4", {})[_key] = _d_mm
    # ★★★ 2026-10-02（第十一轮 · 用户「图A往右移动一点，与图C居中」）：
    #   把 (a) 的**绘图区中心**搬到 (c) 的绘图区中心（色条 `cb.ax` 随行，保持 `cb_pad_mm` 净距）。
    #   ★ 必须在 `align_ytick_label_left()` **之后**：该函数会平移 (a)/(c)/(e) 去对齐刻度文字，
    #     (c) 的最终 x 只有在这之后才确定（第十轮实测 (c) 中心 = 55.851 mm）。
    #   ★★ 与「A/C/E 字母左缘共线」的冲突与化解：字母是 transAxes 文字、会随面板搬走，而
    #     `align_letter_left_edges()` 只会把三者统一到 **min(各标准探针 x)** ⇒ (a) 右移
    #     13.10 mm 后该 min 会被顶到 A 的新位置（= 整列字母跟着右移，与用户裁定
    #     「**字母A不动**，只有图动」相悖）。
    #     ⇒ 先在**平移前**算出「不平移时会被选中的公共左缘 `_chosen_ace_px`」，
    #       平移后用 `x_px=_chosen_ace_px` 显式传回 ⇒ **A/C/E 三个字母纹丝不动**。
    #   ★ 已向用户披露的副作用：A 的 7 个 rs* y 刻度标签**属于绘图区**、随面板右移，
    #     故 A 排刻度文字不再与 C/E 的 9.063 mm 共线 —— 这是「只有图动」的直接结果。
    _plants_ace = [(axa, "a"), (axc, "c"), (axe, "e")]
    fig.canvas.draw()
    _fd_mm = float(fig.dpi) / 25.4
    _pwx = float(fig.get_size_inches()[0])
    _probe_ace = []
    for _ax, _lb in _plants_ace:
        _bb = _ax.get_position().bounds
        _ddx, _ = FS.letter_offset_for(_ax, size=8, label=_lb, anchor="word")
        _probe_ace.append(_bb[0] * _pwx * float(fig.dpi)
                          + _ddx * _bb[2] * _pwx * float(fig.dpi))
    _chosen_ace_px = min(_probe_ace)
    _pa = axa.get_position(); _pcx = axc.get_position()
    _a_cen_before = _x2mm(_pa.x0) + _a_mm / 2.0
    _c_cen = (_x2mm(_pcx.x0) + _x2mm(_pcx.x1)) / 2.0
    _d_center_mm = _c_cen - _a_cen_before
    _dfr = _mm2x(_d_center_mm)
    for _ax in (axa, cb.ax):
        _pp = _ax.get_position()
        _ax.set_position([_pp.x0 + _dfr, _pp.y0, _pp.width, _pp.height])
    ACENTER_RESID["fig4"] = dict(
        a_cen_before_mm=_a_cen_before, c_cen_mm=_c_cen, shift_mm=_d_center_mm,
        a_cen_after_mm=_x2mm(axa.get_position().x0) + _a_mm / 2.0,
        a0_after_mm=_x2mm(axa.get_position().x0),
        a_right_after_mm=_x2mm(axa.get_position().x0) + _a_mm,
        cb0_after_mm=_x2mm(cb.ax.get_position().x0),
        chosen_ace_px=_chosen_ace_px / _fd_mm)
    # ★★★ 2026-10-02（第十二轮 · 用户「图E的图再往左移动，字母E不动，图往左移动5mm」）：
    #   把 (e) 的**绘图区**整体左移 `e_shift_mm`（= −5.0 mm）。
    #   ★ **字母 E 不动** —— 复用第十一轮的「钉住对齐基准 + 移动被对齐对象」手法：
    #     `_chosen_ace_px`（A/C/E 三个标准探针的公共左缘）已在**本轮任何平移之前**算出，
    #     紧接着的 `align_letter_left_edges(_plants_ace, x_px=_chosen_ace_px)` 会把三个字母
    #     **显式钉回**该值 ⇒ (e) 面板左移、而其字母 E 绝对位置不变（A、C 同样不动）。
    #   ★ 副作用（与第十一轮 A 排完全同口径，已按用户口径接受）：(e) 的 y 刻度文字**属于
    #     绘图区**、随面板一起左移 5 mm（`_s240` 实测 9.059 → 约 4.059 mm）⇒
    #     **E 排刻度文字不再与 C 排共线**。判据同步更新（见 `_s240`）。
    #   ★ 只平移 x ⇒ 不影响 `set_row_ink_gaps`（该函数只做纵向平移），也不改变幅面左界
    #     （左界仍由 A 字母墨迹 0.508 mm 决定，E 刻度文字移到 4.059 mm 仍在其右）。
    _e_shift_mm = float(_p.get("e_shift_mm", 0.0))
    if abs(_e_shift_mm) > 1e-9:
        _efr = _mm2x(_e_shift_mm)
        _ep = axe.get_position()
        axe.set_position([_ep.x0 + _efr, _ep.y0, _ep.width, _ep.height])
    ESHIFT_RESID["fig4"] = dict(
        shift_mm=_e_shift_mm,
        e0_after_mm=_x2mm(axe.get_position().x0),
        e_right_after_mm=_x2mm(axe.get_position().x0) + _x2mm(axe.get_position().width),
        e_letter_pinned_px=_chosen_ace_px / _fd_mm)
    # ★ (f) 本轮 y 轴内缩 + 图例下移的实测量（同一 draw 周期内测，避免多余渲染）。
    fig.canvas.draw()
    _f_ax = axf.get_window_extent(fig.canvas.get_renderer())
    _f_leg = axf.get_legend()
    if _f_leg is not None:
        _lb = _f_leg.get_window_extent(fig.canvas.get_renderer())
        _inv = axf.transData.inverted()
        _p0 = _inv.transform((_lb.x0, _lb.y0)); _p1 = _inv.transform((_lb.x1, _lb.y1))
        FYLIM_RESID["fig4"] = dict(
            ymax=_fymax, ymax_prev=float(_p.get("f_ymax_prev", 17.5)),
            yticks=list(_fticks), max_bar=12.0,
            row2_h_mm=_row2_y1 - _row2_y0,
            f_h_mm=_y2mm(axf.get_position().height),
            f_h_ratio=(_y2mm(axf.get_position().height) / (_row2_y1 - _row2_y0)),
            f_top_mm=_y2mm(axf.get_position().y0) + _y2mm(axf.get_position().height),
            e_top_mm=_row2_y1,
            leg_frac_top=(_lb.y1 - _f_ax.y0) / _f_ax.height,
            leg_frac_bot=(_lb.y0 - _f_ax.y0) / _f_ax.height,
            leg_data_top=_p1[1], leg_data_bot=_p0[1],
            leg_inside=(_p0[1] >= 0.0 and _p1[1] <= _fymax))
    # ★★ 全局标准（2026-10-01）：面板字母距离 = Fig2 标准（Δx 6.41 / Δy 2.42 mm）。
    # ★★★ 第九轮新增（用户「三个字母也要最左端对齐」）：A/C/E 三个字母**左缘共线**。
    #   取三者的标准锚点里**最靠左**的那条（= 最宽的那个字母的标准位置）⇒ 至少对最宽的字母
    #   仍精确满足 6.41 mm，其余只会离标签更远 ⇒ 不会因对齐而压字。
    #   ★ 字母可能越过画布左缘（本图原版即如此），统一由 tight 裁剪兜住；
    #     故这里**不做** `_nudge_letter_clear`（对齐后微移会破坏共线）。
    #   ★★★ 第十一轮：`x_px=` 显式钉住公共左缘 ⇒ (a) 右移与字母列**解耦**（见上）。
    _offs = FS.align_letter_left_edges(_plants_ace, x_px=_chosen_ace_px)
    # ★★ 第 7 轮起的既有硬要求「A/B、C/D、E/F **同水平线**」：B/D/F 直接**沿用同行左侧
    #   字母的 dy**。为什么必须显式沿用：标准 dy 是「字母下缘 = **本面板** Y 轴最上面那个词
    #   的上缘 + 2.42」，而各面板"最上面那个词"在面板内的高度位置不同（(a) 的 rs24732930
    #   在 92.9% 高处、(d) 的 WNT4 在 78.6%、(e) 的 CDC42 在 78.6%）⇒ 各用标准 dy
    #   实测会差 A/B 0.41、C/D 3.00、E/F −1.78 mm。
    #   ★ 它们各自的 dx：第十轮之前用各自的**标准锚点**（右列不参与共线）；第十一轮起统一
    #     改用 corner 锚（见下）⇒ 右列三字母**左缘共线**。★ 行 2 自第十二轮起 (f) 比 (e) 矮，
    #     「同一 dy 即同一高度」**只对行 0/1 成立**（A/B、C/D 同高）；行 2 的 F 需另行反解（见下 ③）。
    #   ★★★ 2026-10-02（第十一轮 · 用户「D、F都在图的左上角。B、D、F对齐」）：
    #     ① **三者统一改用 `anchor="corner"`**（锚 = 本面板绘图区左上角；Δx 仍 = 全局标准
    #        6.41 mm）。因 (b)(d)(f) 左缘已同为 `_L1` ⇒ 三字母**右缘**天然全等（108.703）。
    #        ★ 这是对全局标准「有数值刻度就用 `anchor="word"`」的**用户点名例外**
    #          （(d) 的标签是两行基因名、(f) 是 `0.0…17.5` 数值刻度，按标准本应走 "word"）。
    #     ② **左缘共线**：corner 锚只保证右缘等距，而字形宽度不同 —— `D` 宽 2.032 mm、
    #        `F` 宽 1.693 mm ⇒ 首轮实测左缘 106.671 / 106.671 / **107.010**（差 0.339 mm）。
    #        用户要的「对齐」与第九轮 A/C/E 的「最左端对齐」同口径 ⇒ 统一到**最靠左**者
    #        （= B/D 的 106.671 mm），F 再左移 0.339 mm。★ 与 A/C/E 的处理完全同口径。
    #     ③ 纵向：A/B、C/D 沿用同行左侧字母的 dy ⇒ 共线。
    #   ★★★ 2026-10-02（第十二轮澄清 · 用户「**F字母与E字母顶端对齐**」）：
    #     (f) 面板本轮物理变矮（38.536 → 33.030 mm）⇒ 若仍沿用 E 的 **dy 分数**，F 会随
    #     面板顶下降 4.95 mm，与 E 错开。用户要求两者**顶端对齐** ⇒ **F 的 dy 改由绝对 mm 反解**：
    #         dy_f = (E 字母下缘的绝对 mm − (f) 的 y0) / (f) 的高度
    #     因 E/F 同为 8 pt 大写字母（字形上下缘完全相同）⇒ **下缘对齐 ⇔ 顶端对齐**（差 < 0.01 mm）。
    #     ★ 结果 dy_f = **1.04895 > 1** ⇒ 字母落在 (f) 绘图区**上沿之外**约 1.6 mm（浮在
    #       第 2 排与第 1 排之间的行带里）—— 这是「(f) 变矮 + 字母共线」的必然取舍，已向用户披露。
    #     ★ 不改动幅面/行带：(f) 的 tightbbox 顶 = F 字母顶 44.98 mm < (e) 的顶脊 46.21 mm
    #       ⇒ 第 2 排墨迹顶界仍由 (e) 决定 ⇒ `set_row_ink_gaps` 与幅高不变。
    _e_letter_bot_mm = (_y2mm(axe.get_position().y0)
                        + _offs["E"][1] * _y2mm(axe.get_position().height))
    _fpos = axf.get_position()
    _dy_f = (_e_letter_bot_mm - _y2mm(_fpos.y0)) / _y2mm(_fpos.height)
    _plants_bdf = [(axb, "b", _offs["A"][1]),
                   (axd, "d", _offs["C"][1]),
                   (axf, "f", _dy_f)]
    _bdf = []
    for _ax, _lb, _dy in _plants_bdf:
        _dx, _ = FS.letter_offset_for(_ax, label=_lb, anchor="corner")
        _bp = _ax.get_position()
        _bdf.append(dict(ax=_ax, lb=_lb, dx=_dx, dy=_dy, w=_bp.width,
                         left_std=_x2mm(_bp.x0 + _dx * _bp.width)))
    _bdf_tgt = min(_r["left_std"] for _r in _bdf)          # 最靠左者 = 公共左缘 (mm)
    for _r in _bdf:
        _r["dx_final"] = _r["dx"] + _mm2x(_bdf_tgt - _r["left_std"]) / _r["w"]
        FS.panel_label(_r["ax"], _r["lb"], dx=_r["dx_final"], dy=_r["dy"])
    LETTER_ANCHOR_OVR["fig4"] = dict(
        anchor="corner", common_left_mm=_bdf_tgt,
        per_letter={_r["lb"].upper(): dict(
            dx_std=_r["dx"], dx_final=_r["dx_final"],
            left_std_mm=_r["left_std"], left_final_mm=_bdf_tgt) for _r in _bdf},
        f_top_align=dict(
            e_letter_bot_mm=_e_letter_bot_mm, dy_f=_dy_f,
            dy_f_std=_offs["E"][1], dy_exceeds_panel=(_dy_f > 1.0)))
    # ★★ 全局标准（2026-10-01）：行带可见空白 = Fig2 标准 4.74 mm（逐边界测量 + 平移）。
    #   ★ 本图**改前两处行带都是负值**（−1.79 / +1.21 mm，即行1 与行2 的墨迹已经相碰）
    #     ⇒ 必须先把间隙撑开，故本轮的幅高会上升，需靠压 `LAY["fig4"]["h"]` 找回。
    #   ★★ 2026-10-01（逐图第八轮 · Fig4，用户裁定「守住170mm，放大图A，缩小图之间的距离」）：
    #     为给放大的 (a) 腾出高度，行带由全局标准 4.74 → `LAY["fig4"]["row_gap"]` = **3.00 mm**
    #     （每边界省 1.74 mm，两边界共省 3.48 mm 幅高）。改用**参数读取**而非硬编码，
    #     以后调 `row_gap` 一处即可，不会出现「改 LAY 却不生效」的静默不一致。
    ROW_GAP_RESID["fig4"] = FS.set_row_ink_gaps(
        fig, [[axa, cb.ax, axb], [axc, axd], [axe, axf]],
        [_p["row_gap"], _p["row_gap"]])
    if SHOW_NOTES:
        fig.text(0.5, 0.005,
                 "Panel (e), cis-SMR + HEIDI: each point is one cell type; bars are 95% CIs. CDC42 (all positive) and\n"
                 "LINC00339 (all negative) are significant, with opposite signs; WNT4 in GTEx v8 whole blood is a\n"
                 "low-power layer (median TPM 0.228; no variant below 5e-8) and yields no SMR evidence - 'not\n"
                 "detected', not 'refuted'. Panel (f): both instruments regulate CDC42 AND LINC00339 in the same\n"
                 "cell types, and neither regulates WNT4 in any layer -> the instruments are REGION tags, so the MR\n"
                 "signal cannot be assigned to a single gene by eQTL evidence alone. WNT4 is absent from the\n"
                 "OneK1K panel (0 of 14 cell types), which is a lack of evidence, not a negative result.",
                 ha="center", va="bottom", fontsize=5.2, color="#555555")
    return fig


# ============================================================ Fig5
def fig5(scale=1.0):
    ot = load("41d_opentargets_diseases_full.csv")
    dom = load("42b_phewas_safety_domains.csv")
    up = load("42_phewas_safety_assessment.csv")
    pw = load("50c_task35_phewas_power_by_domain.csv")     # task 3.5 power by domain
    pc = load("50d_task35_power_class_by_domain_968.csv")
    cd = load("44c_phewas_3p2_full_domains.csv")           # task 3.2 conditional PheWAS

    _p = LAY["fig5"]
    fig = plt.figure(figsize=(7.2 * scale, _p["h"] * scale))
    # 顶行保持 P0 版式（(a) | (b)/(c) 叠放）：wspace 加大 + left 加大，让 (b) 的长 y 轴标签
    # 落在列间距内、(a) 的标签落在左边距内，从而 tight-bbox ≈ 画布（否则标签外溢会把 scale
    # 压小、字号相对放大 ~18%，与 Fig1-4 不一致）。新增行放 (d)(e)(f)。
    # 下排由「3 列并排」改为「2x2 嵌套」：(d)(e) 占上排两列，(f) 独占整行。
    # 原因：3 列并排时 (f) 面板仅约 27 mm 宽却要容纳 5 个域 x 3 个系列的 x 轴标签，
    # 无论旋转与否都会重叠成不可读色块；独占整行（148 mm）后 5 组标签可横排完整显示。
    outer = gridspec.GridSpec(2, 1, height_ratios=_p["outer_hr"], hspace=_p["outer_hspace"])
    top = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=outer[0], wspace=1.05,
                                           width_ratios=[1.0, 1.14])
    axa = fig.add_subplot(top[0])
    gsr = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=top[1],
                                           height_ratios=_p["gsr_hr"], hspace=_p["gsr_hspace"])
    axb = fig.add_subplot(gsr[0]); axc = fig.add_subplot(gsr[1])
    bot = gridspec.GridSpecFromSubplotSpec(2, 2, subplot_spec=outer[1],
                                           hspace=_p["bot_hspace"], wspace=0.95)
    axd = fig.add_subplot(bot[0, 0]); axe = fig.add_subplot(bot[0, 1])
    axf = fig.add_subplot(bot[1, :])

    # ---------- (a) Open Targets association spectrum ----------
    ot = sorted(ot, key=lambda r: -num(r["overall_score"]))[:18][::-1]
    y = np.arange(len(ot))
    axa.barh(y, [num(r["overall_score"]) for r in ot],
             color=[GENE_C.get(r["gene"], "#999999") for r in ot],
             height=0.74, edgecolor="#4D4D4D", lw=0.3)
    labs = []
    for r in ot:
        n = r["disease_name"]
        labs.append(n if len(n) <= 32 else n[:29] + "...")
    axa.set_yticks(y); axa.set_yticklabels(labs, fontsize=5.6)
    axa.set_xlabel("Open Targets association score"); axa.set_xlim(0, 0.95)
    axa.spines["left"].set_visible(False); axa.tick_params(axis="y", length=0)
    # ★ 全局「点的大小与风格」标准（_fig_style.py）：面积 ×2 + 同色系深 2 度描边。
    #   ★ 图例里的 `*` 是**显式 Line2D 代理**（ms 即绝对尺寸）⇒ 保持 4.2 不变，
    #     即"图例标记的绝对尺寸不变"这条口径。
    _msstar = FS.pt_up(4.2)
    for i, r in enumerate(ot):
        ta = r["therapeutic_areas"].lower()
        if any(k in ta for k in ["reproductive", "hematologic", "immune"]):
            axa.plot(0.905, i, marker="*", ms=_msstar, color=FS.C["region"],
                     mec=FS.edge_color(FS.C["region"]), mew=FS.mew_for(_msstar),
                     clip_on=False)
    hs = [plt.Rectangle((0, 0), 1, 1, fc=GENE_C[g]) for g in ["WNT4", "CDC42", "LINC00339"]]
    hs.append(plt.Line2D([], [], marker="*", ls="none", ms=4.2, color=FS.C["region"]))
    # ★★★ 2026-10-02（Fig5 第一轮）用户：「图A的图标往下和左移动，图标最右端不要超过图的最右端。」
    #   原状：`loc="lower left", bbox_to_anchor=(0.0, 1.015)` ⇒ 图例**左对齐**、宽 50.154 mm，
    #         **右缘越出绘图区右缘 8.601 mm**（实测），且整块浮在轴框之上 8.0 mm。
    #   修法：① 改 `loc="upper right"` ⇒ **右缘恒等于绘图区右缘**（约束由构造保证，与宽度无关）；
    #         ② 字号 5.7 → 5.2（与其余面板一致，仍是 P0 下限）、`columnspacing` 0.9→0.55、
    #            `handletextpad` 0.4→0.30、`handlelength` 2.0→1.4、`borderpad` 0.22
    #            ⇒ 整体变窄，右侧不再外溢、左侧少侵入；
    #         ③ `a_leg_anchor_y` 由 1.12 降到 1.095 ⇒ 整块**下移**（底缘贴近轴顶）。
    axa.legend(hs, ["WNT4", "CDC42", "LINC00339", "reproductive / hematologic / immune"],
               loc="upper right", bbox_to_anchor=(1.0, float(_p.get("a_leg_anchor_y", 1.095))),
               ncol=2, fontsize=5.2, columnspacing=0.55, handletextpad=0.30,
               handlelength=1.4, borderpad=0.22)
    FS.finalize(axa)

    # ---------- (b) FinnGen R12 PheWAS: 12 risk-increasing signals by lineage ----------
    LG = OrderedDict([
        ("A. 盆底支持结构/结缔组织", (FS.C["sig"], "A pelvic-floor")),
        ("B. 胎盘附着/剥离异常与产科出血", (FS.C["lineage_b"], "B obstetric")),
        ("C. 早产", (FS.C["lineage_c"], "C preterm")),
        ("D. 关节软骨退变", (FS.C["gene_ym"], "D knee OA")),
        # ★★ 2026-10-02（指令 S34）：E / F 由硬编码灰改为语义键。
        #    用户原话「图5B灰色的点线换成彩色的」。
        #    GB12 侧 = #9ED17B / #BAD2E1（见 _gb12_palette.GB12_C 的 S34 块）；
        #    npg 侧仍是原灰 #8C8C8C / #C7C7C7 ⇒ FIGPAL=npg 逐字节不变。
        #    ★ 已知代价：六色同框在 GB12 内不可达（28 组合 0 通过），本轮取主轴 D 最优
        #      minD 9.16 -> 15.55（过 D 轴）、minT 12.31 -> 4.46（tritan 轴下降）。
        # ★★ 2026-10-02（M5 落地，第二轮）：图例由未完成状态词 "E pending" 改为真实谱系名。
        #    M5 原话「补出 E 组实际谱系名或合并删除」；正文图注已写「E 软组织残留异物」。
        #    ★ 本字典的 **key 必须逐字符等于 42_phewas_safety_assessment.csv 的
        #      lineage_group 值**（`LG.get(r["lineage_group"])` 的 join key）——
        #      该值由 s36 的 LINEAGE 生成，而 s36 依赖外置盘 E:\（本轮不在线）⇒
        #      **key 保持原样不动**，本轮只改**显示标签**。key 的正名另案处理。
        ("E. 可疑/待复核", (FS.C["lineage_e"], "E foreign body")),
        ("F. 非疾病生育表型(不计为风险)", (FS.C["lineage_f"], "F non-disease")),
    ])
    inc = [r for r in up if r["risk_if_cdc42_inhibited"] == "increased"]
    inc = sorted(inc, key=lambda r: -(-num(r["beta_MR"])))
    yy = np.arange(len(inc))[::-1]
    # ★ 全局「点的大小与风格」标准（_fig_style.py）：面积 ×2 + 同色系深 2 度描边
    _msb = FS.pt_up(3.2)
    for i, r in enumerate(inc):
        col = LG.get(r["lineage_group"], ("#666666", ""))[0]
        m = -num(r["beta_MR"]); e = 1.96 * num(r["se_MR"])
        axb.errorbar([m], [yy[i]], xerr=[[e], [e]], fmt="o", ms=_msb,
                     mec=FS.edge_color(col), mew=FS.mew_for(_msb),
                     color=col, elinewidth=0.6, capsize=1.3)
    hi = max(-num(r["beta_MR"]) + 1.96 * num(r["se_MR"]) for r in inc)
    axb.axvline(0, color="#4D4D4D", lw=0.6)
    axb.set_yticks(yy)
    axb.set_yticklabels(["%s | %s" % (r["phenotype"][:30], r["cell_type"]) for r in inc],
                        fontsize=5.4)
    # P0：标题改为「区域相关风险升高信号」；归因区域而非 CDC42，并标注工具 LD
    axb.set_title("Region-related risk-increasing signals", fontsize=6.8, pad=3.0)
    # 2026-09-28 裁定：图内文字与图注统一为 β 符号体系（原为 ASCII 拼写 "beta_MR"）
    axb.set_xlabel("-β_MR   (positive = higher risk; 95% CI)")
    axb.set_xlim(-0.46, hi * 1.08)
    axb.spines["left"].set_visible(False); axb.tick_params(axis="y", length=0)
    # P0：归因区域而非 CDC42。板内说明放**左上空白带**（该带 12 行数据全部 x>0），
    # 图例放**左下空白带**，两者纵向分居，彻底消除「图例压注记 / 压 y 轴标签」。
    # 图例文字用缩写（全称见 Fig5 图注），否则图例宽过左侧空白带会骑到数据点上。
    axb.text(0.015, 0.985,
             "region-attributed,\nnot CDC42 (see note)",
             transform=axb.transAxes, ha="left", va="top",
             fontsize=5.2, color=FS.C["region"], linespacing=1.35)
    hh = [plt.Line2D([], [], marker="o", ls="none", ms=3.2, color=c) for c, _t in LG.values()]
    axb.legend(hh, [t for _c, t in LG.values()], loc="lower left", fontsize=5.2,
               handletextpad=0.3, labelspacing=0.25, borderpad=0.25)
    FS.finalize(axb)

    # ---------- (c) FinnGen R12 PheWAS: domain rollup ----------
    dm = {r["domain"]: r for r in dom}
    order = ["免疫/血液/肿瘤", "女性生殖/不孕", "妊娠/分娩/产褥", "肌肉骨骼/结缔组织", "其他"]
    en = {"免疫/血液/肿瘤": "immune / blood / tumour", "女性生殖/不孕": "female reproductive",
          "妊娠/分娩/产褥": "pregnancy / puerperium", "肌肉骨骼/结缔组织": "musculoskeletal",
          "其他": "all other"}
    yc = np.arange(len(order))[::-1]
    axc.barh(yc, [int(dm[d]["n_risk_up"]) for d in order], color=FS.C["risk_up"], height=0.55,
             edgecolor="#4D4D4D", lw=0.3, label="risk-increasing")
    axc.barh(yc, [-int(dm[d]["n_risk_down"]) for d in order], color=FS.C["risk_down"], height=0.55,
             edgecolor="#4D4D4D", lw=0.3, label="risk-decreasing")
    axc.axvline(0, color="#4D4D4D", lw=0.6)
    axc.set_yticks(yc)
    lab_c = []
    for d in order:
        if d == "免疫/血液/肿瘤":
            lab_c.append("immune / blood / tumour (968; 0 up)")
        else:
            lab_c.append("%s (%s)" % (en[d], dm[d]["n_tests"]))
    axc.set_yticklabels(lab_c, fontsize=5.2)
    axc.set_xlabel("FDR < 0.05 signals (up + / down -)")
    axc.set_xlim(-34, 34)
    axc.spines["left"].set_visible(False); axc.tick_params(axis="y", length=0)
    axc.legend(fontsize=5.2, loc="lower left", handletextpad=0.4)
    # P0 修正（2026-10-01）：此处**原缺** panel_label —— Fig5 面板字母一直是 a,b,d,e,f，
    #   没有 "c"。P0 验收项「面板字母个数正确」实测 5/6 ⇒ 补上。
    FS.finalize(axc)

    # ---------- (d) statistical power (task 3.5): median minimum detectable OR by domain ----------
    dmP = {r["domain"]: r for r in pw}
    dlab = {"免疫/血液/肿瘤": "immune / blood\n/ tumour", "其他": "all other",
            "肌肉骨骼/结缔组织": "musculo-\nskeletal", "女性生殖/不孕": "female\nreprod.",
            "妊娠/分娩/产褥": "pregnancy"}
    dorder = ["免疫/血液/肿瘤", "其他", "肌肉骨骼/结缔组织", "女性生殖/不孕", "妊娠/分娩/产褥"]
    dmed = [num(dmP[d]["min_detectable_OR_median"]) for d in dorder]
    yq = np.arange(len(dorder))[::-1]
    # ★★★ 2026-10-02（Fig5 第二轮）用户：「图D中灰色的柱状图可以换成彩色的吗？」
    #   口径（AskUserQuestion 已确认）= **复用 (e) 的 5 级效力等级配色**：
    #     颜色 = 该域「中位最小可检出 OR」所落入的等级，与 (e) 面板的 5 个 class 同色同义。
    #     等级阈值与 (e) 的 `keys5` 逐字一致：<=1.2 / 1.2-1.5 / 1.5-2.0 / 2.0-3.0 / >3.0。
    #   实测落位：免疫/血液/肿瘤 3.15 -> >3.0（robust_strong #367DB0，与 (e) 的 500 根同色）；
    #             其他 2.22 / 肌肉骨骼 2.09 -> 2.0-3.0（suggest_strong #519D78）；
    #             女性生殖 1.74 / 妊娠 1.71 -> 1.5-2.0（amber #C4E9CA）。
    #   ★ 因免疫域柱改为 #367DB0，**两条阈值虚线原用的 risk_down(#367DB0) / region(#3D9F3C)
    #     会与柱撞色**，故一并改为中性灰 #4D4D4D（= 全图 spine/edge 色）；图注只写
    #     "dashed lines, OR = 1.5 and OR = 3.0"、**无颜色词** ⇒ 换色不需改图注。
    def _pw_class_color(v):
        if v <= 1.2:
            return FS.C["green_dark"]
        if v <= 1.5:
            return FS.C["cyan"]
        if v <= 2.0:
            return FS.C["amber"]
        if v <= 3.0:
            return FS.C["suggest_strong"]
        return FS.C["robust_strong"]

    dcol = [_pw_class_color(v) for v in dmed]
    axd.barh(yq, dmed, color=dcol, height=0.62, edgecolor="#4D4D4D", lw=0.3)
    for i, d in enumerate(dorder):
        axd.text(dmed[i] + 0.08, yq[i], "%.2f" % dmed[i], va="center", fontsize=5.2,
                 color=(FS.C["region"] if d == "免疫/血液/肿瘤" else "#222222"),
                 fontweight=("bold" if d == "免疫/血液/肿瘤" else "normal"))
    axd.axvline(1.5, color="#4D4D4D", ls="--", lw=0.6)
    axd.axvline(3.0, color="#4D4D4D", ls=":", lw=0.6)
    axd.set_ylim(-0.7, 5.4)
    axd.text(1.5, 4.9, "OR 1.5", fontsize=5.2, color="#4D4D4D", ha="center")
    axd.text(3.0, 4.9, "OR 3.0", fontsize=5.2, color="#4D4D4D", ha="center")
    axd.set_yticks(yq)
    # ★★★ 2026-10-02（Fig5 第一轮）用户：「图D的Y轴标签重叠了，可适当拉宽Y轴长度，不要使文字重叠。」
    #   原状：标签写成 `"%s\n(%s tests)"` ⇒ **3 行**（域名本身已折行）⇒ 单个标签高 6.047 mm，
    #         而相邻类别中心间距只有 26.056/6.1 = 4.270 mm ⇒ 实测 4 处相邻净空全为负
    #         （-0.768 / -0.768 / -1.776 / -0.758 mm）＝**肉眼可见的重叠**。
    #   修法：把 "(N tests)" **并到域名末行**（文字一字不改，只是不再单独起行）⇒ 每个标签恒为
    #         **2 行**（高 ≈4.031 mm），并把 `linespacing` 收到 1.15；配合本轮 scale 因字母对齐
    #         上抬而带来的轴高增长，相邻净空转正。轴高本身**不改**（避免牵动 (e)/(f) 行）。
    axd.set_yticklabels(["%s (%s tests)" % (dlab[d], "{:,}".format(int(dmP[d]["n_tests"])))
                         for d in dorder], fontsize=5.2,
                        linespacing=float(_p.get("d_lab_linespacing", 1.15)))
    axd.set_xlabel("median minimum detectable OR\n(80% power, domain Bonferroni)")
    axd.set_xlim(0, 3.95)
    axd.spines["left"].set_visible(False); axd.tick_params(axis="y", length=0)
    FS.finalize(axd)

    # ---------- (e) statistical power: class distribution of the 968 tests ----------
    cls = pc[0]
    keys5 = [("1_adequate_OR<=1.2", "<=1.2"), ("2_adequate_1.2<OR<=1.5", "1.2-1.5"),
             ("3_limited_1.5<OR<=2.0", "1.5-2.0"), ("4_inadequate_2.0<OR<=3.0", "2.0-3.0"),
             ("5_severe_OR>3.0", ">3.0")]
    cl5 = [FS.C["green_dark"], FS.C["cyan"], FS.C["amber"], FS.C["suggest_strong"], FS.C["robust_strong"]]
    vals5 = [int(cls[k]) for k, _l in keys5]
    xp5 = np.arange(len(keys5))
    axe.bar(xp5, vals5, color=cl5, width=0.66, edgecolor="#4D4D4D", lw=0.3)
    for i, v in enumerate(vals5):
        axe.text(xp5[i], v + 14, "%d" % v, ha="center", va="bottom", fontsize=5.2)
    axe.set_xticks(xp5)
    axe.set_xticklabels([l for _k, l in keys5], fontsize=5.2, rotation=25, ha="right",
                        rotation_mode="anchor")
    axe.set_xlabel("smallest detectable OR (80% power)")
    axe.set_ylabel("number of tests\n(of 968)")
    axe.set_ylim(0, 620)
    axe.text(0.03, 0.98, "the immune / blood /\ntumour domain", transform=axe.transAxes,
             ha="left", va="top", fontsize=5.2, color=FS.C["region"], linespacing=1.3)
    FS.finalize(axe)

    # ---------- (f) conditional PheWAS (task 3.2) ----------
    fam = {}
    for r in cd:
        fam.setdefault(r["domain"], {})[r["family"]] = int(r["n_sig"])
    dord2 = ["免疫/血液/肿瘤", "女性生殖/不孕", "妊娠/分娩/产褥", "肌肉骨骼/结缔组织", "其他"]
    dlab2 = ["immune / blood / tumour", "female reproductive\n/ infertility",
             "pregnancy / childbirth\n/ puerperium", "musculoskeletal\n/ connective tissue",
             "all other"]
    xq = np.arange(len(dord2)); wq = 0.26
    for k, (fk, lab, cl) in enumerate([("uncond", "unconditional", FS.C["risk_up"]),
                                       ("cond_eqtl_only", "conditioned on eQTL", FS.C["cyan"]),
                                       ("cond_both", "conditioned on both", "#9E9E9E")]):
        vs = [fam.get(d, {}).get(fk, 0) for d in dord2]
        axf.bar(xq + (k - 1) * wq, vs, width=wq, color=cl, edgecolor="#4D4D4D", lw=0.3, label=lab)
        for i, v in enumerate(vs):
            # 条件组在全部 5 个域恒为 0 -> 显式标 0，避免读者误判为「数据缺失」。
            # 注意：(f) 的 x 轴刻度为 direction="in"（长 2.2 pt ≈ 0.78 mm），零值标注
            # 若放在 y=0.9 会被刻度线穿透 -> 零值单独抬高到 y=2.6。
            axf.text(xq[i] + (k - 1) * wq, v + (0.9 if v > 0 else 2.6), "%d" % v,
                     ha="center", va="bottom", fontsize=5.2,
                     color=cl if v > 0 else "#4D4D4D")
    axf.set_xticks(xq)
    axf.set_xticklabels(dlab2, fontsize=5.2, linespacing=1.35)
    axf.set_ylabel("FDR < 0.05 signals")
    axf.set_xlabel("disease domain (Bonferroni within domain)")
    axf.set_ylim(0, 46)
    # ★★★ 2026-10-02（Fig5 第一轮）用户：「图F右上角的图标适当往下移动。」
    #   原状：`loc="upper right"` 无 bbox ⇒ 图例紧贴绘图区右上角（上缘距轴顶 0.917 mm、右缘距轴右 0.917 mm）。
    #   修法：`bbox_to_anchor=(1.0, f_leg_anchor_y=0.900)` ⇒ 整块下移约 0.10×轴高 ≈ 2.6 mm，右缘保持齐平。
    axf.legend(fontsize=5.2, loc="upper right",
               bbox_to_anchor=(1.0, float(_p.get("f_leg_anchor_y", 0.900))),
               ncol=1, handletextpad=0.35, labelspacing=0.28)
    axf.text(0.012, 0.99, "conditioning removes\nevery hit (52 -> 0)", transform=axf.transAxes,
             ha="left", va="top", fontsize=5.2, color=FS.C["region"], linespacing=1.3)
    FS.finalize(axf)

    # 为脚注预留底部空间（14 行脚注 + 各行 x 轴标签 -> bottom 提到 0.20）
    fig.subplots_adjust(bottom=_p["bottom"], top=_p["top"], left=_p["left"], right=_p["right"])
    # ★★ 全局标准（2026-10-01）：「各行面板组居中于整幅」（须在 subplots_adjust 之后）。
    #   本图是三段式：(a)|(b)/(c) 的顶行 → (d)(e) → (f) 整行；顶行里的 (a) 与右列 (b)(c)
    #   属同一"行带"，必须一起平移，否则 (c) 会被单独拉到中间（它按设计应在右列）。
    _G5 = [[axa, axb, axc], [axd, axe], [axf]]
    CENTER_RESID["fig5"] = FS.center_row_groups(fig, _G5)
    # ★★★ 2026-10-02（Fig5 第三轮）用户：「图F往右移动，使F的Y轴跟D的Y轴对齐。字母F不动」
    #   实测：(d) 绘图区左缘 = **24.040 mm**、(f) 绘图区左缘 = **18.280 mm** ⇒ (f) 需右移 **5.760 mm**。
    #   ★ 时机硬性（沿用 Fig4 第十一轮的先例）：必须在 `center_row_groups()` **之后**
    #     （居中会逐行独立平移，提前搬会被搬回去）、且在**字母锚定之前**
    #     （`panel_label_anchored` 读的是 axes 的最终位置）。
    #   ★ 这是「各行组居中」标准的**显式例外**：第 3 行（(f) 独占一行）的墨迹中心不再等于
    #     整幅中心 ⇒ `CENTER_RESID["fig5"]` 会保留该残差（登记，不视为缺陷）。
    #   ★ 副作用（预期，且视觉更整齐）：(f) 右缘 166.772 + 5.760 = **172.532 mm** = (e) 右缘
    #     ⇒ (f) 恰好与上一行 (d)+(e) 整块的外边界 `24.040 .. 172.532` 完全重合。
    #   ★ 「字母F不动」由既有机制自动保证：A/F 的 x 由「对齐到 D」（`letter_align_members`）
    #     钉到 **D 的文字左缘**，而 D 本身未移动 ⇒ F 的绝对 x 不变；
    #     字母 y 由本面板 y 轴最上面那个词决定，**水平平移不影响 y**。
    AXALIGN5_RESID["fig5"] = dict(shift_mm=list(FS.align_axes_left(axd, [axf], fig=fig,
                                                                  verbose=False)))
    # ★ 如实记录「各行组居中」被打破的量：第 3 行（(f)）墨迹中心相对整幅中心的偏移。
    fig.canvas.draw()
    _rrA = fig.canvas.get_renderer()
    _pmA = float(fig.dpi) / 25.4
    _fb = axf.get_tightbbox(_rrA)
    _gb = fig.get_tightbbox(_rrA)
    AXALIGN5_RESID["fig5"].update(
        d_x0_mm=axd.get_window_extent(_rrA).x0 / _pmA,
        f_x0_mm=axf.get_window_extent(_rrA).x0 / _pmA,
        f_x1_mm=axf.get_window_extent(_rrA).x1 / _pmA,
        e_x1_mm=axe.get_window_extent(_rrA).x1 / _pmA,
        row3_ink_center_mm=0.5 * (_fb.x0 + _fb.x1) / _pmA,
        # ★ 口径陷阱：`Figure.get_tightbbox()` 返回**英寸**，`Axes.get_tightbbox()` 返回**像素**
        #   ⇒ 必须 ×25.4 归一（曾误按 px/mm 除，得 0.590 的荒谬值）。
        fig_ink_center_mm=0.5 * (_gb.x0 + _gb.x1) * 25.4,
        row3_center_offset_mm=0.5 * ((_fb.x0 + _fb.x1) / _pmA - (_gb.x0 + _gb.x1) * 25.4),
        note="row3 不再居中 = 「各行组居中」标准的显式例外（用户点名要求 (f) 与 (d) 轴对齐）")
    # ★★ 全局标准（2026-10-01）：面板字母距离 = Fig2 标准（Δx 6.41 / Δy 2.42 mm）
    for _ax, _lb in ((axa, "a"), (axb, "b"), (axc, "c"), (axd, "d"), (axe, "e"), (axf, "f")):
        FS.panel_label_anchored(_ax, _lb)
    # ★★★ 2026-10-02（Fig5 第一轮）用户：「以字母D为标准，字母A、F和D纵向对齐。」
    #   为什么必须**事后平移**而不是改锚：`panel_label_anchored` 是「字母右缘 ← 本面板 Y 轴
    #   最上面那个词左缘 − Δx(6.41 mm)」⇒ 各面板最上面那个词的宽度不同，字母左缘自然错开
    #   （实测 A 0.000 / D 13.646 / F 17.841 mm）。要对齐「文字左缘」，只能在锚定之后按
    #   **实测 bbox** 做一次水平搬家（与 Fig4 第十一轮「统一到最靠左者」同一手法，只是基准换成 D）。
    #   ★ 位移换算：字母用 `transform=ax.transAxes` 定位，故 Δdx = Δpx / 本面板轴宽(px)。
    _ALIGN_TO = str(_p.get("letter_align_to", "d"))
    # ★★ 只动用户点名的两个字母（A、F）；B/C/E 保持各面板原位 —— 首版误写成「除基准外全搬」，
    #    实测把 B/C/E 也搬到了 D 上（shift −506.5 / −521.5 / −640.9 px），已修正。
    _ALIGN_MEMBERS = tuple(_p.get("letter_align_members", ("a", "f")))
    _ltx = {}
    for _ax, _lb in ((axa, "a"), (axb, "b"), (axc, "c"), (axd, "d"), (axe, "e"), (axf, "f")):
        for _t in _ax.texts:
            if _t.get_text() == _lb.upper() and abs(float(_t.get_fontsize()) - 8.0) < 1e-6:
                _ltx[_lb] = _t
    fig.canvas.draw()
    _rr = fig.canvas.get_renderer()
    if _ALIGN_TO in _ltx:
        _ref = _ltx[_ALIGN_TO].get_window_extent(_rr).x0
        _rec = {}
        for _k, _t in _ltx.items():
            _ax_k = _t.axes
            _b = _t.get_window_extent(_rr)
            _aw = _ax_k.get_window_extent(_rr).width
            _dx = _t.get_position()[0]
            _rec[_k.upper()] = dict(left_before_px=_b.x0, dx_before=_dx)
            if _k in _ALIGN_MEMBERS and _k != _ALIGN_TO:
                _shift_px = _ref - _b.x0
                _t.set_position((_dx + _shift_px / _aw, _t.get_position()[1]))
                _rec[_k.upper()]["shift_px"] = _shift_px
        fig.canvas.draw()
        _rr = fig.canvas.get_renderer()
        for _k, _t in _ltx.items():
            _rec[_k.upper()]["left_after_px"] = _t.get_window_extent(_rr).x0
        _moved = [_rec[k.upper()]["left_after_px"] for k in _ALIGN_MEMBERS if k.upper() in _rec]
        LETTER5_ALIGN["fig5"] = dict(
            reference=_ALIGN_TO.upper(), members=[m.upper() for m in _ALIGN_MEMBERS],
            per_letter=_rec,
            moved_spread_px=max(_moved) - min(_moved) if _moved else float("nan"),
            ref_left_after_px=_rec[_ALIGN_TO.upper()]["left_after_px"],
            others_untouched=[k for k in _rec if k not in
                              [m.upper() for m in _ALIGN_MEMBERS] + [_ALIGN_TO.upper()]])
    # ★★ 全局标准（2026-10-01）：行带可见空白 = Fig2 标准 4.74 mm（逐边界测量 + 平移）。
    #   ★ 本图改前行带为 +17.89 / −2.05 mm（上宽下重叠）⇒ 净效果是把行 2/3 **上提**，
    #     幅高会**变矮**，故本图无需额外压 `h`。
    ROW_GAP_RESID["fig5"] = FS.set_row_ink_gaps(fig, _G5, [4.74, 4.74])
    # ★★★ 2026-10-02（Fig5 第二轮）用户：「字母E往左上角移动5mm。」
    #   口径（AskUserQuestion 已确认）= 水平、垂直**各** 5 mm（负 dx = 左，正 dy = 上）。
    #   ★★ 必须放在 `set_row_ink_gaps()` **之后**：字母抬升会抬高本行（row1）的墨迹顶界，
    #      若在行带调整之前施加，行带机制会为维持 4.74 mm 净空而把整行下压（≈3.2 mm），
    #      字母的**页面绝对位移**被抵消到只剩 ≈1.4 mm。放在之后 ⇒ 位移恒为精确 5 mm。
    #   ★ 换算：字母用 `transAxes` 定位 ⇒ Δf = Δmm × (dpi/25.4) / 本面板该方向的像素尺寸。
    _SHIFT_MM = dict(_p.get("letter_shift_mm", {}))
    SHIFT5_RESID["fig5"] = {}
    if _SHIFT_MM and _ltx:
        fig.canvas.draw()
        _rr2 = fig.canvas.get_renderer()
        _pm2 = float(fig.dpi) / 25.4
        for _k, (_dxmm, _dymm) in _SHIFT_MM.items():
            _t = _ltx.get(_k)
            if _t is None:
                continue
            _b0 = _t.get_window_extent(_rr2)
            _abk = _t.axes.get_window_extent(_rr2)
            _fx, _fy = _t.get_position()
            _t.set_position((_fx + _dxmm * _pm2 / _abk.width,
                             _fy + _dymm * _pm2 / _abk.height))
            fig.canvas.draw()
            _b1 = _t.get_window_extent(_rr2)
            SHIFT5_RESID["fig5"][_k.upper()] = dict(
                dx_mm=_dxmm, dy_mm=_dymm,
                x0_before_mm=_b0.x0 / _pm2, y0_before_mm=_b0.y0 / _pm2,
                x0_after_mm=_b1.x0 / _pm2, y0_after_mm=_b1.y0 / _pm2,
                dx_real_mm=(_b1.x0 - _b0.x0) / _pm2,
                dy_real_mm=(_b1.y0 - _b0.y0) / _pm2)
    # ---- 本轮四处改动的**成品口径**留痕（最后统一 draw 一次，避免多余渲染）----
    fig.canvas.draw()
    _rr = fig.canvas.get_renderer()
    _pxmm = float(fig.dpi) / 25.4

    def _mmb(_bbox):
        return tuple(v / _pxmm for v in (_bbox.x0, _bbox.x1, _bbox.y0, _bbox.y1))

    _ab = axa.get_window_extent(_rr)
    _lg = axa.get_legend()
    if _lg is not None:
        _l = _lg.get_window_extent(_rr)
        A5LEG_RESID["fig5"] = dict(
            leg_w_mm=(_l.x1 - _l.x0) / _pxmm,
            leg_x0_mm=_l.x0 / _pxmm, leg_x1_mm=_l.x1 / _pxmm,
            ax_x0_mm=_ab.x0 / _pxmm, ax_x1_mm=_ab.x1 / _pxmm,
            right_over_mm=(_l.x1 - _ab.x1) / _pxmm,
            top_over_mm=(_l.y1 - _ab.y1) / _pxmm,
            gap_to_axes_top_mm=(_ab.y1 - _l.y0) / _pxmm)
    _dbb = [(t.get_text().replace("\n", "|"), _mmb(t.get_window_extent(_rr)))
            for t in axd.get_yticklabels() if t.get_text().strip()]
    _gaps = [_dbb[i][1][2] - _dbb[i + 1][1][3] for i in range(len(_dbb) - 1)]
    D5LAB_RESID["fig5"] = dict(
        labels=[dict(text=s, x0=b[0], x1=b[1], y0=b[2], y1=b[3]) for s, b in _dbb],
        height_mm=[b[3] - b[2] for _s, b in _dbb],
        gaps_mm=_gaps, min_gap_mm=min(_gaps) if _gaps else float("nan"),
        overlap_pairs=sum(1 for g in _gaps if g <= 0.0),
        ax_h_mm=(axd.get_window_extent(_rr).height) / _pxmm,
        linespacing=float(_p.get("d_lab_linespacing", 1.15)))
    D5BAR_RESID["fig5"] = dict(
        domains=list(dorder), median_or=dmed, bar_color=dcol,
        class_color_map={"<=1.2": FS.C["green_dark"], "1.2-1.5": FS.C["cyan"],
                         "1.5-2.0": FS.C["amber"], "2.0-3.0": FS.C["suggest_strong"],
                         ">3.0": FS.C["robust_strong"]},
        vline_color="#4D4D4D",
        n_distinct_colors=len(set(dcol)))
    _flg = axf.get_legend()
    if _flg is not None:
        _l = _flg.get_window_extent(_rr)
        _fb = axf.get_window_extent(_rr)
        F5LEG_RESID["fig5"] = dict(
            leg_x1_mm=_l.x1 / _pxmm, ax_x1_mm=_fb.x1 / _pxmm,
            right_over_mm=(_l.x1 - _fb.x1) / _pxmm,
            leg_top_frac=(_l.y1 - _fb.y0) / _fb.height,
            leg_bot_frac=(_l.y0 - _fb.y0) / _fb.height,
            top_gap_mm=(_fb.y1 - _l.y1) / _pxmm,
            anchor_y=float(_p.get("f_leg_anchor_y", 0.900)))
    if SHOW_NOTES:
        fig.text(0.5, 0.006,
                 "FDR < 0.05: 52 signals = 12 risk-increasing (b) + 40 risk-decreasing; 4,938 tests in total.\n"
                 "Panel (b): risk-increasing signals are attributed to the chr1p36.12 REGION, not to CDC42 - both\n"
                 "instruments are in high LD with the WNT4 credible-set lead rs56318008 (r2 = 0.899 / 0.847).\n"
                 "Panel (c): within the immune / blood / tumour domain (968 tests) NO risk-increasing signal was\n"
                 "detected within the available statistical power - this is not equivalent to established safety.\n"
                 "Panels (d, e): at 80% power only 23 of the 968 tests could detect OR <= 1.2 and 138 could detect\n"
                 "OR <= 1.5 (median detectable OR 3.15). The immune / blood / tumour domain is the WEAKEST of the\n"
                 "five domains, whereas the better-powered reproductive (1.74) and pregnancy (1.71) domains are\n"
                 "where the risk-increasing signals were found. Detecting OR = 1.5 needs a median of 4,202 extra\n"
                 "cases per endpoint (OR = 1.2: 20,842). Panel (f): conditioning the whole PheWAS on the regional\n"
                 "eQTL signal removes every FDR < 0.05 hit (52 -> 0) - the profile is region-driven, not gene-driven.\n"
                 "HLH and pancytopenia: no FinnGen R12 endpoint -> NOT assessable (a lack of evidence, "
                 "not evidence of 'no risk').\n"
                 "Panel (b) lineage groups: A pelvic-floor / connective tissue ; B placental / obstetric "
                 "haemorrhage ; C preterm birth ;\n"
                 "D knee osteoarthritis ; E residual foreign body in soft tissue ; "
                 "F non-disease fertility phenotype (not counted as a risk signal).",
                 ha="center", va="bottom", fontsize=5.2, color="#555555")
    return fig


if __name__ == "__main__":
    import json, traceback
    JOBS = [(fig1, "Fig1_study_design"), (fig2, "Fig2_discovery_replication"),
            (fig3, "Fig3_coloc_sensitivity"), (fig4, "Fig4_chr1p36_12_finemap"),
            (fig5, "Fig5_opentargets_phewas_safety")]
    ok = 0
    for fn, name in JOBS:
        try:
            _, wmm, hmm = fit_save(fn, name)
            print("   [OK] %-32s scale=%.4f  %.1f x %.1f mm" % (name, CAL[name], wmm, hmm))
            ok += 1
        except Exception as e:
            print("   [!! FAILED] %s: %s" % (name, e))
            traceback.print_exc()
    with open(os.path.join(FS.FIGDIR, "_fig_scale.json"), "w", encoding="utf-8") as fh:
        json.dump({"target_mm": TARGET_MM, "scales": CAL}, fh, indent=1)
    print("DONE  %d/%d" % (ok, len(JOBS)))
