#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""s39_delivery_manifest.py -- 生成最终交付清单（含大小与 SHA256）

只读；输出 27_最终交付清单_暂停点11.md
"""
import os, hashlib, datetime, csv, io

ROOT = r"D:\endometriosis_project\11_sc_eqtl_mr_project"
PRJ = r"D:\endometriosis_project"


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def human(n):
    for u in ("B", "KB", "MB", "GB"):
        if n < 1024:
            return "%.1f %s" % (n, u)
        n /= 1024.0
    return "%.1f TB" % n


FIG = [("1", "Fig1_study_design"), ("2", "Fig2_discovery_replication"),
       ("3", "Fig3_coloc_sensitivity"), ("4", "Fig4_chr1p36_12_finemap"),
       ("5", "Fig5_opentargets_phewas_safety")]

rows_a = []   # 主交付物
def add(group, name, rel, note=""):
    p = os.path.join(ROOT, rel)
    if os.path.exists(p):
        rows_a.append([group, name, rel, os.path.getsize(p), sha(p), note])
    else:
        rows_a.append([group, name, rel, -1, "MISSING", note])


add("正文", "正文框架与限制段定稿（含 Part A/B/C/D）", r"26_论文正文框架与限制段定稿_暂停点10.md", "暂停点 10 定稿 + 裁定修正")
add("正文", "限制段最终整合（源头）", r"22_论文限制段最终整合_暂停点8.md", "A–F 七组边界")
add("补充数据", "补充数据压缩包（15 条目）", r"Supplementary_Data.zip", "14 CSV + manifest")
for i, nm in FIG:
    add("图件", "Fig%s PDF" % i, r"figures\%s.pdf" % nm, "@600 dpi")
    add("图件", "Fig%s PNG" % i, r"figures\%s.png" % nm, "@600 dpi")
add("脚本", "绘图脚本", r"scripts\s37_final_figures.py", "Fig1–5 组装")
add("脚本", "绘图样式库", r"scripts\_fig_style.py", "Nature 规格 + 字形核验")
add("脚本", "补充数据导出脚本", r"scripts\s38_export_supplementary_data.py", "14 面板 CSV")
add("元数据", "图宽标定记录", r"figures\_fig_scale.json", "5 图 scale")
add("元数据", "补充数据 manifest", r"supplementary_data\_manifest.csv", "14 文件 SHA256")

# 补充数据逐个
supp = os.path.join(ROOT, "supplementary_data")
rows_s = []
if os.path.isdir(supp):
    for f in sorted(os.listdir(supp)):
        if f.endswith(".csv"):
            p = os.path.join(supp, f)
            rows_s.append([f, os.path.getsize(p), sha(p)])

# PheWAS 归档骨架
pw = os.path.join(ROOT, "phewas_delivery_20260924")
rows_p = []
for dp, dn, fn in os.walk(pw):
    for f in sorted(fn):
        p = os.path.join(dp, f)
        rows_p.append([os.path.relpath(p, pw).replace("\\", "/"),
                       os.path.getsize(p), sha(p)])

# 清理审计
cl = r"D:\_transfer_logs\_cleanup_backup_20260927"
rows_c = []
if os.path.isdir(cl):
    for f in sorted(os.listdir(cl)):
        p = os.path.join(cl, f)
        if os.path.isfile(p):
            rows_c.append([f, os.path.getsize(p), sha(p)])

L = []
L.append("# 27 · 最终交付清单（暂停点 11）")
L.append("")
L.append("**生成时间**：%s" % datetime.datetime.now().isoformat(timespec="seconds"))
L.append("**性质**：定稿收尾批处理完成 → **全部打包完毕，等待最后审查**")
L.append("**上游**：`26_论文正文框架与限制段定稿_暂停点10.md`（本清单为其裁定修正与交付封板）")
L.append("")
L.append("> ★**两项必须报告的偏差（不得静默处理）**：")
L.append("> ① 零引用文件「269 个」系 2026-09-24 过期值；重跑扫描实为 **19 个零引用文件 / 52.90 KB / 0 目录**，保留 2 项工具与审计后**实删 17 项**（详见 §3）。")
L.append("> ② PheWAS 交付包归档**受阻**：`Get-Volume` 权威枚举确认本机**无 E 盘、无任何可移动卷**，C/D 盘 0 副本；已建归档骨架 + 一键恢复脚本（详见 §4）。")
L.append("")
L.append("---")
L.append("")
L.append("## 1. 主交付物")
L.append("")
L.append("| 组 | 交付物 | 路径（相对 `11_sc_eqtl_mr_project/`） | 大小 | SHA256（前 16） | 备注 |")
L.append("|---|---|---|---|---|---|")
for g, n, rel, sz, h, note in rows_a:
    L.append("| %s | %s | `%s` | %s | `%s` | %s |" % (g, n, rel, human(sz) if sz >= 0 else "**MISSING**", h[:16], note))
L.append("")
L.append("## 2. 补充数据（`supplementary_data/`，逐文件）")
L.append("")
L.append("| # | 文件 | 大小 | SHA256（前 16） |")
L.append("|---|---|---|---|")
for i, (f, sz, h) in enumerate(rows_s, 1):
    L.append("| %d | `%s` | %s | `%s` |" % (i, f, human(sz), h[:16]))
L.append("")
L.append("**面板 → 文件映射**：Fig1 → `Fig1_study_design_parameters.csv`；Fig2(a/b/c) → `Fig2a_manhattan_source_8612tests.csv` / `Fig2b_forest_discovery_15pairs_plus_IVW.csv` / `Fig2c_replication_comparison.csv`；Fig3(a/b/c) → `Fig3a_coloc_pph4_by_pair.csv` / `Fig3b_sensitivity_matrix_5loci_x9settings.csv` / `Fig3c_pph4_vs_p12_mean_over_windows.csv`；Fig4(a/b/c/d) → `Fig4a_LD_matrix_7x7.csv` / `Fig4b_gene_variant_track.csv` / `Fig4c_finemap_pip_outcome_only.csv` / `Fig4d_sumpip_nearest_gene.csv`；Fig5(a) → `Fig5a_opentargets_top18.csv`；Fig5(b,c) → `Fig5bc_phewas_risk_increasing_12.csv` / `Fig5bc_phewas_domain_rollup.csv`。")
L.append("")
L.append("> 导出脚本 `scripts/s38_export_supplementary_data.py` 复用 `s37_final_figures.py` 的**同一加载/筛选/排序逻辑**，保证 CSV 逐行对应图上所绘元素；**只读源表，未修改任何 `tables/` 文件**。")
L.append("")
L.append("## 3. 零引用文件清理（Task #8 · 已完成）")
L.append("")
L.append("| 项 | 值 |")
L.append("|---|---|")
L.append("| 裁定文数字 | 「269 个」（2026-09-24 `_cleanup_refs.txt` 口径） |")
L.append("| **重跑实测** | **19 个零引用文件 / 52.90 KB / 0 个目录**（2026-09-27 20:11:32） |")
L.append("| 独立宽口径交叉验证 | 19/19 仍零引用（`_cleanup_verify19.py`） |")
L.append("| 保留未删 | `_cleanup_refs_v3.py`（扫描工具）、`_del_log.txt`（历史审计） |")
L.append("| **实删** | **17 项**（备份 17/17 sha256 一致 → 删除 0 失败 → 批后核验 0 残留） |")
L.append("| 备份目录 | `D:\\_transfer_logs\\_cleanup_backup_20260927\\`（可完整还原） |")
L.append("")
L.append("备份区内容（`%d` 文件）：" % len(rows_c))
L.append("")
L.append("| 文件 | 大小 | SHA256（前 16） |")
L.append("|---|---|---|")
for f, sz, h in rows_c:
    L.append("| `%s` | %s | `%s` |" % (f, human(sz), h[:16]))
L.append("")
L.append("## 4. PheWAS 交付包归档（Task #7 · **受阻**）")
L.append("")
L.append("**阻塞结论**：源包仅存于 `E:\\_handover_phewas_20260924\\`（外置盘）；`Get-Volume` 枚举确认本机仅 5 个固定卷（D/Onekey/C/WINPE/WinRE），**无 E 盘、无 Removable 卷**；C/D 盘全盘检索 `20_phewas_all_all.csv`、`package_meta.json`、`SHA256_checksums.txt` **均 0 命中**。")
L.append("")
L.append("已在 D 盘落地归档骨架（`%d` 文件）：" % len(rows_p))
L.append("")
L.append("| 文件 | 大小 | SHA256（前 16） |")
L.append("|---|---|---|")
for f, sz, h in rows_p:
    L.append("| `%s` | %s | `%s` |" % (f, human(sz), h[:16]))
L.append("")
L.append("**一键恢复命令 + 已核验事实表**见 `phewas_delivery_20260924/README_PENDING.md`。")
L.append("**唯一实质缺口**：补充材料 Table S9 所需 4,938 行全量表 `20_phewas_all_all.csv`（论文主结论/限制段/图件**均不受影响**）。")
L.append("")
L.append("## 5. 裁定执行核对（2026-09-27）")
L.append("")
L.append("| 裁定 | 要求 | 状态 |")
L.append("|---|---|---|")
L.append("| 1 | Fig5 双面板（左 a / 右 b+c） | ✅ 图已产出，图注已写明 |")
L.append("| 2 | 不补分子对接，Future Directions 一句（1CEE/1GRN） | ✅ §4.4；PDB 已向 RCSB 实查 |")
L.append("| 3 | 大内存机不再启动 | ✅ 本机完成全部任务 |")
L.append("| 4 | 结局来源更正（`GCST90483463`） | ✅ §0.1 已裁定；§2.1 口径统一 |")
L.append("| 5 | 不补 §3.9，仅 Supplementary Methods 一句 | ✅ §7 已写入 |")
L.append("| 6 | B1 补入 D2 等位异质性 | ✅ B1/B2 末句；β_exp 实查 −0.5939/+0.4135 |")
L.append("| 7 | Fig5 图注加 phenome-wide safety panel | ✅ 首句 + (c) + 脚注 |")
L.append("| 8 | 生成补充数据 CSV + zip | ✅ 14 CSV + manifest + zip（sha256 `10cf03f4…d338a7`） |")
L.append("| 9 | 删除零引用文件 | ⚠️ **实删 17 项**（「269」为过期值，见 §3） |")
L.append("| 10 | PheWAS 交付包归档 D 盘 | ⚠️ **受阻**（E 盘未连接，见 §4） |")
L.append("")
L.append("---")
L.append("")
L.append("**状态**：✅ 全部可交付物已打包完毕（除受 E 盘阻塞的 PheWAS 完整归档外），**等待最后审查**。")

outp = os.path.join(ROOT, "27_最终交付清单_暂停点11.md")
with io.open(outp, "w", encoding="utf-8") as fh:
    fh.write("\n".join(L) + "\n")
print("WROTE", outp, os.path.getsize(outp), "B")
