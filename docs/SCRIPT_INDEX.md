# 脚本索引（148 个 `s*` 脚本，按文件名排序）

> 「说明」列取自脚本自身的 docstring（`.py`）或首个有效注释行（`.R`），**机械抽取、未改写**。
> 生成脚本：`_s554_make_script_index.py`。

| # | 脚本 | 说明（脚本内自述） |
|---|---|---|
| 1 | `s01_inspect_top_eqtl.py` | （无自述） |
| 2 | `s02_onek1k_power_probe.py` | （无自述） |
| 3 | `s03_extract_instruments.py` | 从 OneK1K top cis-eQTL 表抽取 P<5e-8 工具变量，输出 MR 输入表。 |
| 4 | `s03b_extract_instruments_qval05.py` | s03b_extract_instruments_qval05.py —— 敏感性分析工具变量集：TensorQTL 置换 qval < 0.05 |
| 5 | `s04_positive_control_screen.py` | 在 OneK1K P<5e-8 工具变量表中筛选候选阳性对照基因，并统计其跨细胞类型强度。 |
| 6 | `s05_check_opengwas.py` | （无自述） |
| 7 | `s05_liftover_instruments.py` | s05_liftover_instruments.py —— 通用工具变量位置 liftover GRCh37 -> GRCh38 |
| 8 | `s06_mr_extract.py` | s06_mr_extract.py —— 通用 cis-eQTL → 结局 的 Wald ratio MR（流式抽取 + harmonise + FDR） |
| 9 | `s07_descriptive_locus_panel.py` | 描述性位点面板：对已知内异症/不孕症 GWAS 基因 + 本项目关注基因， |
| 10 | `s08_resume_ab_test.py` | 对照实验：Range（断点续传）请求是否导致 EBI / molgenis 挂起。 |
| 11 | `s09_gating_report.py` | s09_gating_report.py —— 阳性对照门控的形式化判定 |
| 12 | `s107_fig1_overview.py` | s107_fig1_overview.py —— Fig1「研究设计总览」重绘（分带式 study-design overview） |
| 13 | `s108_fig1_overview_v2.py` | s108_fig1_overview_v2.py —— Fig1 **V2「图形化研究设计总览」**（graphical-abstract 版） |
| 14 | `s109_make_v8_fig1v2.py` | s109_make_v8_fig1v2.py —— 中文版论文 v8：把 Fig1 换成【V2 图形化版】 |
| 15 | `s10_discovery_aggregate.py` | s10_discovery_aggregate.py —— 把 s06 的逐变异 MR 结果聚合为 gene x cell_type 结果表 |
| 16 | `s110_make_v9_gb12.py` | s110_make_v9_gb12.py —— v8 → v9：把论文里**全部 9 张图**换成 GB12／α=0.35 版。 |
| 17 | `s111_fig_audit_nature.py` | s111_fig_audit_nature.py —— 对 figures/ 全部图件做「Nature 规范符合性」实测盘点。 |
| 18 | `s112_font_vector_audit.py` | s112_font_vector_audit.py —— 只读审计：figures/ 各矢量 PDF 的 |
| 19 | `s113_raster_dpi_check.py` | s113_raster_dpi_check.py —— 只读：figures/ 各 PDF 内嵌位图的"有效分辨率"。 |
| 20 | `s114_p0_shrink_height.py` | s114_p0_shrink_height.py —— P0：把 Fig3 / Fig4 / Fig5 压到幅高 <= 170 mm。 |
| 21 | `s115_p1_panel_labels.py` | s115_p1_panel_labels.py —— P1：重出 s37 的 Fig2/Fig3/Fig4/Fig5。 |
| 22 | `s116_fig2_diagnose.py` | s116_fig2_diagnose.py —— Fig2 布局诊断（只读，不 save）。 |
| 23 | `s117_fig2_bgcolor_audit.py` | s117_fig2_bgcolor_audit.py —— Fig2(a) 背景浅色系的客观审计（只读）。 |
| 24 | `s117b_select_bg6.py` | s117b_select_bg6.py —— 为 Fig2(a) 选 6 色「背景浅色系」（只看不改）。 |
| 25 | `s118_fig2_revise.py` | s118_fig2_revise.py —— Fig2 排版与配色专改（图件细改 · Fig2 · 问题1+2+3）。 |
| 26 | `s119_regen_panels_upper.py` | s119_regen_panels_upper.py —— 面板字母大写后重出 s37 的 Fig3/4/5。 |
| 27 | `s11_outcome_loci_scan.py` | s11_outcome_loci_scan.py —— 结局文件自身的全基因组显著位点扫描（QC + 位点归属） |
| 28 | `s120_fig2_preview.py` | s120_fig2_preview.py —— 生成 Fig2 重绘后的预览图（只读 figures/，不重绘）。 |
| 29 | `s121_update_deliverables.py` | s121_update_deliverables.py —— 更新交付清单 + 核验源数据完整性。 |
| 30 | `s122_bg6_yellow_swap.py` | s122_bg6_yellow_swap.py —— 给 Fig2(a) 的「黄色」换色（只算不改）。 |
| 31 | `s12_discovery_consolidate.py` | s12_discovery_consolidate.py —— Discovery 阶段结果整合 |
| 32 | `s130_letter_candidates.py` | s130_letter_candidates.py —— 生成 Fig2 面板字母 B/C 的**位置候选对比图**（只读正本）。 |
| 33 | `s13_discovery_figures.py` | s13_discovery_figures.py —— Discovery 阶段诊断图 |
| 34 | `s149_make_v10.py` | s149_make_v10.py —— v9 → v10：把论文里的 Fig2–Fig5 / FigS1–S4 换成**点风格统一后**的版本。 |
| 35 | `s14_discovery_locus_and_power.py` | Discovery 小任务交付物： |
| 36 | `s15_chr1_ld_data.py` | chr1 基因座 LD 结构数据构建 |
| 37 | `s166_make_v11.py` | s166_make_v11.py —— v10 → v11：把论文里的 Fig2–Fig5 / FigS1–S4 换成**三项全局标准统一后**的版本。 |
| 38 | `s16_chr1_ld_figure.py` | FigD5 —— chr1p36.12 基因座 LD 结构图 |
| 39 | `s17_replication_mr.py` | s17 —— Replication MR（1M-scBloodNL 作为暴露，结局沿用 GCST90483463） |
| 40 | `s18_flag_weak_instrument.py` | s18 —— 复现表的「弱工具效应量」标记（Yan 2026-09-23 裁定） |
| 41 | `s19_p2_raw_verify.py` | s19_p2_raw_verify.py —— P2 OneK1K raw tensorQTL summary tar.gz 格式核验（暂停点 2） |
| 42 | `s20_p2_extract_and_profile.py` | s20_p2_extract_and_profile.py —— 暂停点 2：从 tar.gz 流式抽取样本成员并画像 |
| 43 | `s21_p2_member_profile.py` | s21_p2_member_profile.py —— 逐成员（cell_type × chr）画像，带断点续跑 |
| 44 | `s22_p2_link_discovery.py` | s22_p2_link_discovery.py —— raw(tar.gz) × Discovery 工具表 连接性实测 |
| 45 | `s23_raw_instruments.py` | s23_raw_instruments.py —— 从 OneK1K raw parquet 提取三座工具变量（P<5e-8） |
| 46 | `s24_clump_raw_instruments.py` | s24_clump_raw_instruments.py —— 对三座 raw 工具变量做 PLINK clumping（per gene×cell_type） |
| 47 | `s25_prep_ivw_inputs.py` | s25_prep_ivw_inputs.py —— 生成 s06 / s10 可直接消费的工具表 |
| 48 | `s26_coloc_prep.py` | s26_coloc_prep.py —— 共定位数据准备（Python 侧） |
| 49 | `s27_ivw_and_merge.py` | s27_ivw_and_merge.py —— IVW 汇总 + 与 Discovery 合并（用户裁定 (b)：并列报告） |
| 50 | `s28_coloc.R` | -*- coding: utf-8 -*- |
| 51 | `s29_coloc_summary.py` | s29_coloc_summary.py -- collapse tables/33_coloc_results.csv (per-row: 1 abf + K susie pairs) |
| 52 | `s30_finemap_chr1.R` | -*- coding: utf-8 -*- |
| 53 | `s31_sensitivity_outcome_mr.py` | s31_sensitivity_outcome_mr.py —— Task A：P3 敏感性结局（GCST90483469）MR |
| 54 | `s32_steiger_egger.py` | s32_steiger_egger.py —— 任务二/四：Steiger 方向性检验 + MR-Egger |
| 55 | `s33_coloc_sensitivity_prep.py` | s33_coloc_sensitivity_prep.py —— 为 coloc 敏感性分析准备 3 窗口 × 输入数据 |
| 56 | `s34_coloc_sensitivity.R` | -*- coding: utf-8 -*- |
| 57 | `s35_opentargets.py` | s35_opentargets.py —— 下游：Open Targets 药物靶点 + 简化版 PheWAS（关联疾病谱） |
| 58 | `s35b_opentargets_ext.py` | s35b_opentargets_ext.py —— 补充：化学探针 + 通路 + 大疾病列表（生殖/自身免疫筛选） |
| 59 | `s36_phewas_safety_integrate.py` | s36_phewas_safety_integrate.py |
| 60 | `s37_final_figures.py` | s37_final_figures.py -- 论文最终图表 Fig1-Fig5 组装 |
| 61 | `s38_export_supplementary_data.py` | s38_export_supplementary_data.py -- 补充数据导出（P3 版） |
| 62 | `s39_delivery_manifest.py` | s39_delivery_manifest.py -- 生成最终交付清单（含大小与 SHA256） |
| 63 | `s40a_recon_3p1.py` | s40a_recon_3p1.py -- 任务 3.1 条件共定位：输入盘点 |
| 64 | `s40b_prep_cond_inputs.py` | s40b_prep_cond_inputs.py -- 任务 3.1 步骤 1：条件共定位输入准备 |
| 65 | `s40c_cond_coloc.R` | -*- coding: utf-8 -*- |
| 66 | `s40d_cond_coloc_negctl.R` | -*- coding: utf-8 -*- |
| 67 | `s40d_prep_negctl_ld.py` | s40d_prep_negctl_ld.py -- 任务 3.1 补丁：为「程序性阴性对照」计算**配对正确的带符号 r** |
| 68 | `s41a_fetch_pheweb.py` | s41a_fetch_pheweb.py -- 任务 3.2 步骤1：抓取 FinnGen R12 三个 SNP 的全端点关联结果 |
| 69 | `s41b_phewas_3p2.py` | s41b_phewas_3p2.py -- 任务 3.2：条件 PheWAS（a）与 区域工具 PheWAS（b） |
| 70 | `s42a_tabix_probe.py` | s42a_tabix_probe.py -- HTTP Range + tabix 单点取数（可行性验证，含逐块 BGZF 解压） |
| 71 | `s42b_phewas_tabix_fetch.py` | s42b_phewas_tabix_fetch.py -- 任务 3.2 数据层：用 HTTP Range + tabix 取**全端点**面板 |
| 72 | `s42c_phewas_3p2_full.py` | s42c_phewas_3p2_full.py -- 任务 3.2 主分析（**完整 2,469 端点面板**，来自 FinnGen R12 公开桶 tabix） |
| 73 | `s42d_phewas_region_raw_supp.py` | s42d_phewas_region_raw_supp.py -- 任务 3.2 (b) 的**假设自由**佐证 |
| 74 | `s43a_build_smr_inputs.py` | s43a_build_smr_inputs.py -- 任务 3.3 SMR 输入构建 |
| 75 | `s43b_build_besd.py` | s43b_build_besd.py — 任务三 3.3 步骤 1：把 14 个 per-cell eQTL query 文件转成 SMR BESD 格式。 |
| 76 | `s43c_run_smr.py` | s43c_run_smr.py — 任务三 3.3 步骤 2/3：SMR + HEIDI 主分析 + 敏感性，并汇总。 |
| 77 | `s43d_smr_summary.py` | s43d_smr_summary.py — 任务三 3.3 步骤 4：SMR/HEIDI 结果汇总与判读 |
| 78 | `s43e_fetch_gtex8_wb_stream.py` | s43e_fetch_gtex8_wb_stream.py — 选项 A：流式抓取 eQTL Catalogue GTEx v8 全血（区域过滤，提前停） |
| 79 | `s43f_build_gtex8_wnt4.py` | s43f_build_gtex8_wnt4.py — 由 GTEx v8 全血区域数据构建 SMR 输入（BESD） |
| 80 | `s43g_run_gtex8_smr.py` | s43g_run_gtex8_smr.py — 选项 A 收尾：GTEx v8 全血（bulk）3 基因 SMR + HEIDI |
| 81 | `s43h_fix_query_freq.py` | s43h_fix_query_freq.py — 修正 GTEx v8 全血 SMR query 的 Freq 取向（选项 A 阻塞点修复） |
| 82 | `s43i_run_gtex8_smr_fix.py` | s43i_run_gtex8_smr_fix.py — 选项 A 正式运行：GTEx v8 全血（bulk）3 基因 SMR + HEIDI |
| 83 | `s44_task34_tool_neighbour_genes.py` | s44_task34_tool_neighbour_genes.py — 任务 3.4：工具邻近基因 cis 效应检查 |
| 84 | `s45_task35_phewas_power.py` | s45_task35_phewas_power.py — 任务三 3.5：PheWAS 统计功效分析 |
| 85 | `s51a2_coord_forensics.py` | s51a2_coord_forensics.py -- 判定 parquet variant_id 的坐标体系，并定位真正的工具变量 |
| 86 | `s51a_verify_L2L3_instruments.py` | s51a_verify_L2L3_instruments.py -- 任务A 步骤1/2 前置核验（修正版） |
| 87 | `s51b_coloc_sens_prep_L2L3.py` | s51b_coloc_sens_prep_L2L3.py -- 任务A：为 L2 (YME1L1/CD4_NC) 生成 coloc 敏感性输入 |
| 88 | `s51c_coloc_sens_L2L3.R` | -*- coding: utf-8 -*- |
| 89 | `s52a_prep_cond_L2L3.py` | s52a_prep_cond_L2L3.py -- 任务A：L2/L3 条件共定位输入准备（两阶段） |
| 90 | `s52b0_diag_L3_allele.py` | s52b0_diag_L3_allele.py -- 诊断 L3 判据A 偏差来源（只读） |
| 91 | `s52b1_diag_L3_N.py` | s52b1_diag_L3_N.py -- 决定性检验：parquet af 的分母是否为 980 供者 |
| 92 | `s52b2_scan_celltype_N.py` | s52b2_scan_celltype_N.py -- 扫描全部细胞类型 parquet 的有效样本量 |
| 93 | `s52b_cond_coloc_L2L3.R` | -*- coding: utf-8 -*- |
| 94 | `s52c_inspect_coloc.R` | -*- coding: utf-8 -*- |
| 95 | `s52d_N_sensitivity.R` | -*- coding: utf-8 -*- |
| 96 | `s52e_L3_sens_N690.R` | -*- coding: utf-8 -*- |
| 97 | `s53_finemap_L2L3.R` | -*- coding: utf-8 -*- |
| 98 | `s54a_build_smr_inputs_L2L3.py` | s54a_build_smr_inputs_L2L3.py -- 任务A：L2 (chr10 YME1L1/CD4_NC) 与 L3 (chr2 ANXA4/Mono_NC) |
| 99 | `s54b_build_besd_L2L3.py` | s54b_build_besd_L2L3.py — 任务A：L2/L3 SMR BESD 构建（复刻 s43b） |
| 100 | `s54c_run_smr_L2L3.py` | s54c_run_smr_L2L3.py — 任务A：L2/L3 的 cis-SMR + HEIDI（复刻 s43c） |
| 101 | `s55a_fetch_gtex8_wb_L2L3.py` | s55a_fetch_gtex8_wb_L2L3.py -- 裁定1：从 EBI 定向抽取 GTEx v8 全血的 L2(YME1L1) / L3(ANXA4) 区域 |
| 102 | `s55b_build_gtex8_query_L2L3.py` | s55b_build_gtex8_query_L2L3.py -- 裁定1：由 GTEx v8 区域数据构建 SMR query（L2 YME1L1 / L3 ANXA4） |
| 103 | `s55c_run_gtex8_smr_L2L3.py` | s55c_run_gtex8_smr_L2L3.py -- 裁定1：GTEx v8 全血（bulk, n=670）L2/L3 规范 cis-SMR + HEIDI |
| 104 | `s57_build_3locus_table.py` | s57_build_3locus_table.py -- 任务A 步骤4：三位点（chr1 / chr10 / chr2）比较表 |
| 105 | `s60b_build_mvmr_inputs.py` | s60b_build_mvmr_inputs.py -- 任务B 第1步：构建 CDC42 + LINC00339 的 MVMR 输入（逐细胞类型） |
| 106 | `s60c_run_mvmr.R` | -*- coding: utf-8 -*- |
| 107 | `s60d_mvmr_decomposition.py` | s60d_mvmr_decomposition.py -- 任务B 第3步：效应分解对照（同一工具集下 单暴露 IVW vs 双暴露 MVMR） |
| 108 | `s61_apply_v4_insertions.py` | s61_apply_v4_insertions.py |
| 109 | `s62_build_techdoc_v8.py` | s62_build_techdoc_v8.py |
| 110 | `s63_build_v4_from_docx.py` | s63_build_v4_from_docx.py |
| 111 | `s63_precheck.py` | s63_precheck.py —— 只读预检（不写任何产物） |
| 112 | `s63b_zipdiff.py` | s63b_zipdiff.py —— 对底本 docx 与输出 docx 做 zip 部件级比对（只读） |
| 113 | `s63c_probe_diag.py` | s63c_probe_diag.py —— 定位 verify 中唯一 miss 的探针在 PDF 文本层的实际形态 |
| 114 | `s63d_lnnum_diag.py` | s63d_lnnum_diag.py —— 核查 Word 底本是否启用了行号（w:lnNumType），并 dump 相关页原始文本 |
| 115 | `s63e_locate_inserts.py` | s63e_locate_inserts.py —— 定位 4 处插入与 5 张图在新 v4 预览 PDF 中的页码 |
| 116 | `s63f_tail_dump.py` | （无自述） |
| 117 | `s64a_recon.py` | s64a_recon.py —— 只读侦察 |
| 118 | `s64b_make_v5.py` | s64b_make_v5.py —— 由 v4（docx 血统）产出 v5 |
| 119 | `s64c_refscan.py` | s64c_refscan.py —— 引用扫描（修正版：按**目录名分量**排除，不再用子串 '_'） |
| 120 | `s64d_rename.py` | s64d_rename.py —— 裁定 6 的档务重命名（★纯重命名，被重命名文件的内容零改动） |
| 121 | `s65_figure_overview.py` | s65_figure_overview.py —— 一页式「图件总览」拼版图（Nature 体例） |
| 122 | `s66_edrive_figure_scan.py` | s66_edrive_figure_scan.py —— E 盘图件资产检索 + 本项目快照差异比对（只读） |
| 123 | `s67_make_supp_figs.py` | s67_make_supp_figs.py -- 补充图 Fig S1 / Fig S2 重做 |
| 124 | `s67a_probe_supp_data.py` | Probe candidate data sources for supplementary Fig S1 / S2. |
| 125 | `s67b_prep_supp.py` | Prep + verify data for supplementary Fig S1 / S2 (read-only). |
| 126 | `s67c_verify_supp_figs.py` | s67c_verify_supp_figs.py — Fig S1 / Fig S2 验收核验（独立判据，不依赖出图脚本的内建日志） |
| 127 | `s67d_crop_check.py` | s67d_crop_check.py — 裁剪 FigS2 中缝区域，用于目视排查「文字压文字」 |
| 128 | `s67e_diff_manifest.py` | s67e_diff_manifest.py — 对比 s38 重跑前后的 manifest，确认旧行逐位不变、新行已加入 |
| 129 | `s67f_check_overview.py` | s67f_check_overview.py — 总览图 PDF 字形/几何核验 |
| 130 | `s67g_crop_s4.py` | s67g_crop_s4.py — 原分辨率裁切 Fig S4 局部，用于「文字压点」目视判据 |
| 131 | `s69_patch_manuscript.py` | s69_patch_manuscript.py — 按用户裁定（2026-09-29）修改当前主线稿件 [v2: 幂等 + 行尾自适应] |
| 132 | `s69b_verify_manuscript.py` | s69b_verify_manuscript.py — 复核 s69 补丁结果（裁定 2/3/5 的判据） |
| 133 | `s69c_hash_report.py` | s69c_hash_report.py — 汇总当前交付物 sha256 与登记信息（供 65_ 交付报告引用） |
| 134 | `s70_patch_figure_inventory.py` | s70 · 把 63_figure_inventory.csv 中 5 张 D 系列诊断图的「状态」由笼统的 |
| 135 | `s71_recon_v4.py` | s71_recon_v4.py —— v4 docx 结构侦察（为 v6 构建定位锚点） |
| 136 | `s72_make_v6.py` | s72_make_v6.py —— 由 v4（docx 血统）产出 44_论文初稿_带图_v6 |
| 137 | `s72a_recon_body.py` | s72a_recon_body.py —— 侦察 v4 docx body 的 XML 子元素序列（含空段落）与关键段落的 pPr/style |
| 138 | `s72b_probe_api.py` | s72b_probe_api.py —— 探查 python-docx 版本与图片插入 API |
| 139 | `s72c_crosscheck.py` | s72c_crosscheck.py —— v6 与 v5（docx + md）的交叉核验 |
| 140 | `s72d_render_pages.py` | s72d_render_pages.py —— 渲染 v6 预览 PDF 的关键页为 PNG（目视复核用） |
| 141 | `s73_patch_inventory_v6.py` | s73_patch_inventory_v6.py —— 给 tables/63_figure_inventory.csv 增列「在v6预览PDF页」 |
| 142 | `s74_reorder_inventory_col.py` | s74 · 修正 tables/63_figure_inventory.csv 的列序 |
| 143 | `s75a_add_notes_switch.py` | s75a · 给 s37_final_figures.py 增加 FIG_NOTES 开关（默认 =1，保持原行为） |
| 144 | `s75b_regen_figs_no_notes.py` | s75b · FIG_NOTES=0 重出 Fig2–Fig5（**不重出 Fig1**） |
| 145 | `s75c_verify_no_notes.py` | s75c · 核验「仅脚注消失、其余逐一同」 |
| 146 | `s76a_recon_v6_images.py` | s76a · 侦察 v6 docx 的图片结构：media 文件 / rId / 每个 inline 的 extent |
| 147 | `s76b_make_v7.py` | s76b · 构建 45_论文初稿_带图_v7.docx |
| 148 | `s77_patch_inventory_v7.py` | s77 · 更新 tables/63_figure_inventory.csv 到 v7 |
