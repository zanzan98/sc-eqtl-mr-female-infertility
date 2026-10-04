# 运行顺序（RUN ORDER）

> 脚本编号即**推荐运行顺序**（同号内字母后缀按 a→b→c…）。每一步的输入/输出见 `SCRIPT_INDEX.md` 的自述列。
> ★ 全部脚本**硬编码 `D:/` 绝对路径**，移植前见 `README.md` §4。
> ★ 重跑会**覆盖**同名产物；正式图件重出须 `FIGPAL=gb12 FIG_NOTES=0`，且先用 `FIG_OUTDIR` 指向临时目录做回归。

---

## 阶段总览

| 阶段 | 编号区间 | 内容 | 关键产物 |
|---|---|---|---|
| **L0** 环境与数据准备 | `s01`–`s02`、`s05_check_opengwas`、`s08` | 资源盘点、功效探针、连通性/断点续传对照实验 | `logs/` 探针日志 |
| **L1a** 工具变量与门控 | `s03`–`s09`、`s19`–`s25` | 工具抽取、liftover、阳性对照门控、raw 核验、clumping | 工具表（`tables/`） |
| **L1b** Discovery | `s10`–`s16`、`s18` | 逐变异 MR 聚合、结局位点扫描、诊断图、chr1 LD | 8,612 个基因×细胞类型检验 |
| **L1c** 独立复核 | `s17` | 1M-scBloodNL 作为暴露的复核 MR | `tables/` 复核表 |
| **L1d** 共定位与精细定位 | `s26`–`s30` | coloc 输入、coloc.abf（R）、susie 汇总、chr1 精细定位（R） | PP.H4 / 可信集 |
| **L1e** 敏感性分析 | `s31`–`s34` | 敏感性结局 MR、Steiger/MR-Egger、coloc 三窗口（R） | 敏感性矩阵 |
| **L1f** 下游注释 | `s35`–`s36` | Open Targets 关联谱、PheWAS 安全性整合 | 关联/安全性表 |
| **L2/L3** 三层级扩展 | `s51a`–`s55c` | chr10 `YME1L1`(L2) 与 chr2 `ANXA4`(L3) 的完整平行分析 | 三位点可比输出 |
| **区域比较** | `s57` | 三位点（chr1 / chr10 / chr2）比较表 | 比较表 |
| **MVMR** | `s60b`–`s60d` | `CDC42` + `LINC00339` 双暴露 MVMR 与效应分解 | 分解表 |
| **PheWAS** | `s40a`–`s40d`、`s41a`–`s42d` | 条件共定位（含阴性对照）、FinnGen R12 全端点扫描 | 2,469 端点面板（表 **S9a**） |
| **cis-SMR** | `s43a`–`s43i`、`s54a`–`s54c`、`s55a`–`s55c` | BESD 构建 + SMR/HEIDI（单细胞面板 + GTEx v8 全血替代层） | SMR/HEIDI 输出 |
| **任务 3.4/3.5** | `s44`–`s45` | 工具邻近基因 cis 效应、PheWAS 功效分析 | 功效边界表 |
| **补充数据导出** | `s38`–`s39` | 补充数据打包、交付清单（含 SHA256） | `Supplementary_Data.zip`、清单 |
| **图件** | `s37`（主图）、`s67`/`s67a`–`s67g`（补充图）、`s65`（总览） | Fig1–Fig5 与 FigS5–FigS8 出图与验收 | `figures/` |
| **图件审计与迭代** | `s111`–`s122`、`s130` | Nature 规范审计、幅高/字号/字母/配色迭代 | 图件定稿参数 |
| **稿件构建** | `s61`–`s77`、`s109`、`s110`、`s149`、`s166` | docx 版本构建、逐节替换、页级核验 | 各版 docx + 预览 PDF |

---

## 1 L0 准备

| 脚本 | 要点 |
|---|---|
| `s01_inspect_top_eqtl.py` | 盘点 OneK1K top cis-eQTL 表结构 |
| `s02_onek1k_power_probe.py` | 各细胞类型工具数 / 功效探针 |
| `s05_check_opengwas.py` | OpenGWAS 连通性与端点可用性 |
| `s08_resume_ab_test.py` | **对照实验**：HTTP Range（断点续传）是否导致 EBI / molgenis 挂起 ⇒ 决定后续下载策略 |

## 2 L1a 工具变量与门控

| 脚本 | 要点 |
|---|---|
| `s03` / `s03b` | 抽 `P<5e-8` 工具变量（`s03b` 为 `qval<0.05` 敏感性集） |
| `s04` | 筛选候选**阳性对照**基因并统计跨细胞类型强度 |
| `s05_liftover_instruments.py` | GRCh37 → GRCh38（暴露/结局坐标系不同，按 rsID 匹配） |
| `s06_mr_extract.py` | 通用 cis-eQTL→结局 Wald ratio MR（流式抽取 + harmonise + FDR）；同时用于**阳性对照门控**（`GCST90132222`，类风湿关节炎）与主分析 |
| `s07` | 描述性位点面板 |
| `s09` | 阳性对照门控的**形式化判定**（不通过则不得进入主分析） |
| `s19`–`s22` | 暂停点 2：raw `tar.gz` 格式核验、成员画像、与 Discovery 工具表的连接性实测 |
| `s23`–`s25` | 从 raw parquet 抽三座工具变量 → PLINK clumping（per gene × cell_type）→ 生成 IVW 输入 |
| `s18` | 复现表弱工具效应量标记 |

## 3 L1b–L1c Discovery 与复核

| 脚本 | 要点 |
|---|---|
| `s10_discovery_aggregate.py` | 逐变异 MR → gene × cell_type 结果表（**8,612 个检验**） |
| `s11` | 结局文件自身全基因组显著位点扫描（QC + 位点归属） |
| `s12` / `s13` | Discovery 整合与诊断图 |
| `s14` | Discovery 小任务交付（位点 + 功效） |
| `s15` / `s16` | chr1 基因座 LD 数据构建 / LD 结构图 |
| `s17` | **独立复核**：1M-scBloodNL 为暴露，结局 `GCST90483463` |
| `s18` | 弱工具标记 |

★ **复核口径**（红线）：复核面板为**合并注释**（`B_MEM`/`B_IN`→B、`Mono_C`/`Mono_NC`→mono）⇒ 只能称
`direction-consistent nominal replication`，**禁称"亚型级独立复现"**。

## 4 L1d 共定位与精细定位

| 脚本 | 要点 |
|---|---|
| `s26` → `s28_coloc.R` → `s29` | coloc 输入 → `coloc.abf` + `coloc.susie`（R）→ 汇总（每行 1 个 abf + K 个 susie 对） |
| `s30_finemap_chr1.R` | chr1p36.12 精细定位（R） |
| `s27_ivw_and_merge.py` | IVW 汇总并与 Discovery 合并（用户裁定：**并列报告**） |

## 5 L1e 敏感性

| 脚本 | 要点 |
|---|---|
| `s31` | 敏感性结局 `GCST90483469` 的 MR |
| `s32` | Steiger 方向性检验 + MR-Egger |
| `s33` → `s34_coloc_sensitivity.R` | 3 窗口 × 设定下的 coloc 敏感性 |

## 6 L2/L3 三层级扩展（`s51a`–`s55c`）

对 **L2 = chr10 `YME1L1` / `CD4_NC`** 与 **L3 = chr2 `ANXA4` / `Mono_NC`** 做与 chr1 **完全平行**的分析：

```
s51a / s51a2   工具变量核验 + 坐标体系取证
s51b → s51c.R  coloc 敏感性输入 → coloc（R）
s52a → s52b0/s52b1/s52b2（诊断）→ s52b.R → s52c.R / s52d.R / s52e.R
               条件共定位输入 → 等位/样本量诊断 → 条件共定位 → 检查/样本量敏感性（R）
s53.R          精细定位（R）
s54a → s54b → s54c   SMR 输入 → BESD → cis-SMR + HEIDI
s55a → s55b → s55c   GTEx v8 全血（bulk, n=670）区域抽取 → SMR query → 规范 cis-SMR + HEIDI
s57             三位点比较表
```

> `s52b0`/`s52b1`/`s52b2` 为**只读诊断**脚本，仅在结果异常时排查用，不参与正式产物。

## 7 MVMR（`s60b`–`s60d`）

| 脚本 | 要点 |
|---|---|
| `s60b_build_mvmr_inputs.py` | `CDC42` + `LINC00339` 双暴露 MVMR 输入（逐细胞类型） |
| `s60c_run_mvmr.R` | 运行 MVMR（R） |
| `s60d_mvmr_decomposition.py` | 效应分解对照：同一工具集下**单暴露 IVW vs 双暴露 MVMR** |

## 8 PheWAS（FinnGen R12）

| 路线 | 脚本 | 要点 |
|---|---|---|
| 条件共定位（任务 3.1） | `s40a` → `s40b` → `s40c.R` / `s40d.R`（+ `s40d_prep_negctl_ld.py`） | 条件共定位输入 → 条件共定位 → **程序性阴性对照**（配对正确的带符号 r） |
| PheWeb API | `s41a` → `s41b` | 抓取三个 SNP 全端点关联 → 条件 PheWAS（a）+ 区域工具 PheWAS（b） |
| tabix 全端点 | `s42a`（可行性）→ `s42b`（HTTP Range + tabix 取数）→ `s42c`（**完整 2,469 端点主分析**）→ `s42d`（假设自由佐证） | 表 **S9a** 来源 |

★ **口径红线**：PheWAS 升高一律归因「**chr1p36.12 区域**」；「0 条升高」写作「**现有功效下未检出**」并给出功效边界。

## 9 cis-SMR + HEIDI

| 脚本 | 要点 |
|---|---|
| `s43a` → `s43b` → `s43c`（SMR 1.3.1 + PLINK LD）→ `s43d` | 单细胞面板 BESD → SMR/HEIDI 主分析 + 敏感性 → 汇总判读 |
| `s43e` → `s43f` → `s43g` → `s43h` → `s43i` | GTEx v8 全血**替代层**（`WNT4` 在单细胞面板缺失时的替代）；`s43h` 修 Freq 取向 |
| `s54a`–`s54c` | L2/L3 单细胞面板 SMR |
| `s55a`–`s55c` | L2/L3 GTEx v8 全血 SMR |

## 10 任务 3.4 / 3.5

| 脚本 | 要点 |
|---|---|
| `s44` | 工具邻近基因 cis 效应检查 |
| `s45` | PheWAS 统计功效分析（产出功效边界；`3.3A WNT4` 阴性为**低功效**，禁写"因果被否证"） |

## 11 补充数据与清单

| 脚本 | 要点 |
|---|---|
| `s38_export_supplementary_data.py` | 补充数据导出（21 面板级 + 3 补充图源 + 25 表 S9–S15） |
| `s39_delivery_manifest.py` | 交付清单（含大小与 SHA256） |

## 12 图件

| 脚本 | 要点 |
|---|---|
| `s37_final_figures.py` | **主图 Fig1–Fig5 组装**（`FIGPAL` 选调色板、`FIG_NOTES=0` 去脚注、`FIG_OUTDIR` 改输出目录） |
| `s67_make_supp_figs.py`（+ `s67a`–`s67g`） | **补充图 FigS5–FigS8** 出图与独立验收 |
| `s65_figure_overview.py` | 一页式图件总览（Nature 体例） |
| `s111` / `s112` / `s113` | **只读审计**：Nature 规范符合性、矢量 PDF 字体、位图有效分辨率 |
| `s114`–`s122`、`s130` | Fig2–Fig5 幅高 / 字号 / 字母位置 / 配色迭代（每轮均为"只读诊断 → 定稿重出"） |

★ **重出图件的三条硬约束**（见 `README.md` §6）：`FIGPAL=gb12 FIG_NOTES=0`；`_gb12_palette.py` 是配色唯一真源；
`_fig_style.save()` 的比例迭代在不可达目标处会振荡 ⇒ 正式图件**钉住一个窗内 scale 直接落盘**，且**一个 scale 一个进程**。

## 13 稿件构建

| 脚本 | 要点 |
|---|---|
| `s61`–`s64` | v4/v5 构建：段落级插入、引用扫描、纯重命名档务 |
| `s69`–`s77` | 稿件补丁（幂等）、图件清单同步、docx 图片结构侦察与 v6/v7 构建 |
| `s109` / `s110` / `s149` / `s166` | v8（Fig1 换 V2 图形化版）/ v9（全图换 GB12）/ v10（点风格统一）/ v11（三项全局标准统一） |
| `s63b`–`s63f`、`s72d`、`s75c` | zip 部件级比对、PDF 文本层定位、页渲染目视核验、脚注消失核验 |

> 稿件构建遵循「**字符级重建**」范式：只替换 `word/document.xml`，其余部件逐字节保留；
> 变更后必须做**段数 / 部件数 / `word/media` 逐字节 / `sectPr` 全子元素**四项守恒核验。

---

## 附：结果表与图件的落点

| 产物 | 目录 | 说明 |
|---|---|---|
| 结果表 | `tables/` | **150 个 CSV** |
| 图件 | `figures/` | 主图 `Fig1`–`Fig5`、补充图 `FigS5`–`FigS8`（各 PDF + PNG） |
| 补充数据包 | 工程根 | `Supplementary_Data.zip`（50 条） |
| 日志 / 中间产物 | `logs/`、`logs/` 与工程外 `D:/_transfer_logs/` | 探针与验收日志 |

★ 磁盘上的解包目录名为 `supplementary_data/`（全小写），而 **zip 内条目前缀为 `Supplementary_Data/`**（首字母大写）；
两者由构建脚本保证内容一一对应。`_manifest.csv` 的 `source` 列指向内部溯源件 `tables/64a_/64b_/64d_`。
