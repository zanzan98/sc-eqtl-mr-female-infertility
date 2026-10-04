# Code release — Single-cell eQTL Mendelian randomization at chr1p36.12 in female infertility

**单细胞 eQTL 孟德尔随机化揭示女性不孕症的免疫细胞类型特异性因果关联区域：chr1p36.12 区域关联与 WNT4 优先候选证据**

本仓库为论文的**分析代码与复现材料**（DAS / Data Availability Statement 对应件）。原始数据不在本仓库，请按
[`docs/DATA_SOURCES.md`](docs/DATA_SOURCES.md) 从各公共来源获取。

| 项 | 值 |
|---|---|
| 主稿版本 | `59_中文版论文_带图_v32_中文投稿模板.docx`（**176 段 / 30 页 / 9 内嵌图 / 0 表格**） |
| 分析脚本 | **148 个**（136 Python + 12 R），见 [`docs/SCRIPT_INDEX.md`](docs/SCRIPT_INDEX.md) |
| 结果表 | **150 个 CSV**（`tables/`，未随本仓库分发，见"结果表"一节） |
| 图件 | 主图 `Fig1`–`Fig5`、补充图 `FigS5`–`FigS8`，各 PDF + PNG |
| 补充数据包 | `Supplementary_Data.zip`，**50 条**（21 面板级源数据 + 3 补充图源数据 + 25 补充表 S9–S15 + 1 清单） |
| Python | **3.13**（`requirements.txt`，89 个 pinned 依赖） |
| R | **4.6.1**（`R_package_versions.tsv`，287 个已装包） |
| 外部可执行程序 | PLINK（`plink.exe`）、SMR **1.3.1**（Windows x86_64）、`tabix` / `bgzip` |

---

## 1 研究一句话

以 **OneK1K 单细胞 cis-eQTL**（980 名供者 / 14 种免疫细胞类型）为暴露、**女性不孕症 GWAS**
（`GCST90483463`，40,024 例 / 665,658 例对照）为结局，用两样本孟德尔随机化框架，
结合共定位、条件分析、精细定位、cis-SMR/HEIDI 与全表型扫描，识别承载因果信号的**免疫细胞类型与区域**。

**三条主结论**（口径见主稿）：

1. **区域层级** —— `chr1p36.12` 与女性不孕症存在因果关联证据（15 个 FDR 显著检验对中 13 个归于该区域）。
2. **细胞类型层级** —— `CDC42` 提供的是**细胞类型特异性的区域标签**，不代表独立因果基因。
3. **基因层级** —— 精细定位指向 `WNT4` 为**优先候选**，但 `WNT4`、`CDC42` **均从未承载结局可信集**；
   受强 LD、单细胞面板缺失与 LD 参考不匹配三重限制，**不能确立因果基因**。

---

## 2 环境

### 2.1 Python

```bash
python -m venv .venv
# Windows: .venv\Scripts\activate
pip install -r requirements.txt
```

`requirements.txt` 由 `pip freeze` 直接生成（Python **3.13.12**）。核心版本：

```
numpy 2.5.3        pandas 2.3.3       scipy 1.18.1      statsmodels 0.15.0
matplotlib 3.11.1  seaborn 0.13.2     lxml 5.4.0        pymupdf 1.28.2
requests 2.34.2    pyarrow 25.0.1     openpyxl 3.1.5    twosamplemr 0.1.3
```

### 2.2 R

R **4.6.1**（本机 `C:\Program Files\R\R-4.6.1`）。包版本全表见 `R_package_versions.tsv`（287 项）。分析直接相关：

```
coloc 5.2.3        susieR 0.14.2      TwoSampleMR 0.7.9  MVMR 0.4.8
MRPRESSO 1.0       RadialMR 1.2.4     MRMix 0.1.0        ieugwasr 1.1.0
data.table 1.18.6.1  ggplot2 4.0.3    qvalue 2.44.0
```

★ **R 调用纪律**（本机实测）：R 脚本必须用 PowerShell 调用；`Rscript` 退出码不可作为成败判据
（加载 `rlang`/`cli` 后退出时段错误，但**不影响计算与产物**）——**判据是产物与日志中的 `Error` 行**。

### 2.3 外部可执行程序

| 程序 | 用途 | 来源 |
|---|---|---|
| `plink.exe` | LD clumping、LD 矩阵、`--extract` 子集 | https://www.cog-genomics.org/plink/ |
| `smr-1.3.1-win.exe` | cis-SMR + HEIDI | https://yanglab.westlake.edu.cn/software/smr/ |
| `tabix` / `bgzip` | FinnGen R12 tabix 端点扫描 | htslib |

---

## 3 目录结构

```
code_release/
├── README.md                      本文件
├── requirements.txt               Python 依赖（pip freeze，89 项）
├── environment.yml                conda 环境（python 3.13 + pip -r requirements.txt）
├── R_package_versions.tsv         R 包版本全表（287 项）
├── LICENSE                        MIT
├── CITATION.cff                   引用信息（含 DOI 与仓库地址）
├── .gitignore
├── repath.py                      路径重定位工具（移植到新机器后运行一次）
├── scripts/                       148 个分析脚本 + 2 个共享模块
│   ├── s01_… – s45_…              主分析链（Python + 5 个 R）
│   ├── s51a_… – s55c_…            三层级扩展分析（Python + 6 个 R）
│   ├── s57_ / s60b_–s60d_         三位点比较表、多效性分解（Python + 1 个 R）
│   ├── s61_… – s166_…             表 / 图件 / 稿件的构建与核验（Python；含 s107–s166 的图件与稿件链）
│   ├── _fig_style.py              **图件全局标准**（字号、字母锚定、轴对齐、间距、点风格、配色）
│   └── _gb12_palette.py           GB12 配色映射的**唯一真源**
└── docs/
    ├── DATA_SOURCES.md            数据集清单（accession / URL / 本地目录 / 体积）
    ├── RUN_ORDER.md               分阶段运行顺序与各脚本输入→输出
    └── SCRIPT_INDEX.md            148 个脚本的自述索引（机械抽取）
```

---

## 4 快速上手

1. **取数**：按 [`docs/DATA_SOURCES.md`](docs/DATA_SOURCES.md) 下载数据集，放入工程根同名目录
   （`00_data_raw/onek1k`、`00_data_raw/scbloodnl`、`00_data_raw/gwas`、`00_data_raw/finngen_r12_*`…）。
2. **重定位路径**：见下方"★ 路径重定位"（一条命令）。
3. **跑主链**：按 [`docs/RUN_ORDER.md`](docs/RUN_ORDER.md) 的 **L1 → L2 → L3 → 敏感性 → PheWAS → 图件** 顺序执行。
4. **出图**：`FIGPAL=gb12 FIG_NOTES=0 python s37_final_figures.py`（重出正式图件时的既定口径）。

### ★ 路径重定位（移植前必须做一次）

脚本按**作者本机的 Windows 绝对路径**编写，为**保证与出图/出表时的脚本逐字节一致**（可复现性），
仓库内不预先改写。克隆后用仓库自带的 `repath.py` 一次性重定位：

```bash
# 1) 先 dry-run（默认不写盘）——原值写死在脚本里，此处给出你机器上的新位置
python repath.py \
    --project  "D:/my/11_sc_eqtl_mr_project" \
    --external "D:/my/endometriosis_project"

# 2) 确认后落盘（自动备份到 _repath_backup_<时间戳>/）
python repath.py \
    --project  "D:/my/11_sc_eqtl_mr_project" \
    --external "D:/my/endometriosis_project" --apply

# 3) 若也想改 docs/ 与根级 *.md 里的示例路径，加 --include-docs
```

**它替换什么**

| 规则 | 原前缀 | 说明 |
|---|---|---|
| 1 `--project` | `D:/endometriosis_project/11_sc_eqtl_mr_project` | 工程根（等价于 `scripts/_fig_style.py` 里的 `ROOT`） |
| 2 `--external` | `D:/endometriosis_project` | 外置数据根，其下含 `_panel`（22 GB）、`_onek1k_plink`（1.4 GB）、`_chr1_ld`（29 MB）、`_tools`（PLINK / SMR） |

**它的保证**（任一不成立即整体中止、**不写任何文件**）

1. **字节级替换**：按 bytes 读写，不编解码、不改换行符 —— 除被替换的前缀外逐字节不变；
2. 命中旧前缀的字符串字面量**全部分隔符统一为 `/`**（避免 `X:/a/b\c\d` 在非 Windows 上失效）；
3. **写盘前自证**：逐文件「整体字节增量 == 已记录替换的实测增量之和」，
   且替换后不得残留任何形态的旧根前缀；
4. `--apply` 先备份到 `_repath_backup_<时间戳>/`（保留相对路径）；
5. 末尾列出**仍未处理**的绝对路径（如 `C:/Users/...`、`E:/...` 之类）及出处文件，供人工确认。

作者本机实测（2026-10-05）：150 个脚本中 **132 个**含硬编码路径；命中 **136（工程根）+ 112（外置根）**；
总字节增量 **−2012**；`py_compile` 138/138 通过；重跑**幂等（0 变更）**。

> **注意**：`_panel` 是外置数据根下的目录名，脚本内实际出现的路径前缀为
> `D:\endometriosis_project\_panel` / `_onek1k_plink` / `_chr1_ld` / `_tools`。
> 重定位后仍需按 [`docs/DATA_SOURCES.md`](docs/DATA_SOURCES.md) 把数据放到对应位置。
>
> 若你的新工程根与外置根**恰好与作者相同**（同盘符同名），可跳过本节；此时脚本可直接运行。

---

## 5 结果表与图件

- **结果表（150 个 CSV）**：中间与最终结果表体积约 51 MB，未随仓库分发。
  表内包含全部发现/复核检验对、共定位后验概率、敏感性矩阵、SMR/HEIDI 输出、PheWAS 全表扫描结果等。
- **图件**：主图 `Fig1`–`Fig5`、补充图 `FigS5`–`FigS8`（`FigS5`=功效门控、`FigS6`=敏感性结局对照、
  `FigS7`=分细胞类型曼哈顿、`FigS8`=发现期效应量分布），各提供 PDF 与 PNG。
- **补充数据包 `Supplementary_Data.zip`（50 条）**：
  21 个面板级源数据 + 3 个补充图源数据 + 25 个补充表 S9–S15 + 1 个 `_manifest.csv`。

---

## 6 复现要点与已知边界

| # | 事项 | 说明 |
|---|---|---|
| 1 | **图件配色** | GB12 为默认调色板（`_gb12_palette.py` 为唯一真源）。重出须 `FIGPAL=gb12 FIG_NOTES=0`；`FIG_OUTDIR` 可把输出重定向到临时目录，避免污染 `figures/`。 |
| 2 | **图件落盘比例** | `_fig_style.save()` 的比例迭代在不可达目标处会振荡 ⇒ 正式图件采用**钉住一个窗内 scale 直接落盘**，且**一个 scale 一个进程**。 |
| 3 | **内存** | 机器可用内存约 3.6 GB；`_fig_overlap_check` 类核验脚本一次只跑单图。 |
| 4 | **R 退出码** | 见 §2.2；判据是产物与日志 `Error` 行。 |
| 5 | **GWAS Catalog API** | 旧 REST `/gwas/rest/api/studies/…` 已弃用并限速（HTTP 429）⇒ 使用 **V2** `/gwas/rest/api/v2/studies/<GCST>`。 |
| 6 | **LD 参考** | OneK1K 980 名供者面板（Zenodo `18870747`）；`chr1p36.12` 区域 r² 矩阵 **354,903 对**。 |
| 7 | **报告口径** | 「未检验 / 不可评估 / 已检验但功效不足」三分互不替代；PheWAS 升高归因**区域**而非基因；「0 条升高」写作「现有功效下未检出」。详见主稿与补充材料。 |

---

## 7 引用

见 [`CITATION.cff`](CITATION.cff)。

- 仓库：<https://github.com/zanzan98/sc-eqtl-mr-female-infertility>
- 归档（**版本 DOI**，引用此版本）：**`10.5281/zenodo.23141009`**（归档版 `v1.0.0`）
- 概念 DOI（恒指最新版）：`10.5281/zenodo.23141008`

> 建议引文（Force 11 数据引用；`[dataset]` 前缀仅供参考文献解析，**出版时移除**）：
> Zan, Y., Lu, Y., Xia, H., Yuan, X., Xia, L., & Xia, Y. (2026). *sc-eqtl-mr-female-infertility* (v1.0.0) [Computer software]. Zenodo. https://doi.org/10.5281/zenodo.23141009

## 8 许可

[MIT](LICENSE)。
