# 数据集清单（DATA SOURCES）

> 全部数据均为**公开汇总级数据**或**公开单细胞资源**；本研究**未访问任何个体级数据、未招募新受试者**。
> 下表 URL / accession 均自脚本内常量与实际下载记录实查；「本地目录」为作者机器上的落地路径（**按需重建**）。

---

## 0 一览

| # | 数据集 | 版本 / accession | 用途 | 本地目录 | 体积 |
|---|---|---|---|---|---|
| 1 | OneK1K 单细胞 cis-eQTL | TensorQTL 重算原始汇总；980 名供者 × 14 种免疫细胞类型（GRCh37） | **暴露**（主分析） | `00_data_raw/onek1k/` | **14 GB** |
| 2 | 1M-scBloodNL | `eqtls_20201106_genome_wide`（UT 细胞类型） | 暴露（血液细胞类型复核） | `00_data_raw/scbloodnl/` | **15 GB** |
| 3 | 女性不孕症 GWAS（主结局） | **`GCST90483463`** | **结局**（发现期） | `00_data_raw/gwas/GCST90483463.h.tsv.gz` | 0.95 GB |
| 4 | 女性不孕症 GWAS（敏感性结局） | **`GCST90483469`** | 结局（敏感性） | `00_data_raw/gwas/GCST90483469.h.tsv.gz` | 1.13 GB |
| 5 | 类风湿关节炎 GWAS（阳性对照） | **`GCST90132222`** | **阳性对照门控**（`s06_mr_extract.py`） | `00_data_raw/gwas/GCST90132222.h.tsv.gz` | 0.60 GB |
| 6 | GTEx v8 全血 eQTL | `GTEx_Analysis_v8_eQTL_all_associations`（Whole_Blood） | `WNT4` 的 **cis-SMR 替代层**（血液面板缺失时的替代） | `00_data_raw/smr_3p3/gtex8*` | — |
| 7 | GTEx v10 单组织 cis-QTL | `GTEx_Analysis_v10_eQTL_all_associations` | 多组织 eQTL 注释 | `00_data_raw/smr_L2L3/gtex8*` | — |
| 8 | eQTLGen cis-eQTL | `2019-12-11-cis-eQTLsFDR0.05-ProbeLevel`（血液） | 血液 cis-eQTL 对照注释 | 按需 | ~1 GB |
| 9 | FinnGen R12（DF12） | 数据冻结 12，2024-11-04 公开；PheWeb + tabix 汇总统计 | **PheWAS 端点库**（2,469 端点） | `00_data_raw/finngen_r12_pheweb/`、`finngen_r12_tabix/` | 3.9 MB + 1.2 MB |
| 10 | Open Targets Platform | GraphQL **v4** | `chr1p36.12` 关联谱 / 安全性谱 | （API 缓存） | — |
| 11 | OneK1K PLINK 基因型面板 | Zenodo **`18870747`**，980 名供者 | LD 参考 / clumping / r² | `_onek1k_plink/`、`_panel/` | **1.4 GB + 22 GB** |
| 12 | chr1 区域 r² 矩阵 | 由面板算出，**354,903 对** | 精细定位 / LD 核验 | `_chr1_ld/` | 29 MB |

---

## 1 暴露：单细胞 cis-eQTL

### 1.1 OneK1K（主分析）

- **内容**：980 名供者的 14 种免疫细胞类型 cis-eQTL；本分析使用 TensorQTL 重算版本。
- **落地件**：`00_data_raw/onek1k/OneK1K_TensorQTL_raw_eQTL_summary.tar.gz`（10.3 GB 压缩包）。
- **引用**：Yazar S, et al. *Science* 2022（OneK1K 资源）。
- **下游工作目录**（脚本自动建立）：`coloc_prep/`、`coloc_work/`、`cond_coloc/`、`finemap_work/`
  及其 L2/L3 版本 `coloc_sens_prep/`、`cond_coloc_L2L3/`、`finemap_work_L2L3/`。

### 1.2 1M-scBloodNL（复核）

- **下载**：`https://molgenis26.gcc.rug.nl/downloads/1m-scbloodnl/eqtls_20201106_genome_wide.tar.gz`
  （13.6 GB；另有较小的 `eqtls_genome_wide.tar.gz` 43 MB）。
- **用途**：UT 细胞类型的独立血液细胞类型对照。

---

## 2 结局：女性不孕症 GWAS

### 2.1 `GCST90483463`（主结局）

| 项 | 值 |
|---|---|
| 性状 | Female infertility (all causes) |
| 样本 | **40,024 例欧洲血统女性病例 / 665,658 例对照** |
| PMID | **40229599** |
| 变异数 | 27,311,379 |
| 许可 | **CC0** |
| MONDO | MONDO_0021124 |
| FTP | `https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90483001-GCST90484000/GCST90483463/harmonised/` |
| API（V2） | `https://www.ebi.ac.uk/gwas/rest/api/v2/studies/GCST90483463` |

### 2.2 `GCST90483469`（敏感性结局）

- Female infertility (all causes)；**42,629 例多血统病例 / 740,619 例对照**；同一 PMID。
- 本地另有 `.tbi` 索引（1.8 MB），供区域抽取。

### 2.3 `GCST90132222`（阳性对照）

- **类风湿关节炎**；用于 `s06_mr_extract.py` 的**阳性对照门控**
  （传入已知阳性基因 `GLIPR1, XBP1, IKZF1, IL2RA, …` 验证 MR 管道灵敏度），**不进入主结局分析**。

> ★ **GWAS Catalog API 纪律**：旧 REST `/gwas/rest/api/studies/…` 已弃用并限速（HTTP 429）
> ⇒ 一律使用 **V2** `/gwas/rest/api/v2/studies/<GCST>` 或网页版 `https://www.ebi.ac.uk/gwas/studies/<GCST>`。

---

## 3 组织 eQTL（SMR 层）

| 资源 | 下载地址 | 用途 |
|---|---|---|
| GTEx v8 全血 | `https://storage.googleapis.com/adult-gtex/bulk-qtl/v8/single-tissue-cis-qtl/GTEx_Analysis_v8_eQTL_all_associations.tar` | `WNT4` cis-SMR **替代层**（`s43e`–`s43i`、`s55a`–`s55c`） |
| GTEx v8（EBI eQTL Catalogue） | `https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/imported/GTEx_V8/ge/Whole_Blood.tsv.gz` | 同上（流式抽取） |
| GTEx v10 | `https://storage.googleapis.com/adult-gtex/bulk-qtl/v10/single-tissue-cis-qtl/GTEx_Analysis_v10_eQTL_all_associations.tar` | 多组织 eQTL 注释 |
| GTEx Portal API | `https://gtexportal.org/api/v2/association/singleTissueEqtl`、`/metadata/dataset`、`/metadata/tissueInfo` | 单组织 eQTL 查询 |
| eQTLGen | `https://download.eqtlgen.org/cis-eqtl/2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz` | 血液 cis-eQTL 对照 |
| EBI eQTL Catalogue API | `https://www.ebi.ac.uk/eqtl/api/v2/datasets/?quant_method=ge` | 数据集元数据 |

**SMR 输入构建链**：`s43a_build_smr_inputs` → `s43b_build_besd`（BESD/ESI/EPI）→
`s43c_run_smr`（SMR 1.3.1 + PLINK LD）→ `s43d_smr_summary`；GTEx v8 替代层为 `s43e`–`s43i`。

---

## 4 PheWAS 端点库：FinnGen R12

| 项 | 值 |
|---|---|
| 版本 | **DF12 = R12**，公开于 **2024-11-04** |
| 规模 | 500,348 样本 / 2,502 端点（**本文使用 2,469 个疾病端点**） |
| 浏览器 | `https://r12.finngen.fi/` |
| API | `https://r12.finngen.fi/api/phenos`、`/api/pheno/<PHENO>`、`/api/variant/<CHR>-<POS>-<REF>-<ALT>`、`/api/region/?chrom=&start=&end=` |
| 汇总统计桶 | `https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release/` |
| 单端点 gz | `https://r12.finngen.fi/pheno_gz/<PHENO>.gz` |
| 推荐引用 | Kurki MI, et al. *Nature* 2023;**613**(7944):508–18. doi:10.1038/s41586-022-05473-8 |
| 致谢文本 | *We want to acknowledge the participants and investigators of the FinnGen study.* |

- 两条抓取路线：`s41a_fetch_pheweb`（PheWeb API）与 `s42b_phewas_tabix_fetch`（**tabix 全端点扫描**，产出 `phewas_3snp_all_endpoints.csv`）。
- 3 个工具的**全部端点 × 2 个工具 × 4 个分析臂**结果集已完整收入 `Supplementary_Data.zip`（表 **S9a**）。

---

## 5 LD 参考与基因型面板

| 项 | 值 |
|---|---|
| 来源 | OneK1K 980 名供者合并基因型面板 |
| 下载 | `https://zenodo.org/records/18870747/files/plink_genotype_merged_980_donors.zip?download=1` |
| 落地 | `_onek1k_plink/`（1.4 GB）、`_panel/`（22 GB，含分染色体面板） |
| 区域 r² | `_chr1_ld/chr1_r2.ld`（**354,904 行 = 1 表头 + 354,903 对**） |
| 工具 | `_tools/plink.exe`（LD clumping / `--extract` / r² 计算）、`_tools/smr-1.3.1-win-x86_64/`（SMR 二进制） |

★ **LD 实查值**（正文引用口径）：
`r²(rs10917151, rs56318008) = 0.84732`；`r²(rs12037376, rs56318008) = 0.899359`；`0.934748` 属 `rs61768001`。

---

## 6 其他在线服务

| 服务 | 端点 | 用途 |
|---|---|---|
| Open Targets Platform | `https://api.platform.opentargets.org/api/v4/graphql` | 区域关联谱 / 安全性谱（`s35`、`s35b`） |
| Ensembl REST | `https://rest.ensembl.org/lookup/symbol/homo_sapiens/<gene>`、`/variation/human/<rsid>` | 基因符号 / 变异坐标与链向 |
| Europe PMC | `https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=` | 文献核验 |
| Crossref | `https://api.crossref.org/works/<doi>` | 文献元数据核验 |

★ **文献纪律**：正文与参考文献中的文献条目**逐条联网实查**；无法实查者以占位符保留，**严禁编造**。

---

## 7 未随仓库分发的内容

| 项 | 原因 |
|---|---|
| `00_data_raw/`（**31 GB**） | 可从上述公共来源重新下载 |
| `_panel/` + `_onek1k_plink/` + `_chr1_ld/`（**约 23.4 GB**） | 同上（Zenodo 面板 + 本地算出） |
| `tables/` 结果表 **150 个 CSV**（约 51 MB） | 可由脚本链重算；随 `Supplementary_Data.zip` 提供投稿所需子集 |
| `figures/` 全部图件（约 171 MB，含历史版本） | 正式图件可由 `s37_final_figures.py` + `s67_make_supp_figs.py` 重出 |
