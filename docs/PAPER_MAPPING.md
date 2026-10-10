# 论文图表与代码输出

本文档对应2026-10-09整理的论文1。计算输出的图号是模块编号，投稿图号是排版编号，二者不自动等同。

## 核心结果来源

| 论文内容 | 本仓库计算结果 |
|---|---|
| 人类临床状态系数、SE及阈值跨越 | `human_results.csv.gz`，仅 `clinical_status` 对比 |
| 配对系数变化与HC3敏感性 | `human_paired_changes.csv.gz` |
| 检验家族及分母选择 | `correction_scope.csv`、`human_correction_scope_summary.csv` |
| 小鼠五类模型结果 | `mouse_results.csv.gz`、`mouse_background_fits.csv.gz` |
| 笼位标准误分解 | `mouse_precision.csv.gz` |
| 笼方差、ICC及边界拟合 | `mouse_variance.csv.gz`、`mouse_dream_diagnostics.csv` |
| 四模型设计模拟 | `simulation_summary.csv`、`simulation_failures.csv` |
| 配对差值检验条件校准 | `calibration/calibration_summary.csv` |

表内路径相对于本次工作目录的 `results/`。完整列说明见 `OUTPUTS.md`。

## 图号对应

| 代码生成的图 | 内容 | 投稿稿件中的主要对应内容 |
|---|---|---|
| Figure1_human_effect_precision | 临床系数、q值及配对比值 | 主图1，部分内容用于主图7 |
| Figure2_paired_changes | 所选临床记录的系数差及区间 | 主图1及配对检验补充表 |
| Figure3_mouse_precision | 六种设定的SE因子分布 | 主图4 |
| Figure4_cage_variance | ICC和模型可用性 | 主图3及补充图4 |
| Figure5_design_simulation | 四布局模拟性能 | 主图5、6及补充图5 |
| Figure6_paired_test_calibration | 检验校准、SE/SD及功效 | 主图8及配对校准补充图 |
| Figure7_mouse_testing_families | 小鼠显著集合及校正family | 主图3及相关源表 |
| Figure8_SE_implementation | 直接系数、对比及经验贝叶斯设置 | SE实现核查数据 |

确切输出文件名见 `pipeline/figures.py`。这些图从本次统计输出生成；它们不自动复刻Word中的版式，也不自动重排全部补充图。不能按编号直接覆盖投稿目录图片。

## 当前运行的核对点

下列数值来自最后通过核查的运行，用于帮助读者确认读表范围，不作为选择基因或调整参数的条件：

- 人类主家族每模型臂17,064个计划位置，其中11,232个可检验。
- 临床状态对比有3,744个可检验位置，联合校正下基础40条显著、质量调整后0条。
- 这40条记录的系数绝对值配对比中位数约0.928，SE配对比中位数约1.019。
- 小鼠共同背景为13,312基因，六种设定共79,872条SE分解记录。
- 当前小鼠背景基因SE比的各设定中位数范围为1.495–1.745，不是1.572–1.948。
- 两套dream权重有159,744条方差组分记录；最后运行无未解决dream警告。

`pipeline/audit.py` 内的 `manuscript_comparison.csv` 仍是早期稿件锚点对照，六个小鼠SE参照值不等于当前稿件值。保留这一分析源码是为了与实际完成运行严格一致；锚点仅用于附加比较，不影响计算、校正或退出码的统计判定。

## 不由此入口重算的内容

补充信息中其他队列的独立探索结果、额外的广义细胞类别混合实验，以及投稿版本的全部表格命名和图件排版，不由本入口完整生成。该范围应随论文补充方法和数据说明一起阅读，不能将本包称为所有补充研究的单一重算入口。

重新下载公开资源时，文件或同源注释可能随提供方更新；本地记录URL、下载时间、大小和SHA-256。已有数据完整且带收据时可以复用，不需要为了重跑模型重新下载。
