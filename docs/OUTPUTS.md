# 输出与论文逐项核对

## 主要结果表

| 文件 | 一行表示什么 | 主要列 |
|---|---|---|
| panel.csv | 本次额叶选中基因 | gene_symbol、beta、p、q、n |
| ortholog_panel.csv | 人类面板端点及其小鼠对应 | gene_symbol、mouse_gene |
| human_results.csv.gz | 一个基因×亚类×对比×模型臂 | beta、SE、df、p、q_joint、区间、status |
| human_paired_changes.csv.gz | 同一端点的质量调整差值 | delta、SE_delta、p、q_joint、SE_HC3、q_HC3_joint |
| human_designs.csv | 一个亚类×对比 | n、设计列数、可估计状态 |
| human_gene_membership.csv | 一个亚类×面板基因 | 输入是否有该符号、过滤后是否保留、匹配规则 |
| correction_scope.csv | 一个明确校正family×分母策略 | planned、available、denominator、significant |
| human_correction_scope_summary.csv | 模型臂×范围×分母策略 | 临床状态可用数、显著数、基础临床集合保留数 |
| mouse_results.csv.gz | 一个面板端点×规格×模型 | beta、SE、df、p、q_base_fixed_joint、q_three_model_joint、q_dream_joint、rho |
| mouse_precision.csv.gz | 一个背景基因×规格 | SE_ratio、design_factor、posterior_factor、reweight_factor |
| mouse_variance.csv.gz | 一个背景基因×规格×权重策略 | 笼方差、残差尺度、ICC、同笼动物对平均相关、边界与收敛 |
| mouse_SE_implementation_audit.csv.gz | 一个背景基因×规格×主模型 | 精确QR、直接系数、精确对比、contrastAsCoef、近似对比及两类收缩SE/p |
| mouse_background_fits.csv.gz | 一个背景基因×规格×模型 | 全背景效应、区间和可用状态 |
| mouse_model_status.csv | 一个规格×模型 | 计划数、可用数、失败数 |
| mouse_fit_failures.csv | 一个失败的背景拟合 | 基因、规格、模型和状态 |
| mouse_model_messages.csv | 一条模型错误或警告 | 规格、模型、原始消息 |
| mouse_dream_diagnostics.csv | 一个问题基因或一条未定位警告 | gene_symbol、spec、model、status、message |
| simulation_failures.csv | 一条数值失败或结构性不可估计 | 完整情境、重复编号、原因 |
| result_audit.json | 整体结果核查 | 分母检查、失败数、旧稿差异数 |
| manuscript_comparison.csv | 一个旧稿关键数值 | 旧值、新值、比较容差和是否一致 |
| mouse_designs.csv | 小鼠六种设定 | 动物数、笼数、混合笼数、信息比例 |
| simulation_summary.csv | 布局×相关×偏移×效应×模型 | 偏差、SD、SE、覆盖、拒绝、宽度、失败率 |
| calibration/calibration_summary.csv | 配对检验校准情境 | 同上及质量项F拒绝率 |

`simulation_replicates.csv.gz`保留逐次结果便于审计，不装进审阅ZIP。阶段计数和下载文件保留在本地，GitHub仓库不携带它们。

`mouse_dream_diagnostics.csv`只有表头表示未捕获需定位的数值警告。基因名为空表示警告尚未定位，不等于不存在问题。拟合失败或数值警告都不能作为“阴性基因”计数。

## 图表

1. Figure1：临床状态端点系数；该对比基础显著集合的q值变化；效应比和SE比分布。
2. Figure2：同一临床状态显著集合的配对系数差与95%区间，不混入条件临床对比。
3. Figure3：小鼠六种设定的SE三因子经验分布。
4. Figure4：共享/原生两套权重的可用基因ICC分布；CSV同时保留可用、失败和边界计数。
5. Figure5：不同布局和模型的偏差、拒绝率、宽度、覆盖与失败。
6. Figure6：配对检验经验SD与估计SE、零假设拒绝率及两SD备择功效。
7. Figure7：小鼠方法内与两模型联合校正计数，分开展示校正family。
8. Figure8：相同输入下，对比变换和经验贝叶斯设置各自带来的SE比。

每幅图同名CSV提供实际绘图及其相关明细；大表压缩保存在results，发表时可据需要只提供图示汇总。这里的图号为新代码图号，**不是自动覆盖当前Word中的图号**。

## 重新分析完成后如何核对论文

先核对来源和单位，再看p值：

1. 新面板基因数及基因身份是否与稿件一致。
2. SEA-AD供体、亚类、各模型完整病例数量是否一致；质量变量定义是否相同。
3. 小鼠纳入/排除原因、动物ID、cage、slide及皮层标签是否一致。
4. 共同背景基因数、TMM、voom、经验贝叶斯设置与family是否一致。
5. 最后比较系数、SE、q值、区间及论文计数。

发现不一致时应先输出差异解释，不能把新结果覆盖成旧数值。尤其不能用“都不显著”代替数值复现，也不能因需要40条记录就改BH范围。

## 不在本版中的内容

其他探索队列、旧版18个规格的候选图、全部旧补充图的版式及历史布局模拟不在本次原样重绘范围。本包覆盖新的核心分析链，并不声称与现有Word每个补充数字一一对应。研究者完成核对后，才能据新结果决定是否更新论文和补充材料。
