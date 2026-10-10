# 文件与函数说明

## 唯一入口：run.py

`main`解析参数、检查Linux环境、取得目录锁并按依赖顺序执行。`fingerprint`记录代码与参数；`stage`核验完整阶段的输出后恢复；`rrun`启动独立R进程；`bundle`打包审阅所需的小型结果。

`prepare_human`调用SEA-AD聚合后连接供体质量表；`prepare_mouse`重建动物身份并读取Seurat counts。它们不会接受旧结果路径作为备用输入。

## Python模块

| 文件 | 主要函数 | 输入和输出 |
|---|---|---|
| pipeline/io.py | download、sha、dump、command | HTTP文件、完整性记录、子进程日志 |
| pipeline/sources.py | fetch_small | 官方URL到当前raw目录 |
| 同上 | mouse_manifest | GEO与S2到动物身份、资格表 |
| 同上 | unwrap_mouse | 多层gzip到R可读取对象 |
| 同上 | quality_join | 官方供体表到伪批量临床/质量元数据 |
| 同上 | read_geo_matrix、prepare_panel | GEO样本特征和芯片表到供体×脑区表达 |
| 同上 | map_orthologs | MGI与新面板到一对一映射及缺失记录 |
| pipeline/seaad.py | sea_inventory、prepare_seaad | 指定S3发行版到供体亚类聚合 |
| 同上 | h5_column、count_sum_sparse | HDF5分类字段解码、稀疏计数求和 |
| 同上 | collapse_genes、combine | 重复符号合并、供体矩阵连接 |
| pipeline/calibration.py | geometry、evaluate | 嵌套设计到差值、SE、F检验 |
| 同上 | mean_profiles、run | 已知均值条件模拟到性能汇总 |
| pipeline/figures.py | generate、save、ecdf | 新统计结果到图示CSV、PDF、PNG |
| pipeline/reuse.py | collect、import_inputs | 收据及SHA-256核验，复用已有计数，不导入旧模型 |
| pipeline/simulation_design.py | synthetic、matrices、prepare | 精确布局、设计空间、32情境和逐次种子 |
| pipeline/audit.py | clinical_pairs、assess | 按对比汇总、分母闭合、失败与旧稿数字核查 |

每个函数都有中文说明。`seaad.py`中的检查点只保存本次已下载供体的聚合计算；不附带已有计数或结果。

## R模块

| 文件 | 用途 |
|---|---|
| analysis/preflight.R | 软件安装与独立进程验证；已有依赖完整时不请求Bioconductor在线版本检查 |
| analysis/common.R | 计数合法性、BH位置分母、WLS差值、模型结果和方差组分提取 |
| analysis/diagnostics.R | 精确加权QR、四组合SE核查、dream官方区间核对与失败记录 |
| analysis/panel.R | 芯片供体级模型与额叶面板重建 |
| analysis/prepare_mouse.R | 验证Seurat身份、读取counts、全脑/皮层聚合 |
| analysis/models.R | 人类两臂、配对变化、小鼠五模型、SE分解和笼方差 |
| analysis/simulation.R | 新实际布局与合成布局的Gaussian四模型比较 |

`common.R`不是单独运行入口，其他R文件也由run.py传入参数。直接双击R文件不能替代主程序。

## 测试

`tests/test_calibration.py`使用独立最小二乘表达核对系数、协方差、HC3、退化情况和已知模拟均值。`tests/test_inputs.py`检查样本编号、基因聚合、GEO解析、坏gzip、缓存摘要、部分HTTP响应和说明文档的存在。

Python测试不依赖真实样本，可以先执行：

```bash
python -m unittest discover -s tests -v
```

R测试在`--task check`中执行。若R冒烟测试失败，不能继续实证分析。

`test_reanalysis.py`检查不联网复用、被修改的输入、路径越界、模型参数与准备参数的区分、重复端点、四布局身份和执行代码打包。`test_analysis.R`要求limma和statmod齐备，用正式EList及显式矩阵加权接口交叉检查，再以base R的加权最小二乘核对QR；缺少依赖时不会跳过。它也检查普通列表、零权重、维度错误和秩亏输入。

## 修改约定

新增模型须同时写清估计对象、输入单位、family和失败分母。不要只在绘图函数改数字。修改下载URL时保留资料版本说明，使用新work目录。改变作者细胞注释、年龄处理或动物身份规则属于分析变更，不能作为普通文件路径修复悄悄引入。
## 缓存与dream复核

- `pipeline/cache.py`：完整入口的数据检查与复用选择。先用合格计数，其次用已下载原始文件，最后下载缺失项。
- `run_dream_repair.py`：可选的定点复核入口。核验已有完整运行，只重算有警告的dream模型并合并结果。
- `analysis/repair_dream.R`：重建受影响模型的权重，核对默认结果复现，再执行优化器复核和共同收缩。
- `analysis/diagnostics.R`：新增警告收集、原始估计比较、两优化器复核及背景纳入/排除对照。
- `tests/test_cache.py`、`tests/test_dream_repair.py`、`tests/test_dream_retry.R`：分别验证缓存复用、结果合并和数值判据；最后一个文件还包含服务器正式API冒烟测试。

