# 论文1：从公开输入完整重算

本仓库复现论文 Distinguishing coefficient changes from precision costs in brain transcriptome analysis 的主分析与补充分析。统一入口从公开表达对象和元数据生成分析输入、拟合、模拟、图件和完整数值表。

## 运行

**已完成实验只需重新导出时，不要启动全量复跑。** 使用：

```bash
python paper1_reexport_results.py --run /home/research/git/repository/paper1_complete_from_public_20261022_v2
```

入口只检查原运行记录中的 NumPy、pandas、SciPy、Matplotlib，不要求 openpyxl 或 R。
原实验目录保持只读，新建 `p1_review_年月日_时分秒` 短目录；出现
`COLLECTION_COMPLETE.json` 后将整个新目录带回核对，不生成 ZIP。
导出采用论文定义下的 ICC/dream SE 比、native 精度分解、完整 BH 家族与补图编号。
结果按已完成实验原样导出。详见 [导出说明](docs/REEXPORT.md)。

S4.2 使用全部可分解临床记录，显著子集另计；S9 交互项不作疾病主效应方向比较；
S19.4 分列五模型与三模型的实际布局归属；S15.1 分列绝对差和有符号差。
投稿表以 `reporting/block_sources.csv` 指向的文件为准，原运行副本保持不变。
Windows 建议放在 `C:\p1\` 等短父目录，不要重复嵌套同名结果目录。

使用已有 Linux、Python 3.12、R 4.5.3 环境。启动器不安装、不下载、不升级软件包。

```bash
conda activate paper1
cd /home/research/git/repository/paper1-from-publicV1
python reproduce_from_public.py --check-env
python reproduce_from_public.py --workers 60 \
  --output /home/research/git/repository/paper1_complete_from_public_v1
```

新输出目录必须为空。默认执行主分析和全部补充分析，不需要添加模块开关。只有明确使用 `--primary-only` 才会省略补充分析。

便捷入口为 `bash run_full_reproduction.sh`，默认 60 进程，输出到代码目录同级的
paper1_complete_from_public_20261022，并自动检查同级已有的论文1原始下载缓存。
加 `--check-env` 只检查环境，加 `--resume` 继续同一次实验。
`--check-env` 同时检查 `config.json` 与代码发行版本是否一致；正式运行也在建立输出目录前执行该检查。源码与配置应整套同步，不要单独改版本号绕过不匹配。
可用环境变量 PAPER1_OUTPUT 和 PAPER1_WORKERS 指定输出及进程数。

允许复用的只有已下载的公开原始对象，需有来源 URL、文件大小及 SHA256 下载记录。添加 `--raw-cache /已有下载目录或工作目录` 即可，参数可重复。程序不会从那里导入发现面板、预处理计数、细胞标签、拟合结果或模拟批次。缺失的公开数据会下载。

这里的“从头”从公开计数对象、作者存储的芯片表达矩阵及来源元数据开始，不包含 FASTQ 比对或 CEL 重新归一化。由这些输入重新生成面板、聚合矩阵、标签和全部分析结果。

## 环境

中心分析精确版本见 requirements.txt 和 inputs/environment.json。补充细胞分析还需要 scikit-learn、threadpoolctl、scrublet；版本按实际安装记录并在同一实验中锁定。

可用 `--python` 指定已有 Python，重复 `--r-library` 合并已有 R 库。选定 Python 缺少附加包时，可用 `--python-package-path /已有/site-packages` 暴露已安装的附加包；数值核心版本及导入必须通过检查。复用软件不等于复用分析结果。

启动器逐一检查当前解释器、仓库及其 `paper1-from-public` 子目录中的虚拟环境，以及同级 `paper1_public_evidence` 下的已有源码虚拟环境。Python 符号链接保留虚拟环境内的调用路径，不解析为底层 Conda 解释器。当前解释器显示 `PYTHON NOT SELECTED` 后仍会继续检查其他候选；也可显式指定 `--python "$PWD/.raw-reproduction-venv/bin/python"`。未找到合格环境时停止，不自动安装或放宽版本要求。

60 进程用于独立模拟任务；单细胞数据准备与注释不会盲目同时占用 60 份大型矩阵。工作进程内部限制 BLAS/R 并行，不使用 GPU。完整计数流程可能需要数天，以实际进度为准。

## 查看和续跑

```bash
python run.py --status --output /home/research/git/repository/paper1_complete_from_public_v1
python reproduce_from_public.py --workers 60 --resume \
  --output /home/research/git/repository/paper1_complete_from_public_v1
```

续跑保留首次运行的缓存及环境参数，只接着使用同一次实验中校验通过的完成任务。未完成的补充模块在新的 attempt 目录重算；独立模拟场景有单独的完成记录。改变科学代码、参数或软件后必须另开输出目录，不能强行续跑。进程数可调整。不要删除锁文件强行启动第二个写入进程。

`python run.py --plan` 只列运行步骤；`python run.py --coverage` 只检查代码覆盖范围，不代表数值验证成功。

## 结果与收集

默认流程覆盖：
- 发现面板、人脑配对模型、检验家族、小鼠笼位模型和标准误分解。
- 人脑八套敏感性、小鼠七臂比较、1,600 个完整计数流程数据集。
- 固定设计校准、五模型笼位基准。
- 年龄与质量敏感性、外部队列、认知关联。
- 从公开计数重新注释细胞、doublet 敏感性、六类细胞交互及合成混合。
- 双模型条件模拟、三模型布局模拟和实际动物布局。

结果导出为 reporting/S01 至 S20 下的 CSV，block_sources.csv 对应补充工作簿的112个内容区块。另直接生成 reporting/Supplementary_Tables.xlsx，含20张表、112个区块及S01表头行目录；workbook_manifest.json记录工作簿校验值和区块位置。重复展示共用数值来源。

工作簿纳入S6.6人类情境性别计数、S10.2 CAMERA，以及S20.13–S20.17的数据集级计数结果、配对比较、不可用位置、生成参数和完整模拟任务。S1.2使用公共数据入口与输入版本说明。参数表中的mouse_units.csv指流程生成的输入；机器目录只在工作簿显示中省略，原始协议和下载凭据保留。超长任务数组完整展开至S20.17，不截断。

已归档v2.1.0的默认导出仍为105块。本次源码更新的导出版本为 publication_exports_20261026_v3，可直接生成112块，且兼容对旧105块完整运行重新导出。此次不修改已发布标签或DOI，不更改统计分析。

已有完整结果只需运行 `python paper1_reexport_results.py --run /已有完整结果目录`。脚本读取既有结果并生成新的短路径目录，不重新拟合、不重新下载，也不安装软件。

S18 的模拟保持预处理小鼠表的行顺序，并记录 simulation_animal_order.csv 和 simulation_input_order.json，包括动物与随机抽样位置的对应关系、输入哈希及随机种子。相同种子只有在样本顺序、输入、模型和软件一致时才指定同一模拟。程序不为对齐某组结果重排动物。两个记录文件不改变科学计算；已有完整运行仍可由其预处理元数据和 layout_units.csv 核实。

另外生成八张主结果图和十二张补充数值展示图。它们用于核对当前数据，不保证与手工排版图逐像素一致；不自动改Word。Excel由当前结果自动组装，不读取作者手工整理的投稿工作簿。

跑完后只收集最终表格、图和必要环境记录：

```bash
python collect_results.py \
  --run /home/research/git/repository/paper1_complete_from_public_v1 \
  --output /home/research/git/repository/paper1_complete_results_for_review
```

不打 ZIP，不复制原始计数、细胞级标签或模拟任务目录。只有完整分析和导出均完成后才允许收集。

运行测试：`python -m unittest discover -s tests -v`。合成数据测试通过不能替代服务器全量科学重跑；新结果应按基因、模型、供体和检验家族逐项与论文比较。

[数据定义](docs/DATA.md) | [方法定义](docs/METHODS.md) | [输出说明](docs/OUTPUTS.md)
