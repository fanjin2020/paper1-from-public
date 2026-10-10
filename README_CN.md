# AD转录组模型调整分析

[English](README.md)

Linux分析代码，版本 **1.2.2**。

本仓库比较组织质量调整和笼位建模对转录组效应估计、标准误与显著关联集合的影响。人类分析使用SEA-AD单核计数，小鼠分析使用GSE233208空间表达；GSE118553用于建立候选基因面板，MGI用于一对一同源映射。

完整入口为 `run.py`。代码包含公开数据下载、计数聚合、统计模型、配对系数检验、SE分解、dream收敛复核、条件模拟和结果绘图，不是只读取整理结果的绘图包。函数和类的说明为中文。

## 目录

```text
paper1-from-public/
  run.py                    完整分析入口
  run_dream_repair.py        已有运行中dream警告的定点复核
  config.json               样本条件、统计设置和随机种子
  environment.yml           Linux Conda环境
  requirements.txt          Python依赖
  analysis/                 R数据处理、模型和模拟
  pipeline/                 Python下载、聚合、校准、核查和绘图
  tests/                    Python与R回归测试
  docs/                     数据、方法、运行、输出和发布说明
  CITATION.cff              软件引用信息
  release_manifest.json     本包文件摘要
```

仓库不包含原始数据、聚合矩阵、运行结果、安装后的R包、日志或临时文件。

## 第一次运行

在Linux中解压后进入含 `run.py` 的目录：

```bash
conda env create -f environment.yml
conda activate paper1-public
set -euo pipefail

python -u run.py \
  --task install-r \
  --work "$PWD/work_public" \
  --rscript "$CONDA_PREFIX/bin/Rscript" \
  --library "$PWD/R_library"

python -u run.py \
  --task all \
  --work "$PWD/work_public" \
  --rscript "$CONDA_PREFIX/bin/Rscript" \
  --library "$PWD/R_library" \
  --workers 1 --resume \
  2>&1 | tee run.log
```

依赖安装失败时，`set -euo pipefail` 会停止后续分析。已有Python环境可以按 `requirements.txt` 检查依赖，无需重复创建环境。R需为4.5.x，Bioconductor系列为3.22；运行时记录实际包版本。

首次下载需要较大的磁盘和内存。这里的起点是作者公开的UMI/counts层与芯片处理矩阵，不包含FASTQ比对或芯片原始强度处理。

## 已有数据，不重新下载

使用原分析的完整工作目录，不是删去大文件后的 `review_bundle`。先把 `DATA_WORK` 和 `R_LIB` 改为本机真实路径：

```bash
conda activate paper1
set -euo pipefail

DATA_WORK="/absolute/path/to/existing/work"
R_LIB="/absolute/path/to/existing/R_library"

python -u run.py \
  --task all \
  --data-work "$DATA_WORK" \
  --work "$PWD/work_public_new" \
  --rscript "$CONDA_PREFIX/bin/Rscript" \
  --library "$R_LIB" \
  --workers 1 --resume \
  2>&1 | tee run.log
```

数据复用顺序：

1. 已有聚合计数及其元数据通过SHA-256和来源记录核验，直接进入模型。
2. 聚合输入不齐但原始缓存完整，读取缓存重新聚合。
3. 只下载缺失的原始文件。

不覆盖已有数据。文件存在但缺少完整性记录，或摘要不一致时会停止核查，不会静默采用。严格禁止联网下载的重新分析可将 `--task all` 改为 `--task reanalyse`，此时四个准备阶段必须完整。

新结果目录与已有数据目录不能相同，也不能互相嵌套。中断后使用相同代码、参数及工作目录重跑并保留 `--resume`。不要删除来源收据或阶段记录。

## 检查与结果

不读取真实数据的Python测试：

```bash
python -m unittest discover -s tests -v
```

完整运行先检查依赖，再执行R数值回归和dream官方API小样本测试。分析期间不自动安装软件。

| 输出 | 用途 |
|---|---|
| `prepared/` | 聚合计数、样本纳入清单和设计矩阵 |
| `results/` | 效应、SE、p/q值、方差组分、校准和模拟CSV |
| `figures/` | 八组分析图，含PDF、600 dpi PNG及同名CSV |
| `status.json` | 请求的任务状态和结果核查摘要 |
| `review_bundle.zip` | 便于审阅的小型结果包，排除大计数对象 |

退出码0表示所请求任务完成，2表示输入、依赖或计算检查失败，3表示流程结束但仍有拟合失败或未解决警告。纯基因型笼中固定效应不可估计单独记录，不与数值失败混合。警告经两优化器复核后通过的拟合可保留。

## 与论文的对应

本发布包的32个分析、配置和测试文件与最后完成运行时保存的执行副本逐字节一致。此次只整理文档与GitHub打包，没有改动拟合、筛选、校正、种子或优化器规则。

代码图号按分析模块命名，与投稿排版图号不同，见 [图表对应说明](docs/PAPER_MAPPING.md)。本仓库覆盖论文的核心公开数据分析；补充材料中的其他队列及其独立分析不在这一入口内。

`results/manuscript_comparison.csv` 是保留在原执行代码中的早期稿件数值对照，不能用它判断2026-10-09稿件是否匹配。尤其其中六个小鼠SE参照值已不适用于当前稿件。当前结果应按图表对应说明直接读取计算结果；这些参照常数不参与拟合或显著性判定。

## 文档

- [数据来源、下载地址和处理规则](docs/DATA_SOURCES.md)
- [模型、检验家族与模拟定义](docs/METHODS.md)
- [文件和函数说明](docs/CODE_GUIDE.md)
- [结果表及列含义](docs/OUTPUTS.md)
- [论文图表对应与复现范围](docs/PAPER_MAPPING.md)
- [dream定点复核](docs/DREAM_REPAIR.md)
- [版本验证](docs/VALIDATION.md)
- [上传GitHub](docs/GITHUB_UPLOAD.md)

数据分别遵循原提供方使用条款。软件许可由作者确定；未指定许可证前，公开代码不等于授予任意再分发或修改许可。本包没有虚构仓库地址、DOI或论文发表信息。
