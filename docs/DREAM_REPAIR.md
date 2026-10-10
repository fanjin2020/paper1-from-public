# v1.2.2：dream收敛警告的定点复核

## 作用范围

v1.2.2只替换通过两优化器复核的原始拟合行，正常基因保持原始估计。负特征值警告须定位并核查，不能直接忽略。

定点入口从原结果表读取有警告的任务，不预先指定基因。目前只接受v1.2.0完整工作目录，不接受缺少计数的审阅ZIP。新分析直接使用 `run.py` 即可，其中已包含同一套复核逻辑。已完成的v1.2.2结果无需再用这个旧版本迁移入口。

## 复核顺序

1. 核验原结果、计数和动物清单的SHA-256；新旧工作目录分开，来源文件不改动。
2. 检查R依赖及版本。先执行回归测试和官方dream接口冒烟测试，成功后才分析真实计数。
3. 按原背景、样本顺序、设计、TMM和voom设置重建有警告模型的输入。
4. 重现默认拟合。原可用估计和缺失位置必须在规定容差内一致，否则停止合并。
5. 对已定位基因分别使用 `bobyqa` 和 `nloptwrap`。两个优化器均不许有未解决警告或后台错误；目标系数、原始SE、自由度、似然、残差尺度和笼位方差均须一致。
6. 通过后，只用预先指定的bobyqa替换该基因的原始拟合行，包括系数、残差尺度、自由度、协方差及方差组分。合并时逐项核对对象结构、设计和基因标识。未知字段、已收缩对象或无法对应的行均不允许合并。
7. 对整个有效背景统一执行经验贝叶斯收缩、区间计算及原计划范围的BH校正。不把单基因收缩后的p值拼到旧表。
8. 在同一最终未收缩拟合上，再计算一次排除恢复基因的结果，量化它对其它基因SE和p值的影响。

若两优化器不一致、任一个仍有数值警告，或原始对象不能安全合并，就不恢复该基因。不改变模型公式、不降低收敛检查、不删笼位项，也不按显著性选择优化器。

## 数值规则

比较使用 `abs(a-b) <= tolerance * (1 + max(abs(a), abs(b)))`。系数、原始SE、自由度、似然、残差尺度、笼位方差的容差依次为 `1e-5, 1e-4, 1e-3, 1e-6, 1e-4, 1e-4`。这些数值只是优化解一致性的检查精度，不是生物学等效界值。零笼位方差允许存在；非有限值、非正SE或自由度不允许通过。

原输出复现采用 `max(abs(old-new)/(1+abs(old))) < 1e-6`，并同时核对全部缺失位置及状态。原R包版本必须一致，防止把升级依赖引起的变化混入本次诊断。

## 运行方法

把ZIP解压到新的代码目录。不要覆盖正在运行的代码或原结果。以下路径用自己的完整工作目录替换：

```bash
python -u run_dream_repair.py \
  --source-work /path/to/v1.2.0/work_reanalysis_v1_2_0 \
  --work /path/to/new/work_dream_repair_v1_2_2 \
  --rscript "$CONDA_PREFIX/bin/Rscript" \
  --library /path/to/existing/R_library \
  --workers 1
```

不需要重新下载数据。来源目录必须保留，修复目录通过符号链接只读访问两个小鼠计数文件。结果表复制到新目录后修改。定点入口目前不支持中途续跑；若被中断，请指定另一个全新输出目录，避免混入半次复核结果。

`run.py`的完整分析入口采用相同的收敛复核逻辑；只有处理已有v1.2.0结果时才需要上面的定点入口。

## 主要输出

| 文件 | 用途 |
| --- | --- |
| `results/mouse_dream_retry_attempts.csv` | 两优化器的原始估计、警告、是否通过数值一致性检查 |
| `results/dream_repair_default_reproduction.csv` | 默认拟合与v1.2.0结果的逐指标复现检查 |
| `results/dream_repair_before_after.csv.gz` | 受影响模型所有背景基因的前后结果 |
| `results/mouse_dream_background_impact.csv.gz` | 同一最终拟合中，纳入或排除恢复基因的共同收缩结果；未恢复时为空表 |
| `results/mouse_dream_diagnostics.csv` | 最终警告定位及恢复状态 |
| `results/mouse_fit_failures.csv` | 仍不可用的记录；不会把失败当作不显著 |
| `repair_changes.json` | 哪些原文件更新、哪些文件未改动 |
| `review_bundle.zip` | 用于回传核对的结果、日志及本次执行代码 |

`accepted`表示两个单基因优化解通过检查；是否最终恢复，以诊断表的 `recovered_optimizer_agreement` 和完整结果的 `available` 为准。

原 `mouse_variance.csv.gz` 是独立的方差组分提取结果，本入口不改写它，也不把它称为重试后dream的方差估计。重试基因的笼位方差直接见 `mouse_dream_retry_attempts.csv`。

退出码0：完整性检查通过，所有拟合位置可用。

退出码3：计算结束，但仍有不可用拟合或未解决警告。需要阅读诊断明细，不要自行删除警告记录。

退出码2：执行、环境或来源复现检查失败；不能把新目录用于论文定稿。

## 保持不变的分析

人类结果、固定效应/无笼项/相关性模型、笼位SE分解、模拟和配对检验校准全部复用并检查文件摘要。修正后的精确SE计算不回退，旧稿中的SE差异与此次收敛复核是两件事。

## 实现依据与验证范围

- [lme4收敛检查说明](https://lme4.github.io/lme4/reference/convergence.html)：比较不同优化器、加强停止精度，而不是直接屏蔽警告。
- [variancePartition的dream接口](https://gabrielhoffman.github.io/variancePartition/reference/dream-method.html)：通过公开的 `control` 参数传入优化设置，保留REML及Satterthwaite自由度。

发布核验见VALIDATION.md。服务器在分析前仍会执行自身版本的冒烟测试；真实计数中的收敛结论由本次拟合与复核输出决定。
