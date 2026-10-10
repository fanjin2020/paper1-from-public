# 数据来源与处理

## SEA-AD：人类死后脑单核计数

官方说明：https://brain-map.org/consortia/sea-ad/our-data

S3根地址：https://sea-ad-single-cell-profiling.s3.amazonaws.com/

程序列出`PFC/RNAseq/donors_objects/`，选取文件名以`_SEAAD_DFC_RNAseq_final-nuclei.2026-06-22.h5ad`结尾的80个对象。发行日期写在config中，不自动使用未来版本。清单保存对象名、大小和ETag。

供体工作簿来自官方页面的`sea-ad_cohort_donor_metadata_072524.xlsx`链接；完整URL在`pipeline/sources.py`的`URLS`中。工作簿实际列名须包含Donor ID、Cognitive Status、PMI、RIN和Brain pH。

处理顺序：

1. 检查H5AD内的供体编号与文件名一致，脑区为DFC，细胞条码唯一。
2. 保留作者`Used in analysis`为真且不是`Neurotypical reference`的细胞核。Braak、年龄及性别必须能解释，不能从其他标签猜诊断。
3. 只读取`layers/UMIs`原始稀疏计数。`X`或标准化表达不能替代这一层。
4. 按供体×作者Subclass汇总计数；重复基因符号相加。未提供的基因保持缺失，不填零。不同供体基因空间不一致时，计数模型停止并要求核查。
5. 以Donor ID连接官方质量信息，与H5AD临床状态核对。
6. 临床比较仅纳入Dementia和No dementia，使用质量变量完整的同一批供体。每个亚类至少20个回收核；每个临床组至少5位供体；残差自由度至少5。

PMI、RIN及pH的取材定义不相同：供体工作簿中的RIN并不等于每个测序文库的质量测量。论文应沿用来源定义，不能把这些值写成逐文库质量控制指标。作者注释在本流程中被使用，不被重新宣称为独立认证。

## GSE233208：小鼠空间表达及动物身份

官方记录：https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE233208

论文：https://doi.org/10.1038/s41588-024-01961-x

下载GEO公开Seurat对象`GSE233208_5XFAD_seurat_processed_annotated.rds.gz`以及六份样本表：pilot、8samples和四个日期的sample_meta。实际下载地址在代码中逐项列出。

动物身份来自论文补充资料压缩包的`SupplementaryTables.xlsx`，读取`Supp. Table 2`。GEO表负责确认分析资格，S2提供动物ID、笼号和耳号。二者按SAMPLE精确连接，核对年龄、性别和基因型。

关键处理：

- 只用对象中的`SAMPLE`，不采用不完整的`SampleID`，不从spot条码或文件名提取动物号。
- 性别缺失的pilot样本不静默补入；所有来源样本的资格记录保留。
- 一个动物若对应多个SAMPLE，本版本明确停止，要求先核查重复切片关系。
- 笼号及耳号组合须唯一；不把每个空间spot当作独立动物。
- 使用Spatial assay的counts层，兼容SeuratObject 5的`layer`接口。
- 全脑聚合所有合格spot；皮层仅聚合`ctx-deep-layers`和`ctx-upper-layers`，不包含嗅皮层。
- 皮层至少100个spot，两个区域及所有模型共享满足条件的动物。
- 最终动物数、笼数和混合基因型笼数由本次资料计算。不会从历史74动物清单选择样本。

公开RDS可能有两层gzip。Python只去掉多余外层，R读取剩下一层。压缩包内容不是通用R对象或Spatial counts层缺失时停止，不改读处理矩阵。

## GSE118553：面板来源

官方记录：https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE118553

采用GEO series matrix和GPL10558官方平台注释。处理矩阵是作者已处理的芯片表达，不是芯片原始探针强度。本流程不擅自追加log转换。

1. 读取每个GSM的individual、tissue、disease state、gender和age。
2. 仅保留额叶、内嗅及颞叶皮层。供体×脑区中的技术重复先在探针层取中位数。
3. 探针只接受唯一、合法的基因符号；同基因多探针再取中位数。
4. 额叶和内嗅分别拟合AD对control，调整中心化年龄与性别，女性为参照；完整病例每组至少10。
5. 每个脑区在全部有效基因上BH校正。额叶q<0.05形成本次面板，不由SEA-AD或小鼠结果筛选。

历史论文使用237个基因。本代码不随包携带这237个名字；如果新下载和本次流程得到不同数量，需要核对平台注释、处理矩阵和模型定义，再决定论文数字是否更新。**在未核对前，不能称为历史237面板完全复现。**

## MGI：人鼠基因对应

官方报告：https://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt

按DB Class Key查找一个人类符号与一个小鼠符号的同源类，只保留全表唯一的一对一对应。缺失或歧义映射不从旧清单补齐；端点位置仍保留，并标记不可检验。报告随MGI更新，本地保存下载日期、URL和SHA-256，便于解释版本差异。

## 文件完整性

下载保留TLS验证，支持ETag或Last-Modified保护的续传，核对响应范围与完整长度。完成后记录SHA-256。已有文件只有收据、URL、长度和摘要全部一致才可使用。

SHA-256记录确认本地文件后续未变；出版方未提供原始摘要时，不能把本地摘要称为出版方身份认证。官方元数据如有矛盾，程序报错，不以“合理猜测”修补。
