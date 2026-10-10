# 上传GitHub

## 上传哪些文件

解压 `Paper1_GitHub_v1.2.2_20261009.zip`，进入 `paper1-from-public`。把这个目录里的文件放在仓库根目录，打开仓库首页应直接看到 `README.md`、`run.py`、`analysis` 和 `pipeline`。

不要只把ZIP作为一个文件上传到仓库首页。ZIP可以另外作为GitHub Release附件。

本包已经排除原始数据、供体级数据文件、聚合计数、R库、运行缓存、日志、临时测试文件和旧版压缩包。运行后，`.gitignore` 会忽略默认工作目录和输出文件。使用自定义工作目录时，建议放在仓库外。

## 网页上传

1. 在GitHub新建空仓库，例如 `paper1-model-adjustment`。
2. 选择上传文件，将解压目录中的文件和子目录加入仓库根目录。
3. 确认 `.gitignore`、`.gitattributes`、配置和环境文件也已包含。
4. 检查页面上目录结构与README一致，再提交。

网页拖拽可能漏掉隐藏文件。Linux终端上传更容易完整保留文件。

## Linux终端上传

在含 `run.py` 的目录执行。先在GitHub创建空仓库，再将最后一段URL替换为自己的仓库地址：

```bash
git init -b main
git add .
git diff --cached --stat
git status --short
```

确认暂存区只有代码、配置、测试和文档，没有 `raw`、`prepared`、`results` 或 `R_library` 后继续：

```bash
git commit -m "Release Paper1 analysis v1.2.2"
git remote add origin https://github.com/YOUR_ACCOUNT/YOUR_REPOSITORY.git
git push -u origin main
```

如仓库已经有提交，不要强制推送；在现有仓库中复制本包内容，检查差异后正常提交。

## 发布版本

在新建仓库且确认版本号未占用后，可创建 `v1.2.2` tag及Release。`CITATION.cff` 已提供作者和软件版本信息，尚未填写仓库URL或DOI。

如需Zenodo DOI，完成实际归档后再把真实地址补入引用信息和论文。不要用GitHub网页链接替代尚未生成的DOI。

作者还需决定软件许可证；数据许可与软件许可分别处理。此包不代为选择许可证，也没有自动建立或公开任何远程仓库。

## 文件摘要

`release_manifest.json` 列出所有发布文件的SHA-256，本文件自身不进入自己的摘要清单。分析源码保留已验证字节形式，`.gitattributes` 对这些文件关闭自动换行转换，以便ZIP和Git检出版本都能核对同一摘要。

编辑任意源码后应更新版本与发布清单。不要在完成记录仍指向旧代码的工作目录内悄悄替换源码后续跑。
