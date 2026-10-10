# 安装和验证分开。已有依赖齐备时，不联系Bioconductor版本服务器。
a <- commandArgs(TRUE)
mode <- a[1]; lib <- normalizePath(a[2],mustWork=FALSE)
dir.create(lib,recursive=TRUE,showWarnings=FALSE); .libPaths(c(lib,.libPaths()))
if(getRversion()<'4.5.0'||getRversion()>='4.6.0')stop('需要R 4.5.x；请先建立对应环境')
packages <- c('jsonlite','Matrix','edgeR','limma','SeuratObject','variancePartition',
              'BiocParallel','snow','lme4','lmerTest','pbkrtest','statmod')
missing <- packages[!vapply(packages,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing)&&mode=='install') {
  # 固定Bioconductor发行系列，不使用不验证TLS的镜像或curl参数。
  if(!requireNamespace('BiocManager',quietly=TRUE))
    install.packages('BiocManager',lib=lib,repos='https://cloud.r-project.org')
  need_version <- getRversion()>='4.5.0' && getRversion()<'4.6.0'
  if(!need_version)stop('本流程要求R 4.5.x与Bioconductor 3.22，请先调整环境')
  BiocManager::install(missing,lib=lib,version='3.22',ask=FALSE,update=FALSE)
}
missing <- packages[!vapply(packages,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing))stop(paste('缺少R依赖：',paste(missing,collapse=', ')))
for(p in packages)cat(p,as.character(packageVersion(p)),'\n')
if(packageVersion('SeuratObject')<'5.0.0')stop('需要SeuratObject 5或更高版本的layer接口')
if(!('adaptive.span'%in%names(formals(limma::voom))))stop('limma::voom版本不符合要求')
if(mode=='check') {
  source(file.path(a[3],'analysis','common.R'))
  source(file.path(a[3],'analysis','diagnostics.R'))
  source(file.path(a[3],'tests','test_analysis.R'))
  cat('预检查：EList数据类、逐观测权重和SE核查\n')
  run_analysis_regression()
  source(file.path(a[3],'tests','test_preflight.R'))
  run_preflight_regression()
  source(file.path(a[3],'tests','test_dream_retry.R'))
  run_dream_retry_unit_tests()
  run_dream_merge_regression()
  run_dream_retry_smoke(if(length(a)>=4L)a[4] else NULL)
  # 真正启动一个SOCK工作进程，发现仅requireNamespace无法暴露的snow问题。
  cat('预检查：SOCK工作进程\n')
  param <- BiocParallel::SnowParam(workers=1,type='SOCK',fallback=FALSE)
  got <- tryCatch(BiocParallel::bplapply(1:2,function(x)x+1,BPPARAM=param),
                  finally=BiocParallel::bpstop(param))
  values <- unlist(got,use.names=FALSE)
  need(length(got)==2L && all(vapply(got,function(x)is.numeric(x)&&length(x)==1L,logical(1))) &&
       same_numeric_values(values,2:3), 'SOCK返回的结果数量或数值错误；期望两个标量2和3')
  cat('SOCK结果数值核对通过；返回存储类型：',typeof(values),'\n')
  set.seed(17); md<-data.frame(genotype=rep(c(0,1),12),age=rnorm(24),cage=factor(rep(1:6,each=4)))
  rownames(md)<-paste0('s',1:24)
  y<-matrix(rnorm(120*24),120,24,dimnames=list(paste0('g',1:120),rownames(md)))
  cat('预检查：dream拟合、官方区间与p值\n')
  f<-variancePartition::dream(y,~genotype+age+(1|cage),md,ddf='Satterthwaite',
    BPPARAM=BiocParallel::SerialParam(),hideErrorsInBackend=FALSE)
  stopifnot('genotype'%in%colnames(f$coefficients),all(is.finite(f$coefficients[,'genotype'])))
  check<-dream_rows(f,rownames(y),'smoke','dream',robust=FALSE)$table
  stopifnot(all(check$status=='available'))
  raw<-round(exp(y/3+6));dge<-edgeR::calcNormFactors(edgeR::DGEList(raw))
  cat('预检查：dream原生权重\n')
  native<-variancePartition::voomWithDreamWeights(dge,~genotype+age+(1|cage),md,BPPARAM=BiocParallel::SerialParam())
  stopifnot(all(is.finite(native$weights)),all(native$weights>0))
  cat('预检查：Seurat counts层\n')
  object<-SeuratObject::CreateSeuratObject(counts=raw,assay='Spatial')
  recovered<-SeuratObject::GetAssayData(object,assay='Spatial',layer='counts')
  observed <- as.matrix(recovered)
  need(identical(dimnames(observed),dimnames(raw)) && same_numeric_values(observed,raw),
       'Seurat counts的行列标识、维度或数值不一致')
  cat('预检查：方差组分与配对WLS\n')
  vf<-variancePartition::fitVarPartModel(y[1:4,,drop=FALSE],~genotype+age+(1|cage),md,
       REML=TRUE,fxn=variance_extract,BPPARAM=BiocParallel::SerialParam())
  stopifnot(length(vf)==4,all(vapply(vf,is.data.frame,logical(1))),
            setequal(names(vf),rownames(y)[1:4]))
  x0<-model.matrix(~genotype+age,md);x1<-model.matrix(~genotype+age+cage,md)
  w<-exp(rnorm(24,0,.2));yy<-y[1,]
  pair<-paired_wls(yy,w,x0,x1)
  b0<-coef(lm(yy~0+x0,weights=w))[2];b1<-coef(lm(yy~0+x1,weights=w))[2]
  stopifnot(abs(pair$delta-(b1-b0))<1e-9)
  aa<-solve(crossprod(x0,w*x0),t(x0)*rep(w,each=ncol(x0)))
  cv<-diag(ncol(x1))[2,]-as.numeric(aa[2,]%*%x1)
  full<-lm(yy~0+x1,weights=w)
  expected<-summary(full)$sigma^2*as.numeric(cv%*%solve(crossprod(x1,w*x1))%*%cv)
  stopifnot(abs(pair$SE_delta^2-expected)<1e-9)
  cat('独立R进程、SOCK和dream官方API检查通过\n')
}
sessionInfo()
