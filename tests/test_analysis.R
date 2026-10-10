# 精确WLS与两种limma参数化的数值回归测试；由依赖预检查调用。
run_analysis_regression <- function() {
  packages<-c('limma','statmod')
  missing<-packages[!vapply(packages,requireNamespace,logical(1),quietly=TRUE)]
  need(!length(missing),paste('R模型回归测试缺少依赖：',paste(missing,collapse=', ')))
  set.seed(71);n<-36
  d<-data.frame(genotype=rep(c(0,1),18),age=rnorm(n),sex=rep(c(0,0,1),12))
  x<-model.matrix(~genotype+age+sex,d)
  y<-matrix(rnorm(120*n),120,n,dimnames=list(paste0('g',1:120),paste0('s',1:n)))
  w<-matrix(exp(rnorm(120*n,sd=.7)),120,n)
  dimnames(w)<-dimnames(y)
  exact<-exact_wls_target(y,w,x,'genotype')
  ref<-lm(y[1,]~genotype+age+sex,data=d,weights=w[1,])
  stopifnot(abs(exact$beta[1]-coef(ref)['genotype'])<1e-10)
  unscaled<-sqrt(vcov(ref)['genotype','genotype'])/summary(ref)$sigma
  stopifnot(abs(exact$unscaled[1]-unscaled)<1e-10)
  # voom返回EList，普通list并非lmFit接受的表达对象；测试使用相同的正式数据类。
  v<-methods::new('EList',list(E=y,weights=w))
  stopifnot(inherits(v,'EList'),identical(v$E,y),identical(v$weights,w))
  z<-audit_se(v,x,'synthetic','baseline')
  stopifnot(nrow(z)==120L,max(abs(z$unscaled_exact/z$unscaled_direct-1))<1e-8)
  stopifnot(max(abs(z$unscaled_contrast_exact/z$unscaled_direct-1))<1e-8,
            max(abs(z$unscaled_contrast_refit/z$unscaled_direct-1))<1e-8)
  # 年龄平移不改变精确目标SE；近似对比在相同模型空间内可能改变。
  xc<-x;xc[,'age']<-x[,'age']-mean(x[,'age']);attr(xc,'raw_age')<-x[,'age']+20
  zc<-audit_se(v,xc,'synthetic','baseline')
  stopifnot(max(abs(zc$unscaled_direct-z$unscaled_direct))<1e-10)
  stopifnot(max(abs(zc$SE_contrast_uncentered_ordinary-zc$SE_contrast_ordinary))>1e-4)
  stopifnot(max(abs(z$unscaled_contrast/z$unscaled_exact-1))>1e-4)
  # 零/常数权重差异不存在时，对比近似恢复到精确解。
  flat<-v;flat$weights[,]<-1
  flat_audit<-audit_se(flat,x,'unweighted','baseline')
  stopifnot(max(abs(flat_audit$unscaled_contrast/flat_audit$unscaled_exact-1))<1e-10)
  # 测试警告定位逻辑，不以安装了variancePartition为前提。
  fake<-function(ids)if('g3'%in%ids)warning('Model failed to converge: negative eigenvalue')
  located<-locate_dream_warnings(paste0('g',1:9),fake,4L)
  stopifnot(identical(located$gene_symbol,'g3'),located$status=='numerical_warning')
  absent<-locate_dream_warnings(c('a','b'),function(ids)NULL,2L)
  stopifnot(absent$status=='warning_not_reproduced',is.na(absent$gene_symbol))
  # 一块已定位不代表其他块也已解决；保留只能在整块出现的警告。
  mixed_warning<-function(ids) {
    if('a'%in%ids)warning('negative eigenvalue')
    if(length(ids)>1&&'c'%in%ids)warning('Model failed to converge')
  }
  mixed<-locate_dream_warnings(letters[1:4],mixed_warning,2L)
  stopifnot(nrow(mixed)==2L,mixed$gene_symbol[1]=='a',is.na(mixed$gene_symbol[2]),
    mixed$status[2]=='block_warning_unresolved')
  # 与显式矩阵加权接口核对，确认EList中的逐观测权重确实进入了拟合。
  matrix_fit<-limma::lmFit(y,x,weights=w)
  stopifnot(max(abs(z$beta_direct-matrix_fit$coefficients[,'genotype']))<1e-10,
            max(abs(z$unscaled_direct-matrix_fit$stdev.unscaled[,'genotype']))<1e-10)
  # 普通与稳健收缩都必须返回；不能只保存更接近旧数字的那一个。
  stopifnot(all(is.finite(z$SE_direct_ordinary)),all(is.finite(z$SE_direct_robust)),
            all(is.finite(z$SE_contrast_ordinary)),all(is.finite(z$SE_contrast_robust)))
  # 错误对象和无效权重在进入limma之前给出具体错误，不能自动丢弃权重。
  wrong<-try(audit_se(list(E=y,weights=w),x,'synthetic','baseline'),silent=TRUE)
  stopifnot(inherits(wrong,'try-error'),grepl('EList',as.character(wrong),fixed=TRUE))
  invalid<-v;invalid$weights[1,1]<-0
  stopifnot(inherits(try(audit_se(invalid,x,'synthetic','baseline'),silent=TRUE),'try-error'))
  invalid<-v;invalid$weights<-w[,-1,drop=FALSE]
  stopifnot(inherits(try(audit_se(invalid,x,'synthetic','baseline'),silent=TRUE),'try-error'))
  bad<-x;bad[,4]<-bad[,3]
  stopifnot(inherits(try(exact_wls_target(y,w,bad,'genotype'),silent=TRUE),'try-error'))
  cat('EList接口、精确加权系数、四组合标准误及无效输入回归检查通过\n')
}
