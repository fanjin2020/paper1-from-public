# 不依赖真实数据，检验数值门槛、警告保留和两优化器一致性判据。
run_dream_retry_unit_tests <- function() {
  a<-c(beta=.1,SE_raw=.2,df_raw=30,logLik=-100,sigma=1,cage_variance=0)
  need(dream_metrics_agree(a,a),'相同估计未通过')
  b<-a;b['beta']<-.2;need(!dream_metrics_agree(a,b),'系数明显不同时错误通过')
  b<-a;b['df_raw']<-0;need(!dream_metrics_agree(a,b),'零自由度错误通过')
  b<-a;b['SE_raw']<-NA_real_;need(!dream_metrics_agree(a,b),'缺失SE错误通过')
  b<-a;b['cage_variance']<-1;need(!dream_metrics_agree(a,b),'方差明显不同时错误通过')
  z<-capture_dream_call(function(){warning('negative eigenvalue');list(ok=TRUE)})
  need(!z$failed&&length(z$messages)==1L,'警告丢失')
  z<-capture_dream_call(function()stop('failed'))
  need(z$failed&&is.null(z$fit),'错误被当作成功')
  fake<-list(coefficients=matrix(.1,1,1,dimnames=list('g','genotype')),
    stdev.unscaled=matrix(.2,1,1),df.residual=matrix(30,1,1),logLik=-100,sigma=1)
  attr(fake,'varComp')<-data.frame(cage=0,row.names='g')
  controls<-list(bobyqa='one',nloptwrap='two')
  z<-probe_dream_retry('g',function(ids,control)fake,controls)
  need(z$accepted&&all(z$table$clean),'零方差且一致的有限估计未通过')
  z<-probe_dream_retry('g',function(ids,control){
    if(control=='two')warning('negative eigenvalue')
    fake
  },controls)
  need(!z$accepted,'带数值警告的重试被接受')
  z<-probe_dream_retry('g',function(ids,control){
    if(control=='two')attr(fake,'errors')<-c(g='backend failed')
    fake
  },controls)
  need(!z$accepted,'后台错误被遗漏')
  z<-probe_dream_retry('g',function(ids,control){
    if(control=='two')fake$logLik<--99
    fake
  },controls)
  need(!z$accepted,'似然不一致的重试被接受')
  need(is.null(z$fit),'不一致的原始拟合被保留为可合并结果')
  z<-probe_dream_retry('g',function(ids,control)fake,controls)
  need(identical(z$fit,fake),'通过检查后没有保留原始拟合')
  cat('dream复核判据：13项检查通过\n')
}

# 用正式拟合对象检查行合并，不允许顺手修改正常基因或已经收缩的结果。
run_dream_merge_regression <- function() {
  set.seed(419)
  md<-data.frame(genotype=rep(c(0,1),24),age=rep(seq(-1,1,length.out=12),each=4),
    cage=factor(rep(seq_len(12),each=4)))
  rownames(md)<-paste0('s',seq_len(nrow(md)))
  y<-matrix(rnorm(24*48),24,48)+matrix(rnorm(24*12,sd=2),24,12)[,md$cage]
  dimnames(y)<-list(paste0('g',seq_len(nrow(y))),rownames(md))
  v<-methods::new('EList',list(E=y,weights=matrix(exp(rnorm(length(y),0,.15)),nrow(y))))
  formula<-~genotype+age+(1|cage)
  full<-variancePartition::dream(v,formula,md,ddf='Satterthwaite',BPPARAM=BiocParallel::SerialParam())
  one<-variancePartition::dream(v['g1',,drop=FALSE],formula,md,ddf='Satterthwaite',
    BPPARAM=BiocParallel::SerialParam())
  merged<-replace_dream_raw_row(full,one,'g1')
  ids<-rownames(full$coefficients);unchanged<-ids!='g1'
  # 不依赖S3子集方法识别字段，逐项检查真实存储内容。
  for(k in names(full)) {
    a<-full[[k]];b<-merged[[k]]
    if(is.matrix(a)&&identical(rownames(a),ids))
      need(identical(a[unchanged,,drop=FALSE],b[unchanged,,drop=FALSE]),paste('正常矩阵行改变：',k))
    else if(is.atomic(a)&&length(a)==length(ids))
      need(identical(a[unchanged],b[unchanged]),paste('正常逐基因量改变：',k))
    else if(k=='cov.coefficients.list')
      need(identical(a[unchanged],b[unchanged]),'正常基因协方差改变')
  }
  need(identical(attr(full,'varComp')[unchanged,,drop=FALSE],attr(merged,'varComp')[unchanged,,drop=FALSE]),
       '正常基因方差组分改变')
  expect_error<-function(a,b,g)inherits(tryCatch(replace_dream_raw_row(a,b,g),error=function(e)e),'error')
  bad<-one;bad$unknown_gene_statistic<-1
  need(expect_error(full,bad,'g1'),'未知字段未被拦截')
  bad<-one;bad$s2.post<-1
  need(expect_error(full,bad,'g1'),'收缩后结果未被拦截')
  bad<-one;colnames(bad$coefficients)[1]<-'other'
  need(expect_error(full,bad,'g1'),'系数标识错误未被拦截')
  need(expect_error(full,one,'unknown'),'不存在的基因未被拦截')
  bad<-one;bad$design[1,1]<-2
  need(expect_error(full,bad,'g1'),'不一致设计未被拦截')
  bad<-one;attr(bad,'varComp')<-NULL
  need(expect_error(full,bad,'g1'),'缺失方差组分未被拦截')
  cat('dream原始合并：正常行不变及6项错误输入检查通过\n')
}

# Save numerical diagnostics before checking the smoke-test result.
write_dream_smoke_report <- function(ans,report_dir) {
  if(!is.null(report_dir)) {
    dir.create(report_dir,recursive=TRUE,showWarnings=FALSE)
    for(k in c('table','diagnostics','retry_attempts','background_impact'))
      write.csv(ans[[k]],file.path(report_dir,paste0('dream_retry_smoke_',k,'.csv')),row.names=FALSE)
    writeLines(enc2utf8(ans$messages),file.path(report_dir,'dream_retry_smoke_messages.txt'),useBytes=TRUE)
  }
  if(any(ans$table$status!='available')) {
    cat('dream重试冒烟测试详细信息：\n')
    print(table(ans$table$status));print(ans$diagnostics);print(ans$retry_attempts)
    cat(paste(ans$messages,collapse='\n'),'\n')
  }
}

# 用官方API走重试、原始行合并和共同收缩；不因测试困难而降低收敛门槛。
run_dream_retry_smoke <- function(report_dir=NULL) {
  set.seed(821)
  md<-data.frame(genotype=rep(c(0,1),24),cage=factor(rep(seq_len(12),each=4)))
  rownames(md)<-paste0('s',seq_len(nrow(md)))
  noise<-matrix(rnorm(120*48),120,48)
  re<-matrix(rnorm(120*12,sd=2),120,12)
  y<-noise+re[,as.integer(md$cage)]
  dimnames(y)<-list(paste0('g',1:120),rownames(md))
  v<-methods::new('EList',list(E=y,weights=matrix(1,120,48)))
  ans<-checked_dream(v,~genotype+(1|cage),md,rownames(y),'smoke','dream_native',
    BiocParallel::SerialParam(),known_excluded='g1')
  write_dream_smoke_report(ans,report_dir)
  need(all(ans$table$status=='available'),'官方重试冒烟测试未全部可用')
  need(any(ans$diagnostics$status=='recovered_optimizer_agreement'),'未走到重试恢复分支')
  need(nrow(ans$background_impact)==120,'没有共同背景的纳入/排除对照')
  need(all(ans$retry_attempts$clean)&&all(ans$retry_attempts$accepted),'两个优化器未一致通过')
  z<-ans$background_impact
  need(is.na(z$p_without[z$gene_symbol=='g1']),'排除对照仍包含测试基因')
  need(all(is.finite(z$p_with)),'恢复后存在无效p值')
  cat('dream官方API：两优化器复核、定点原始行合并及共同收缩通过\n')
}
