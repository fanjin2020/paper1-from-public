# Gaussian条件模拟：四模型使用同一响应，正式调用相关模型和dream。
a<-commandArgs(TRUE);work<-a[1];.libPaths(c(a[2],.libPaths()));root<-a[3]
source(file.path(root,'analysis','common.R'))
source(file.path(root,'analysis','diagnostics.R'))
suppressPackageStartupMessages({library(limma);library(variancePartition);library(BiocParallel);library(jsonlite)})
cfg<-fromJSON(file.path(work,'run_config.json'))
tasks<-rd(file.path(work,'results','simulation_scenarios.csv'))
folder<-file.path(work,'prepared','simulation');cache<-file.path(work,'simulation_cache')
dir.create(cache,showWarnings=FALSE)

# 未收缩的点态t区间；dream额外用官方p值复核内部SE与自由度。
extract_simulation <- function(f,dream=FALSE) {
  j<-match('genotype',colnames(f$coefficients));need(!is.na(j),'目标系数未返回')
  b<-as.numeric(f$coefficients[1,j]);se<-as.numeric(f$stdev.unscaled[1,j]*f$sigma[1])
  df<-as.numeric(if(is.matrix(f$df.residual))f$df.residual[1,j] else f$df.residual[1])
  need(all(is.finite(c(b,se,df)))&&se>0&&df>0,'拟合未产生有效区间')
  p<-2*pt(-abs(b/se),df)
  if(dream)need(is.finite(f$p.value[1,j])&&abs(p-f$p.value[1,j])<1e-7,'dream未收缩区间与官方p不一致')
  c(beta=b,SE=se,df=df,p=p,low=b-qt(.975,df)*se,high=b+qt(.975,df)*se)
}

parts<-list()
for(index in seq_len(nrow(tasks))) {
  task<-tasks[index,];name<-task$layout
  file<-file.path(cache,paste0(task$scenario,'.csv.gz'));receipt<-paste0(file,'.json')
  if(file.exists(file)&&file.exists(receipt)) {
    old<-fromJSON(receipt)
    need(identical(unname(tools::md5sum(file)),old$md5),'模拟缓存校验失败')
    z<-rd(file)
    need(nrow(z)==task$planned*4&&!anyDuplicated(z[,c('replicate','model')]),'模拟缓存位置不完整')
    parts[[index]]<-z;cat('已核验模拟情境',task$scenario,'\n');next
  }
  d<-rd(file.path(folder,paste0(name,'_units.csv')),colClasses=c(cage='character'))
  # Python已输出连续整数编码。显式按整数level排序，避免0,1,10,11,2的字符顺序。
  codes<-suppressWarnings(as.integer(d$cage))
  need(!anyNA(codes)&&identical(sort(unique(codes)),seq.int(0L,length(unique(codes))-1L)),
       '模拟笼位必须是连续整数编码')
  d$cage<-factor(codes,levels=sort(unique(codes)));rownames(d)<-paste0('s',seq_len(nrow(d)))
  x<-as.matrix(rd(file.path(folder,paste0(name,'_base.csv'))))
  xf<-as.matrix(rd(file.path(folder,paste0(name,'_fixed.csv'))))
  fixed_ok<-qr(xf)$rank==ncol(xf)&&nrow(xf)-ncol(xf)>1
  # 干扰效应真值为零，但拟合仍保留真实布局的年龄和性别。
  nuisance<-x[,setdiff(colnames(x),c('Intercept','genotype')),drop=FALSE]
  md<-data.frame(genotype=d$genotype,cage=d$cage,row.names=rownames(d))
  if(ncol(nuisance))for(j in seq_len(ncol(nuisance)))md[[paste0('x',j)]]<-nuisance[,j]
  terms<-c('genotype',if(ncol(nuisance))paste0('x',seq_len(ncol(nuisance))),'(1|cage)')
  form<-as.formula(paste('~',paste(terms,collapse='+')))
  cage_shift<-ave(d$genotype,d$cage,FUN=mean)-mean(d$genotype)
  betas<-rep(c(0,.3),length.out=64);betas[1]<-task$truth
  signal<-outer(betas,d$genotype)+matrix(rep(task$shift*cage_shift,each=64),64,nrow(d))
  records<-list()
  for(r in seq_len(task$planned)) {
    set.seed(as.integer(task$seed+r))
    y<-signal+sqrt(1-task$rho_truth)*matrix(rnorm(64*nrow(d)),64,nrow(d))+
      sqrt(task$rho_truth)*matrix(rnorm(64*nlevels(d$cage)),64,nlevels(d$cage))[,as.integer(d$cage)]
    dimnames(y)<-list(paste0('g',seq_len(64)),rownames(d))
    for(model in c('baseline','cage_fixed','cage_correlation','dream')) {
      record<-data.frame(scenario=task$scenario,layout=name,rho_truth=task$rho_truth,
        shift=task$shift,truth=task$truth,replicate=r,model=model,beta=NA_real_,SE=NA_real_,
        df=NA_real_,p=NA_real_,low=NA_real_,high=NA_real_,rho_estimated=NA_real_,boundary=NA,
        status='failed',reason='',warning='')
      warnings<-character()
      if(model=='cage_fixed'&&!fixed_ok) {
        record$status<-'not_estimable';record$reason<-'基因型与笼位不能分离'
      } else tryCatch(withCallingHandlers({
        if(model=='baseline')f<-lmFit(y,x)
        if(model=='cage_fixed')f<-lmFit(y,xf)
        if(model=='cage_correlation') {
          rho<-duplicateCorrelation(y,x,block=d$cage)$consensus.correlation
          lower<- -1/(max(table(d$cage))-1)
          need(is.finite(rho)&&rho>lower&&rho<1,'相关参数不可估计或不满足正定性')
          record$rho_estimated<-rho;record$boundary<-min(rho-lower,1-rho)<1e-6
          f<-lmFit(y,x,block=d$cage,correlation=rho)
        }
        if(model=='dream') {
          f<-variancePartition::dream(y[1,,drop=FALSE],form,md,REML=TRUE,ddf='Satterthwaite',
            BPPARAM=SerialParam(),hideErrorsInBackend=FALSE)
          vv<-variancePartition::fitVarPartModel(y[1,,drop=FALSE],form,md,REML=TRUE,
            fxn=variance_extract,BPPARAM=SerialParam())
          need(length(vv)==1&&is.data.frame(vv[[1]]),'dream方差诊断未返回')
          need(isTRUE(vv[[1]]$converged),'dream随机模型未收敛')
          record$rho_estimated<-vv[[1]]$ICC;record$boundary<-vv[[1]]$singular
        }
        z<-extract_simulation(f,model=='dream')
        # 只有一个目标基因，数值警告可以直接归属本次检验，不用另一轮拟合的收敛状态代替。
        if(model=='dream')need(!length(numerical_dream_messages(warnings)),
          paste('dream数值警告：',paste(numerical_dream_messages(warnings),collapse=' | ')))
        record[1,names(z)]<-as.list(z);record$status<-'available'
      },warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')}),
        error=function(e)record$reason<<-conditionMessage(e))
      record$warning<-paste(unique(warnings),collapse=' | ');records[[length(records)+1]]<-record
    }
    if(r%%100==0||r==task$planned)cat(sprintf('模拟 %s %s rho %.2f shift %.2f truth %.2f %d / %d\n',
      task$scenario,name,task$rho_truth,task$shift,task$truth,r,task$planned))
  }
  z<-do.call(rbind,records);wr(z,file)
  write_json(list(md5=unname(tools::md5sum(file)),planned=nrow(z)),receipt,auto_unbox=TRUE)
  parts[[index]]<-z
}
all<-do.call(rbind,parts)
wr(all,file.path(work,'results','simulation_replicates.csv.gz'))
wr(all[all$status!='available',],file.path(work,'results','simulation_failures.csv'))
key<-interaction(all$scenario,all$model,drop=TRUE)
summary<-lapply(split(all,key),function(z) {
  d<-z[z$status=='available',];n<-nrow(d);r<-z[1,c('scenario','layout','rho_truth','shift','truth','model')]
  r$planned<-nrow(z);r$available<-n;r$numerical_failed<-sum(z$status=='failed')
  r$not_estimable<-sum(z$status=='not_estimable');r$failure_rate<-r$numerical_failed/nrow(z)
  r$unavailable_rate<-1-n/nrow(z);r$not_estimable_rate<-r$not_estimable/nrow(z)
  r$bias<-if(n)mean(d$beta-d$truth) else NA_real_;r$empirical_SD<-if(n>1)sd(d$beta) else NA_real_
  r$mean_SE<-if(n)mean(d$SE) else NA_real_;r$RMS_SE<-if(n)sqrt(mean(d$SE^2)) else NA_real_
  r$coverage<-if(n)mean(d$low<=d$truth & d$high>=d$truth) else NA_real_
  r$rejection<-if(n)mean(d$p<.05) else NA_real_
  r$coverage_all_planned<-sum(d$low<=d$truth & d$high>=d$truth)/nrow(z)
  r$rejection_all_planned<-sum(d$p<.05)/nrow(z)
  r$width<-if(n)mean(d$high-d$low) else NA_real_
  r$rejection_MCSE<-if(n)sqrt(r$rejection*(1-r$rejection)/n) else NA_real_
  r$rho_available<-sum(is.finite(d$rho_estimated))
  r$rho_mean<-if(r$rho_available)mean(d$rho_estimated,na.rm=TRUE) else NA_real_
  r$boundary_evaluable<-sum(!is.na(d$boundary))
  r$boundary_rate<-if(r$boundary_evaluable)mean(d$boundary,na.rm=TRUE) else NA_real_
  r
})
wr(do.call(rbind,summary),file.path(work,'results','simulation_summary.csv'))
