# 只重建有警告的dream模型；其它模型、模拟和校准由Python按摘要核验后复用。
a<-commandArgs(TRUE);work<-a[1];.libPaths(c(a[2],.libPaths()));root<-a[3]
source(file.path(root,'analysis','common.R'))
source(file.path(root,'analysis','diagnostics.R'))
suppressPackageStartupMessages({library(edgeR);library(limma);library(variancePartition);library(BiocParallel);library(jsonlite)})
cfg<-fromJSON(file.path(work,'run_config.json'));out<-file.path(work,'results')
plan<-fromJSON(file.path(work,'repair_plan.json'))
session_old<-paste(readLines(file.path(out,'mouse_R_session.txt'),warn=FALSE),collapse='\n')
versions<-vapply(c('edgeR','limma','variancePartition','lme4','lmerTest','pbkrtest'),function(p) {
  version<-as.character(packageVersion(p))
  need(grepl(paste0(p,'_',version),gsub('-','.',session_old,fixed=TRUE),fixed=TRUE),paste('R包版本与原运行不一致：',p,version))
  version
},character(1))
d<-rd(file.path(work,'prepared','mouse_units.csv'),colClasses=c(animal_id='character',cage='character',slide='character'))
d$genotype<-as.numeric(d$genotype=='5X');d$sex<-factor(d$sex,levels=c('F','M'))
d$age_factor<-factor(d$age);d$cage<-factor(d$cage);d$slide<-factor(d$slide);rownames(d)<-d$animal_id
bg<-rd(file.path(out,'mouse_background.csv'))$gene_symbol
before<-rd(file.path(out,'mouse_background_fits.csv.gz'))
all_retry<-list();all_impact<-list();all_problems<-list();all_messages<-list();all_tables<-list();comparisons<-list()
param<-backend(cfg$workers)

# 核对旧输出的每个可用估计，确保本轮仍在同一响应、权重及模型上复核。
check_original <- function(old,new,spec,model) {
  need(identical(old$gene_symbol,new$gene_symbol),'复现结果的基因顺序改变')
  need(identical(old$status,new$status),'默认拟合未复现原状态；先核查环境和权重')
  rows<-lapply(c('beta','SE','df','p','low','high'),function(name) {
    x<-old[[name]];y<-new[[name]]
    need(identical(is.na(x),is.na(y)),paste('原结果缺失位置改变',name))
    ok<-is.finite(x)&is.finite(y)
    err<-if(any(ok))max(abs(x[ok]-y[ok])/(1+abs(x[ok]))) else 0
    data.frame(spec=spec,model=model,metric=name,max_scaled_difference=err,passed=err<1e-6)
  })
  z<-do.call(rbind,rows)
  need(all(z$passed),'默认拟合未能在1e-6容差内复现，不能直接合并新旧输出')
  z
}

for(k in seq_len(nrow(plan$targets))) {
  spec<-plan$targets$spec[k];model<-plan$targets$model[k]
  region<-if(startsWith(spec,'whole_'))'whole' else 'cortex'
  kind<-substring(spec,nchar(region)+2)
  need(kind%in%c('linear','age_factor','age_factor_slide'),'未知小鼠设计')
  dd<-d;dd$age<-dd$age-mean(dd$age)
  terms<-c('genotype',if(kind=='linear')'age' else 'age_factor','sex')
  if(kind=='age_factor_slide')terms<-c(terms,'slide')
  x<-model.matrix(reformulate(terms),dd)
  formula<-as.formula(paste('~',paste(c(terms,'(1|cage)'),collapse='+')))
  y<-counts(file.path(work,'prepared',paste0('mouse_',region,'_counts.csv.gz')))
  need(identical(colnames(y),rownames(dd))&&all(bg%in%rownames(y)),'计数与原背景或样本身份不一致')
  g<-calcNormFactors(DGEList(y[bg,,drop=FALSE]))
  cat('重建权重',spec,model,length(bg),'基因\n')
  weights_run<-capture_dream_call(function() {
    if(model=='dream_native')variancePartition::voomWithDreamWeights(g,formula,dd,BPPARAM=param) else
      voom(g,x,plot=FALSE,adaptive.span=TRUE)
  })
  need(!weights_run$failed&&!length(numerical_dream_messages(weights_run$messages)),
    paste('权重重建出现错误或数值警告：',paste(weights_run$messages,collapse=' | ')))
  v<-weights_run$fit
  known<-plan$failures$gene_symbol[plan$failures$spec==spec&plan$failures$model==model]
  cat('复现及复核',spec,model,'\n')
  ans<-checked_dream(v,formula,dd,bg,spec,model,param,cfg$dream_ebayes_robust,known_excluded=known)
  wr(ans$retry_attempts,file.path(work,'repair_details',paste0(spec,'_',model,'_attempts.csv')))
  wr(ans$diagnostics,file.path(work,'repair_details',paste0(spec,'_',model,'_diagnostics.csv')))
  need(!is.null(ans$default_table),paste('默认拟合未完成：',paste(ans$messages,collapse=' | ')))
  old<-before[before$spec==spec&before$model==model,,drop=FALSE]
  old<-old[match(bg,old$gene_symbol),,drop=FALSE]
  comparisons[[k]]<-check_original(old,ans$default_table,spec,model)
  all_tables[[k]]<-ans$table;all_retry[[k]]<-ans$retry_attempts;all_impact[[k]]<-ans$background_impact
  all_problems[[k]]<-ans$diagnostics
  all_messages[[k]]<-data.frame(spec=rep(spec,length(ans$messages)),model=rep(model,length(ans$messages)),message=ans$messages)
}
wr(do.call(rbind,comparisons),file.path(out,'dream_repair_default_reproduction.csv'))
wr(do.call(rbind,all_tables),file.path(out,'dream_repair_replacement_rows.csv.gz'))
wr(do.call(rbind,all_retry),file.path(out,'mouse_dream_retry_attempts.csv'))
wr(do.call(rbind,all_impact),file.path(out,'mouse_dream_background_impact.csv.gz'))
wr(do.call(rbind,all_problems),file.path(out,'dream_repair_diagnostics.csv'))
wr(do.call(rbind,all_messages),file.path(out,'dream_repair_messages.csv'))
capture.output(sessionInfo(),file=file.path(out,'dream_repair_R_session.txt'))
BiocParallel::bpstop(param)
