# 新计数上的人类匹配模型、小鼠笼位模型和标准误分解。
a<-commandArgs(TRUE);mode<-a[1];work<-a[2];.libPaths(c(a[3],.libPaths()));root<-a[4]
source(file.path(root,'analysis','common.R'))
source(file.path(root,'analysis','diagnostics.R'))
suppressPackageStartupMessages({library(edgeR);library(limma);library(variancePartition);library(BiocParallel);library(jsonlite)})
cfg<-fromJSON(file.path(work,'run_config.json'));out<-file.path(work,'results')
param<-backend(cfg$workers)

# 在共同voom响应及权重下拟合一个设计，不因面板结果改变背景。
fit_voom <- function(v,x,cage=NULL,rho=NULL) moderate(lmFit(v,x,block=cage,correlation=rho),cfg$ebayes_robust)

if(mode=='human') {
  md<-rd(file.path(work,'prepared','human','metadata.csv'))
  y<-counts(file.path(work,'prepared','human','expression.csv.gz'))
  need(identical(colnames(y),md$unit_id),'人类表达与单位对应错误')
  panel<-rd(file.path(out,'panel.csv'))$gene_symbol
  contexts<-sort(unique(md$context)); contrasts<-c('clinical_status','clinical_status_adjusted_for_braak','braak_same_clinical_population')
  result<-list();paired<-list();status<-list();family<-list();membership<-list();i<-0
  for(context in contexts) {
    selected<-md$context==context & md$n_cells>=cfg$minimum_nuclei & md$clinical_group%in%c('Dementia','No_dementia')
    selected<-selected & complete.cases(md[,c('age','sex','pathology','PMI','RIN','pH')])
    d<-md[selected,,drop=FALSE];need(!anyDuplicated(d$donor_id),'类群内供体重复')
    xfilter<-human_design(d,'clinical_status')
    possible<-nrow(d)>=10 && all(table(factor(d$clinical_group,levels=c('No_dementia','Dementia')))>=5) && estimable(xfilter)
    v<-NULL
    if(possible) {
      counts0<-y[,selected,drop=FALSE];g<-DGEList(counts0)
      keep<-filterByExpr(g,design=xfilter);g<-g[keep,,keep.lib.sizes=FALSE]
      if(nrow(g)>=100) {
        g<-calcNormFactors(g);v<-voom(g,xfilter,plot=FALSE,adaptive.span=TRUE)
      }
    }
    membership[[context]]<-data.frame(context=context,gene_symbol=panel,
      present_in_counts=panel%in%rownames(y),
      present_after_filter=if(is.null(v))FALSE else panel%in%rownames(v$E),
      match_policy='source_symbol_case_normalized',
      note='计数中的符号已统一大小写；不按旧稿可用数量删基因')
    for(contrast in contrasts) {
      i<-i+1;spec<-paste(context,contrast,sep='::')
      x0<-human_design(d,contrast);x1<-human_design(d,contrast,TRUE)
      ok<-!is.null(v)&&estimable(x0)&&estimable(x1)
      status[[i]]<-data.frame(spec=spec,n=nrow(d),p0=ncol(x0),p1=ncol(x1),status=if(ok)'available' else 'not_estimable')
      for(arm in c('base','quality')) {
        if(ok) {
          f<-fit_voom(v,if(arm=='base')x0 else x1)
          z<-fit_rows(f,panel,'target',spec,arm)
        } else z<-data.frame(gene_symbol=panel,spec=spec,model=arm,beta=NA_real_,SE=NA_real_,df=NA_real_,p=NA_real_,low=NA_real_,high=NA_real_,status='not_estimable')
        z$context<-context;z$contrast<-contrast;z$n<-nrow(d)
        result[[length(result)+1]]<-z
      }
      if(ok) {
        common<-intersect(panel,rownames(v$E))
        changes<-lapply(common,function(gene) {
          k<-match(gene,rownames(v$E));z<-paired_wls(v$E[k,],v$weights[k,],x0,x1)
          cbind(gene_symbol=gene,spec=spec,z)
        })
        if(length(changes))paired[[length(paired)+1]]<-do.call(rbind,changes)
        # 设计只来自本次完整病例，用于后续条件校准和可识别性核查。
        wr(data.frame(unit_id=d$unit_id,x0,check.names=FALSE),file.path(work,'prepared','designs',paste0(i,'_base.csv')))
        wr(data.frame(unit_id=d$unit_id,x1,check.names=FALSE),file.path(work,'prepared','designs',paste0(i,'_quality.csv')))
      }
      cat('人类模型',i,'/',length(contexts)*3,spec,if(ok)'完成' else '不可估计','\n')
    }
  }
  all<-do.call(rbind,result)
  all$q_joint<-ave(all$p,all$model,FUN=bh)
  wr(all,file.path(out,'human_results.csv.gz'))
  wr(do.call(rbind,status),file.path(out,'human_designs.csv'))
  wr(do.call(rbind,membership),file.path(out,'human_gene_membership.csv'))
  grid<-expand.grid(gene_symbol=panel,spec=unique(all$spec),stringsAsFactors=FALSE)
  changes<-if(length(paired))do.call(rbind,paired) else data.frame(gene_symbol=character(),spec=character(),p=numeric(),p_HC3=numeric())
  changes<-merge(grid,changes,by=c('gene_symbol','spec'),all.x=TRUE,sort=FALSE)
  changes$q_joint<-bh(changes$p);changes$q_HC3_joint<-bh(changes$p_HC3)
  wr(changes,file.path(out,'human_paired_changes.csv.gz'))
  # 显式枚举校正范围，不把不同family计算出的q值混用。
  scope_summary<-list()
  for(arm in unique(all$model))for(scope in c('all','contrast','context','context_contrast'))for(denominator in c('planned','available_only')) {
    z<-all[all$model==arm,];key<-switch(scope,all=rep('all',nrow(z)),contrast=z$contrast,context=z$context,context_contrast=z$spec)
    adjust<-if(denominator=='planned')bh else function(p) {
      q<-rep(NA_real_,length(p));ok<-is.finite(p);q[ok]<-p.adjust(p[ok],'BH');q
    }
    z$q_scope<-ave(z$p,key,FUN=adjust)
    selected<-all[all$model=='base' & all$contrast=='clinical_status' & !is.na(all$q_joint) & all$q_joint<.05,]
    clinical<-z[z$contrast=='clinical_status',]
    selected_key<-paste(selected$spec,selected$gene_symbol)
    scope_summary[[length(scope_summary)+1]]<-data.frame(model=arm,scope=scope,denominator=denominator,
      planned=nrow(clinical),available=sum(is.finite(clinical$p)),
      significant=sum(clinical$q_scope<.05,na.rm=TRUE),
      original40_retained=sum(clinical$q_scope<.05 & paste(clinical$spec,clinical$gene_symbol)%in%selected_key,na.rm=TRUE))
    for(k in unique(key)) {
      zz<-z[key==k,];family[[length(family)+1]]<-data.frame(model=arm,scope=scope,family=k,
        denominator=denominator,
        planned=nrow(zz),available=sum(is.finite(zz$p)),significant=sum(zz$q_scope<.05,na.rm=TRUE))
    }
  }
  wr(do.call(rbind,family),file.path(out,'correction_scope.csv'))
  wr(do.call(rbind,scope_summary),file.path(out,'human_correction_scope_summary.csv'))
}

# 六种小鼠设定共享动物和全脑/皮层共同的表达背景。
mouse_design <- function(d,kind) {
  raw_age<-d$age;d$age<-d$age-mean(d$age)
  terms<-c('genotype',if(kind=='linear')'age' else 'age_factor','sex')
  if(kind=='age_factor_slide')terms<-c(terms,'slide')
  form<-reformulate(terms)
  x<-model.matrix(form,d);fixed<-model.matrix(reformulate(c(terms,'cage')),d)
  attr(x,'raw_age')<-raw_age;attr(fixed,'raw_age')<-raw_age
  list(x=x,fixed=fixed,data=d,
       dream=as.formula(paste('~',paste(c(terms,'(1|cage)'),collapse='+'))))
}

if(mode=='mouse') {
  d<-rd(file.path(work,'prepared','mouse_units.csv'),colClasses=c(animal_id='character',cage='character',slide='character'))
  d$genotype<-as.numeric(d$genotype=='5X');d$sex<-factor(d$sex,levels=c('F','M'))
  d$age_factor<-factor(d$age);d$cage<-factor(d$cage);d$slide<-factor(d$slide);rownames(d)<-d$animal_id
  sources<-lapply(c('whole','cortex'),function(r)counts(file.path(work,'prepared',paste0('mouse_',r,'_counts.csv.gz'))))
  names(sources)<-c('whole','cortex')
  for(y in sources)need(identical(colnames(y),d$animal_id),'小鼠表达与身份不匹配')
  linear<-mouse_design(d,'linear')$x
  background<-Reduce(intersect,lapply(sources,function(y)rownames(y)[filterByExpr(DGEList(y),design=linear)]))
  need(length(background)>=100,'小鼠共同表达背景不足')
  wr(data.frame(gene_symbol=background),file.path(out,'mouse_background.csv'))
  panel<-rd(file.path(out,'ortholog_panel.csv'));panel$mouse_gene[is.na(panel$mouse_gene)]<-paste0('__UNMAPPED__',panel$gene_symbol[is.na(panel$mouse_gene)])
  results<-list();decomposition<-list();variance<-list();diagnostics<-list()
  audits<-list();background_rows<-list();messages<-list();model_status<-list();dream_diagnostics<-list()
  retry_attempts<-list();background_impact<-list()
  for(region in names(sources))for(kind in c('linear','age_factor','age_factor_slide')) {
    spec<-paste(region,kind,sep='_');design<-mouse_design(d,kind)
    need(estimable(design$x)&&estimable(design$fixed),paste('小鼠设计不可估计',spec))
    g<-DGEList(sources[[region]][background,,drop=FALSE]);g<-calcNormFactors(g)
    v<-voom(g,design$x,plot=FALSE,adaptive.span=TRUE)
    vf<-voom(g,design$fixed,plot=FALSE,adaptive.span=TRUE)
    base<-fit_voom(v,design$x);fixed_shared<-fit_voom(v,design$fixed);fixed<-fit_voom(vf,design$fixed)
    # 先用初始权重估计相关，再更新voom权重，第二次估计相关并作最终GLS。
    rho1<-duplicateCorrelation(v,design$x,block=d$cage)$consensus.correlation
    need(is.finite(rho1),paste('第一遍相关参数不可用',spec))
    vcorr<-voom(g,design$x,block=d$cage,correlation=rho1,plot=FALSE,adaptive.span=TRUE)
    rho<-duplicateCorrelation(vcorr,design$x,block=d$cage)$consensus.correlation
    need(is.finite(rho),paste('相关参数不可用',spec))
    corr<-fit_voom(vcorr,design$x,d$cage,rho)
    audits[[paste0(spec,'_base')]]<-audit_se(v,design$x,spec,'baseline')
    audits[[paste0(spec,'_fixed')]]<-audit_se(vf,design$fixed,spec,'cage_fixed')
    cat('运行',spec,'dream，共',length(background),'基因\n')
    dream_shared<-checked_dream(v,design$dream,design$data,background,spec,'dream_shared',param,cfg$dream_ebayes_robust)
    vn<-variancePartition::voomWithDreamWeights(g,design$dream,design$data,BPPARAM=param)
    dream_native<-checked_dream(vn,design$dream,design$data,background,spec,'dream_native',param,cfg$dream_ebayes_robust)
    fits<-list(baseline=base,cage_fixed=fixed,cage_correlation=corr)
    tables<-lapply(names(fits),function(model)fit_rows(fits[[model]],background,'genotype',spec,model))
    names(tables)<-names(fits)
    tables$dream_shared<-dream_shared$table;tables$dream_native<-dream_native$table
    for(model in names(tables)) {
      tab<-tables[[model]];background_rows[[paste(spec,model)]]<-tab
      model_status[[paste(spec,model)]]<-data.frame(spec=spec,model=model,planned=nrow(tab),
        available=sum(tab$status=='available'),failed=sum(tab$status!='available'))
      ix<-match(panel$mouse_gene,tab$gene_symbol);z<-tab[ix,,drop=FALSE]
      z$gene_symbol<-panel$mouse_gene;z$spec<-spec;z$model<-model
      z$status[is.na(ix)]<-'outside_background_or_unmapped'
      z$mouse_gene<-z$gene_symbol;z$gene_symbol<-panel$gene_symbol
      z$n<-nrow(d);z$cages<-nlevels(d$cage);z$rho<-if(model=='cage_correlation')rho else NA_real_
      results[[length(results)+1]]<-z
    }
    for(model in c('dream_shared','dream_native')) {
      msg<-get(model)$messages
      dream_diagnostics[[paste(spec,model)]]<-get(model)$diagnostics
      retry_attempts[[paste(spec,model)]]<-get(model)$retry_attempts
      background_impact[[paste(spec,model)]]<-get(model)$background_impact
      if(length(msg))messages[[length(messages)+1]]<-data.frame(spec=spec,model=model,message=msg)
    }
    common<-Reduce(intersect,lapply(list(base,fixed_shared,fixed),function(f)rownames(f$coefficients)))
    b<-base[common,];s<-fixed_shared[common,];f<-fixed[common,];j<-match('genotype',colnames(b$coefficients))
    design_factor<-s$stdev.unscaled[,j]/b$stdev.unscaled[,j]
    posterior_factor<-sqrt(s$s2.post/b$s2.post)
    weight_factor<-(f$stdev.unscaled[,j]*sqrt(f$s2.post))/(s$stdev.unscaled[,j]*sqrt(s$s2.post))
    observed<-(f$stdev.unscaled[,j]*sqrt(f$s2.post))/(b$stdev.unscaled[,j]*sqrt(b$s2.post))
    need(max(abs(observed-design_factor*posterior_factor*weight_factor))<1e-8,'SE分解恒等式不成立')
    # 无权设计信息是前瞻近似；不能把精确加权恒等式当作独立验证。
    u0<-sqrt(solve(crossprod(design$x))[j,j]);u1<-sqrt(solve(crossprod(design$fixed))[j,j])
    decomposition[[length(decomposition)+1]]<-data.frame(gene_symbol=common,spec=spec,SE_ratio=observed,
      design_factor=design_factor,posterior_factor=posterior_factor,reweight_factor=weight_factor,
      unweighted_design_prediction=u1/u0,df_base=nrow(d)-ncol(design$x),df_fixed=nrow(d)-ncol(design$fixed))
    cat('提取',spec,'笼位方差组分\n')
    for(policy in c('shared','native')) {
      obj<-if(policy=='shared')v else vn
      vs<-checked_variance(obj,design$dream,design$data,background,spec,policy,param)
      variance[[length(variance)+1]]<-vs$table
      if(length(vs$messages))messages[[length(messages)+1]]<-data.frame(spec=spec,
        model=paste0('variance_',policy),message=vs$messages)
    }
    diagnostics[[length(diagnostics)+1]]<-data.frame(spec=spec,animals=nrow(d),cages=nlevels(d$cage),
      mixed_cages=sum(tapply(d$genotype,d$cage,function(g)length(unique(g)))==2),rho_first=rho1,rho=rho,
      information_fraction=(u0/u1)^2,background=length(background))
  }
  all<-do.call(rbind,results)
  # 两模型主比较与三模型比较分开保存；不能把每模型较小family的q冒充联合q。
  all$q_within_model<-ave(all$p,all$model,FUN=bh)
  all$q_base_fixed_joint<-NA_real_;all$q_three_model_joint<-NA_real_;all$q_dream_joint<-NA_real_
  k<-all$model%in%c('baseline','cage_fixed');all$q_base_fixed_joint[k]<-bh(all$p[k])
  k<-all$model%in%c('baseline','cage_fixed','cage_correlation');all$q_three_model_joint[k]<-bh(all$p[k])
  k<-all$model%in%c('dream_shared','dream_native');all$q_dream_joint[k]<-bh(all$p[k])
  wr(all,file.path(out,'mouse_results.csv.gz'))
  wr(do.call(rbind,decomposition),file.path(out,'mouse_precision.csv.gz'))
  wr(do.call(rbind,variance),file.path(out,'mouse_variance.csv.gz'))
  wr(do.call(rbind,diagnostics),file.path(out,'mouse_designs.csv'))
  wr(do.call(rbind,audits),file.path(out,'mouse_SE_implementation_audit.csv.gz'))
  bg<-do.call(rbind,background_rows)
  wr(bg,file.path(out,'mouse_background_fits.csv.gz'))
  wr(bg[bg$status!='available',],file.path(out,'mouse_fit_failures.csv'))
  wr(do.call(rbind,model_status),file.path(out,'mouse_model_status.csv'))
  msgs<-if(length(messages))do.call(rbind,messages) else data.frame(spec=character(),model=character(),message=character())
  wr(msgs,file.path(out,'mouse_model_messages.csv'))
  wr(do.call(rbind,dream_diagnostics),file.path(out,'mouse_dream_diagnostics.csv'))
  wr(do.call(rbind,retry_attempts),file.path(out,'mouse_dream_retry_attempts.csv'))
  wr(do.call(rbind,background_impact),file.path(out,'mouse_dream_background_impact.csv.gz'))
}
BiocParallel::bpstop(param)
capture.output(sessionInfo(),file=file.path(out,paste0(mode,'_R_session.txt')))
