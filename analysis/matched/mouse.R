fit_mouse_models <- function(g,d,kind,spec,workers=1L,empirical=TRUE) {
  genes<-rownames(g);design<-mouse_designs(d,kind)
  models<-c('baseline','cage_fixed','cage_correlation',
             if(empirical)'dream_shared','dream_native',if(empirical)'mundlak_shared','mundlak_native')
  tables<-list();messages<-list();checks<-list();retry<-list()
  param<-backend(workers);on.exit(BiocParallel::bpstop(param),add=TRUE)
  if(!design$base$estimable)return(list(table=bind_rows(lapply(models,function(m)
    unavailable(genes,m,spec,'target_not_estimable'))),messages=data.frame(),diagnostics=data.frame(),
    retries=data.frame(),design=design,vref=NULL,fixed_shared=NULL))
  v<-voom(g,design$base$x,plot=FALSE,adaptive.span=TRUE)
  tables$baseline<-fit_linear(v,design$base$x,genes,'baseline',spec)
  fixed_shared<-NULL
  if(design$fixed$estimable) {
    fixed_shared<-fit_linear(v,design$fixed$x,genes,'cage_fixed_shared_check',spec)
    vf<-voom(g,design$fixed$x,plot=FALSE,adaptive.span=TRUE)
    tables$cage_fixed<-fit_linear(vf,design$fixed$x,genes,'cage_fixed',spec)
  } else tables$cage_fixed<-unavailable(genes,'cage_fixed',spec,'target_not_estimable')
  tables$cage_correlation<-tryCatch({
    rho1<-duplicateCorrelation(v,design$base$x,block=d$cage)$consensus.correlation
    need(is.finite(rho1),'Non-finite first consensus correlation')
    vc<-voom(g,design$base$x,block=d$cage,correlation=rho1,plot=FALSE,adaptive.span=TRUE)
    rho<-duplicateCorrelation(vc,design$base$x,block=d$cage)$consensus.correlation
    need(is.finite(rho),'Non-finite second consensus correlation')
    z<-fit_linear(vc,design$base$x,genes,'cage_correlation',spec,block=d$cage,rho=rho)
    z$rho_estimate<-rho;z
  },error=function(e){
    messages[[length(messages)+1]]<<-data.frame(model='cage_correlation',message=conditionMessage(e))
    unavailable(genes,'cage_correlation',spec,'fit_failed')
  })
  for(type in c('dream','mundlak'))for(policy in if(empirical)c('shared','native') else 'native') {
    model<-paste(type,policy,sep='_');ds<-if(type=='dream')design$base else design$mundlak
    if(!ds$estimable){tables[[model]]<-unavailable(genes,model,spec,'target_not_estimable');next}
    inp<-mixed_input(ds,d)
    weight_messages<-character()
    ans<-tryCatch({
      vv<-if(policy=='shared')v else withCallingHandlers(
        variancePartition::voomWithDreamWeights(g,inp$formula,inp$data,BPPARAM=param),
        warning=function(w){weight_messages<<-c(weight_messages,conditionMessage(w));invokeRestart('muffleWarning')})
      fit_mixed(vv,inp,genes,model,spec,param)
    },error=function(e)list(table=unavailable(genes,model,spec,'fit_failed'),messages=conditionMessage(e),
                            diagnostics=data.frame(),retry_attempts=data.frame()))
    tables[[model]]<-ans$table
    if(length(weight_messages))ans$messages<-unique(c(ans$messages,paste0('weight_estimation: ',weight_messages)))
    if(length(ans$messages))messages[[length(messages)+1]]<-data.frame(model=model,message=ans$messages)
    checks[[length(checks)+1]]<-ans$diagnostics
    retry[[length(retry)+1]]<-ans$retry_attempts
    cat('Mouse fit',spec,model,sum(ans$table$status=='available'),'/',length(genes),'\n')
  }
  list(table=bind_rows(tables),messages=bind_rows(messages),diagnostics=bind_rows(checks),
       retries=bind_rows(retry),design=design,vref=v,fixed_shared=fixed_shared)
}

run_mouse <- function(job,cfg) {
  Sys.setenv(PAPER1_DIAGNOSTICS_DIR=file.path(job$output,'failed_raw_fits'))
  obj<-readRDS(file.path(cfg$out,'cache',paste0('mouse_',job$spec,'.rds')))
  ans<-fit_mouse_models(obj$dge,obj$metadata,obj$kind,job$spec,cfg$inner_workers,TRUE)
  z<-ans$table
  z$q_background_spec<-ave(z$p,z$model,FUN=bh)
  wr(z,file.path(job$output,'background_fits.csv.gz'))
  panel<-obj$panel
  mapped<-lapply(unique(z$model),function(model) {
    part<-z[z$model==model,];i<-match(panel$mouse_gene,part$gene_symbol)
    zz<-part[i,,drop=FALSE];zz$mouse_gene<-panel$mouse_gene;zz$gene_symbol<-panel$gene_symbol
    zz$spec<-job$spec;zz$model<-model;zz$status[is.na(i)]<-'outside_background_or_unmapped'
    zz
  })
  wr(bind_rows(mapped),file.path(job$output,'panel_fits.csv.gz'))
  identity<-if(!is.null(ans$vref))weighted_center_identity(ans$vref,ans$design,obj$metadata) else data.frame()
  wr(identity,file.path(job$output,'weighted_centering_identity.csv.gz'))
  wr(ans$fixed_shared,file.path(job$output,'fixed_shared_check.csv.gz'))
  geometry<-lapply(c('base','fixed','mundlak'),function(model) {
    ds<-ans$design[[model]];x<-ds$x
    nuisance<-x[,colnames(x)!='genotype',drop=FALSE]
    info<-sum(qr.resid(qr(nuisance),x[,'genotype'])^2)
    data.frame(spec=job$spec,model=model,n=nrow(x),columns=ncol(x),rank=qr(x)$rank,
      df_residual=nrow(x)-qr(x)$rank,estimable=ds$estimable,unweighted_target_information=info,
      dropped_redundant_nuisance=paste(ds$dropped,collapse=';'),
      mean_policy=if(model=='mundlak')'arithmetic_means_of_all_baseline_covariates' else 'not_applicable')
  })
  wr(bind_rows(geometry),file.path(job$output,'designs.csv'))
  wr(ans$messages,file.path(job$output,'messages.csv'))
  wr(ans$diagnostics,file.path(job$output,'diagnostics.csv'))
  wr(ans$retries,file.path(job$output,'retry_attempts.csv'))
  raw_files<-list.files(file.path(job$output,'failed_raw_fits'),recursive=TRUE,full.names=FALSE)
  usable<-vapply(split(z,z$model),function(part)any(part$status=='available' & is.finite(part$p)),logical(1))
  need(length(usable)==7L && all(usable),paste('Whole empirical mouse model failed:',
    paste(names(usable)[!usable],collapse=', '),'; inspect messages.csv'))
  finish_job(job,c('background_fits.csv.gz','panel_fits.csv.gz','weighted_centering_identity.csv.gz',
    'fixed_shared_check.csv.gz','designs.csv','messages.csv','diagnostics.csv','retry_attempts.csv',
    if(length(raw_files))file.path('failed_raw_fits',raw_files)))
}
