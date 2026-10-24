save_json <- function(x,path) {
  dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
  temp<-paste0(path,'.tmp')
  jsonlite::write_json(x,temp,auto_unbox=TRUE,pretty=TRUE,na='null')
  need(file.rename(temp,path),paste('Cannot publish',path))
}

save_rds <- function(x,path) {
  dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
  temp<-paste0(path,'.tmp');saveRDS(x,temp,compress=FALSE)
  need(file.rename(temp,path),paste('Cannot publish',path))
}

finish_job <- function(job,files)save_json(list(id=job$id,files=files),file.path(job$output,'complete.json'))

unavailable <- function(genes,model,spec,reason) {
  data.frame(gene_symbol=genes,spec=spec,model=model,beta=NA_real_,SE=NA_real_,df=NA_real_,
             p=NA_real_,low=NA_real_,high=NA_real_,status=reason)
}

bind_rows <- function(x) {
  x<-Filter(function(z)!is.null(z)&&nrow(z)>0,x)
  if(!length(x))return(data.frame())
  names0<-unique(unlist(lapply(x,names)))
  x<-lapply(x,function(z){for(n in setdiff(names0,names(z)))z[[n]]<-NA;z[,names0,drop=FALSE]})
  do.call(rbind,x)
}

# Drop only redundant nuisance columns; never replace the target with another contrast.
target_basis <- function(raw,target='genotype') {
  need(target%in%colnames(raw),'Target column missing')
  colnames(raw)<-sub('^\\(Intercept\\)$','Intercept',colnames(raw))
  nuisance<-raw[,setdiff(colnames(raw),target),drop=FALSE]
  chosen<-matrix(numeric(),nrow(raw),0)
  keep<-character();rank0<-0L
  for(n in colnames(nuisance)) {
    test<-cbind(chosen,nuisance[,n])
    rank1<-qr(test,tol=1e-9)$rank
    if(rank1>rank0){chosen<-test;keep<-c(keep,n);rank0<-rank1}
  }
  colnames(chosen)<-keep
  augmented<-cbind(chosen,raw[,target]);full<-qr(augmented,tol=1e-9)$rank==ncol(augmented)
  x<-cbind(Intercept=1,genotype=raw[,target],chosen[,setdiff(keep,'Intercept'),drop=FALSE])
  if(target!='genotype')colnames(x)[2]<-target
  list(x=x,estimable=full&&nrow(x)-ncol(x)>=5,
       dropped=setdiff(colnames(nuisance),keep),rank=qr(augmented)$rank)
}

mouse_designs <- function(d,kind) {
  d$age<-d$age-mean(d$age)
  terms<-c('genotype',if(kind=='linear')'age' else 'age_factor','sex')
  if(kind=='age_factor_slide')terms<-c(terms,'slide')
  base_raw<-model.matrix(reformulate(terms),d)
  colnames(base_raw)<-sub('^\\(Intercept\\)$','Intercept',colnames(base_raw))
  fixed_raw<-model.matrix(reformulate(c(terms,'cage')),d)
  means<-sapply(setdiff(colnames(base_raw),'Intercept'),function(k)ave(base_raw[,k],d$cage,FUN=mean))
  if(is.null(dim(means)))means<-matrix(means,ncol=1)
  colnames(means)<-paste0('cage_mean_',make.names(setdiff(colnames(base_raw),'Intercept')))
  list(base=target_basis(base_raw),fixed=target_basis(fixed_raw),
       mundlak=target_basis(cbind(base_raw,means)),raw=base_raw,data=d)
}

mixed_input <- function(design,d) {
  frame<-data.frame(design$x,check.names=FALSE)
  # Slide labels can contain hyphens; formula syntax must not interpret them as subtraction.
  names(frame)<-make.names(names(frame),unique=TRUE)
  fixed_names<-names(frame)
  frame$cage<-d$cage;rownames(frame)<-d$animal_id
  formula<-as.formula(paste('~ 0 +',paste(fixed_names,collapse=' + '),'+ (1|cage)'))
  list(data=frame,formula=formula)
}

fit_linear <- function(v,x,genes,model,spec,robust=FALSE,block=NULL,rho=NULL,target='genotype') {
  tryCatch(fit_rows(moderate(limma::lmFit(v,x,block=block,correlation=rho),robust),
                    genes,target,spec,model),error=function(e){
    z<-unavailable(genes,model,spec,'fit_failed');attr(z,'message')<-conditionMessage(e);z
  })
}

fit_mixed <- function(v,input,genes,model,spec,param) {
  ans<-checked_dream(v,input$formula,input$data,genes,spec,model,param,FALSE)
  unresolved<-nrow(ans$diagnostics)>0 && any(is.na(ans$diagnostics$gene_symbol))
  if(unresolved) {
    ans$table$status[ans$table$status=='available']<-'unresolved_numerical_warning'
    ans$table[,c('p','low','high')]<-NA_real_
  }
  ans
}

paired_weights <- function(y,w0,w1,x0,x1) {
  need(estimable(x0)&&estimable(x1),'Non-estimable paired design')
  need(qr(cbind(x1,x0))$rank==ncol(x1),'Non-nested mean spaces')
  a0<-solve(crossprod(x0,w0*x0),t(x0)*rep(w0,each=ncol(x0)))
  a1<-solve(crossprod(x1,w1*x1),t(x1)*rep(w1,each=ncol(x1)))
  d<-a1[2,]-a0[2,];delta<-sum(d*y);df<-length(y)-ncol(x1)
  e<-as.numeric(y-x1%*%(a1%*%y));h<-rowSums(x1*t(a1))
  unit<-sum(d*d/w1);se<-sqrt(sum(w1*e*e)/df*unit)
  hc3<-if(max(h)<1-1e-10)sqrt(sum(d*d*(e/(1-h))^2)) else NA_real_
  valid<-is.finite(se)&&se>1e-12&&unit>1e-14
  data.frame(delta=delta,SE_delta=if(valid)se else NA_real_,
    p=if(valid)2*pt(-abs(delta/se),df) else NA_real_,
    low=if(valid)delta-qt(.975,df)*se else NA_real_,
    high=if(valid)delta+qt(.975,df)*se else NA_real_,
    SE_HC3=hc3,p_HC3=if(is.finite(hc3)&&hc3>1e-12)2*pt(-abs(delta/hc3),df) else NA_real_,
    low_HC3=if(is.finite(hc3)&&hc3>1e-12)delta-qt(.975,df)*hc3 else NA_real_,
    high_HC3=if(is.finite(hc3)&&hc3>1e-12)delta+qt(.975,df)*hc3 else NA_real_,
    df=df,max_leverage=max(h),status=if(valid)'available' else 'degenerate_difference')
}

weighted_center_identity <- function(v,design,d) {
  if(!design$fixed$estimable)return(data.frame())
  x<-design$raw;fx<-design$fixed$x;mx<-design$mundlak$x
  ids<-seq_len(nrow(v$E))
  records<-lapply(ids,function(i) {
    w<-v$weights[i,];y<-v$E[i,]
    wm<-sapply(setdiff(colnames(x),'Intercept'),function(k)
      ave(w*x[,k],d$cage,FUN=sum)/ave(w,d$cage,FUN=sum))
    colnames(wm)<-paste0('weighted_mean_',make.names(setdiff(colnames(x),'Intercept')))
    aug<-target_basis(cbind(x,wm))
    bfix<-unname(lm.wfit(fx,y,w)$coefficients['genotype'])
    bweighted<-if(aug$estimable)unname(lm.wfit(aug$x,y,w)$coefficients['genotype']) else NA_real_
    barithmetic<-if(design$mundlak$estimable)unname(lm.wfit(mx,y,w)$coefficients['genotype']) else NA_real_
    data.frame(gene_symbol=rownames(v$E)[i],beta_FE=bfix,beta_weighted_augmented_WLS=bweighted,
      beta_arithmetic_augmented_WLS=barithmetic,difference_weighted=bweighted-bfix,
      difference_arithmetic=barithmetic-bfix)
  })
  z<-do.call(rbind,records)
  need(all(is.finite(z$difference_weighted))&&max(abs(z$difference_weighted))<1e-6,
       'Weighted-centering WLS identity failed')
  z
}
