extra <- Sys.getenv("R_LIBS", "")
if(nzchar(extra)) .libPaths(c(strsplit(extra,.Platform$path.sep,fixed=TRUE)[[1]],.libPaths()))
if(!('adaptive.span' %in% names(formals(limma::voom))))stop('limma::voom lacks adaptive.span; use the recorded compatible environment')
suppressPackageStartupMessages(library(limma))
options(stringsAsFactors=FALSE, digits=17, warn=1)
a<-commandArgs(TRUE)
decode_path<-function(s) {v<-rawToChar(as.raw(strtoi(substring(s,seq(1,nchar(s),2),seq(2,nchar(s),2)),16L)));Encoding(v)<-'UTF-8';v}
out<-decode_path(a[1]); source<-decode_path(a[2])
need<-function(ok,msg) if(!isTRUE(ok)) stop(msg)
wr<-function(x,n) write.csv(x,file.path(out,n),row.names=FALSE,na='')
read_expression<-function(p) {
 e<-read.csv(p,check.names=FALSE); y<-as.matrix(e[,-1,drop=FALSE]); storage.mode(y)<-'double'
 rownames(y)<-e[[1]]; need(!anyDuplicated(rownames(y)),'Duplicate gene symbols'); y
}
meta<-read.csv(file.path(out,'metadata.csv'),check.names=FALSE,colClasses=c(donor_id='character'))
ep<-read.csv(file.path(out,'endpoints.csv'),check.names=FALSE)
ep<-ep[ep$discovery_region=='frontal_cortex',,drop=FALSE]
need(nrow(ep)>0 && !anyDuplicated(ep$gene_symbol),'Invalid frontal endpoints')
fixed_results<-list(); signature_results<-list(); design_results<-list(); sample_results<-list()
gene_connection<-gzfile(file.path(out,'all_gene_models.csv.gz'),'wt'); first_gene<-TRUE
adjust_family<-function(d,field='p') {
 valid<-is.finite(d[[field]]); p<-d[[field]]; p[!valid]<-1
 d$joint_q<-p.adjust(p,'BH'); d$joint_q[!valid]<-NA_real_; d$family_size<-nrow(d); d
}
make_design<-function(d,effect,extra=FALSE,quadratic=FALSE) {
 x<-cbind(intercept=rep(1,nrow(d)),effect=effect)
 age<-d$age-mean(d$age)
 if(length(unique(d$age))>1) x<-cbind(x,age=age)
 if(length(unique(d$sex))>1) x<-cbind(x,male=as.numeric(d$sex=='male'))
 if(extra) x<-cbind(x,braak=d$pathology-mean(d$pathology))
 if(quadratic) x<-cbind(x,age_squared=age^2)
 x
}
fit_model<-function(y,d,effect,raw=FALSE,extra=FALSE,quadratic=FALSE,group_gate=TRUE) {
 need(identical(as.character(colnames(y)),as.character(d$unit_id)),'Matrix/sample order mismatch')
 need(!anyDuplicated(d$donor_id),'Repeated donor inside a model')
 need(all(is.finite(d$age)) && all(d$sex %in% c('male','female')),'Invalid age/sex')
 x<-make_design(d,effect,extra,quadratic)
 rank<-if(nrow(x)) qr(x)$rank else 0L
 ready<-nrow(x)>=10 && rank==ncol(x) && nrow(x)-ncol(x)>=5 && group_gate
 status<-if(ready) 'design_ready' else 'not_testable_design'
 gate<-data.frame(n=nrow(x),rank=rank,columns=ncol(x),residual_df=nrow(x)-rank,
    condition_number=if(ready) kappa(x) else NA_real_,status=status)
 g<-data.frame(gene_symbol=rownames(y),beta=NA_real_,SE=NA_real_,CI_low=NA_real_,CI_high=NA_real_,p=NA_real_,status='not_testable_design')
 input<-NULL
 if(ready) {
  complete<-rowSums(is.finite(y))==ncol(y); keep<-rep(FALSE,nrow(y))
  if(raw) {
   need(all(y>=0,na.rm=TRUE) && all(abs(y-round(y))<1e-8,na.rm=TRUE),'Noninteger counts')
   need(all(colSums(y[complete,,drop=FALSE])>0),'Empty count library')
   dg<-edgeR::DGEList(y[complete,,drop=FALSE]); keep[complete]<-edgeR::filterByExpr(dg,design=x)
   g$status<-'low_count_or_missing'
   if(sum(keep)>=100) {
    dg<-edgeR::calcNormFactors(dg[keep[complete],,keep.lib.sizes=FALSE]); input<-voom(dg,x,plot=FALSE,span=0.5,adaptive.span=TRUE)
   }
  } else {
   keep[complete]<-apply(y[complete,,drop=FALSE],1,var)>0; g$status<-'constant_or_missing'
   if(sum(keep)>=100) input<-y[keep,,drop=FALSE]
  }
  if(!is.null(input)) {
   lf<-lmFit(input,x); good<-is.finite(lf$sigma) & lf$sigma>1e-10 & lf$df.residual>=5
   g$status[keep]<-'unassessable_residual'
   if(sum(good)>=100) {
    ff<-eBayes(lf[good,],trend=FALSE,robust=FALSE)
    se<-ff$stdev.unscaled[,'effect']*sqrt(ff$s2.post); beta<-ff$coefficients[,'effect']
    ci<-qt(.975,ff$df.total)*se; p<-ff$p.value[,'effect']; pos<-which(keep)[good]
    g[pos,c('beta','SE','CI_low','CI_high','p')]<-cbind(beta,se,beta-ci,beta+ci,p)
    g$status[pos]<-'tested'; input<-input[good,,drop=FALSE]
   } else input<-NULL
  }
 }
 list(genes=g,input=input,design=x,gate=gate)
}
record_model<-function(fit,d,id,context,model,scale) {
 info<-data.frame(test_id=id,context=context,model=model,effect_scale=scale)
 decorate<-function(t) cbind(info[rep(1,nrow(t)),,drop=FALSE],t,row.names=NULL)
 design_results[[id]]<<-decorate(fit$gate)
 sample_results[[id]]<<-decorate(d)
 g<-fit$genes; q<-adjust_family(g); g$within_model_background_q<-q$joint_q
 write.table(decorate(g),gene_connection,sep=',',row.names=FALSE,col.names=first_gene,na='',quote=TRUE)
 first_gene<<-FALSE
 pos<-match(ep$gene_symbol,g$gene_symbol); f<-g[pos,,drop=FALSE]; rownames(f)<-NULL
 f$gene_symbol<-ep$gene_symbol; f$status[is.na(pos)]<-'unmapped_endpoint'
 f$discovery_direction<-ep$discovery_direction
 f$direction_relation<-ifelse(!is.finite(f$beta),'not_assessable',ifelse(f$beta==0,'zero',ifelse(sign(f$beta)==ifelse(ep$discovery_direction=='Up',1,-1),'same','opposite')))
 fixed_results[[id]]<<-decorate(f)
 for(direction in c('Up','Down')) {
  symbols<-ep$gene_symbol[ep$discovery_direction==direction]; tested<-g$gene_symbol[g$status=='tested']
  ix<-which(tested %in% symbols)
  sr<-data.frame(signature=direction,fixed_size=length(symbols),measured_size=length(ix),background_size=length(tested),low_background_review=length(tested)<200,Direction=NA_character_,Correlation=NA_real_,p=NA_real_,status='not_testable')
  if(!is.null(fit$input) && length(ix)>=15 && length(tested)-length(ix)>=15) {
   need(identical(rownames(fit$input),tested),'CAMERA background alignment failure')
   cc<-camera(fit$input,list(signature=ix),design=fit$design,contrast=2,use.ranks=FALSE,inter.gene.cor=NA,allow.neg.cor=FALSE)
   sr$Direction<-cc$Direction; sr$Correlation<-cc$Correlation; sr$p<-cc$PValue; sr$status<-'tested'
  }
  signature_results[[paste(id,direction)]]<<-decorate(sr)
 }
}
finish_models<-function() {
 close(gene_connection)
 fr<-adjust_family(do.call(rbind,fixed_results)); sr<-adjust_family(do.call(rbind,signature_results))
 wr(fr,'all_fixed_endpoint_results.csv');wr(sr,'all_signature_results.csv')
 wr(do.call(rbind,design_results),'design_status.csv');wr(do.call(rbind,sample_results),'model_sample_manifest.csv')
 wr(do.call(rbind,lapply(split(fr,fr$test_id),function(z) data.frame(test_id=z$test_id[1],context=z$context[1],model=z$model[1],planned=nrow(z),tested=sum(z$status=='tested'),significant=sum(z$joint_q<=.05,na.rm=TRUE),same=sum(z$direction_relation=='same'),opposite=sum(z$direction_relation=='opposite')))),'endpoint_summary.csv')
 writeLines(capture.output(sessionInfo()),file.path(out,'R_session.txt'))
 invisible(fr)
}



sensitivity_scaled_kappa <- function(X) {
 if(is.null(X) || !nrow(X) || any(!is.finite(X)) || qr(X)$rank<ncol(X)) return(Inf)
 Z<-X
 for(j in seq_len(ncol(Z))) {
  if(length(unique(Z[,j]))>1) Z[,j]<-Z[,j]-mean(Z[,j])
  norm<-sqrt(sum(Z[,j]^2))
  if(norm==0) return(Inf)
  Z[,j]<-Z[,j]/norm
 }
 kappa(Z,exact=TRUE)
}
sensitivity_fixed_fit <- function(input,X,coefficient,endpoints,robust=FALSE,
                            block=NULL,correlation=NULL,min_genes=100,min_df=5) {
 stopifnot(length(endpoints)>0,!anyDuplicated(endpoints))
 if(is.null(X)) X<-matrix(numeric(),nrow=0,ncol=1,dimnames=list(NULL,coefficient))
 tab<-data.frame(gene_symbol=endpoints,beta=NA_real_,SE=NA_real_,
                 CI_low=NA_real_,CI_high=NA_real_,p=NA_real_,status='not_testable_design',
                 error='',stringsAsFactors=FALSE)
 rank<-if(nrow(X) && all(is.finite(X))) qr(X)$rank else 0L
 tab$n<-nrow(X);tab$columns<-ncol(X);tab$rank<-rank
 tab$residual_df<-nrow(X)-rank;tab$scaled_condition_number<-sensitivity_scaled_kappa(X)
 tab$condition_review_required<-tab$scaled_condition_number>30
 tab$robust_ebayes<-robust;tab$ci_scope<-'pointwise_moderated_not_simultaneous'
 if(is.null(input) || rank<ncol(X) || nrow(X)-rank<min_df || !(coefficient%in%colnames(X))) return(tab)
 tryCatch({
  lf<-if(is.null(block)) limma::lmFit(input,X) else
      limma::lmFit(input,X,block=block,correlation=correlation)
  good<-is.finite(lf$sigma) & lf$sigma>1e-10 & lf$df.residual>=min_df
  tab$status<-'insufficient_residual_background'
  if(sum(good)<min_genes) return(tab)
  ff<-limma::eBayes(lf[good,],trend=FALSE,robust=robust)
  b<-ff$coefficients[,coefficient];se<-ff$stdev.unscaled[,coefficient]*sqrt(ff$s2.post)
  h<-qt(.975,ff$df.total)*se
  pos<-match(endpoints,rownames(ff));ok<-!is.na(pos)
  tab$status<-'unmapped_or_unassessable_endpoint'
  tab[ok,c('beta','SE','CI_low','CI_high','p')]<-cbind(b,se,b-h,b+h,ff$p.value[,coefficient])[pos[ok],,drop=FALSE]
  tab$status[ok]<-'tested'
  tab
 },error=function(e){tab$status<-'sensitivity_fit_error';tab$error<-conditionMessage(e);tab})
}
sensitivity_write_family <- function(rows,path,family_fields) {
 if(!length(rows)) return(invisible(NULL))
 z<-do.call(rbind,rows);z$sensitivity_q<-NA_real_;z$sensitivity_family_size<-NA_integer_
 key<-do.call(paste,c(z[family_fields],sep='::'))
 for(id in unique(key)) {
  ix<-which(key==id);p<-z$p[ix];valid<-is.finite(p)
  stopifnot(all(p[valid]>=0 & p[valid]<=1))
  p[!valid]<-1;q<-p.adjust(p,'BH');q[!valid]<-NA_real_
  z$sensitivity_q[ix]<-q;z$sensitivity_family_size[ix]<-length(ix)
 }
 z$analysis_role<-'quality_sensitivity'
 write.csv(z,path,row.names=FALSE,na='')
 invisible(z)
}

if(!('adaptive.span' %in% names(formals(limma::voom))))stop('limma::voom lacks adaptive.span; use the recorded compatible environment')
suppressPackageStartupMessages(library(edgeR))
y<-read_expression(file.path(source,'expression.csv.gz'))
need(identical(colnames(y),meta$unit_id),'Source matrix order mismatch')
base_design<-make_design
make_design<-function(d,effect,extra=FALSE,quadratic=FALSE) {
 x<-base_design(d,effect,extra,quadratic)
 if(extended) for(nm in c('PMI','RIN','pH')) x<-cbind(x,setNames(data.frame(d[[nm]]-mean(d[[nm]])),nm))
 as.matrix(x)
}
old_fit<-fit_model
members<-list(); backgrounds<-list(); conditions<-list(); sensitivity_rows<-list()
for(context in sort(unique(meta$context))) {
 s<-meta$context==context & meta$n_cells>=20 & meta$clinical_group %in% c('Dementia','No_dementia')
 cc<-complete.cases(meta[,c('age','sex','pathology','PMI','RIN','pH')])
 members[[context]]<-data.frame(unit_id=meta$unit_id[s],context=context,
  clinical_group=meta$clinical_group[s],age=meta$age[s],
  base_covariates_complete=complete.cases(meta[s,c('age','sex','pathology')]),
  missing_PMI=is.na(meta$PMI[s]),missing_RIN=is.na(meta$RIN[s]),missing_pH=is.na(meta$pH[s]),
  included=cc[s],reason=ifelse(cc[s],'included','missing_covariate'))
 d<-meta[s & cc,,drop=FALSE]; z<-y[,s & cc,drop=FALSE]
 gate<-nrow(d)>=10 && all(table(factor(d$clinical_group,levels=c('Dementia','No_dementia')))>=5)
 keep<-rep(FALSE,nrow(z))
 if(gate) {
  extended<-FALSE
  x0<-make_design(d,as.numeric(d$clinical_group=='Dementia'))
  complete<-rowSums(is.finite(z))==ncol(z)
  need(all(z[complete,]>=0) && all(abs(z[complete,]-round(z[complete,]))<1e-8),'Expected raw integer pseudobulk counts')
  keep[complete]<-filterByExpr(DGEList(z[complete,,drop=FALSE]),design=x0)
 }
 backgrounds[[context]]<-data.frame(context=context,gene_symbol=rownames(z),included=keep)
 input<-NULL
 if(sum(keep)>=100 && gate) {
  dg<-calcNormFactors(DGEList(z[keep,,drop=FALSE]))
  need(all(dg$samples$lib.size>0),'Empty library')
  # Fixed voom weights and filtering isolate coefficient adjustment sensitivity.
  input<-voom(dg,x0,plot=FALSE,span=0.5,adaptive.span=TRUE)
 }
 for(ext in c(FALSE,TRUE)) for(model in c('clinical_status','braak_same_clinical_population','clinical_status_adjusted_for_braak')) {
  extended<-ext
  effect<-if(model=='braak_same_clinical_population') d$pathology else as.numeric(d$clinical_group=='Dementia')
  x<-make_design(d,effect,extra=model=='clinical_status_adjusted_for_braak')
  rank<-if(nrow(x)) qr(x)$rank else 0
  ready<-gate && !is.null(input) && rank==ncol(x) && nrow(x)-rank>=5 && (model!='braak_same_clinical_population' || length(unique(effect))>=3)
  g<-data.frame(gene_symbol=rownames(z),beta=NA_real_,SE=NA_real_,CI_low=NA_real_,CI_high=NA_real_,p=NA_real_,status=ifelse(keep,'not_testable_design','filtered_background'))
  fitted_input<-NULL
  if(ready) {
   lf<-lmFit(input,x)
   good<-is.finite(lf$sigma) & lf$sigma>1e-10 & lf$df.residual>=5
   if(sum(good)>=100) {
    ff<-eBayes(lf[good,],trend=FALSE,robust=FALSE)
    se<-ff$stdev.unscaled[,'effect']*sqrt(ff$s2.post); beta<-ff$coefficients[,'effect']; ci<-qt(.975,ff$df.total)*se
    pos<-which(keep)[good]
    g[pos,c('beta','SE','CI_low','CI_high','p')]<-cbind(beta,se,beta-ci,beta+ci,ff$p.value[,'effect'])
    g$status[keep]<-'unassessable_residual';g$status[pos]<-'tested';fitted_input<-input[good,,drop=FALSE]
   }
  }
  fit<-list(genes=g,input=fitted_input,design=x,gate=data.frame(n=nrow(d),rank=rank,columns=ncol(x),residual_df=nrow(x)-rank,condition_number=if(ready) kappa(x) else NA_real_,scaled_condition_number=sensitivity_scaled_kappa(x),condition_review_required=sensitivity_scaled_kappa(x)>30,status=if(ready) 'design_ready' else 'not_testable_design'))
  arm<-if(ext) 'quality_extended' else 'matched_base'
  record_model(fit,d,paste(context,arm,model,sep='__'),context,paste(arm,model,sep='::'),if(model=='braak_same_clinical_population') 'log2_expression_per_Braak_stage' else 'log2_Dementia_minus_No_dementia')
  wr(as.data.frame(x),paste0('design_',match(context,sort(unique(meta$context))),'_',arm,'_',model,'.csv'))
  # Each variant and adjustment arm has its own reserved family.
  for(variant in c('design_specific_voom','fixed_weights_robust_ebayes')) {
   alt_input<-input; weight_error<-''
   if(variant=='design_specific_voom' && ready) {
    alt_input<-tryCatch(voom(dg,x,plot=FALSE,span=0.5,adaptive.span=TRUE),
      error=function(e){weight_error<<-conditionMessage(e);NULL})
   }
   a<-sensitivity_fixed_fit(if(ready) alt_input else NULL,x,'effect',ep$gene_symbol,
                      robust=variant=='fixed_weights_robust_ebayes')
   if(nzchar(weight_error)) {a$status<-'sensitivity_weight_error';a$error<-weight_error}
   a$context<-context;a$arm<-arm;a$model<-model;a$variant<-variant
   a$reference_beta<-g$beta[match(a$gene_symbol,g$gene_symbol)]
   a$reference_SE<-g$SE[match(a$gene_symbol,g$gene_symbol)]
   a$beta_change<-a$beta-a$reference_beta
   a$SE_ratio<-ifelse(a$reference_SE>0,a$SE/a$reference_SE,NA_real_)
   sensitivity_rows[[length(sensitivity_rows)+1]]<-a
  }
  if(model=='clinical_status_adjusted_for_braak') {
   xf<-x[,colnames(x)!='braak',drop=FALSE]
   if(nrow(d) && length(unique(d$pathology))>1) {
    bf<-model.matrix(~factor(d$pathology))[,-1,drop=FALSE]
    colnames(bf)<-paste0('braak_category_',seq_len(ncol(bf)));xf<-cbind(xf,bf)
   }
   a<-sensitivity_fixed_fit(if(ready) input else NULL,xf,'effect',ep$gene_symbol)
   a$context<-context;a$arm<-arm;a$model<-model;a$variant<-'braak_factor_fixed_weights'
   a$reference_beta<-g$beta[match(a$gene_symbol,g$gene_symbol)]
   a$reference_SE<-g$SE[match(a$gene_symbol,g$gene_symbol)]
   a$beta_change<-a$beta-a$reference_beta
   a$SE_ratio<-ifelse(a$reference_SE>0,a$SE/a$reference_SE,NA_real_)
   sensitivity_rows[[length(sensitivity_rows)+1]]<-a
  }
 }
}
finish_models()
wr(do.call(rbind,members),'217_sample_membership.csv')
wr(do.call(rbind,backgrounds),'217_common_background.csv')
sensitivity_write_family(sensitivity_rows,file.path(out,'quality_variants.csv'),c('variant','arm'))
membership<-do.call(rbind,members)
flow<-lapply(split(membership,list(membership$context,membership$clinical_group),drop=TRUE),function(z)
 data.frame(context=z$context[1],clinical_group=z$clinical_group[1],eligible_before_quality=nrow(z),
 base_covariates_complete=sum(z$base_covariates_complete),complete_case_included=sum(z$included),
 missing_PMI=sum(z$missing_PMI),missing_RIN=sum(z$missing_RIN),missing_pH=sum(z$missing_pH),
 mean_age_included=if(any(z$included)) mean(z$age[z$included],na.rm=TRUE) else NA_real_,
 mean_age_excluded=if(any(!z$included)) mean(z$age[!z$included],na.rm=TRUE) else NA_real_))
wr(do.call(rbind,flow),'217_quality_selection_flow.csv')
# Same-sized planned BH families, one family for each adjustment arm.
for(file in c('all_fixed_endpoint_results.csv','all_signature_results.csv')) {
 t<-read.csv(file.path(out,file),check.names=FALSE)
 t$arm<-sub('::.*','',t$model);t$model<-sub('^[^:]+::','',t$model)
 for(arm in unique(t$arm)) {
  ix<-which(t$arm==arm);q<-adjust_family(t[ix,,drop=FALSE]);t$joint_q[ix]<-q$joint_q;t$family_size[ix]<-length(ix)
 }
 wr(t,file)
}
