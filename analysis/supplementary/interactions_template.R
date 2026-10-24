extra <- Sys.getenv("R_LIBS", "")
if(nzchar(extra)) .libPaths(c(strsplit(extra,.Platform$path.sep,fixed=TRUE)[[1]],.libPaths()))

suppressPackageStartupMessages({library(edgeR);library(limma)})
a<-commandArgs(TRUE); d<-a[1]; mode<-a[2]
need<-function(ok,msg)if(!isTRUE(ok))stop(msg)
m<-read.csv(file.path(d,'metadata.csv'),colClasses=c(donor='character',unit_id='character'))
ep<-read.csv(file.path(d,'endpoints.csv'),colClasses='character')
need(!anyNA(m$donor)&&all(nzchar(trimws(m$donor)))&&!anyDuplicated(m$unit_id),'Invalid donor/unit identities')
need(nrow(ep)>0&&!anyDuplicated(ep$gene_symbol),'Invalid endpoint family')
need(all(m$case %in% c(0,1))&&length(unique(m$case))==2&&all(m$sex %in% c('F','M')),'Invalid case/sex coding')
need(all(is.finite(as.matrix(m[,c('age','pmi','ph_cerebellum')]))),'Missing/nonfinite quality covariate')
for(field in c('case','sex','age','pmi','ph_cerebellum')){
 consistent<-vapply(split(m[[field]],m$donor),function(z){
  if(is.numeric(z))all(abs(z-z[1])<=1e-12*(1+abs(z[1]))) else length(unique(z))==1
 },logical(1))
 need(all(consistent),paste('Conflicting paired covariates:',field))
}
u<-m[!duplicated(m$donor),,drop=FALSE]
need(min(table(factor(u$case,levels=c(0,1))))>=8,'Insufficient donors per group')
# Reconstruct the design independently from metadata, not from exported coefficients.
B<-cbind(Intercept=rep(1,nrow(m)),case=m$case)
for(field in c('age','pmi','ph_cerebellum')){
 sd0<-sd(u[[field]])
 if(sd0>0){B<-cbind(B,(m[[field]]-mean(u[[field]]))/sd0);colnames(B)[ncol(B)]<-field}
}
if(length(unique(m$sex))>1)B<-cbind(B,male=as.numeric(m$sex=='M'))
readcount<-function(name){
 z<-read.csv(gzfile(file.path(d,name)),check.names=FALSE,colClasses=c(gene_symbol='character'))
 Y<-as.matrix(z[,-1,drop=FALSE]);rownames(Y)<-z[[1]];storage.mode(Y)<-'double'
 need(!anyDuplicated(rownames(Y))&&!anyDuplicated(colnames(Y))&&all(is.finite(Y))&&all(Y>=0)&&all(Y==floor(Y))&&all(colSums(Y)>0),'Invalid counts')
 Y
}
dd<-read.csv(file.path(d,'design.csv'),check.names=FALSE,colClasses=c(unit_id='character'))
X<-as.matrix(dd[,-1,drop=FALSE]);rownames(X)<-dd$unit_id
need(identical(dd$unit_id,m$unit_id)&&qr(X)$rank==ncol(X)&&nrow(X)-ncol(X)>=8,'Invalid design')
expected<-B
if(mode=='interaction'){
 need(setequal(unique(m$cell_type),c('Ast','Oli')),'Requires exactly Ast and Oli labels')
 Zast<-B*as.numeric(m$cell_type=='Ast');colnames(Zast)<-paste0('ast_',colnames(B))
 expected<-cbind(B,Zast)
}
need(identical(colnames(X),colnames(expected))&&isTRUE(all.equal(unname(X),unname(expected),tolerance=1e-10)),'Design differs from metadata coding')
exportfit<-function(fit,coef,arm){
 tt<-topTable(fit,coef=coef,number=Inf,sort.by='none',confint=.95)
 se<-fit$stdev.unscaled[,coef]*sqrt(fit$s2.post)
 k<-match(ep$analysis_symbol,rownames(tt));k[ep$mapping_status!='mapped'|is.na(ep$mapping_status)]<-NA_integer_
 z<-data.frame(arm=arm,gene_symbol=ep$gene_symbol,testable=!is.na(k),logFC=tt$logFC[k],
 CI_low=tt$CI.L[k],CI_high=tt$CI.R[k],P_Value=tt$P.Value[k],
 SE_moderated=se[k],df_total=rep_len(fit$df.total,nrow(tt))[k],n_donors=length(unique(m$donor)))
 write.csv(data.frame(gene_symbol=rownames(tt),tt),gzfile(file.path(d,paste0('all_genes_',arm,'.csv.gz'))),row.names=FALSE)
 z
}
norm<-function(Y,keep){calcNormFactors(DGEList(Y[keep,,drop=FALSE]),method='TMM')}
results<-list()
if(mode=='annotation'){
 need(!anyDuplicated(m$donor),'Repeated donor')
 A<-readcount('original.csv.gz'); B<-readcount('consensus.csv.gz')
 need(identical(dimnames(A),dimnames(B))&&identical(colnames(A),m$unit_id),'Counts not paired')
 ka<-filterByExpr(DGEList(A),group=m$case);kb<-filterByExpr(DGEList(B),group=m$case);keep<-ka&kb
 write.csv(data.frame(gene_symbol=rownames(A),original=ka,consensus=kb,common=keep),file.path(d,'common_background.csv'),row.names=FALSE)
 need(sum(keep)>=500,'Fewer than 500 shared genes')
 for(label in c('original','consensus')){
  dy<-norm(if(label=='original') A else B,keep)
  write.csv(dy$samples,file.path(d,paste0('normalization_',label,'.csv')))
  v<-voom(dy,X,plot=FALSE,span=.5,adaptive.span=TRUE);f<-lmFit(v,X)
  for(rb in c(FALSE,TRUE)){
   arm<-paste0('matched_',label,if(rb)'_robust' else '')
   results[[arm]]<-exportfit(eBayes(f,trend=FALSE,robust=rb),'case',arm)
  }
 }
}else if(mode=='interaction'){
 Y<-readcount('counts.csv.gz');need(identical(colnames(Y),m$unit_id),'Interaction sample mismatch')
 need(all(table(m$donor)==2)&&all(table(m$donor,m$cell_type)==1),'Requires paired Ast and Oli')
 ia<-m$cell_type=='Ast';io<-m$cell_type=='Oli'
 ka<-filterByExpr(DGEList(Y[,ia,drop=FALSE]),group=m$case[ia])
 ko<-filterByExpr(DGEList(Y[,io,drop=FALSE]),group=m$case[io]);keep<-ka&ko
 write.csv(data.frame(gene_symbol=rownames(Y),Ast=ka,Oli=ko,retained=keep),file.path(d,'background.csv'),row.names=FALSE)
 need(sum(keep)>=500,'Fewer than 500 interaction genes')
 dy<-DGEList(Y[keep,,drop=FALSE])
 for(ct in unique(m$cell_type)){
  ix<-which(m$cell_type==ct)
  dy$samples$norm.factors[ix]<-calcNormFactors(dy[,ix],method='TMM')$samples$norm.factors
 }
 write.csv(dy$samples,file.path(d,'normalization_within_celltype.csv'))
 v0<-voom(dy,X,plot=FALSE,span=.5,adaptive.span=TRUE)
 dc<-duplicateCorrelation(v0,X,block=m$donor);rho<-dc$consensus.correlation
 need(is.finite(rho)&&abs(rho)<.99,'Unusable within-donor correlation')
 v<-voom(dy,X,block=m$donor,correlation=rho,plot=FALSE,span=.5,adaptive.span=TRUE)
 fit<-eBayes(lmFit(v,X,block=m$donor,correlation=rho),trend=FALSE,robust=TRUE)
 # Equivalent two-cell-type coefficient parameterization independently checks the contrast.
 ZA<-B*as.numeric(ia);ZO<-B*as.numeric(io)
 colnames(ZA)<-paste0('Ast_',colnames(B));colnames(ZO)<-paste0('Oli_',colnames(B))
 Z<-cbind(ZO,ZA);C<-setNames(rep(0,ncol(Z)),colnames(Z));C['Ast_case']<-1;C['Oli_case']<--1
 # Fit the contrast directly: contrasts.fit can approximate SEs with gene-specific weights.
 Zc<-contrastAsCoef(Z,matrix(C,ncol=1))
 alt<-eBayes(lmFit(v,Zc$design,block=m$donor,correlation=rho),trend=FALSE,robust=TRUE)
 f1<-topTable(fit,coef='ast_case',number=Inf,sort.by='none',confint=.95)
 f2<-topTable(alt,coef=Zc$coef,number=Inf,sort.by='none',confint=.95)
 checks<-lapply(c('logFC','t','P.Value','CI.L','CI.R'),function(field){
  a<-f1[[field]];b<-f2[[field]]
  pass<-all(is.finite(a))&&all(is.finite(b))&&all(abs(a-b)<=1e-9+1e-7*pmax(abs(a),abs(b)))
  data.frame(field=field,pass=pass,max_absolute_difference=max(abs(a-b)))
 })
 checks<-do.call(rbind,checks)
 write.csv(checks,file.path(d,'independent_contrast_verification.csv'),row.names=FALSE)
 need(all(checks$pass),'Independent interaction contrast check failed')
 results[[1]]<-exportfit(fit,'ast_case','Ast_minus_Oli_interaction')
 write.csv(data.frame(rho=rho,donors=length(unique(m$donor)),libraries=nrow(m),
  contrast='AD_minus_CT_Ast_minus_AD_minus_CT_Oli'),file.path(d,'paired_correlation.csv'),row.names=FALSE)
}else stop('Unknown mode')
write.csv(do.call(rbind,results),file.path(d,'fixed_results.csv'),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(d,'R_session.txt'))

