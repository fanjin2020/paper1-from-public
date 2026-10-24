extra <- Sys.getenv("R_LIBS", "")
if(nzchar(extra)) .libPaths(c(strsplit(extra,.Platform$path.sep,fixed=TRUE)[[1]],.libPaths()))

scaled_condition_number <- function(X) {
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

suppressPackageStartupMessages(library(limma))
a<-commandArgs(TRUE); h<-a[1];out<-rawToChar(as.raw(strtoi(substring(h,seq(1,nchar(h),2),seq(2,nchar(h),2)),16L)));Encoding(out)<-'UTF-8'
options(stringsAsFactors=FALSE)
wr<-function(t,n) write.csv(t,file.path(out,n),row.names=FALSE,na='')
m<-read.csv(file.path(out,'230_all_samples.csv'))
y<-as.matrix(read.csv(gzfile(file.path(out,'230_expression_log2.csv.gz')),row.names=1,check.names=FALSE))
ep<-read.csv(file.path(out,'229_fixed_endpoints.csv'))
stopifnot(identical(colnames(y),m$sample_id),all(is.finite(y)))
c0<-function(x) x-mean(x)
x0<-cbind(intercept=1,MMSE_per_5=c0(m$mmse)/5,age=c0(m$age),male=as.numeric(m$sex=='M'),PMI=c0(m$pmi))
designs<-list(demographic=x0,braak_adjusted=cbind(x0,Braak=c0(m$braak)),nft_adjusted=cbind(x0,log1p_NFT=c0(log1p(m$nft))))
safe_cor<-function(a,b) if(length(unique(a))>1 && length(unique(b))>1) cor(a,b) else NA_real_
wr(data.frame(pair=c('MMSE_Braak','MMSE_log1p_NFT'),
 pearson_r=c(safe_cor(m$mmse,m$braak),safe_cor(m$mmse,log1p(m$nft))),
 scope='descriptive_correlation_not_a_causal_or_independence_test'),'231_covariate_correlations.csv')
gs<-ss<-ds<-list()
for(nm in names(designs)) {
 x<-designs[[nm]];rank<-qr(x)$rank;df<-nrow(x)-rank
 ready<-rank==ncol(x) && df>=10 && nrow(x)>=25
 ds[[nm]]<-data.frame(model=nm,n=nrow(x),rank=rank,columns=ncol(x),residual_df=df,
  condition_number=kappa(x,exact=TRUE),scaled_condition_number=scaled_condition_number(x),
  condition_review_required=scaled_condition_number(x)>30,status=if(ready)'ready' else 'not_estimable')
 wr(data.frame(sample_id=m$sample_id,x,check.names=FALSE),paste0('231_design_',nm,'.csv'))
 g<-data.frame(model=nm,gene_symbol=ep$gene_symbol,discovery_direction=ep$discovery_direction,beta=NA_real_,SE=NA_real_,CI_low=NA_real_,CI_high=NA_real_,p=NA_real_,status='unmapped_or_not_estimable')
 z<-NULL
 if(ready) {
  f<-lmFit(y,x);good<-is.finite(f$sigma)&f$sigma>1e-10
  if(sum(good)>=1000) {
   z<-y[good,,drop=FALSE];f<-eBayes(f[good,],trend=TRUE,robust=FALSE)
   se<-f$stdev.unscaled[,2]*sqrt(f$s2.post);b<-f$coefficients[,2];ci<-qt(.975,f$df.total)*se
   ix<-match(ep$gene_symbol,rownames(z));ok<-!is.na(ix)
   g[ok,c('beta','SE','CI_low','CI_high','p')]<-cbind(b,se,b-ci,b+ci,f$p.value[,2])[ix[ok],,drop=FALSE];g$status[ok]<-'tested'
   wr(data.frame(gene_symbol=rownames(z),beta=b,SE=se,p=f$p.value[,2]),paste0('231_all_gene_',nm,'.csv'))
  }
 }
 gs[[nm]]<-g
 for(d in c('Up','Down')) {
  idx<-which(rownames(z) %in% ep$gene_symbol[ep$discovery_direction==d])
  s<-data.frame(model=nm,signature=d,fixed_size=sum(ep$discovery_direction==d),measured_size=length(idx),Direction=NA_character_,Correlation=NA_real_,p=NA_real_,status='not_testable')
  if(!is.null(z) && length(idx)>=15 && nrow(z)-length(idx)>=15) {
   cc<-camera(z,list(fixed=idx),design=x,contrast=2,use.ranks=FALSE,inter.gene.cor=NA,allow.neg.cor=FALSE,trend.var=TRUE)
   s$Direction<-cc$Direction;s$Correlation<-cc$Correlation;s$p<-cc$PValue;s$status<-'tested'
  }
  ss[[paste(nm,d)]]<-s
 }
}
adj<-function(t){ok<-is.finite(t$p);p<-t$p;p[!ok]<-1;t$joint_q<-p.adjust(p,'BH');t$joint_q[!ok]<-NA;t$family_size<-nrow(t);t}
wr(do.call(rbind,ds),'231_design_status.csv')
wr(adj(do.call(rbind,gs)),'231_all_fixed_gene_results.csv')
wr(adj(do.call(rbind,ss)),'232_all_fixed_signature_results.csv')
writeLines(capture.output(sessionInfo()),file.path(out,'R_session.txt'))

