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
if(!('adaptive.span' %in% names(formals(limma::voom))))stop('limma::voom lacks adaptive.span; use the recorded compatible environment')

args <- commandArgs(TRUE); d <- args[1]
suppressPackageStartupMessages({library(edgeR); library(limma); library(statmod)})
writeLines(capture.output(sessionInfo()), file.path(d,"R_session.txt"))
m <- read.csv(file.path(d,"metadata.csv"),check.names=FALSE,stringsAsFactors=FALSE)
ep <- read.csv(file.path(d,"endpoints.csv"),stringsAsFactors=FALSE)
cfg <- readLines(file.path(d,"formulas.txt"))
mode <- cfg[1]; forms <- strsplit(cfg[2],";",fixed=TRUE)[[1]]; coeff <- cfg[3]
empty <- function(reason) {
  tab <- expand.grid(gene_symbol=ep$gene_symbol,model=seq_along(forms),stringsAsFactors=FALSE)
  for(k in c("logFC","CI_low","CI_high","P_Value","all_gene_BH"))tab[[k]]<-NA_real_
  tab$testable <- FALSE; tab$reason <- reason
  write.csv(tab,file.path(d,"endpoint_results.csv"),row.names=FALSE,na="")
  write.csv(data.frame(status="not_testable",reason=reason,n=nrow(m)),file.path(d,"status.csv"),row.names=FALSE)
  quit(status=0)
}
if(nrow(m)<8 || length(unique(m$case))!=2 || min(table(unique(m[,c("donor","case")])$case))<4)empty("fewer_than_four_donors_per_group")
des <- tryCatch(lapply(forms,function(f)model.matrix(as.formula(f),m)),
                error=function(e)empty(paste0("design_unavailable: ",conditionMessage(e))))
for(i in seq_along(des)) {
 X <- des[[i]]
 if(nrow(X)!=nrow(m)||any(!is.finite(X)))empty("missing_or_nonfinite_design_values")
 write.csv(cbind(unit_id=m$unit_id,as.data.frame(X)),file.path(d,paste0("design_",i,".csv")),row.names=FALSE)
 if(qr(X)$rank<ncol(X)||nrow(X)-ncol(X)<4||!coeff %in% colnames(X))empty("rank_or_residual_df_or_coefficient")
}
Y <- as.matrix(read.csv(gzfile(file.path(d,"expression.csv.gz")),row.names=1,check.names=FALSE))
stopifnot(identical(colnames(Y),m$unit_id),!anyDuplicated(rownames(Y)),all(is.finite(Y)),all(Y>=0),all(colSums(Y)>0))
y <- DGEList(Y); keep <- rowSums(cpm(y)>=1)>=4
write.csv(data.frame(gene_symbol=rownames(Y),retained=keep),file.path(d,"common_background.csv"),row.names=FALSE)
if(sum(keep)<100)empty("insufficient_expressed_background")
y <- calcNormFactors(y[keep,,keep.lib.sizes=FALSE]); results <- list(); statuses<-list()
for(i in seq_along(des)) {
 X <- des[[i]]; v <- voom(y,X,plot=FALSE,span=0.5,adaptive.span=TRUE); rho <- rho_second <- NA_real_
 if(mode=="paired_celltype") {
   stopifnot(all(table(m$donor)==2))
   corfit <- duplicateCorrelation(v,X,block=m$donor); rho <- corfit$consensus.correlation
   if(!is.finite(rho))empty("donor_correlation_not_estimable")
   v <- voom(y,X,plot=FALSE,span=0.5,adaptive.span=TRUE,block=m$donor,correlation=rho)
   # The first estimate defines this fit; the second is a sensitivity diagnostic.
   rho_second <- tryCatch(duplicateCorrelation(v,X,block=m$donor)$consensus.correlation,
                          error=function(e)NA_real_)
   fit <- eBayes(lmFit(v,X,block=m$donor,correlation=rho),robust=TRUE)
 } else fit <- eBayes(lmFit(v,X),robust=TRUE)
 j <- match(coeff,colnames(X)); b<-fit$coefficients[,j]; se<-sqrt(fit$s2.post)*fit$stdev.unscaled[,j]
 ci<-qt(.975,fit$df.total)*se
 tab<-data.frame(gene_symbol=rownames(fit),logFC=b,CI_low=b-ci,CI_high=b+ci,P_Value=fit$p.value[,j],all_gene_BH=p.adjust(fit$p.value[,j],"BH"))
 con<-gzfile(file.path(d,paste0("all_gene_model_",i,".csv.gz")),"wt");write.csv(tab,con,row.names=FALSE);close(con)
 z<-merge(ep,tab,by="gene_symbol",all.x=TRUE,sort=FALSE); z$model<-i;z$testable<-!is.na(z$P_Value)
 z$reason<-ifelse(z$testable,"tested","not_in_expressed_background");results[[i]]<-z
 statuses[[i]]<-data.frame(status="fitted",model=i,n=nrow(m),donors=length(unique(m$donor)),rank=ncol(X),residual_df=nrow(X)-ncol(X),rho=rho,rho_second_diagnostic=rho_second,
   correlation_review_required=mode=="paired_celltype" && (!is.finite(rho_second) || abs(rho_second-rho)>.1),
   correlation_fit_policy="first_estimate_fit_second_estimate_diagnostic",
   scaled_condition_number=scaled_condition_number(X),formula=forms[i])
}
write.csv(do.call(rbind,results),file.path(d,"endpoint_results.csv"),row.names=FALSE,na="")
write.csv(do.call(rbind,statuses),file.path(d,"status.csv"),row.names=FALSE)


