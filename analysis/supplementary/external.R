extra <- Sys.getenv("R_LIBS", "")
if(nzchar(extra)) .libPaths(c(strsplit(extra,.Platform$path.sep,fixed=TRUE)[[1]],.libPaths()))
if(!('adaptive.span' %in% names(formals(limma::voom))))stop('limma::voom lacks adaptive.span; use the recorded compatible environment')

args <- commandArgs(TRUE); d <- args[1]; kind <- args[2]; form <- args[3]; coefname <- args[4]
suppressPackageStartupMessages({library(limma);library(edgeR);library(statmod)})
writeLines(capture.output(sessionInfo()),file.path(d,"R_session.txt"))
ep <- read.csv(file.path(d,"endpoints.csv"),check.names=FALSE,stringsAsFactors=FALSE)
m <- read.csv(file.path(d,"metadata.csv"),check.names=FALSE,stringsAsFactors=FALSE,colClasses=c(unit_id="character"))
out <- ep
for(k in c("logFC","CI_low","CI_high","moderated_SE","P_Value","all_gene_BH"))out[[k]] <- NA_real_
out$testable <- FALSE; out$reason <- "not_in_tested_background"
finish <- function(status,reason,rank=NA,df=NA) {
  if(status!="fitted") out$reason <<- reason
  write.csv(out,file.path(d,"fixed_endpoint_results.csv"),row.names=FALSE,na="")
  write.csv(data.frame(status=status,reason=reason,n=nrow(m),rank=rank,residual_df=df,formula=form),
            file.path(d,"design_status.csv"),row.names=FALSE)
}
if(nrow(m)<6 || length(unique(m$case))!=2 || min(table(m$case))<3) {
  finish("not_testable","fewer_than_three_units_per_group");quit(status=0)
}
term_factors <- attr(terms(as.formula(form)), "factors")
has_case_treat <- all(c("case","treat") %in% rownames(term_factors)) &&
  any(term_factors["case",]>0 & term_factors["treat",]>0)
if(has_case_treat && (length(table(m$case,m$treat))!=4 || min(table(m$case,m$treat))<3)) {
  finish("not_testable","fewer_than_three_units_per_factorial_cell");quit(status=0)
}
design <- model.matrix(as.formula(form),m)
write.csv(cbind(unit_id=m$unit_id,as.data.frame(design)),file.path(d,"design.csv"),row.names=FALSE)
r <- qr(design)$rank; df <- nrow(design)-r
if(r<ncol(design) || df<3 || !coefname %in% colnames(design)) {
  finish("not_testable","nonidentifiable_or_insufficient_residual_df",r,df);quit(status=0)
}
Y <- as.matrix(read.csv(gzfile(file.path(d,"expression.csv.gz")),row.names=1,check.names=FALSE))
stopifnot(identical(colnames(Y),m$unit_id),all(is.finite(Y)),!anyDuplicated(rownames(Y)))
if(kind %in% c("counts","estimated_counts")) {
  stopifnot(all(Y>=0),all(colSums(Y)>0))
  if(kind=="counts")stopifnot(all(Y==floor(Y)))
  y <- DGEList(Y);keep <- rowSums(cpm(y)>=1)>=3
  bg <- data.frame(gene_symbol=rownames(Y),retained=keep,rule="CPM>=1_in_at_least_3_units")
  if(sum(keep)<100) {write.csv(bg,file.path(d,"background.csv"),row.names=FALSE);finish("not_testable","insufficient_background",r,df);quit(status=0)}
  y <- calcNormFactors(y[keep,,keep.lib.sizes=FALSE]);v <- voom(y,design,plot=FALSE,span=0.5,adaptive.span=TRUE)
  fit <- eBayes(lmFit(v,design),robust=TRUE)
} else {
  keep <- apply(Y,1,var)>0
  bg <- data.frame(gene_symbol=rownames(Y),retained=keep,rule="finite_nonzero_variance")
  if(sum(keep)<100) {write.csv(bg,file.path(d,"background.csv"),row.names=FALSE);finish("not_testable","insufficient_background",r,df);quit(status=0)}
  fit <- eBayes(lmFit(Y[keep,,drop=FALSE],design),trend=TRUE,robust=TRUE)
}
write.csv(bg,file.path(d,"background.csv"),row.names=FALSE)
j <- match(coefname,colnames(design));beta <- fit$coefficients[,j]
se <- sqrt(fit$s2.post)*fit$stdev.unscaled[,j];crit <- qt(.975,fit$df.total)
tab <- data.frame(gene_symbol=rownames(fit),logFC=beta,CI_low=beta-crit*se,CI_high=beta+crit*se,
                  moderated_SE=se,P_Value=fit$p.value[,j],all_gene_BH=p.adjust(fit$p.value[,j],"BH"))
con <- gzfile(file.path(d,"all_gene_models.csv.gz"),"wt");write.csv(tab,con,row.names=FALSE);close(con)
ix <- match(ep$analysis_gene,tab$gene_symbol);valid <- !is.na(ix)
for(k in c("logFC","CI_low","CI_high","moderated_SE","P_Value","all_gene_BH"))out[[k]][valid] <- tab[[k]][ix[valid]]
out$testable <- valid;out$reason[valid] <- "tested"
finish("fitted","",r,df)
try({
  z <- out[out$human_gene_symbol %in% c("GFAP","GSDMD","IGFBP5"),]
  if(any(z$testable)) {
    pdf(file.path(d,"fixed_candidates.pdf"),width=7,height=4)
    lim <- range(c(0,z$CI_low,z$CI_high),na.rm=TRUE)
    plot(z$logFC,seq_len(nrow(z)),xlim=lim,ylim=c(.5,nrow(z)+.5),yaxt="n",ylab="",xlab="Model coefficient; pointwise 95% CI",pch=19)
    axis(2,seq_len(nrow(z)),z$human_gene_symbol,las=1);abline(v=0,lty=2,col="grey")
    segments(z$CI_low,seq_len(nrow(z)),z$CI_high,seq_len(nrow(z)));dev.off()
  }
},silent=TRUE)

