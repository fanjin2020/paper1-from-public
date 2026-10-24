extra <- Sys.getenv("R_LIBS", "")
if(nzchar(extra)) .libPaths(c(strsplit(extra,.Platform$path.sep,fixed=TRUE)[[1]],.libPaths()))

a <- commandArgs(TRUE); mode <- a[1]
options(timeout=7200)
for(p in c("limma","statmod")) if(!requireNamespace(p,quietly=TRUE)) stop(paste("Missing R package:",p))
need <- function(ok,msg) if(!isTRUE(ok)) stop(msg)
fit_one <- function(y,x,block=NULL,rho=NULL) {
 need(qr(x)$rank==ncol(x),"rank_deficient")
 f <- limma::lmFit(y,x,block=block,correlation=rho)
 j <- match("genotype",colnames(x)); df <- as.numeric(f$df.residual[1])
 # limma may retain the coefficient name on a scalar. Strip it before naming the result.
 b <- as.numeric(f$coefficients[1,j]); se <- as.numeric(f$stdev.unscaled[1,j]*f$sigma[1])
 need(is.finite(b)&&is.finite(se)&&se>0&&is.finite(df)&&df>0,"invalid_fit")
 result <- c(beta=b,se=se,df=df,low=b-qt(.975,df)*se,high=b+qt(.975,df)*se,p=2*pt(-abs(b/se),df))
 need(identical(names(result),c("beta","se","df","low","high","p"))&&all(is.finite(result)),"Invalid fit result schema")
 result
}
if(mode=="preflight") {
 set.seed(325); n<-48; cage<-factor(rep(1:12,each=4)); g<-rep(c(0,0,1,1),12)
 x<-cbind(Intercept=1,genotype=g); xf<-cbind(x,model.matrix(~cage)[,-1,drop=FALSE])
 y<-matrix(rnorm(64*n),64,n)+matrix(rnorm(64*12,sd=.4),64,12)[,as.integer(cage)]
 o<-fit_one(y,x); f<-fit_one(y,xf); z<-fit_one(y,x,cage,0)
 need(max(abs(o-z))<1e-7,"GLS rho=0 differs from OLS")
 ref<-coef(summary(lm(y[1,]~g)))["g",]
 fixed_ref<-coef(summary(lm(y[1,]~g+cage)))["g",]
 checks<-data.frame(check=c("OLS_beta","OLS_SE","fixed_beta","fixed_SE"),
   observed=unname(c(o['beta'],o['se'],f['beta'],f['se'])),
   reference=unname(c(ref[1],ref[2],fixed_ref[1],fixed_ref[2])))
 checks$absolute_difference<-abs(checks$observed-checks$reference)
 checks$passed<-is.finite(checks$absolute_difference)&checks$absolute_difference<1e-8
 print(checks,digits=16)
 write.csv(checks,paste0(a[2],"_preflight_checks.csv"),row.names=FALSE)
 need(all(checks$passed),paste("OLS reference mismatch; maximum absolute difference:",max(checks$absolute_difference)))
 dc<-limma::duplicateCorrelation(y,x,block=cage)
 rr<-dc$consensus.correlation; need(is.finite(rr)&&rr> -1/3&&rr<1,"Invalid consensus correlation")
 fit_one(y,x,cage,rr)
 write.csv(data.frame(package=c("R","limma","statmod"),version=c(as.character(getRversion()),
  as.character(packageVersion("limma")),as.character(packageVersion("statmod")))),a[2],row.names=FALSE)
 capture.output(sessionInfo(),file=a[3])
 capture.output(limma::duplicateCorrelation,file=paste0(a[3],"_duplicateCorrelation.txt"))
 cat("PREFLIGHT_PASS\n");quit(status=0)
}
need(mode=="simulate","Unknown mode")
cfg<-read.csv(a[2],stringsAsFactors=FALSE,check.names=FALSE)
x<-as.matrix(read.csv(cfg$base_design,check.names=FALSE));storage.mode(x)<-"double"
xf<-as.matrix(read.csv(cfg$fixed_design,check.names=FALSE));storage.mode(xf)<-"double"
u<-read.csv(cfg$units,stringsAsFactors=FALSE);cage<-factor(u$cage)
n<-nrow(x); G<-cfg$background_genes; reps<-cfg$reps
g<-x[,"genotype"];cg<-ave(g,cage,FUN=mean); shift<-cg-mean(g)
beta<-rep(c(0,.3),length.out=G);beta[1]<-cfg$truth
signal<-outer(beta,g)+matrix(rep(cfg$shift*shift,each=G),G,n)
lower_pd<- -1/(max(table(cage))-1)
base_ok<-qr(x)$rank==ncol(x)&&n>ncol(x)+1
fixed_ok<-as.character(cfg$fixed_estimable)%in%c("True","TRUE","true","1")&&qr(xf)$rank==ncol(xf)&&n>ncol(xf)+1
rows<-vector("list",3*reps); counter<-0L
for(i in seq_len(reps)) {
 set.seed(as.integer(cfg$seed+i))
 y<-signal+sqrt(1-cfg$rho)*matrix(rnorm(G*n),G,n)+
    sqrt(cfg$rho)*matrix(rnorm(G*nlevels(cage)),G,nlevels(cage))[,as.integer(cage)]
 rho<-NA_real_; finite_fraction<-NA_real_; clip_fraction<-NA_real_; warning_text<-character()
 cor_error<-""
 tryCatch(withCallingHandlers({
   dc<-limma::duplicateCorrelation(y,x,block=cage)
   rho<-dc$consensus.correlation; raw<-tanh(dc$atanh.correlations)
   finite_fraction<-mean(is.finite(raw))
   need(finite_fraction>=.95,"less_than_95_percent_correlations_estimable")
   clip_fraction<-mean(raw[is.finite(raw)]<=lower_pd+.0100001|raw[is.finite(raw)]>=.9899999)
   need(is.finite(rho)&&rho>lower_pd&&rho<1,"correlation_not_positive_definite")
 },warning=function(w){warning_text<<-c(warning_text,conditionMessage(w));invokeRestart("muffleWarning")}),
 error=function(e)cor_error<<-conditionMessage(e))
 for(model in c("no_cage","cage_fixed","cage_correlation")) {
  counter<-counter+1L; status<-"ok"; reason<-""; values<-rep(NA_real_,7)
  names(values)<-c("beta","se","df","low","high","p","unused")
  if(model=="cage_fixed"&&!fixed_ok) {status<-"not_estimable";reason<-"genotype_not_identifiable_or_insufficient_df"}
  if(model=="no_cage"&&!base_ok) {status<-"not_estimable";reason<-"base_design_rank_or_df"}
  if(model=="cage_correlation"&&nzchar(cor_error)) {status<-"failed";reason<-cor_error}
  if(status=="ok") tryCatch(withCallingHandlers({
   values<-if(model=="no_cage") fit_one(y,x) else if(model=="cage_fixed") fit_one(y,xf) else fit_one(y,x,cage,rho)
  },warning=function(w){warning_text<<-c(warning_text,conditionMessage(w));invokeRestart("muffleWarning")}),
  error=function(e){status<<-"failed";reason<<-conditionMessage(e)})
  rows[[counter]]<-data.frame(scenario_id=cfg$scenario_id,replicate=i,model=model,status=status,reason=reason,
   beta=unname(values['beta']),se=unname(values['se']),df=unname(values['df']),
   ci_low=unname(values['low']),ci_high=unname(values['high']),p=unname(values['p']),
   rho_estimate=if(model=="cage_correlation")rho else NA_real_,
   correlation_finite_fraction=if(model=="cage_correlation")finite_fraction else NA_real_,
   gene_correlation_near_clip_fraction=if(model=="cage_correlation")clip_fraction else NA_real_,
   consensus_near_pd_boundary=if(model=="cage_correlation"&&is.finite(rho))min(rho-lower_pd,1-rho)<.02 else NA,
   warning=paste(unique(warning_text),collapse=" | "))
 }
 if(i%%50==0) cat("replicate",i,"of",reps,"\n")
}
z<-do.call(rbind,rows);tmp<-paste0(cfg$output,".tmp")
con<-gzfile(tmp,"wt");write.csv(z,con,row.names=FALSE,na="");close(con)
need(file.rename(tmp,cfg$output),"Could not publish scenario result")
