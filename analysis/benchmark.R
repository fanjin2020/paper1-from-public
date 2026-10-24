
options(stringsAsFactors=FALSE, digits=17, warn=1)
RNGkind('Mersenne-Twister','Inversion','Rejection')
a <- commandArgs(TRUE)
need <- function(x,msg) if(!isTRUE(x)) stop(msg,call.=FALSE)
atomic_json <- function(x,path) {
  tmp<-paste0(path,'.tmp'); jsonlite::write_json(x,tmp,auto_unbox=TRUE,pretty=TRUE,na='null')
  if(file.exists(path)) unlink(path)
  need(file.rename(tmp,path),paste('Cannot publish',path))
}
libraries<-strsplit(Sys.getenv('PAPER1_R_LIBRARIES'),.Platform$path.sep,fixed=TRUE)[[1]]
libraries<-libraries[nzchar(libraries)&dir.exists(libraries)]
.libPaths(unique(c(libraries,.libPaths())))
packages<-c('jsonlite','limma','variancePartition','BiocParallel','lme4','lmerTest')
if(a[1]=='check') {
  setup<-a[3]=='TRUE'; lib<-a[4]; dir.create(lib,recursive=TRUE,showWarnings=FALSE)
  missing<-packages[!vapply(packages,requireNamespace,logical(1),quietly=TRUE)]
  if(length(missing)&&setup) {
    options(repos=c(CRAN=Sys.getenv('PAPER1_CRAN','https://cloud.r-project.org')),timeout=1200)
    cran<-intersect(missing,c('jsonlite','lme4','lmerTest'))
    if(length(cran)) install.packages(cran,lib=lib,Ncpus=2)
    bioc<-setdiff(missing,cran)
    if(length(bioc)) {
      Sys.setenv(BIOCONDUCTOR_ONLINE_VERSION_DIAGNOSIS='FALSE')
      if(!requireNamespace('BiocManager',quietly=TRUE)) install.packages('BiocManager',lib=lib)
      BiocManager::install(bioc,lib=lib,ask=FALSE,update=FALSE,Ncpus=2)
    }
  }
  missing<-packages[!vapply(packages,requireNamespace,logical(1),quietly=TRUE)]
  need(!length(missing),paste('Missing R packages:',paste(missing,collapse=', '),
    '; use --setup or --r-library PATH. No model will be silently skipped.'))
  versions<-setNames(lapply(packages,function(p)list(version=as.character(packageVersion(p)),
    path=find.package(p))),packages)
  jsonlite::write_json(list(R=R.version.string,RNG=RNGkind(),packages=versions),a[2],auto_unbox=TRUE,pretty=TRUE)
  quit(status=0)
}
suppressPackageStartupMessages({library(jsonlite);library(limma);library(variancePartition);library(BiocParallel)})
job<-fromJSON(a[2]); spec<-job$spec
d<-read.csv(job$units,check.names=FALSE)
d$cage<-factor(d$cage,levels=sort(unique(d$cage)))
rownames(d)<-paste0('s',seq_len(nrow(d)))
g<-d$genotype; n<-nrow(d)
base<-data.frame(Intercept=rep(1,n),genotype=g)
if('age'%in%names(d)) {base$age<-d$age-mean(d$age);base$male<-d$male}
x<-as.matrix(base)
basis<-function(x) {
  s<-svd(x,nu=min(dim(x)),nv=0); tol<-max(dim(x))*.Machine$double.eps*s$d[1]
  s$u[,s$d>tol,drop=FALSE]
}
z0<-x[,setdiff(colnames(x),'genotype'),drop=FALSE]
zf<-basis(cbind(z0,model.matrix(~cage-1,d)))
mean_cov<-sapply(setdiff(colnames(x),'Intercept'),function(k)ave(x[,k],d$cage,FUN=mean))
zm<-basis(cbind(z0,mean_cov))
xf<-cbind(genotype=g,zf);colnames(xf)<-c('genotype',paste0('f',seq_len(ncol(zf))))
xm<-cbind(genotype=g,zm);colnames(xm)<-c('genotype',paste0('m',seq_len(ncol(zm))))
ss<-function(z)sum((g-z%*%crossprod(z,g))^2)
fixed_ok<-ss(zf)>1e-9 && n-ncol(xf)>1
mundlak_ok<-ss(zm)>1e-9 && n-ncol(xm)>1
need(fixed_ok==mundlak_ok,'Fixed/Mundlak estimability disagrees')
md<-data.frame(cage=d$cage,row.names=rownames(d))
for(k in setdiff(names(base),'Intercept'))md[[k]]<-base[[k]]
form<-as.formula(paste('~',paste(c(setdiff(names(base),'Intercept'),'(1|cage)'),collapse='+')))
# Use a true formula intercept: a numerically constant SVD column can trigger
# lme4's scale checks despite representing precisely the same nuisance space.
raw_nuisance<-cbind(z0,mean_cov)
centered<-sweep(raw_nuisance,2L,colMeans(raw_nuisance),'-')
center_basis<-basis(centered)*sqrt(n)
mm<-data.frame(cage=d$cage,genotype=g,row.names=rownames(d))
if(ncol(center_basis))for(k in seq_len(ncol(center_basis)))mm[[paste0('m',k)]]<-center_basis[,k]
mform<-as.formula(paste('~',paste(c(setdiff(names(mm),'cage'),'(1|cage)'),collapse='+')))
need(qr(cbind(zm,1,center_basis))$rank==ncol(zm),'Mundlak reparameterization changed nuisance space')
extract<-function(f,official=FALSE) {
  j<-match('genotype',colnames(f$coefficients));need(!is.na(j),'Target coefficient missing')
  beta<-as.numeric(f$coefficients[1,j]);se<-as.numeric(f$stdev.unscaled[1,j]*f$sigma[1])
  df<-as.numeric(if(is.matrix(f$df.residual))f$df.residual[1,j] else f$df.residual[1])
  need(all(is.finite(c(beta,se,df)))&&se>0&&df>0,'Invalid coefficient, SE or df')
  p<-2*pt(-abs(beta/se),df)
  if(official) {
    need(is.finite(f$p.value[1,j])&&abs(p-f$p.value[1,j])<1e-7,'Official dream p-value mismatch')
    errors<-attr(f,'errors'); need(!length(errors),paste('dream errors:',paste(errors,collapse=' | ')))
  }
  c(beta=beta,SE=se,df=df,p=p,low=beta-qt(.975,df)*se,high=beta+qt(.975,df)*se)
}
numeric_messages<-function(w)w[grepl('converg|eigen|Hessian|positive definite',w,ignore.case=TRUE)&
  !grepl('boundary.*singular',w,ignore.case=TRUE)]
shift<-ave(g,d$cage,FUN=mean)-mean(g)
betas<-rep(c(0,.3),length.out=64);betas[1]<-spec$truth
signal<-outer(betas,g)+matrix(rep(spec$shift*shift,each=64),64,n)
records<-list();start_time<-proc.time()[3]
for(rep in seq.int(job$start,job$end-1L)) {
  # Preserve the specified order of animal and cage random draws.
  set.seed(as.integer(spec$seed+rep+1L))
  y<-signal+sqrt(1-spec$rho_truth)*matrix(rnorm(64*n),64,n)+
    sqrt(spec$rho_truth)*matrix(rnorm(64*nlevels(d$cage)),64,nlevels(d$cage))[,as.integer(d$cage)]
  dimnames(y)<-list(paste0('g',seq_len(64)),rownames(d))
  current<-list()
  for(model in c('baseline','cage_fixed','cage_correlation','dream','mundlak_dream')) {
    row<-data.frame(scenario=spec$scenario,layout=spec$layout,rho_truth=spec$rho_truth,
      shift=spec$shift,truth=spec$truth,replicate=rep+1L,seed=spec$seed+rep+1L,model=model,
      beta=NA_real_,SE=NA_real_,df=NA_real_,p=NA_real_,low=NA_real_,high=NA_real_,
      rho_estimated=NA_real_,boundary=NA,status='failed',reason='',warning='',
      mundlak_fixed_beta_difference=NA_real_)
    w<-character()
    if((model=='cage_fixed'&&!fixed_ok)||(model=='mundlak_dream'&&!mundlak_ok)) {
      row$status<-'not_estimable';row$reason<-'No within-cage genotype information'
    } else tryCatch(withCallingHandlers({
      if(model=='baseline')f<-limma::lmFit(y,x)
      if(model=='cage_fixed')f<-limma::lmFit(y,xf)
      if(model=='cage_correlation') {
        rho<-limma::duplicateCorrelation(y,x,block=d$cage)$consensus.correlation
        lower<- -1/(max(table(d$cage))-1)
        need(is.finite(rho)&&rho>lower&&rho<1,'Invalid plug-in cage correlation')
        row$rho_estimated<-rho;row$boundary<-min(rho-lower,1-rho)<1e-6
        f<-limma::lmFit(y,x,block=d$cage,correlation=rho)
      }
      if(model%in%c('dream','mundlak_dream')) {
        f<-variancePartition::dream(y[1,,drop=FALSE],if(model=='dream')form else mform,
          if(model=='dream')md else mm,REML=TRUE,ddf='Satterthwaite',
          BPPARAM=BiocParallel::SerialParam(),hideErrorsInBackend=FALSE)
        need(!length(numeric_messages(w)),paste(numeric_messages(w),collapse=' | '))
        vc<-attr(f,'varComp')
        need(!is.null(vc)&&'cage'%in%colnames(vc),'dream cage variance missing')
        tau<-as.numeric(vc[1,'cage']);sigma2<-as.numeric(f$sigma[1]^2)
        need(is.finite(tau)&&tau>=0&&is.finite(sigma2)&&sigma2>0,'Invalid variance components')
        row$rho_estimated<-tau/(tau+sigma2)
        row$boundary<-tau/(tau+sigma2)<1e-8
      }
      z<-extract(f,model%in%c('dream','mundlak_dream'))
      row[1,names(z)]<-as.list(z);row$status<-'available'
    },warning=function(e){w<<-c(w,conditionMessage(e));invokeRestart('muffleWarning')}),
      error=function(e){row$reason<<-conditionMessage(e)})
    row$warning<-paste(unique(w),collapse=' | ')
    current[[model]]<-row
  }
  if(current$cage_fixed$status=='available'&&current$mundlak_dream$status=='available') {
    difference<-current$mundlak_dream$beta-current$cage_fixed$beta
    current$mundlak_dream$mundlak_fixed_beta_difference<-difference
    need(abs(difference)<1e-6,'Mundlak/fixed coefficient identity failed')
  }
  records<-c(records,current)
  if((rep-job$start+1L)%%10L==0L||rep==job$end-1L) {
    atomic_json(list(job_id=job$id,phase='cage',done=rep-job$start+1L,total=job$end-job$start,
      pid=Sys.getpid(),updated=as.numeric(Sys.time())),job$heartbeat)
    cat(spec$scenario,rep+1L,'/',job$end,'\n');flush.console()
  }
}
all<-do.call(rbind,records)
need(nrow(all)==(job$end-job$start)*5L&&!anyDuplicated(all[,c('replicate','model')]),'Incomplete R batch')
con<-gzfile(paste0(job$result,'.tmp'),'wt');write.csv(all,con,row.names=FALSE,na='');close(con)
if(file.exists(job$result))need(file.remove(job$result),'Cannot replace unreceipted batch file')
need(file.rename(paste0(job$result,'.tmp'),job$result),'Cannot publish R batch')
