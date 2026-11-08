args <- commandArgs(TRUE)
root <- normalizePath(if(length(args))args[1] else '.',winslash='/')
extra <- strsplit(Sys.getenv('PAPER1_CALIBRATION_LIBRARIES'),.Platform$path.sep,fixed=TRUE)[[1]]
.libPaths(c(extra[nzchar(extra)],.libPaths()))
suppressPackageStartupMessages(library(edgeR));suppressPackageStartupMessages(library(limma))
source(file.path(root,'analysis/common.R'),encoding='UTF-8')
source(file.path(root,'analysis/matched/functions.R'),encoding='UTF-8')
source(file.path(root,'paired_count_calibration/engine.R'))
s <- list(n=40,rho=.7,depth='equal',exposure_fraction=.2,design_seed=40040,index=1,id='test')
d <- make_design(s)
stopifnot(abs(cor(d$x,d$x1[,3])-.7)<1e-10,ncol(d$x1)==5)
set.seed(1);w0<-runif(40,.3,3);w1<-runif(40,.3,3)
# A stratum-constant mean belongs to X0: its paired mean contrast is zero for ANY weights.
mu <- 3+.8*d$x
for(k in 1:20) {
  a <- paired_weights(mu,runif(40,.3,3),runif(40,.3,3),d$x0,d$x1)
  stopifnot(abs(a$delta)<1e-10)
}
y<-rnorm(40)+mu
a<-paired_weights(y,w0,w1,d$x0,d$x1)
b0<-lm.wfit(d$x0,y,w0)$coefficients[2];b1<-lm.wfit(d$x1,y,w1)$coefficients[2]
stopifnot(abs(a$delta-(b1-b0))<1e-10)
for(w in list(rep(1,40),w0)) {
  a<-paired_weights(y,w,w,d$x0,d$x1)
  m0<-solve(crossprod(d$x0,w*d$x0),t(d$x0)*rep(w,each=2))
  m1<-solve(crossprod(d$x1,w*d$x1),t(d$x1)*rep(w,each=5))
  e<-y-d$x1%*%(m1%*%y);dd<-m1[2,]-m0[2,]
  manual<-sqrt(sum(w*e^2)/(40-5)*sum(dd^2/w))
  stopifnot(abs(manual-a$SE_delta)<1e-10)
}
p<-list(settings=list(genes=250,seed=20261107),parameter_file=NULL)
par<-generator_parameters(p)
path<-tempfile(fileext='.csv')
write.csv(data.frame(mean_count=seq(1,500),dispersion=rep(.1,500)),path,row.names=FALSE)
p2<-p;p2$parameter_file<-path
empirical<-generator_parameters(p2)
stopifnot(nrow(empirical)==250,empirical$mean_count[1]==1,tail(empirical$mean_count,1)==500)
write.csv(data.frame(mean_count=seq(1,500),dispersion=rep(0,500)),path,row.names=FALSE)
stopifnot(inherits(try(generator_parameters(p2),silent=TRUE),'try-error'))
unlink(path)
c1<-simulate_counts(par,d,s,99);c2<-simulate_counts(par,d,s,99)
stopifnot(identical(c1,c2),all(c1$y>=0),all(c1$y==round(c1$y)),sum(c1$beta!=0)==50)
s2<-s;s2$depth<-'iid_lognormal'
c3<-simulate_counts(par,d,s2,99)
stopifnot(!identical(c3$depth,c1$depth),!identical(c3$y,c1$y))
# Filtering/TMM/logCPM commute with a within-exposure column permutation.
permutation<-c(sample(1:20),sample(21:40))
transform<-function(count) {
  g<-DGEList(count);g<-g[filterByExpr(g,design=d$x0),,keep.lib.sizes=FALSE]
  g<-calcNormFactors(g);voom(g,d$x0,adaptive.span=TRUE,plot=FALSE)$E
}
e1<-transform(c1$y);e2<-transform(c1$y[,permutation])
stopifnot(identical(rownames(e1),rownames(e2)),max(abs(e1[,permutation]-e2))<1e-10)
r<-one_dataset(p,s,1,par)
stopifnot(length(r$metrics)==4,all(vapply(r$checks,function(z)z$lmfit_max_difference<1e-7,TRUE)))
rr<-score_test(c(.01,NA),c(-1,NA),c(1,NA),c(0,NA),c(1,NA),10,2,'shared','HC3')
stopifnot(rr$available==1,rr$unavailable_retained==1,rr$filtered==8,
          rr$discoveries_planned_050==0,rr$discoveries_retained_050==1,rr$coverage95==1)
cat('PASS: design, null maps, two-model coefficients, paired variance, seeds, count generation,\n')
cat('      filtering/TMM permutation equivariance, full pipeline, BH family and availability\n')
