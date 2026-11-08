root<-normalizePath(commandArgs(TRUE)[1],winslash='/')
extra<-strsplit(Sys.getenv('PAPER1_CALIBRATION_LIBRARIES'),.Platform$path.sep,fixed=TRUE)[[1]]
.libPaths(c(extra[nzchar(extra)],.libPaths()))
suppressPackageStartupMessages(library(edgeR));suppressPackageStartupMessages(library(limma))
source(file.path(root,'analysis/common.R'),encoding='UTF-8')
source(file.path(root,'analysis/matched/functions.R'),encoding='UTF-8')
source(file.path(root,'paired_count_calibration/engine.R'))
source(file.path(root,'paired_tail_diagnostics/engine.R'))
p<-list(settings=list(genes=250,seed=120261108),parameter_file=NULL)
s<-list(n=40,rho=.3,depth='iid_lognormal',exposure_fraction=.2,design_seed=40040,index=1,id='s01')
par<-generator_parameters(p);a<-make_tail_inputs(p,s,1,par)
stopifnot(length(a$arms)==6)
for(arm in names(a$arms)[1:5])stopifnot(identical(a$arms[[arm]]$E,a$arms$estimated_shared$E))
stopifnot(!identical(a$arms$gaussian_unit$E,a$arms$estimated_shared$E),
          !identical(a$arms$independent_shared$w0,a$arms$estimated_shared$w0),
          identical(a$arms$independent_shared$w0,a$arms$independent_shared$w1),
          all(a$arms$fixed_unit$w0==1),all(a$arms$gaussian_unit$w1==1))
r<-tail_dataset(p,s,1,par);rr<-tail_dataset(p,s,1,par)
stopifnot(identical(r,rr),length(r$result$metrics)==12)
for(m in r$result$metrics) {
  stopifnot(m$lmfit_max_difference<1e-7,length(m$tails)==8,
            sum(unlist(m$p_histogram))==m$available,
            abs(m$p_histogram[[1]]-m$fpr_050*m$available)<1e-8,
            sum(vapply(m$strata,function(x)x$available,0))==m$available)
  ev<-r$events[r$events$arm==m$weight_mode&r$events$test==m$test,,drop=FALSE]
  stopifnot(nrow(ev)==m$event_candidates-m$events_truncated,
            all(abs(ev$delta-(ev$beta_full-ev$beta_base))<1e-7))
  if(m$events_truncated==0)stopifnot(sum(ev$bh_planned)==m$discoveries_planned_050,
                                   sum(ev$bh_retained)==m$discoveries_retained_050)
}
x<-.05+1e-9
stopifnot(hist(x,breaks=seq(0,1,.05),fuzz=0,plot=FALSE)$counts[1]==0)
cat('PASS: response alignment, independent/control weights, seeded equality, 12 arms/tests,\n')
cat('      lmFit identity, strict tails, strata, event selection, BH counts and histogram boundary\n')
