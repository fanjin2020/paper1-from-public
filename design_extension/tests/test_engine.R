args<-commandArgs(FALSE)
this<-normalizePath(sub('^--file=','',args[grep('^--file=',args)][1]))
root<-dirname(dirname(this))
for(name in c('common.R','diagnostics.R','functions.R','mouse.R','counts.R'))
  source(file.path(root,'engine',name))
invisible(parse(file.path(root,'count_bridge.R')))
set.seed(42)
d<-data.frame(animal_id=paste0('a',1:48),cage=factor(rep(1:12,each=4)),
              genotype=rep(c(0,1),24),age=runif(48,3,12),sex=factor(rep(c('F','F','M','M'),12)),
              slide=factor(rep(1:4,each=12)))
ds<-mouse_designs(d,'linear')
stopifnot(ds$base$estimable,ds$fixed$estimable,ds$mundlak$estimable)
d$genotype<-rep(rep(c(0,1),6),each=4)
ds<-mouse_designs(d,'linear')
stopifnot(ds$base$estimable,!ds$fixed$estimable,!ds$mundlak$estimable)
cal<-list(metadata=d,parameters=data.frame(gene_symbol=paste0('g',1:200),
          mean_count=rep(30,200),dispersion=rep(.2,200)),depth_ratio=rep(1,48))
s<-data.frame(scenario=0,cage_log2_sd=.3,cage_mean_amplitude=.5,de_fraction=.1,
              effect_log2=.5,mean_affected_fraction=.2)
a<-generate_counts(cal,s,3,321);b<-generate_counts(cal,s,3,321)
stopifnot(identical(a,b),all(a$counts>=0),all(a$counts==round(a$counts)),
          sum(a$beta!=0)==20,sum(a$gamma!=0)==40)
tab<-unavailable(cal$parameters$gene_symbol,'baseline','test','fit_failed')
r<-performance_row(tab,a$beta,a$gamma,s,3,'baseline',rep(TRUE,200))
stopifnot(is.na(r$FDP),r$availability==0,r$delivered_coverage==0)
tab$status<-'available';tab$beta<-a$beta;tab$p<-1;tab$low<-a$beta-1;tab$high<-a$beta+1
r<-performance_row(tab,a$beta,a$gamma,s,3,'baseline',rep(TRUE,200))
stopifnot(r$FDP==0,r$coverage==1,r$availability==1,r$power==0,r$bias==0)
stopifnot(isTRUE(all.equal(bh(c(.01,.04,NA))[1:2],c(.03,.06))))
cat('PASS: R syntax, rank/estimability, NB generator, deterministic seeds, planned BH, failure accounting\n')

if(requireNamespace('jsonlite',quietly=TRUE)) {
  local({
    folder<-tempfile('p1_diagnostic_test_');dir.create(folder)
    folder<-normalizePath(folder,winslash='/',mustWork=TRUE)
    temporary_root<-normalizePath(tempdir(),winslash='/',mustWork=TRUE)
    stopifnot(startsWith(folder,paste0(temporary_root,'/')),
              startsWith(basename(folder),'p1_diagnostic_test_'))
    on.exit(unlink(folder,recursive=TRUE),add=TRUE)
    job<-list(design='k18_w1',scenario=list(scenario=3L),start=7L,end=8L,output=folder)
    z<-data.frame(gene_symbol='g1',model='dream_native',status='numerical_warning')
    tagged<-tag_count_rows(z,job,7L)
    stopifnot(tagged$design=='k18_w1',tagged$scenario==3L,tagged$replicate==7L,
              identical(names(z),c('gene_symbol','model','status')),
              is.null(tag_count_rows(data.frame(),job,7L)))
    wr(tagged,file.path(folder,'diagnostics.csv'))
    wr(data.frame(),file.path(folder,'retry_attempts.csv'))
    stopifnot(rd(file.path(folder,'diagnostics.csv'))$gene_symbol=='g1',
              file.exists(file.path(folder,'retry_attempts.csv')))
    Sys.setenv(PAPER1_ATTEMPT_ID='attempt',PAPER1_RUN_ID='run',PAPER1_JOB_HASH='jobhash')
    save_count_progress(job,1L)
    progress<-jsonlite::fromJSON(file.path(folder,'progress.json'))
    stopifnot(progress$done==1L,progress$total==1L,progress$attempt_id=='attempt',
              progress$run_id=='run',progress$job_hash=='jobhash')
  })
  cat('PASS: diagnostic identity, empty-table persistence and attempt-tagged progress\n')
} else cat('SKIP: progress serialization test requires the existing jsonlite package\n')
