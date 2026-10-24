# Numerical identities and omitted-row recovery for the matched analysis.
run_matched_regression <- function(root) {
  source(file.path(root,'analysis','matched','functions.R'))
  source(file.path(root,'analysis','matched','summarise.R'))
  z<-mean_mc(rep(1,200),TRUE)
  stopifnot(z['low']>.98,z['low']<.99,z['high']==1,z['mean']==1)
  z<-mean_mc(rep(0,200),TRUE)
  stopifnot(z['low']==0,z['high']>.01,z['high']<.02)
  stopifnot(is.na(mean_mc(numeric(),TRUE)['low']))
  x<-c(.1,.2,.3,.4)
  stopifnot(mean_mc(x)['mean']==mean(x),mean_mc(x)['MCSE']==sd(x)/2)
  set.seed(19)
  d<-data.frame(genotype=rep(c(0,1,0,1),9),age=rnorm(36),cage=factor(rep(1:12,each=3)))
  rownames(d)<-paste0('s',seq_len(nrow(d)))
  y<-matrix(rnorm(50*36),50,36,dimnames=list(paste0('g',1:50),rownames(d)))
  v<-methods::new('EList',list(E=y,weights=matrix(1,50,36)))
  form<-~genotype+age+(1|cage)
  f<-variancePartition::dream(v,form,d,ddf='Satterthwaite',BPPARAM=BiocParallel::SerialParam())
  gene<-'g3'
  one<-variancePartition::dream(v[gene,,drop=FALSE],form,d,ddf='Satterthwaite',BPPARAM=BiocParallel::SerialParam())
  replacement<-replace_dream_raw_row(f,one,gene)
  without<-f;ids<-rownames(f$coefficients);keep<-ids!=gene
  matrices<-c('coefficients','stdev.unscaled','df.residual','df.total','t','p.value','z.std','residuals','hatvalues','genes')
  vectors<-c('rdf','Amean','sigma','logLik','edf','F','F.p.value')
  for(k in intersect(names(without),matrices))if(!is.null(without[[k]]))without[[k]]<-without[[k]][keep,,drop=FALSE]
  for(k in intersect(names(without),vectors))if(!is.null(without[[k]]))without[[k]]<-without[[k]][keep]
  without$cov.coefficients.list<-without$cov.coefficients.list[keep]
  attr(without,'varComp')<-attr(without,'varComp')[keep,,drop=FALSE]
  attr(without,'edf')<-attr(without,'edf')[keep]
  attr(without,'errors')<-setNames('synthetic backend omission',gene)
  restored<-replace_dream_raw_row(without,one,gene)
  stopifnot(nrow(restored$coefficients)==nrow(f$coefficients),!length(attr(restored,'errors')),
    identical(restored$coefficients[rownames(without$coefficients),,drop=FALSE],without$coefficients))
  actual<-dream_rows(restored,ids,'test','test')$table
  expected<-dream_rows(replacement,ids,'test','test')$table
  stopifnot(identical(actual$status,expected$status),max(abs(actual$p-expected$p))<1e-10,
            max(abs(actual$SE-expected$SE))<1e-10)
  cat('PASS: dataset Monte Carlo intervals, omitted-row restoration, joint moderation and unchanged healthy rows\n')
}
