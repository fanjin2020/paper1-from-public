extra <- Sys.getenv("R_LIBS", "")
if(nzchar(extra)) .libPaths(c(strsplit(extra,.Platform$path.sep,fixed=TRUE)[[1]],.libPaths()))

suppressPackageStartupMessages(library(limma))
options(stringsAsFactors=FALSE,digits=17)
out <- commandArgs(TRUE)[1]
rd <- function(n) read.csv(file.path(out,n),check.names=FALSE,na.strings=c('NA','NaN','null',''))
wr <- function(x,n) write.csv(x,file.path(out,n),row.names=FALSE,na='')
m <- rd('187_selected_samples.csv'); ann <- rd('187_probe_annotation.csv')
ep <- rd('186_fixed_endpoints.csv'); src <- rd('187_selected_probe_matrix.csv.gz')
stopifnot(!anyDuplicated(src$probe_id),!anyDuplicated(ann$probe_id),
 identical(names(src)[-1],m$sample_id),!anyDuplicated(m[c('region','recorded_unit')]))
y <- as.matrix(src[,-1]); storage.mode(y) <- 'double'; rownames(y) <- src$probe_id
aa <- ann[match(rownames(y),ann$probe_id),]
stopifnot(!anyNA(aa$probe_id))
models <- c('overlap_log2_linear','full_age_log2_linear','overlap_log2_quadratic','overlap_native_linear')
regions <- c('entorhinal_cortex','frontal_cortex')
gate <- list(); membership <- list(); cover <- list(); genes <- list(); endpts <- list(); signatures <- list()
profiles <- list(); bins <- list(); ranges <- list()
for(r in regions) {
 d <- m[m$region==r,,drop=FALSE]; ix <- match(d$sample_id,m$sample_id)
 stopifnot(nrow(d)>0,all(c('AD','control') %in% d$group))
 raw <- y[,ix,drop=FALSE]
 positive <- rowSums(!is.finite(raw) | raw<=0)==0
 mapped <- aa$status=='unique_symbol' & !is.na(aa$gene_symbol)
 ok <- positive & mapped
 bg <- data.frame(region=r,probe_id=rownames(raw),gene_symbol=aa$gene_symbol,
   included=ok,reason=ifelse(!mapped,'ambiguous_or_unmapped',ifelse(!positive,'nonpositive_or_nonfinite','included')))
 cover[[r]] <- bg
 stopifnot(sum(ok)>1000)
 split_probe <- split(which(ok),aa$gene_symbol[ok])
 # Probe summaries are fixed before group testing, with no max-effect/probe selection.
 collapse <- function(mat) t(vapply(split_probe,function(i) {
   if(length(i)==1) as.numeric(mat[i,]) else apply(mat[i,,drop=FALSE],2,median)
 },numeric(ncol(mat))))
 logmat <- collapse(log2(raw)); nativemat <- collapse(raw)
 colnames(logmat) <- colnames(nativemat) <- d$sample_id
 write.csv(data.frame(gene_symbol=rownames(logmat),logmat,check.names=FALSE),
   gzfile(file.path(out,paste0('187_',r,'_gene_log2.csv.gz'))),row.names=FALSE)
 low <- max(tapply(d$age,d$group,min)); high <- min(tapply(d$age,d$group,max))
 keep <- d$age>=low & d$age<=high
 membership[[r]] <- data.frame(d,overlap_lower=low,overlap_upper=high,in_age_intersection=keep)
 ranges[[r]] <- data.frame(region=r,lower=low,upper=high,full_AD=sum(d$group=='AD'),
  full_control=sum(d$group=='control'),overlap_AD=sum(keep & d$group=='AD'),overlap_control=sum(keep & d$group=='control'))
 for(gr in c('AD','control')) for(sex in c('all','female','male')) {
  ss <- d[d$group==gr & (sex=='all' | d$sex==sex),,drop=FALSE]
  profiles[[paste(r,gr,sex)]] <- data.frame(region=r,group=gr,sex=sex,n=nrow(ss),
   age_min=if(nrow(ss)) min(ss$age) else NA,age_max=if(nrow(ss)) max(ss$age) else NA,
   age_mean=if(nrow(ss)) mean(ss$age) else NA,age_sd=if(nrow(ss)>1) sd(ss$age) else NA)
 }
 bins[[r]] <- as.data.frame(table(region=rep(r,nrow(d)),group=d$group,sex=d$sex,
   age_bin=cut(d$age,breaks=seq(0,125,5),right=FALSE)))
 fixed <- ep[ep$discovery_region==r,,drop=FALSE]
 for(md in models) {
  take <- if(md=='full_age_log2_linear') rep(TRUE,nrow(d)) else keep
  sub <- d[take,,drop=FALSE]; a <- sub$age-mean(sub$age)
  x <- cbind(intercept=1,AD=as.numeric(sub$group=='AD'),age=a,male=as.numeric(sub$sex=='male'))
  if(md=='overlap_log2_quadratic') x <- cbind(x,age_squared=a^2)
  rank <- if(nrow(x)) qr(x)$rank else 0
  ready <- sum(sub$group=='AD')>=5 && sum(sub$group=='control')>=5 &&
    rank==ncol(x) && nrow(x)-ncol(x)>=5
  id <- paste(r,md)
  gate[[id]] <- data.frame(region=r,model=md,AD=sum(sub$group=='AD'),control=sum(sub$group=='control'),
    rank=rank,n_columns=ncol(x),residual_df=nrow(x)-ncol(x),
    status=ifelse(ready,'tested','not_testable_design'))
  wr(data.frame(sample_id=sub$sample_id,x),paste0('187_',r,'_',md,'_design.csv'))
  z <- (if(md=='overlap_native_linear') nativemat else logmat)[,take,drop=FALSE]
  ans <- data.frame(region=r,model=md,gene_symbol=rownames(z),beta=NA_real_,SE=NA_real_,
    CI_low=NA_real_,CI_high=NA_real_,p=NA_real_,status='not_testable_design',
    scale=ifelse(md=='overlap_native_linear','deposited_normalized_units','log2_of_deposited_normalized_values'))
  numer_ok <- rep(FALSE,nrow(z))
  if(ready) {
    fit <- lmFit(z,x)
    tolerance <- 64*.Machine$double.eps*max(nrow(x),ncol(x))*pmax(1,apply(abs(z),1,max))
    numer_ok <- is.finite(fit$sigma) & fit$sigma>tolerance & fit$df.residual>=5
    ans$status <- 'unassessable_residual_variance'
    if(sum(numer_ok)>=2) {
      fit <- eBayes(fit[numer_ok,],proportion=.01,stdev.coef.lim=c(.1,4),trend=FALSE,robust=FALSE)
      se <- fit$stdev.unscaled[,'AD']*sqrt(fit$s2.post); b <- fit$coefficients[,'AD']
      ci <- qt(.975,fit$df.total)*se
      stopifnot(all(is.finite(b)),all(is.finite(se)),all(se>0),all(is.finite(fit$p.value[,'AD'])))
      ans[numer_ok,c('beta','SE','CI_low','CI_high','p')] <- cbind(b,se,b-ci,b+ci,fit$p.value[,'AD'])
      ans$status[numer_ok] <- 'tested'
    } else numer_ok[] <- FALSE
  }
  genes[[id]] <- ans
  match_gene <- match(fixed$gene_symbol,ans$gene_symbol)
  ff <- ans[match_gene,,drop=FALSE]; rownames(ff) <- NULL
  ff$region <- r; ff$model <- md; ff$gene_symbol <- fixed$gene_symbol
  absent <- is.na(match_gene)
  ff$status[absent] <- 'unmapped_or_no_complete_positive_probe'
  ff$discovery_beta_native <- fixed$discovery_log2_fold_change
  ff$discovery_direction <- fixed$discovery_direction
  ff$direction_relation <- ifelse(is.na(ff$beta),'not_assessable',
   ifelse(ff$beta==0,'zero',ifelse(sign(ff$beta)==ifelse(fixed$discovery_direction=='Up',1,-1),'same','opposite')))
  endpts[[id]] <- ff
  for(direction in c('Up','Down')) {
    set <- fixed$gene_symbol[fixed$discovery_direction==direction]
    eligible <- which(numer_ok & rownames(z) %in% set)
    nbg <- sum(numer_ok); nset <- length(eligible)
    sr <- data.frame(region=r,model=md,signature=direction,fixed_size=length(set),
      measured_size=nset,background_size=nbg,Direction=NA_character_,Correlation=NA_real_,p=NA_real_,
      status=ifelse(ready,'insufficient_members_or_background','not_testable_design'))
    if(ready && nset>=15 && nbg-nset>=15) {
      zz <- z[numer_ok,,drop=FALSE]; index <- which(rownames(zz) %in% set)
      cc <- camera(zz,index=list(fixed_signature=index),design=x,contrast=2,
                   use.ranks=FALSE,inter.gene.cor=NA,allow.neg.cor=FALSE)
      stopifnot(nrow(cc)==1,is.finite(cc$PValue),cc$PValue>=0,cc$PValue<=1)
      sr$Direction <- cc$Direction; sr$Correlation <- cc$Correlation
      sr$p <- cc$PValue; sr$status <- 'tested'
    }
    signatures[[paste(id,direction)]] <- sr
  }
 }
}
wr(do.call(rbind,cover),'187_all_probe_background_decisions.csv')
wr(do.call(rbind,membership),'187_analysis_sample_membership.csv')
wr(do.call(rbind,ranges),'187_age_intersections.csv')
wr(do.call(rbind,profiles),'187_age_sex_summary.csv')
wr(do.call(rbind,bins),'187_age_sex_bins.csv')
wr(do.call(rbind,gate),'188_design_status.csv')
all <- do.call(rbind,genes); endpoints <- do.call(rbind,endpts); sig <- do.call(rbind,signatures)
adjust <- function(t) {
 t$joint_q <- NA_real_; t$family_size <- 0L
 for(md in unique(t$model)) {
  k <- which(t$model==md); pp <- t$p[k]; valid <- is.finite(pp); pp[!valid]<-1
  q <- p.adjust(pp,'BH'); q[!valid] <- NA_real_
  t$joint_q[k] <- q; t$family_size[k] <- length(k)
 }
 t
}
all <- adjust(all); endpoints <- adjust(endpoints); sig <- adjust(sig)
endpoints$all_background_joint_q <- all$joint_q[match(paste(endpoints$region,endpoints$model,endpoints$gene_symbol),
 paste(all$region,all$model,all$gene_symbol))]
sig$direction_relation <- ifelse(is.na(sig$Direction),'not_assessable',
 ifelse(sig$Direction==sig$signature,'same','opposite'))
write.csv(all,gzfile(file.path(out,'188_all_gene_models.csv.gz')),row.names=FALSE,na='')
wr(endpoints,'188_all_fixed_endpoint_results.csv')
wr(sig,'189_all_fixed_signature_camera.csv')
summary <- list()
for(r in regions) for(md in models) {
 a <- endpoints[endpoints$region==r & endpoints$model==md,]
 summary[[paste(r,md)]] <- data.frame(region=r,model=md,fixed_n=nrow(a),
  tested=sum(a$status=='tested'),unassessable=sum(a$status!='tested'),same=sum(a$direction_relation=='same'),
  opposite=sum(a$direction_relation=='opposite'),
  significant_same=sum(a$joint_q<=.05 & a$direction_relation=='same',na.rm=TRUE),
  significant_opposite=sum(a$joint_q<=.05 & a$direction_relation=='opposite',na.rm=TRUE),
  min_q=if(any(is.finite(a$joint_q))) min(a$joint_q,na.rm=TRUE) else NA_real_)
}
wr(do.call(rbind,summary),'190_endpoint_summary.csv')
primary <- endpoints[endpoints$model=='overlap_log2_linear',]
paired <- list()
for(md in setdiff(models,'overlap_log2_linear')) {
 s <- endpoints[endpoints$model==md,]; s <- s[match(paste(primary$region,primary$gene_symbol),paste(s$region,s$gene_symbol)),]
 paired[[md]] <- data.frame(region=primary$region,gene_symbol=primary$gene_symbol,sensitivity_model=md,
 primary_beta=primary$beta,sensitivity_beta=s$beta,primary_q=primary$joint_q,sensitivity_q=s$joint_q,
 direction_agrees=sign(primary$beta)==sign(s$beta),
 beta_difference=if(md=='overlap_native_linear') NA_real_ else s$beta-primary$beta,
 note=if(md=='overlap_native_linear') 'different_scales_no_magnitude_subtraction' else 'dependent_estimates_descriptive_difference_only')
}
wr(do.call(rbind,paired),'190_all_endpoint_sensitivity.csv')
writeLines(capture.output(sessionInfo()),file.path(out,'190_R_session.txt'))
writeLines(c('# GSE48350 external association analysis',
 '', 'Primary: age-intersection, log2 of deposited normalized data, age/sex adjusted.',
 'Full-age, quadratic-age and native-scale models are sensitivities, not alternative primary winners.',
 '188_all_fixed_endpoint_results.csv retains all original endpoints, including unmapped and unassessable.',
 'joint_q on the endpoint table uses the fixed endpoint family; all_background_joint_q is a different family.',
 'BH families span both regions within each model, with p=1 placeholders for untestable planned endpoints.',
 '189_all_fixed_signature_camera.csv uses four planned signature-region tests per model.',
 'CAMERA correlation is residual inter-gene correlation, not a disease effect size.',
 'Model intervals are marginal, not multiplicity-adjusted; discovery and external scales are not pooled.',
 'Direction, q values and effect estimates must be reported separately. Non-significance is inconclusive.',
 'No disease-by-region interaction was fitted; regional significance differences do not prove regional specificity.',
 'This analysis does not certify donor independence, production-stage equivalence or causal mechanisms.',
 'Source-unit labels are metadata labels; clinical/donor identity is not independently authenticated.',
 'Control reuse from GSE11882 and repeated regional source units must be disclosed.',
 'Do not publish private/local source-unit keys without reviewing provenance and access conditions.',
 '', '## Fixed-endpoint results', capture.output(print(do.call(rbind,summary),row.names=FALSE)),
 '', '## All signature tests', capture.output(print(sig,row.names=FALSE))),file.path(out,'190_results_handoff.md'))
cat('Models, endpoint coverage, signatures and handoff completed.\n')
