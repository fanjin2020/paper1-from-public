extra <- Sys.getenv("R_LIBS", "")
if(nzchar(extra)) .libPaths(c(strsplit(extra,.Platform$path.sep,fixed=TRUE)[[1]],.libPaths()))

suppressPackageStartupMessages(library(limma))
options(stringsAsFactors=FALSE,digits=17)
a <- commandArgs(TRUE); inp <- a[1]; out <- a[2]
wr <- function(x,n) write.csv(x,file.path(out,n),row.names=FALSE,na='')
m <- read.csv(file.path(inp,'metadata.csv'),check.names=FALSE)
g <- read.csv(file.path(inp,'models.csv'),check.names=FALSE)
f <- read.csv(file.path(inp,'fixed.csv'),check.names=FALSE)
e <- read.csv(gzfile(file.path(inp,'expression.csv.gz')),check.names=FALSE)
stopifnot(all(c('analysis_unit_id','parsed_individual_id','parsed_brain_region',
 'parsed_group','parsed_age','parsed_sex') %in% names(m)),
 all(c('region','model','gene_symbol','beta_input_scale','BH_q') %in% names(g)),
 all(c('discovery_region','gene_symbol') %in% names(f)))
y <- as.matrix(e[,-1]); storage.mode(y) <- 'double'; rownames(y) <- e[[1]]
stopifnot(!anyDuplicated(rownames(y)),!anyDuplicated(colnames(y)),
 identical(colnames(y),m$analysis_unit_id),all(is.finite(y)),
 !anyDuplicated(g[c('region','model','gene_symbol')]),
 !anyDuplicated(f[c('discovery_region','gene_symbol')]))
design <- function(d,quadratic=FALSE) {
 age <- d$parsed_age-mean(d$parsed_age)
 x <- cbind(intercept=1,AD=as.numeric(d$parsed_group=='AD'),age=age,
            male=as.numeric(d$parsed_sex=='male'))
 if(quadratic) x <- cbind(x,age_squared=age^2)
 x
}
ready <- function(d,x) nrow(d)>0 && sum(d$parsed_group=='AD')>=10 &&
 sum(d$parsed_group=='control')>=10 && qr(x)$rank==ncol(x) && nrow(x)-ncol(x)>=10
fit <- function(d,x) {
 z <- y[,d$analysis_unit_id,drop=FALSE]
 obj <- eBayes(lmFit(z,x),proportion=.01,stdev.coef.lim=c(.1,4),trend=FALSE,robust=FALSE)
 tol <- 64*.Machine$double.eps*max(ncol(z),ncol(x))*pmax(1,apply(abs(z),1,max))
 stopifnot(all(is.finite(obj$sigma)),all(obj$sigma>tol))
 b <- obj$coefficients[,'AD']; se <- obj$stdev.unscaled[,'AD']*sqrt(obj$s2.post)
 p <- obj$p.value[,'AD']; crit <- qt(.975,df=obj$df.total)
 stopifnot(all(is.finite(b)),all(is.finite(se)),all(se>0),all(is.finite(p)))
 data.frame(gene_symbol=rownames(z),beta=b,moderated_SE=se,
            CI_low=b-crit*se,CI_high=b+crit*se,p=p,q=p.adjust(p,'BH'))
}
manifest <- list(); gate <- list(); overlap <- list(); fits <- list(); reproduction <- list()
for(r in c('entorhinal_cortex','frontal_cortex')) {
 d <- m[m$parsed_brain_region==r & m$parsed_group %in% c('AD','control') &
 is.finite(m$parsed_age) & m$parsed_sex %in% c('female','male'),,drop=FALSE]
 stopifnot(!anyNA(d$parsed_individual_id),all(nzchar(d$parsed_individual_id)),
           !anyDuplicated(d$parsed_individual_id))
 d <- d[order(d$parsed_individual_id),,drop=FALSE]
 x <- design(d); stopifnot(ready(d,x))
 b <- fit(d,x)
 ref <- g[g$region==r & g$model=='linear_age',,drop=FALSE]
 stopifnot(setequal(ref$gene_symbol,b$gene_symbol))
 ref <- ref[match(b$gene_symbol,ref$gene_symbol),]
 bd <- max(abs(b$beta-ref$beta_input_scale)); qd <- max(abs(b$q-ref$BH_q))
 stopifnot(bd<1e-8,qd<1e-8)
 stopifnot(setequal(f$gene_symbol[f$discovery_region==r],b$gene_symbol[b$q<=.05]))
 reproduction[[r]] <- data.frame(region=r,max_beta_difference=bd,max_q_difference=qd)
 agesAD <- d$parsed_age[d$parsed_group=='AD']; agesC <- d$parsed_age[d$parsed_group=='control']
 low <- max(min(agesAD),min(agesC)); high <- min(max(agesAD),max(agesC))
 keep <- d$parsed_age>=low & d$parsed_age<=high
 manifest[[r]] <- data.frame(region=r,d,overlap_lower=low,overlap_upper=high,
 retained_overlap=keep,reason=ifelse(keep,'inside_inclusive_range','outside_inclusive_range'))
 overlap[[r]] <- data.frame(region=r,lower=low,upper=high,
 original_AD=length(agesAD),original_control=length(agesC),
 retained_AD=sum(keep & d$parsed_group=='AD'),retained_control=sum(keep & d$parsed_group=='control'))
 b$region <- r; b$model <- 'full_linear'; fits[[paste(r,'full')]] <- b
 for(model in c('overlap_linear','overlap_quadratic')) {
   sub <- d[keep,,drop=FALSE]; dx <- design(sub,model=='overlap_quadratic')
   ok <- ready(sub,dx)
   gate[[paste(r,model)]] <- data.frame(region=r,model=model,n=nrow(sub),
    rank=if(nrow(sub)) qr(dx)$rank else 0,columns=ncol(dx),
    residual_df=nrow(sub)-ncol(dx),status=ifelse(ok,'tested','not_testable_design'))
   if(!ok) next
   wr(data.frame(analysis_unit_id=sub$analysis_unit_id,dx),paste0('182_',r,'_',model,'_design.csv'))
   z <- fit(sub,dx); z$region <- r; z$model <- model; fits[[paste(r,model)]] <- z
 }
}
wr(do.call(rbind,manifest),'182_all_sample_membership.csv')
wr(do.call(rbind,overlap),'182_age_intersections.csv')
wr(do.call(rbind,gate),'182_design_status.csv')
wr(do.call(rbind,reproduction),'183_baseline_reproduction.csv')
all <- do.call(rbind,fits); rownames(all)<-NULL
all$joint_q <- NA_real_
for(md in unique(all$model)) {
 ix <- which(all$model==md); all$joint_q[ix] <- p.adjust(all$p[ix],'BH')
}
wr(all,'183_all_gene_models.csv')
comparisons <- list(); summary <- list()
for(r in unique(all$region)) for(md in c('overlap_linear','overlap_quadratic')) {
 b <- all[all$region==r & all$model=='full_linear',]
 z <- all[all$region==r & all$model==md,]
 if(!nrow(z)) next
 z <- z[match(b$gene_symbol,z$gene_symbol),]
 selected <- b$gene_symbol %in% f$gene_symbol[f$discovery_region==r]
 c <- data.frame(region=r,model=md,gene_symbol=b$gene_symbol,fixed_endpoint=selected,
 baseline_beta=b$beta,baseline_CI_low=b$CI_low,baseline_CI_high=b$CI_high,
 overlap_beta=z$beta,overlap_CI_low=z$CI_low,overlap_CI_high=z$CI_high,
 beta_change=z$beta-b$beta,overlap_q_region=z$q,overlap_q_joint=z$joint_q,
 direction=ifelse(b$beta*z$beta>0,'same',ifelse(b$beta*z$beta<0,'opposite','zero')))
 comparisons[[paste(r,md)]] <- c
 summary[[paste(r,md)]] <- data.frame(region=r,model=md,fixed_n=sum(selected),
 same=sum(c$direction[selected]=='same'),opposite=sum(c$direction[selected]=='opposite'),
 retained_region_q=sum(z$q[selected]<=.05),retained_joint_q=sum(z$joint_q[selected]<=.05),
 median_absolute_beta_change=median(abs(c$beta_change[selected])))
}
if(length(comparisons)) wr(do.call(rbind,comparisons),'183_all_gene_comparisons.csv')
if(length(summary)) wr(do.call(rbind,summary),'183_fixed_endpoint_summary.csv')
writeLines(capture.output(sessionInfo()),file.path(out,'183_R_session.txt'))
writeLines(c('Post hoc sensitivity, not independent replication.',
 'Intervals are model-based marginal 95% moderated-t intervals, not simultaneous intervals.',
 'Subset and full-sample estimates are dependent; no independence-based difference test.',
 'Range intersection does not establish adequate conditional age overlap or remove confounding.',
 'Native expression scale retained; do not label beta log2FC without source-scale evidence.',
 'Region BH and joint-across-region BH are distinct; joint_q is the sensitivity reporting family.',
 'Models are separate sensitivity families; do not select the model with smaller q.',
 'No genes selected anew for external endpoints.'),file.path(out,'183_interpretation.txt'))
cat('182-183 complete: age range sensitivity, not age-confounding resolution.\n')
