# Each gene has a zero paired mean contrast, including genes with an exposure effect.
make_design <- function(s) {
  set.seed(s$design_seed)
  n <- s$n
  x <- rep(c(0,1),each=n/2)
  x0 <- cbind(Intercept=1,target=x)
  u <- qr.Q(qr(qr.resid(qr(x0),matrix(rnorm(n*3),n,3))))*sqrt(n-1)
  z <- s$rho*as.numeric(scale(x))+sqrt(1-s$rho^2)*u[,1]
  x1 <- cbind(x0,quality1=z,quality2=u[,2],quality3=u[,3])
  need(estimable(x0)&&estimable(x1),'Rank-deficient design')
  need(abs(cor(x,z)-s$rho)<1e-10,'Incorrect design correlation')
  list(x0=x0,x1=x1,x=x)
}

generator_parameters <- function(p) {
  g <- p$settings$genes
  if(is.null(p$parameter_file)) {
    m <- exp(seq(log(2),log(2000),length.out=g))
    phi <- .02+.3/sqrt(m/10)
    return(data.frame(gene_id=sprintf('gene%05d',seq_len(g)),mean_count=m,dispersion=phi))
  }
  a <- read.csv(p$parameter_file,check.names=FALSE)
  need(all(c('mean_count','dispersion')%in%names(a)),'Parameter columns missing')
  need(nrow(a)>=g,'Parameter file has fewer genes than requested')
  need(all(is.finite(a$mean_count)&a$mean_count>0&is.finite(a$dispersion)&a$dispersion>0),
       'Invalid NB means/dispersion; no silent dropping or flooring')
  a <- a[order(a$mean_count,seq_len(nrow(a))),,drop=FALSE]
  a <- a[unique(round(seq(1,nrow(a),length.out=g))),,drop=FALSE]
  need(nrow(a)==g,'Parameter subsampling lost genes')
  data.frame(gene_id=sprintf('gene%05d',seq_len(g)),mean_count=a$mean_count,dispersion=a$dispersion)
}

simulate_counts <- function(par,d,s,seed) {
  set.seed(s$design_seed+19L)
  g <- nrow(par)
  affected <- sample.int(g,2L*floor(g*s$exposure_fraction/2))
  b <- numeric(g)
  b[affected] <- rep(c(-.75,.75),length.out=length(affected))
  set.seed(seed)
  # IID depths are resampled within each dataset, not fixed correlated offsets.
  depth <- if(s$depth=='equal')rep(1,s$n) else exp(rnorm(s$n,sd=.6)-.6^2/2)
  mu <- (par$mean_count%o%depth)*2^(b%o%d$x)
  y <- matrix(rnbinom(length(mu),mu=as.numeric(mu),size=rep(1/par$dispersion,s$n)),g,s$n)
  dimnames(y) <- list(par$gene_id,sprintf('sample%03d',seq_len(s$n)))
  need(all(colSums(y)>0),'Empty count library')
  list(y=y,depth=depth,beta=b)
}

score_test <- function(p,lo,hi,delta,se,g,retained,mode,test) {
  available <- is.finite(p)&is.finite(lo)&is.finite(hi)&is.finite(se)&se>0
  need(all(p[available]>=0&p[available]<=1),'Invalid p-value')
  pp <- rep(1,retained); pp[available] <- p[available]
  np <- sum(available)
  bin <- hist(p[available],breaks=seq(0,1,.05),plot=FALSE)$counts
  q0 <- p.adjust(pp,'BH',n=g); q1 <- p.adjust(pp,'BH',n=retained)
  row <- list(weight_mode=mode,test=test,genes=g,retained=retained,available=np,
    unavailable_retained=retained-np,filtered=g-retained,
    paired_truth=0,mean_delta=if(np)mean(delta[available]) else NA_real_,
    rms_delta=if(np)sqrt(mean(delta[available]^2)) else NA_real_,
    mean_se=if(np)mean(se[available]) else NA_real_,
    coverage95=if(np)mean(lo[available]<=0&hi[available]>=0) else NA_real_,
    delivered_coverage95=sum(lo[available]<=0&hi[available]>=0)/g,
    p_histogram=as.list(as.integer(bin)))
  for(alpha in c(.01,.05,.10)) {
    suffix <- sprintf('%03d',round(1000*alpha))
    row[[paste0('fpr_',suffix)]] <- if(np)mean(p[available]<=alpha) else NA_real_
    row[[paste0('discoveries_planned_',suffix)]] <- sum(q0[available]<=alpha)
    row[[paste0('discoveries_retained_',suffix)]] <- sum(q1[available]<=alpha)
    # Under the global paired-contrast null, FDP is one iff any rejection occurs.
    row[[paste0('fdp_planned_',suffix)]] <- if(np)as.numeric(any(q0[available]<=alpha)) else NA_real_
    row[[paste0('fdp_retained_',suffix)]] <- if(np)as.numeric(any(q1[available]<=alpha)) else NA_real_
  }
  row
}

one_dataset <- function(p,s,replicate,par) {
  d <- make_design(s)
  seed <- p$settings$seed+s$index*100000L+replicate
  sim <- simulate_counts(par,d,s,seed)
  g <- nrow(sim$y)
  # Symmetrize reference-column tie handling without changing the fixed design.
  permutation <- unlist(lapply(split(seq_len(s$n),d$x),sample),use.names=FALSE)
  y <- edgeR::DGEList(sim$y[,permutation,drop=FALSE])
  keep <- edgeR::filterByExpr(y,design=d$x0)
  need(sum(keep)>=100,'Fewer than 100 genes retained: pipeline cannot be calibrated')
  y <- y[keep,,keep.lib.sizes=FALSE]
  y <- edgeR::calcNormFactors(y,method='TMM')
  y <- y[,order(permutation),keep.lib.sizes=TRUE]
  need(identical(colnames(y),colnames(sim$y)),'Failed to restore the fixed sample order')
  v0 <- limma::voom(y,d$x0,plot=FALSE,adaptive.span=TRUE)
  v1 <- limma::voom(y,d$x1,plot=FALSE,adaptive.span=TRUE)
  need(identical(dimnames(v0$E),dimnames(v1$E))&&max(abs(v0$E-v1$E))<1e-10,
       'Two models do not have the same response')
  need(all(is.finite(v0$weights)&v0$weights>0)&all(is.finite(v1$weights)&v1$weights>0),'Invalid voom weights')
  base <- limma::lmFit(v0,d$x0)
  rows <- list(); checks <- list()
  for(mode in c('shared','design_specific')) {
    vf <- if(mode=='shared')v0 else v1
    full <- limma::lmFit(vf,d$x1)
    paired <- do.call(rbind,lapply(seq_len(nrow(v0$E)),function(i)
      paired_weights(v0$E[i,],v0$weights[i,],vf$weights[i,],d$x0,d$x1)))
    discrepancy <- max(abs(paired$delta-(full$coefficients[,2]-base$coefficients[,2])))
    need(is.finite(discrepancy)&&discrepancy<1e-7,'Paired kernel differs from the two lmFit coefficients')
    rows[[length(rows)+1L]] <- score_test(paired$p,paired$low,paired$high,paired$delta,
      paired$SE_delta,g,sum(keep),mode,'model_based')
    rows[[length(rows)+1L]] <- score_test(paired$p_HC3,paired$low_HC3,paired$high_HC3,paired$delta,
      paired$SE_HC3,g,sum(keep),mode,'HC3')
    checks[[mode]] <- list(lmfit_max_difference=discrepancy,max_leverage=max(paired$max_leverage),
      degeneracies=sum(paired$status!='available'))
  }
  list(scenario=s$id,replicate=replicate,seed=seed,metrics=rows,checks=checks,
       depth_range=as.list(range(sim$depth)),normalization_factor_range=as.list(range(y$samples$norm.factors)))
}
