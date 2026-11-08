normalize_selected <- function(count,d,keep=NULL) {
  permutation<-unlist(lapply(split(seq_len(ncol(count)),d$x),sample),use.names=FALSE)
  g<-edgeR::DGEList(count[,permutation,drop=FALSE])
  if(is.null(keep))keep<-edgeR::filterByExpr(g,design=d$x0)
  need(sum(keep)>=100,'Too few retained genes')
  g<-g[keep,,keep.lib.sizes=FALSE]
  g<-edgeR::calcNormFactors(g,method='TMM')
  g<-g[,order(permutation),keep.lib.sizes=TRUE]
  need(identical(colnames(g),colnames(count)),'Sample order mismatch')
  list(g=g,keep=keep)
}

make_tail_inputs <- function(p,s,replicate,par) {
  d<-make_design(s);seed<-p$settings$seed+s$index*100000L+replicate
  sim<-simulate_counts(par,d,s,seed)
  obs<-normalize_selected(sim$y,d)
  v0<-limma::voom(obs$g,d$x0,plot=FALSE,adaptive.span=TRUE)
  v1<-limma::voom(obs$g,d$x1,plot=FALSE,adaptive.span=TRUE)
  need(max(abs(v0$E-v1$E))<1e-10,'Responses differ between count models')
  # Independent NB draws share latent depth, not observed library totals or counts.
  set.seed(seed+50000000L)
  mu<-(par$mean_count%o%sim$depth)*2^(sim$beta%o%d$x)
  aux<-matrix(rnbinom(length(mu),mu=as.numeric(mu),size=rep(1/par$dispersion,s$n)),nrow(par),s$n)
  dimnames(aux)<-dimnames(sim$y)
  independent<-normalize_selected(aux,d,obs$keep)
  iv0<-limma::voom(independent$g,d$x0,plot=FALSE,adaptive.span=TRUE)
  iv1<-limma::voom(independent$g,d$x1,plot=FALSE,adaptive.span=TRUE)
  need(identical(dimnames(iv0$E),dimnames(v0$E)),'Independent weight alignment mismatch')
  unit<-matrix(1,nrow(v0$E),s$n,dimnames=dimnames(v0$E))
  set.seed(seed+100000000L)
  gaussian<-matrix(rnorm(length(v0$E)),nrow(v0$E),s$n)+matrix(3+.2*d$x,nrow(v0$E),s$n,byrow=TRUE)
  dimnames(gaussian)<-dimnames(v0$E)
  arms<-list(
    estimated_shared=list(E=v0$E,w0=v0$weights,w1=v0$weights),
    estimated_design=list(E=v0$E,w0=v0$weights,w1=v1$weights),
    independent_shared=list(E=v0$E,w0=iv0$weights,w1=iv0$weights),
    independent_design=list(E=v0$E,w0=iv0$weights,w1=iv1$weights),
    fixed_unit=list(E=v0$E,w0=unit,w1=unit),
    gaussian_unit=list(E=gaussian,w0=unit,w1=unit))
  list(d=d,seed=seed,arms=arms,keep=obs$keep,sim=sim)
}

tail_dataset <- function(p,s,replicate,par) {
  input<-make_tail_inputs(p,s,replicate,par)
  d<-input$d;keep<-input$keep;g<-nrow(par);k<-sum(keep)
  gene<-par[keep,,drop=FALSE];raw<-input$sim$y[keep,,drop=FALSE]
  cutoffs<-c(p1e6=1e-6,p1e5=1e-5,bh_first=.05/g,p5e5=5e-5,p1e4=1e-4,p1e3=1e-3,p01=.01,p05=.05)
  strata<-cut(gene$mean_count,breaks=c(-Inf,10,50,Inf),right=FALSE,labels=c('mean_lt10','mean10_lt50','mean_ge50'))
  metrics<-events<-list()
  for(mode in names(input$arms)) {
    a<-input$arms[[mode]]
    need(all(is.finite(a$w0)&a$w0>0)&all(is.finite(a$w1)&a$w1>0),'Invalid control weights')
    f0<-limma::lmFit(a$E,d$x0,weights=a$w0);f1<-limma::lmFit(a$E,d$x1,weights=a$w1)
    fit<-do.call(rbind,lapply(seq_len(k),function(i)paired_weights(a$E[i,],a$w0[i,],a$w1[i,],d$x0,d$x1)))
    discrepancy<-max(abs(fit$delta-(f1$coefficients[,2]-f0$coefficients[,2])))
    need(is.finite(discrepancy)&&discrepancy<1e-7,'Two-model coefficient identity failed')
    for(test in c('model_based','HC3')) {
      pp<-if(test=='HC3')fit$p_HC3 else fit$p
      se<-if(test=='HC3')fit$SE_HC3 else fit$SE_delta
      lo<-if(test=='HC3')fit$low_HC3 else fit$low
      hi<-if(test=='HC3')fit$high_HC3 else fit$high
      valid<-is.finite(pp)&is.finite(se)&se>0&is.finite(lo)&is.finite(hi)
      row<-score_test(pp,lo,hi,fit$delta,se,g,k,mode,test)
      row$p_histogram<-as.list(as.integer(hist(pp[valid],breaks=seq(0,1,.05),fuzz=0,plot=FALSE)$counts))
      row$lmfit_max_difference<-discrepancy
      row$min_p<-if(any(valid))min(pp[valid]) else NA_real_
      row$max_leverage<-max(fit$max_leverage)
      row$response<-if(mode=='gaussian_unit')'independent_gaussian' else 'observed_logCPM'
      tails<-lapply(names(cutoffs),function(n)list(label=n,cutoff=cutoffs[[n]],count=sum(pp[valid]<=cutoffs[[n]])))
      row$tails<-tails
      adj<-rep(1,k);adj[valid]<-pp[valid]
      qp<-p.adjust(adj,'BH',n=g);qr<-p.adjust(adj,'BH',n=k)
      row$strata<-lapply(levels(strata),function(label) {
        idx<-which(strata==label&valid)
        list(stratum=label,retained=sum(strata==label),available=length(idx),
             raw05=sum(pp[idx]<=.05),tail1e4=sum(pp[idx]<=1e-4),
             discoveries_planned=sum(qp[idx]<=.05),discoveries_retained=sum(qr[idx]<=.05))
      })
      selected<-which(valid&(pp<=1e-4|qp<=.05|qr<=.05))
      if(any(valid))selected<-unique(c(selected,which.min(ifelse(valid,pp,Inf))))
      selected<-selected[order(pp[selected])]
      row$event_candidates<-length(selected);row$events_truncated<-max(0,length(selected)-200)
      selected<-head(selected,200)
      if(length(selected))events[[length(events)+1L]]<-data.frame(scenario=s$id,replicate=replicate,
        arm=mode,test=test,gene_id=gene$gene_id[selected],response=row$response,
        generator_mean=gene$mean_count[selected],generator_dispersion=gene$dispersion[selected],
        observed_count_mean=rowMeans(raw)[selected],observed_zero_fraction=rowMeans(raw==0)[selected],
        beta_base=f0$coefficients[selected,2],beta_full=f1$coefficients[selected,2],
        delta=fit$delta[selected],SE_model=fit$SE_delta[selected],SE_HC3=fit$SE_HC3[selected],
        t_value=fit$delta[selected]/se[selected],df=fit$df[selected],p=pp[selected],
        q_planned=qp[selected],q_retained=qr[selected],max_leverage=fit$max_leverage[selected],
        full_weight_ratio=apply(a$w1,1,function(w)max(w)/min(w))[selected],
        bh_planned=qp[selected]<=.05,bh_retained=qr[selected]<=.05)
      metrics[[length(metrics)+1L]]<-row
    }
  }
  list(result=list(scenario=s$id,replicate=replicate,seed=input$seed,
       independent_seed=input$seed+50000000L,gaussian_seed=input$seed+100000000L,metrics=metrics),
       events=do.call(rbind,events))
}
