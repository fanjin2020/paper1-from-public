human_variant_design <- function(d,contrast,quality,categorical) {
  x<-human_design(d,contrast,quality)
  if(categorical && contrast=='clinical_status_adjusted_for_braak') {
    x<-x[,setdiff(colnames(x),'Braak'),drop=FALSE]
    braak<-if(length(unique(d$pathology))>1)model.matrix(~factor(pathology),d)[,-1,drop=FALSE] else matrix(numeric(),nrow(d),0)
    colnames(braak)<-paste0('Braak_level_',seq_len(ncol(braak)))
    # The target remains clinical status; the Braak-slope contrast remains ordinal.
    x<-cbind(x,braak)
  }
  x
}

run_human <- function(job,cfg) {
  obj<-readRDS(file.path(cfg$out,'cache',sprintf('human_%03d.rds',job$index)))
  d<-obj$metadata;genes<-obj$panel
  contrasts<-c('clinical_status','clinical_status_adjusted_for_braak','braak_same_clinical_population')
  rows<-pairs<-designs<-list()
  for(categorical in c(FALSE,TRUE))for(weight in c('shared','design'))for(contrast in contrasts) {
    spec<-paste(obj$context,contrast,sep='::')
    x0<-human_variant_design(d,contrast,FALSE,categorical)
    x1<-human_variant_design(d,contrast,TRUE,categorical)
    viable<-obj$viable&&estimable(x0)&&estimable(x1)
    v0<-v1<-NULL
    if(viable) {
      v0<-if(weight=='shared')obj$vref else voom(obj$dge,x0,plot=FALSE,adaptive.span=TRUE)
      v1<-if(weight=='shared')obj$vref else voom(obj$dge,x1,plot=FALSE,adaptive.span=TRUE)
      need(identical(rownames(v0$E),rownames(obj$vref$E))&&identical(rownames(v1$E),rownames(obj$vref$E)),
           'Sensitivity changed human gene background')
      need(max(abs(v0$E-v1$E))<1e-10,'Weight sensitivity unexpectedly changed the response matrix')
    }
    braak_label<-if(categorical)'categorical_nuisance' else 'ordinal_nuisance'
    key<-paste(weight,braak_label,sep='__')
    designs[[length(designs)+1]]<-data.frame(context=obj$context,contrast=contrast,weight=weight,
      braak=braak_label,n=nrow(d),p_base=ncol(x0),p_quality=ncol(x1),available=viable,
      background=if(obj$viable)nrow(obj$dge) else 0,
      changed_estimand=FALSE,categorical_applies=contrast=='clinical_status_adjusted_for_braak'&&categorical)
    for(robust in c(FALSE,TRUE)) {
      variant<-paste(key,if(robust)'robust' else 'ordinary',sep='__')
      for(arm in c('base','quality')) {
        z<-if(viable)fit_linear(if(arm=='base')v0 else v1,if(arm=='base')x0 else x1,
          genes,arm,spec,robust,target='target') else unavailable(genes,arm,spec,'not_estimable')
        z$variant<-variant;z$context<-obj$context;z$contrast<-contrast;z$n<-nrow(d)
        z$weight<-weight;z$braak<-braak_label;z$robust<-robust
        rows[[length(rows)+1]]<-z
      }
    }
    ptab<-data.frame(gene_symbol=genes,spec=spec,context=obj$context,contrast=contrast,
      paired_variant=key,weight=weight,braak=braak_label,status='not_expressed_or_not_estimable')
    if(viable) {
      available<-intersect(genes,rownames(v0$E))
      values<-lapply(available,function(gene) {
        i<-match(gene,rownames(v0$E))
        cbind(gene_symbol=gene,paired_weights(v0$E[i,],v0$weights[i,],v1$weights[i,],x0,x1))
      })
      v<-bind_rows(values)
      if(nrow(v)) {
        ix<-match(ptab$gene_symbol,v$gene_symbol)
        for(k in setdiff(names(v),'gene_symbol'))ptab[[k]]<-v[[k]][ix]
        ptab$status[is.na(ix)]<-'not_expressed_or_not_estimable'
      }
    }
    ptab$variance_scope<-'Conditional on estimated weights; model-based variance uses the extended-model precision; HC3 uses extended residuals'
    pairs[[length(pairs)+1]]<-ptab
    cat('Human sensitivity',job$id,key,contrast,'\n')
  }
  result<-bind_rows(rows);paired<-bind_rows(pairs)
  for(k in c('delta','SE_delta','p','low','high','SE_HC3','p_HC3','low_HC3','high_HC3'))
    if(!k%in%names(paired))paired[[k]]<-NA_real_
  wr(result,file.path(job$output,'fits.csv.gz'))
  wr(paired,file.path(job$output,'paired.csv.gz'))
  wr(bind_rows(designs),file.path(job$output,'designs.csv'))
  finish_job(job,c('fits.csv.gz','paired.csv.gz','designs.csv'))
}
