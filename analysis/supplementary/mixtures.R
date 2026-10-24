extra <- Sys.getenv("R_LIBS", "")
if(nzchar(extra)) .libPaths(c(strsplit(extra,.Platform$path.sep,fixed=TRUE)[[1]],.libPaths()))

args<-commandArgs(TRUE); mode<-args[1]; src<-args[2]; out<-args[3]; aux<-args[4]
suppressPackageStartupMessages({library(edgeR);library(limma)})
need<-function(ok,msg)if(!isTRUE(ok))stop(msg)
rd<-function(p)read.csv(p,check.names=FALSE,stringsAsFactors=FALSE)
wr<-function(d,n)write.csv(d,file.path(out,n),row.names=FALSE,na='')
readmat<-function(p){d<-rd(p);need(!anyDuplicated(d[[1]]),'Duplicate genes'); m<-as.matrix(d[,-1,drop=FALSE]);rownames(m)<-d[[1]];need(is.numeric(m)&&all(is.finite(m))&&all(m>=0)&&all(m==round(m)),'Invalid raw counts');m}
bh<-function(p){good<-is.finite(p);q<-p.adjust(ifelse(good,p,1),'BH');q[!good]<-NA_real_;q}
capture.output(sessionInfo(),file=file.path(out,paste0(mode,'_R_session.txt')))
fitcontrast<-function(v,d,c){eBayes(contrasts.fit(lmFit(v,d),c))}
extract<-function(f){data.frame(gene_symbol=rownames(f),beta=f$coefficients[,1],
 SE=f$stdev.unscaled[,1]*sqrt(f$s2.post),unscaled=f$stdev.unscaled[,1],s2post=f$s2.post,
 residual_variance=f$sigma^2,df=rep_len(f$df.total,nrow(f)),p=f$p.value[,1])}
if(mode=='cell'){
 m<-rd(file.path(aux,'sea_metadata.csv'));m$unit_id<-as.character(m$unit_id);m$donor_id<-as.character(m$donor_id)
 need(!anyDuplicated(m$unit_id)&&!anyDuplicated(m[c('donor_id','context')]),'Duplicated cell units')
 contexts<-sort(unique(m$context));need(length(contexts)==24,'Changed annotation context universe')
 cols<-c('age','sex','clinical_group','pathology','PMI','RIN','pH','region')
 for(k in cols) need(all(tapply(m[[k]],m$donor_id,function(x)length(unique(x)))==1),paste('Donor metadata conflict',k))
 need(all(is.finite(m$n_cells))&&all(m$n_cells>=0)&&all(m$n_cells==round(m$n_cells)),'Invalid cell counts')
 need(length(unique(m$region))==1,'Multiple regions require separate models')
 classmap<-c(setNames(rep('Excitatory',9),c('L2/3 IT','L4 IT','L5 ET','L5 IT','L5/6 NP','L6 CT','L6 IT','L6 IT Car3','L6b')),
  setNames(rep('Inhibitory',9),c('Chandelier','Lamp5','Lamp5 Lhx6','Pax6','Pvalb','Sncg','Sst','Sst Chodl','Vip')),
  Astrocyte='Astrocyte',Immune='Immune',OPC='OPC',Oligodendrocyte='Oligodendrocyte',Endothelial='Vascular','VLMC & Perivascular'='Vascular')
 need(setequal(contexts,names(classmap)),'Source subclasses differ from the existing 203 broad-class mapping')
 wr(data.frame(source_context=names(classmap),broad_class=unname(classmap)),'323_subclass_mapping.csv')
 rawcounts<-readmat(src);need(all(m$unit_id%in%colnames(rawcounts)),'Count/metadata IDs differ')
 expected_totals<-rowSums(rawcounts[,m$unit_id,drop=FALSE])
 m$source_context<-m$context;m$context<-unname(classmap[m$context])
 m$target_unit<-paste(m$donor_id,m$context,sep='__broad__')
 wr(m[,c('unit_id','donor_id','source_context','context','target_unit','n_cells')],'323_subclass_aggregation_trace.csv')
 pieces<-split(m,m$target_unit)
 counts<-matrix(0,nrow(rawcounts),length(pieces),dimnames=list(rownames(rawcounts),names(pieces)))
 agg<-list()
 for(id in names(pieces)){
  z<-pieces[[id]];counts[,id]<-rowSums(rawcounts[,z$unit_id,drop=FALSE])
  a<-z[1,,drop=FALSE];a$unit_id<-id;a$n_cells<-sum(z$n_cells);a$observed_subclasses<-nrow(z)
  a$source_context<-paste(sort(z$source_context),collapse=';');agg[[id]]<-a
 }
 aggregation_error<-max(abs(rowSums(counts)-expected_totals))
 wr(data.frame(genes=nrow(counts),maximum_count_conservation_error=aggregation_error),'323_aggregation_checks.csv')
 need(aggregation_error==0,'Broad-class aggregation changed total counts')
 m<-do.call(rbind,agg);rm(rawcounts);contexts<-sort(unique(m$context));need(length(contexts)==7,'Broad class count changed')
 wr(m,'323_all_broad_units.csv')
 membership<-do.call(rbind,lapply(split(m,m$donor_id),function(z)data.frame(donor_id=z$donor_id[1],
   available_contexts=nrow(z),all_7_broad_contexts=nrow(z)==7,min_recovered_cells=min(z$n_cells),
   quality_and_demographics_complete=all(complete.cases(z[,cols])),
   eligible=nrow(z)==7&&all(z$n_cells>=20)&&all(complete.cases(z[,cols])))))
 wr(membership,'323_donor_membership.csv')
 eligible<-setNames(membership$eligible,membership$donor_id)
 donors<-names(eligible)[eligible];need(length(donors)>=20,'Fewer than 20 donors with complete 7-broad-class coverage')
 m<-m[m$donor_id%in%donors,];d<-m[match(donors,m$donor_id),]
 need(all(d$clinical_group%in%c('Dementia','No_dementia'))&&min(table(d$clinical_group))>=8,'Clinical group coverage')
 need(all(d$sex%in%c('male','female')),'Sex coding')
 group<-as.numeric(d$clinical_group=='Dementia');age<-d$age-mean(d$age)
 x<-cbind(intercept=1,effect=group,age=age,male=as.numeric(d$sex=='male'),
          PMI=d$PMI-mean(d$PMI),RIN=d$RIN-mean(d$RIN),pH=d$pH-mean(d$pH))
 need(qr(x)$rank==ncol(x)&&nrow(x)-ncol(x)>=10,'Cell-context design not estimable; no automatic covariate removal')
 wr(data.frame(donor_id=donors,x,check.names=FALSE),'323_matched_design.csv');wr(m,'323_matched_units.csv')
 need(all(m$unit_id%in%colnames(counts)),'Aggregated count/metadata IDs differ')
 nc<-sapply(contexts,function(ct)m$n_cells[match(paste(donors,ct),paste(m$donor_id,m$context))])
 fractions<-nc/rowSums(nc);wr(data.frame(donor_id=rep(donors,length(contexts)),context=rep(contexts,each=length(donors)),fraction=as.vector(fractions)),'323_recovered_fractions.csv')
 logexpr<-list();cpmsets<-list();keep<-rep(TRUE,nrow(counts));rownames(x)<-donors
 for(ct in contexts){
  u<-m$unit_id[match(paste(donors,ct),paste(m$donor_id,m$context))]
  y<-DGEList(counts[,u,drop=FALSE]);need(all(y$samples$lib.size>0),'Empty context library')
  k<-filterByExpr(y,design=x);keep<-keep & k
  y<-calcNormFactors(y);cp<-cpm(y,log=FALSE);cpmsets[[ct]]<-cp
 }
 need(sum(keep)>=100,'Insufficient common measurable background')
 genes<-rownames(counts)[keep];wr(data.frame(gene_symbol=rownames(counts),common=keep),'323_common_background.csv')
 observed<-fixed<-matrix(0,length(genes),length(donors),dimnames=list(genes,donors))
 for(i in seq_along(contexts)){
  cp<-cpmsets[[contexts[i]]][keep,,drop=FALSE]
  observed<-observed+sweep(cp,2,fractions[,i],'*')
  fixed<-fixed+cp*mean(fractions[,i])
  logexpr[[contexts[i]]]<-log2(cp+.5)
 }
 logexpr[['recovered_weighted_synthetic_mixture']]<-log2(observed+.5)
 logexpr[['fixed_weight_synthetic_mixture']]<-log2(fixed+.5)
 ep<-sort(unique(rd(file.path(aux,'panel.csv'))$gene_symbol));need(length(ep)>0,'Empty discovery panel')
 res<-list()
 for(name in names(logexpr)){
  y<-logexpr[[name]];colnames(y)<-donors
  fit<-eBayes(lmFit(y,x),trend=TRUE);pos<-match(ep,genes)
  se<-fit$stdev.unscaled[,'effect']*sqrt(fit$s2.post);beta<-fit$coefficients[,'effect'];ci<-qt(.975,fit$df.total)*se
  res[[name]]<-data.frame(context=name,gene_symbol=ep,beta=beta[pos],SE=se[pos],CI_low=(beta-ci)[pos],CI_high=(beta+ci)[pos],p=fit$p.value[pos,'effect'],n=length(donors))
 }
 z<-do.call(rbind,res);z$q_237_by_9<-bh(z$p);wr(z,'323_matched_cell_and_synthetic_mixture_results.csv')
 frac<-list()
 for(i in seq_along(contexts)){
  y<-qlogis((nc[,i]+.5)/(rowSums(nc)+1));fit<-lm.fit(x,y);df<-length(y)-fit$rank
  se<-sqrt(sum(fit$residuals^2)/df*solve(crossprod(x))[2,2]);b<-fit$coefficients[2]
  frac[[i]]<-data.frame(context=contexts[i],beta=b,SE=se,CI_low=b-qt(.975,df)*se,CI_high=b+qt(.975,df)*se,p=2*pt(-abs(b/se),df))
 }
 f<-do.call(rbind,frac);f$q_7<-bh(f$p);wr(f,'323_recovered_composition_associations.csv')
 writeLines(c('Same-donor clinical comparison in one existing SEA-AD source; NOT independent validation.',
 'Recovered nuclei fractions are not tissue cell fractions or RNA contribution fractions.',
 'Synthetic CPM mixtures are not measured bulk expression and do not estimate a causal percentage.',
 'Annotations are inherited from the archived source; this step does not certify their accuracy.',
 'Cross-context contrasts are descriptive, not tests of cell specificity or mediation.',
 'Synthetic-mixture pointwise intervals condition on estimated cohort reference weights.',
 'PMI/RIN/pH remain source-defined proxies; unrecorded batches and source dependencies remain.'),file.path(out,'323_interpretation.txt'))
}
