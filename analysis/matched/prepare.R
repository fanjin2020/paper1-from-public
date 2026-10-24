prepare_inputs <- function(cfg) {
  src<-cfg$source;out<-cfg$out
  for(k in c('cache','results','figures','runtime'))dir.create(file.path(out,k),recursive=TRUE,showWarnings=FALSE)
  panel<-rd(file.path(src,'results','panel.csv'))$gene_symbol
  ortholog<-rd(file.path(src,'results','ortholog_panel.csv'))
  need(!anyDuplicated(panel)&&all(nzchar(panel)),'Duplicate or empty human panel genes')
  need(all(c('gene_symbol','mouse_gene')%in%names(ortholog)),'Ortholog columns missing')
  md<-rd(file.path(src,'prepared','human','metadata.csv'))
  required<-c('unit_id','donor_id','context','n_cells','clinical_group','age','sex','pathology','PMI','RIN','pH')
  need(all(required%in%names(md)),'Human metadata columns missing')
  need(all(na.omit(md$sex)%in%c('male','female')),'Unknown human source sex coding')
  y<-counts(file.path(src,'prepared','human','expression.csv.gz'))
  need(identical(colnames(y),md$unit_id),'Human sample order differs from prepared input')
  contexts<-sort(unique(md$context))
  if(!cfg$smoke) {
    need(length(panel)==237,'Expected 237-gene panel')
    need(length(unique(md$donor_id))==80 && length(contexts)==24,'Expected 80 donors and 24 contexts')
  }
  human_manifest<-list();sex_tables<-list()
  for(i in seq_along(contexts)) {
    context<-contexts[i]
    selected<-md$context==context & md$n_cells>=cfg$minimum_nuclei &
      md$clinical_group%in%c('Dementia','No_dementia') & complete.cases(md[,c('age','sex','pathology','PMI','RIN','pH')])
    d<-md[selected,,drop=FALSE]
    need(!anyDuplicated(d$donor_id),paste('Repeated donor in',context))
    x<-human_design(d,'clinical_status')
    viable<-nrow(d)>=10 && all(table(factor(d$clinical_group,levels=c('No_dementia','Dementia')))>=5)&&estimable(x)
    g<-v<-NULL
    if(viable) {
      g<-DGEList(y[,selected,drop=FALSE]);keep<-filterByExpr(g,design=x)
      g<-g[keep,,keep.lib.sizes=FALSE]
      viable<-nrow(g)>=100
      if(viable){g<-calcNormFactors(g);v<-voom(g,x,plot=FALSE,adaptive.span=TRUE)}
    }
    save_rds(list(context=context,metadata=d,dge=g,vref=v,panel=panel,viable=viable),
             file.path(out,'cache',sprintf('human_%03d.rds',i)))
    human_manifest[[i]]<-data.frame(context=context,n=nrow(d),background=if(viable)nrow(g) else 0,
      available=viable,units=paste(d$unit_id,collapse=';'))
    z<-as.data.frame(table(clinical_group=d$clinical_group,sex=d$sex));z$context<-context
    sex_tables[[i]]<-z
    cat('Prepared human',i,'/',length(contexts),context,'n=',nrow(d),'\n')
  }
  wr(do.call(rbind,human_manifest),file.path(out,'results','human_input_manifest.csv'))
  wr(bind_rows(sex_tables),file.path(out,'results','human_context_sex_counts.csv'))
  # Count donors once, not once per context. Require source metadata to be consistent.
  donor<-unique(md[,c('donor_id','clinical_group','sex')])
  need(!anyDuplicated(donor$donor_id),'Inconsistent sex or clinical labels across contexts')
  wr(as.data.frame(table(clinical_group=donor$clinical_group,sex=donor$sex)),file.path(out,'results','human_donor_sex_counts.csv'))
  rm(y);gc()
  d<-rd(file.path(src,'prepared','mouse_units.csv'),colClasses=c(animal_id='character',cage='character',slide='character'))
  need(all(c('animal_id','genotype','sex','age','cage','slide')%in%names(d)),'Mouse metadata columns missing')
  need(!anyNA(d[,c('animal_id','genotype','sex','age','cage','slide')]),'Missing mouse metadata')
  need(all(d$genotype%in%c('5X','WT'))&&all(d$sex%in%c('F','M')),'Unknown mouse genotype or sex coding')
  d$genotype<-as.numeric(d$genotype=='5X');d$sex<-factor(d$sex,levels=c('F','M'))
  d$age_factor<-factor(d$age);d$cage<-factor(d$cage);d$slide<-factor(d$slide);rownames(d)<-d$animal_id
  need(!anyDuplicated(d$animal_id),'Repeated mouse animal')
  if(!cfg$smoke)need(nrow(d)==74&&nlevels(d$cage)==36,'Expected 74 animals and 36 cages')
  sources<-setNames(lapply(c('whole','cortex'),function(region)
    counts(file.path(src,'prepared',paste0('mouse_',region,'_counts.csv.gz')))),c('whole','cortex'))
  for(y in sources)need(identical(colnames(y),d$animal_id),'Mouse count columns differ from metadata order')
  x<-mouse_designs(d,'linear')$base$x
  background<-Reduce(intersect,lapply(sources,function(y)rownames(y)[filterByExpr(DGEList(y),design=x)]))
  if(!cfg$smoke)need(length(background)==13312,'Shared mouse background does not contain 13,312 genes')
  wr(data.frame(gene_symbol=background),file.path(out,'results','mouse_background.csv'))
  wr(as.data.frame(table(genotype=d$genotype,sex=d$sex)),file.path(out,'results','mouse_sex_by_genotype.csv'))
  mixed<-tapply(d$genotype,d$cage,function(z)length(unique(z))==2)
  geometry<-data.frame(animals=nrow(d),cages=nlevels(d$cage),mixed_cages=sum(mixed),
    animals_in_mixed_cages=sum(d$cage%in%names(mixed)[mixed]),mutants=sum(d$genotype))
  wr(geometry,file.path(out,'results','mouse_layout.csv'))
  specs<-character()
  for(region in names(sources))for(kind in c('linear','age_factor','age_factor_slide')) {
    spec<-paste(region,kind,sep='_');specs<-c(specs,spec)
    g<-DGEList(sources[[region]][background,,drop=FALSE]);g<-calcNormFactors(g)
    save_rds(list(spec=spec,kind=kind,metadata=d,dge=g,panel=ortholog,genes=background),
             file.path(out,'cache',paste0('mouse_',spec,'.rds')))
  }
  g<-DGEList(sources$whole[background,,drop=FALSE]);g<-calcNormFactors(g)
  calibration<-estimateDisp(g,x,robust=TRUE)
  phi<-calibration$tagwise.dispersion
  need(length(phi)==length(background)&&all(is.finite(phi))&&all(phi>=0),'Invalid calibrated NB dispersion')
  raw_depth<-colSums(g$counts)
  mean_count<-rowMeans(t(t(g$counts)/raw_depth))*median(raw_depth)
  need(all(mean_count>0),'Zero calibrated mean count')
  selected<-seq_along(background)
  if(cfg$sim_genes>0 && cfg$sim_genes<length(selected)) {
    ordered<-order(mean_count,background)
    selected<-sort(ordered[unique(round(seq(1,length(ordered),length.out=cfg$sim_genes)))])
  }
  par<-data.frame(gene_symbol=background[selected],mean_count=mean_count[selected],
    dispersion=pmax(phi[selected],1e-8),dispersion_floored=phi[selected]<1e-8)
  wr(par,file.path(out,'results','count_generator_parameters.csv'))
  save_rds(list(parameters=par,metadata=d,depth_ratio=raw_depth/median(raw_depth),
               full_background=length(background)),file.path(out,'cache','count_calibration.rds'))
  scenarios<-expand.grid(cage_log2_sd=c(0,.3),cage_mean_amplitude=c(0,.5),de_fraction=c(0,.1))
  scenarios$scenario<-seq_len(nrow(scenarios))-1L
  scenarios$effect_log2<-.5;scenarios$mean_affected_fraction<-.2
  scenarios$interpretation<-'Known library-offset log2 count-mean coefficient; latent cage SD is not count ICC'
  wr(scenarios,file.path(out,'results','count_scenarios.csv'))
  save_json(list(contexts=as.list(contexts),mouse_specs=as.list(specs),human_donors=nrow(donor),
    mouse_animals=nrow(d),mouse_genes=length(background),simulation_genes=nrow(par),
    simulation_scope=if(cfg$sim_genes==0)'all_common_background' else 'reduced_expression_stratified_subset',
    validation_only=cfg$smoke),file.path(out,'runtime','prepared_manifest.json'))
}
