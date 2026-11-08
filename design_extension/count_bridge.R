args <- commandArgs(TRUE)
this <- normalizePath(sub('^--file=', '', commandArgs(FALSE)[grep('^--file=',commandArgs(FALSE))][1]))
here <- dirname(this)
libs <- strsplit(Sys.getenv('PAPER1_R_LIBRARIES'), .Platform$path.sep, fixed=TRUE)[[1]]
.libPaths(unique(c(libs[nzchar(libs)], .libPaths())))
options(stringsAsFactors=FALSE, warn=1)
packages <- c('jsonlite','edgeR','limma','variancePartition','BiocParallel','lme4','lmerTest','statmod')
missing <- packages[!vapply(packages,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing))stop(paste('Missing or unloadable existing R packages:',paste(missing,collapse=', '),
  '\nChecked libraries:',paste(.libPaths(),collapse='; ')))
if(!'adaptive.span'%in%names(formals(limma::voom)))
  stop('The existing limma version lacks adaptive.span; select the library used by the completed reproduction')
if(args[1]=='check') {
  info <- list(R=R.version.string, packages=setNames(lapply(packages,function(p)
    list(version=as.character(packageVersion(p)),path=find.package(p))),packages))
  jsonlite::write_json(info,args[2],auto_unbox=TRUE,pretty=TRUE)
  quit(status=0)
}
suppressPackageStartupMessages({library(edgeR);library(limma);library(variancePartition);library(BiocParallel)})
source(file.path(here,'engine','common.R'))
source(file.path(here,'engine','diagnostics.R'))
source(file.path(here,'engine','functions.R'))
source(file.path(here,'engine','mouse.R'))
source(file.path(here,'engine','counts.R'))

metadata <- function(path) {
  d<-rd(path,colClasses=c(animal_id='character',cage='character',slide='character'))
  need(all(c('animal_id','genotype','sex','age','cage','slide')%in%names(d)), 'Missing metadata columns')
  need(!anyNA(d[,c('animal_id','genotype','sex','age','cage','slide')])&&!anyDuplicated(d$animal_id), 'Invalid metadata')
  need(all(d$genotype%in%c('5X','WT'))&&all(d$sex%in%c('F','M')), 'Unknown genotype or sex')
  d$genotype<-as.numeric(d$genotype=='5X');d$sex<-factor(d$sex,levels=c('F','M'))
  d$age_factor<-factor(d$age);d$cage<-factor(d$cage);d$slide<-factor(d$slide)
  rownames(d)<-d$animal_id;d
}

if(args[1]=='prepare') {
  cfg<-jsonlite::fromJSON(args[2]);src<-cfg$data_work
  d<-metadata(file.path(src,'prepared','mouse_units.csv'))
  whole<-counts(file.path(src,'prepared','mouse_whole_counts.csv.gz'))
  cortex<-counts(file.path(src,'prepared','mouse_cortex_counts.csv.gz'))
  need(identical(colnames(whole),d$animal_id)&&identical(colnames(cortex),d$animal_id),'Input sample order differs')
  x<-mouse_designs(d,'linear')$base$x
  background<-intersect(rownames(whole)[filterByExpr(DGEList(whole),design=x)],
                        rownames(cortex)[filterByExpr(DGEList(cortex),design=x)])
  need(length(background)>=100,'Too few common eligible genes')
  g<-calcNormFactors(DGEList(whole[background,,drop=FALSE]))
  fitted<-estimateDisp(g,x,robust=TRUE)
  depth<-colSums(g$counts)
  mu<-rowMeans(t(t(g$counts)/depth))*median(depth)
  phi<-fitted$tagwise.dispersion
  need(all(is.finite(mu)&mu>0)&all(is.finite(phi)&phi>=0),'Invalid count calibration')
  selected<-seq_along(background)
  if(cfg$count_genes>0&&cfg$count_genes<length(selected)) {
    ordered<-order(mu,background)
    selected<-sort(ordered[unique(round(seq(1,length(ordered),length.out=cfg$count_genes)))])
  }
  parameters<-data.frame(gene_symbol=background[selected],mean_count=mu[selected],
    dispersion=pmax(phi[selected],1e-8),dispersion_floored=phi[selected]<1e-8)
  save_rds(list(parameters=parameters,depth_ratio=depth/median(depth),full_background=length(background)),
    file.path(cfg$output,'calibration.rds'))
  wr(parameters,file.path(cfg$output,'calibration_parameters.csv'))
  save_json(list(full_background=length(background),used_genes=nrow(parameters),
    scope=if(cfg$count_genes>0)'expression_stratified_subset' else 'complete_common_background',
    provenance='Fresh edgeR calibration from prepared public counts; no fitted-result reuse'),
    file.path(cfg$output,'calibration_info.json'))
  quit(status=0)
}

if(args[1]!='job')stop('Unknown command')
job<-jsonlite::fromJSON(args[2]);cfg<-jsonlite::fromJSON(job$config)
cal<-readRDS(file.path(cfg$output,'calibration.rds'))
d<-metadata(job$design_path);cal$metadata<-d
# Same empirical depth distribution is deterministically resampled to each layout.
set.seed(as.integer(cfg$seed+91L))
cal$depth_ratio<-sample(cal$depth_ratio,nrow(d),replace=TRUE)
models<-c('baseline','cage_fixed','cage_correlation','dream_native','mundlak_native')
genes<-cal$parameters$gene_symbol;records<-failures<-list()
diagnostics<-retries<-messages<-unavailable_rows<-list()
scenario<-as.data.frame(job$scenario)
for(rep in seq.int(job$start,job$end-1L)) {
  started<-proc.time()[3]
  paired_scenario<-scenario;paired_scenario$scenario<-0L
  generated<-generate_counts(cal,paired_scenario,rep,as.integer(cfg$seed+job$design_index*10007L))
  g<-DGEList(generated$counts);x<-mouse_designs(d,'linear')$base$x
  keep<-filterByExpr(g,design=x)
  Sys.setenv(PAPER1_DIAGNOSTICS_DIR=file.path(job$output,'numerical_details',paste0('r',rep)))
  error<-NULL
  ans<-tryCatch({
    need(sum(keep)>=100,'Fewer than 100 filtered genes')
    g<-calcNormFactors(g[keep,,keep.lib.sizes=FALSE])
    fit_mouse_models(g,d,'linear',job$design,1L,FALSE)
  },error=function(e){error<<-conditionMessage(e);NULL})
  for(model in models) {
    tab<-unavailable(genes,model,job$design,'filtered_out');tab$status[keep]<-'fit_failed'
    if(!is.null(ans)) {
      part<-ans$table[ans$table$model==model,,drop=FALSE];ix<-match(part$gene_symbol,genes)
      need(!anyNA(ix)&&!anyDuplicated(ix),'Unexpected gene alignment')
      common<-intersect(names(tab),names(part));tab[ix,common]<-part[,common]
    }
    bad_rows<-tab[!tab$status%in%c('available','target_not_estimable'),
                  c('gene_symbol','model','status'),drop=FALSE]
    unavailable_rows[[length(unavailable_rows)+1]]<-tag_count_rows(bad_rows,job,rep)
    row<-performance_row(tab,generated$beta,generated$gamma,scenario,rep,model,keep)
    row$design<-job$design;row$split<-job$split
    row$seconds_dataset_so_far<-proc.time()[3]-started
    records[[length(records)+1]]<-row
    bad<-table(tab$status[tab$status!='available'])
    if(length(bad))failures[[length(failures)+1]]<-data.frame(design=job$design,scenario=scenario$scenario,
      replicate=rep,model=model,status=names(bad),genes=as.integer(bad))
  }
  if(!is.null(error)) {
    cat('DATASET_FAILURE',rep,error,'\n')
    messages[[length(messages)+1]]<-tag_count_rows(data.frame(model='all',message=error),job,rep)
  }
  if(!is.null(ans)) {
    diagnostics[[length(diagnostics)+1]]<-tag_count_rows(ans$diagnostics,job,rep)
    retries[[length(retries)+1]]<-tag_count_rows(ans$retries,job,rep)
    messages[[length(messages)+1]]<-tag_count_rows(ans$messages,job,rep)
    if(nrow(ans$messages))print(ans$messages)
  }
  save_count_progress(job,rep-job$start+1L)
  cat('COMPLETED_DATASET',rep,'seconds',proc.time()[3]-started,'\n')
}
wr(bind_rows(records),file.path(job$output,'metrics.csv'))
wr(bind_rows(failures),file.path(job$output,'failure_counts.csv'))
wr(bind_rows(diagnostics),file.path(job$output,'diagnostics.csv'))
wr(bind_rows(retries),file.path(job$output,'retry_attempts.csv'))
wr(bind_rows(messages),file.path(job$output,'messages.csv'))
wr(bind_rows(unavailable_rows),file.path(job$output,'unavailable_genes.csv'))
