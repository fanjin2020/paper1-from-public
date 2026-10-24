args <- commandArgs(TRUE)
this <- normalizePath(sub('^--file=', '', commandArgs(FALSE)[grep('^--file=',commandArgs(FALSE))][1]))
here <- dirname(this)
libs <- strsplit(Sys.getenv('PAPER1_R_LIBRARIES'), .Platform$path.sep, fixed=TRUE)[[1]]
.libPaths(unique(c(libs[nzchar(libs)], .libPaths())))
options(stringsAsFactors=FALSE, warn=1, digits=17)
packages <- c('jsonlite','edgeR','limma','variancePartition','BiocParallel','lme4','lmerTest','statmod')
if (args[1]=='check') {
  missing <- packages[!vapply(packages,requireNamespace,logical(1),quietly=TRUE)]
  if(length(missing) && args[4]=='TRUE') {
    local <- args[3]
    options(repos=c(CRAN='https://cloud.r-project.org'),timeout=1200)
    cran <- intersect(missing,c('jsonlite','lme4','lmerTest','statmod'))
    if(length(cran))install.packages(cran,lib=local,Ncpus=2)
    bioc <- setdiff(missing,cran)
    if(length(bioc)) {
      if(!requireNamespace('BiocManager',quietly=TRUE))install.packages('BiocManager',lib=local)
      BiocManager::install(bioc,lib=local,ask=FALSE,update=FALSE,Ncpus=2)
    }
  }
  for(p in packages)if(!requireNamespace(p,quietly=TRUE))stop(paste('Missing R package:',p))
  if(!capabilities('cairo'))stop('This R build needs Cairo PDF support for the vector result figures')
  info <- list(R=R.version.string,packages=setNames(lapply(packages,function(p)
    list(version=as.character(packageVersion(p)),path=find.package(p))),packages),
    loaded_namespaces=setNames(lapply(sort(loadedNamespaces()),function(p)as.character(packageVersion(p))),
                               sort(loadedNamespaces())),cairo=unname(capabilities('cairo')))
  jsonlite::write_json(info,args[2],auto_unbox=TRUE,pretty=TRUE)
  cat('R dependencies passed\n')
  quit(status=0)
}
suppressPackageStartupMessages({library(edgeR);library(limma);library(variancePartition);library(BiocParallel)})
source(file.path(here,'..','common.R'))
source(file.path(here,'..','diagnostics.R'))
source(file.path(here,'functions.R'))
source(file.path(here,'prepare.R'))
source(file.path(here,'human.R'))
source(file.path(here,'mouse.R'))
source(file.path(here,'counts.R'))
source(file.path(here,'summarise.R'))
if(args[1]=='prepare')prepare_inputs(jsonlite::fromJSON(args[2]))
if(args[1]=='job') {
  job<-jsonlite::fromJSON(args[2]); cfg<-jsonlite::fromJSON(job$config)
  switch(job$phase,human=run_human(job,cfg),mouse=run_mouse(job,cfg),counts=run_counts(job,cfg),stop('Unknown phase'))
}
if(args[1]=='summarise')summarise_all(jsonlite::fromJSON(args[2]))
