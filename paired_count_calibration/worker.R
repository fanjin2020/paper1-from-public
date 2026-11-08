args <- commandArgs(TRUE)
root <- normalizePath(args[2],winslash='/',mustWork=TRUE)
extra <- strsplit(Sys.getenv('PAPER1_CALIBRATION_LIBRARIES'),.Platform$path.sep,fixed=TRUE)[[1]]
.libPaths(c(extra[nzchar(extra)],.libPaths()))
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(limma))
source(file.path(root,'analysis','common.R'),encoding='UTF-8')
source(file.path(root,'analysis','matched','functions.R'),encoding='UTF-8')
source(file.path(root,'paired_count_calibration','engine.R'))
options(warn=1)
save_json <- function(x,path) {
  dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
  tmp <- paste0(path,'.tmp')
  jsonlite::write_json(x,tmp,auto_unbox=TRUE,pretty=FALSE,na='null',digits=NA)
  need(file.rename(tmp,path),paste('Cannot publish',path))
}
if(args[1]=='check') {
  need('adaptive.span'%in%names(formals(limma::voom)),'limma voom lacks adaptive.span')
  out <- list(R=R.version.string,packages=setNames(lapply(c('edgeR','limma','jsonlite'),
    function(z)as.character(packageVersion(z))),c('edgeR','limma','jsonlite')),
    libraries=as.list(.libPaths()),platform=R.version$platform)
  save_json(out,args[3])
} else if(args[1]=='run') {
  job <- jsonlite::read_json(args[3],simplifyVector=FALSE)
  p <- jsonlite::read_json(job$protocol,simplifyVector=FALSE)
  s <- p$scenarios[[job$scenario_index]]
  par <- generator_parameters(p)
  for(rep in unlist(job$replicates)) {
    warnings <- character()
    r <- withCallingHandlers(one_dataset(p,s,rep,par),warning=function(w) {
      warnings <<- unique(c(warnings,conditionMessage(w)));invokeRestart('muffleWarning')
    })
    r$warnings <- as.list(warnings);r$run_id <- p$run_id;r$attempt <- job$attempt
    dest <- file.path(p$output,'tasks',s$id,sprintf('r%05d',rep),'result.json')
    save_json(r,dest)
    cat(s$id,rep,'completed\n');flush.console()
  }
} else stop('Unknown worker command')
