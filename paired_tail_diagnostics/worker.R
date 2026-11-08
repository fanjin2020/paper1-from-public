args<-commandArgs(TRUE);root<-normalizePath(args[1],winslash='/',mustWork=TRUE)
extra<-strsplit(Sys.getenv('PAPER1_CALIBRATION_LIBRARIES'),.Platform$path.sep,fixed=TRUE)[[1]]
.libPaths(c(extra[nzchar(extra)],.libPaths()))
suppressPackageStartupMessages(library(jsonlite));suppressPackageStartupMessages(library(edgeR));suppressPackageStartupMessages(library(limma))
source(file.path(root,'analysis/common.R'),encoding='UTF-8')
source(file.path(root,'analysis/matched/functions.R'),encoding='UTF-8')
source(file.path(root,'paired_count_calibration/engine.R'))
source(file.path(root,'paired_tail_diagnostics/engine.R'))
job<-jsonlite::read_json(args[2],simplifyVector=FALSE)
p<-jsonlite::read_json(job$protocol,simplifyVector=FALSE);s<-p$scenarios[[job$scenario_index]]
par<-generator_parameters(p)
for(rep in unlist(job$replicates)) {
  warnings<-character()
  a<-withCallingHandlers(tail_dataset(p,s,rep,par),warning=function(w) {
    warnings<<-unique(c(warnings,conditionMessage(w)));invokeRestart('muffleWarning')
  })
  dest<-file.path(p$output,'tasks',s$id,sprintf('r%05d',rep))
  dir.create(dest,recursive=TRUE,showWarnings=FALSE)
  wr(a$events,file.path(dest,'events.csv.gz'))
  r<-a$result;r$warnings<-as.list(warnings);r$run_id<-p$run_id;r$attempt<-job$attempt
  jsonlite::write_json(r,file.path(dest,'result.json.tmp'),auto_unbox=TRUE,na='null',digits=NA)
  need(file.rename(file.path(dest,'result.json.tmp'),file.path(dest,'result.json')),'Cannot publish result')
  cat(s$id,rep,'complete\n');flush.console()
}
