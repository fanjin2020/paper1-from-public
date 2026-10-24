# Version checks only: this script never installs or downloads packages.
a <- commandArgs(TRUE)
if(length(a)<1L) stop('An environment specification is required')
libs <- if(length(a)>1L) a[-1] else character()
.libPaths(unique(c(libs[dir.exists(libs)], .libPaths())))
if(!requireNamespace('jsonlite',quietly=TRUE)) stop('Missing jsonlite in the supplied R libraries')
expected <- jsonlite::fromJSON(a[1])
errors <- character()
if(as.character(getRversion())!=expected$R)
  errors <- c(errors,paste('R: expected',expected$R,'found',getRversion()))
versions <- paths <- list()
for(p in names(expected$packages)) {
  ok <- tryCatch(requireNamespace(p,quietly=TRUE),error=function(e)FALSE)
  if(!ok) {
    errors <- c(errors,paste('Missing or unloadable R package:',p))
    next
  }
  versions[[p]] <- as.character(utils::packageVersion(p))
  paths[[p]] <- find.package(p)
  if(versions[[p]]!=expected$packages[[p]])
    errors <- c(errors,paste(p,'expected',expected$packages[[p]],'found',versions[[p]]))
}
if(length(errors))stop(paste(errors,collapse='\n'))
namespaces <- sort(loadedNamespaces())
loaded <- setNames(lapply(namespaces,function(p)as.character(utils::packageVersion(p))),namespaces)
session <- utils::sessionInfo()
record <- list(R=as.character(getRversion()),packages=versions,paths=paths,libraries=.libPaths(),
               loaded_namespaces=loaded,platform=R.version$platform,
               BLAS=session$BLAS,LAPACK=session$LAPACK)
cat('PAPER1_ENVIRONMENT\t',jsonlite::toJSON(record,auto_unbox=TRUE), '\n',sep='')
