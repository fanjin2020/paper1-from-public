generate_counts <- function(calibration,scenario,replicate,seed) {
  p<-calibration$parameters;d<-calibration$metadata;G<-nrow(p);N<-nrow(d)
  set.seed(as.integer(seed+100003L*replicate+1009L*scenario$scenario))
  beta<-gamma<-rep(0,G)
  alternatives<-if(scenario$de_fraction>0)sample.int(G,max(2L,round(G*scenario$de_fraction))) else integer()
  if(length(alternatives))beta[alternatives]<-rep(c(-1,1),length.out=length(alternatives))*scenario$effect_log2
  mean_affected<-sample.int(G,max(2L,round(G*scenario$mean_affected_fraction)))
  gamma[mean_affected]<-rep(c(-1,1),length.out=length(mean_affected))*scenario$cage_mean_amplitude
  cage_index<-as.integer(d$cage);K<-nlevels(d$cage)
  cage_mean<-ave(d$genotype,d$cage,FUN=mean)-mean(d$genotype)
  u<-matrix(rnorm(G*K),G,K)*scenario$cage_log2_sd
  eta<-outer(beta,d$genotype)+outer(gamma,cage_mean)+u[,cage_index,drop=FALSE]
  # The lognormal correction keeps E[2^u]=1 across generated cage effects.
  mu<-outer(p$mean_count,calibration$depth_ratio)*2^eta*exp(-.5*(log(2)*scenario$cage_log2_sd)^2)
  y<-matrix(rnbinom(G*N,mu=as.vector(mu),size=rep(1/p$dispersion,N)),G,N)
  rownames(y)<-p$gene_symbol;colnames(y)<-d$animal_id
  list(counts=y,beta=beta,gamma=gamma,alternatives=alternatives)
}

performance_row <- function(tab,truth,gamma,scenario,replicate,model,kept,whole_status='completed') {
  available<-tab$status=='available' & is.finite(tab$p) & is.finite(tab$low) & is.finite(tab$high)
  alternative<-truth!=0;null<-!alternative
  q<-bh(tab$p);discover<-is.finite(q)&q<.05
  q_available<-rep(NA_real_,length(q));ok<-is.finite(tab$p)
  q_available[ok]<-p.adjust(tab$p[ok],'BH')
  R<-sum(discover);V<-sum(discover&null);TP<-sum(discover&alternative)
  n_available<-sum(available);usable<-n_available>0
  covered<-available & tab$low<=truth & tab$high>=truth
  ratio<-function(a,b)if(b>0)a/b else NA_real_
  data.frame(scenario=scenario$scenario,replicate=replicate,model=model,
    cage_log2_sd=scenario$cage_log2_sd,cage_mean_amplitude=scenario$cage_mean_amplitude,
    de_fraction=scenario$de_fraction,genes=length(truth),filtered=sum(!kept),available=n_available,
    null_genes=sum(null),alternative_genes=sum(alternative),discoveries=R,false_discoveries=V,true_discoveries=TP,
    FDP=if(usable)V/max(R,1) else NA_real_,FWER=if(usable)as.numeric(V>0) else NA_real_,
    power=if(usable)ratio(TP,sum(alternative)) else NA_real_,
    delivered_power=ratio(TP,sum(alternative)),
    coverage=ratio(sum(covered),n_available),delivered_coverage=sum(covered)/length(truth),
    null_coverage=ratio(sum(covered&null),sum(available&null)),
    alternative_coverage=ratio(sum(covered&alternative),sum(available&alternative)),
    availability=n_available/length(truth),fit_availability_among_filtered=ratio(n_available,sum(kept)),
    bias=if(usable)mean(tab$beta[available]-truth[available]) else NA_real_,
    RMSE=if(usable)sqrt(mean((tab$beta[available]-truth[available])^2)) else NA_real_,
    null_abs_bias_mean_affected=if(any(available&null&gamma!=0))mean(abs(tab$beta[available&null&gamma!=0])) else NA_real_,
    mean_width=if(usable)mean(tab$high[available]-tab$low[available]) else NA_real_,
    discoveries_available_family=sum(q_available<.05,na.rm=TRUE),
    FDP_available_family=if(usable)sum(q_available[null]<.05,na.rm=TRUE)/max(1,sum(q_available<.05,na.rm=TRUE)) else NA_real_,
    usable_run=usable,run_status=if(usable)whole_status else 'no_available_inference')
}

run_counts <- function(job,cfg) {
  cal<-readRDS(file.path(cfg$out,'cache','count_calibration.rds'))
  scenarios<-rd(file.path(cfg$out,'results','count_scenarios.csv'))
  scenario<-scenarios[scenarios$scenario==job$scenario,,drop=FALSE]
  need(nrow(scenario)==1,'Unknown count scenario')
  genes<-cal$parameters$gene_symbol;models<-c('baseline','cage_fixed','cage_correlation','dream_native','mundlak_native')
  results<-errors<-pair_rows<-checks<-retries<-unavailable_rows<-inventory<-list()
  for(rep in seq.int(job$start,job$end-1L)) {
    Sys.setenv(PAPER1_DIAGNOSTICS_DIR=file.path(job$output,'failed_raw_fits',paste0('rep_',rep)))
    start<-proc.time()[3]
    generated<-generate_counts(cal,scenario,rep,cfg$seed)
    g<-DGEList(generated$counts)
    x<-mouse_designs(cal$metadata,'linear')$base$x
    keep<-filterByExpr(g,design=x)
    message<-NULL
    ans<-tryCatch({
      need(sum(keep)>=100,'Fewer than 100 genes passed filtering')
      g<-g[keep,,keep.lib.sizes=FALSE];g<-calcNormFactors(g)
      fit_mouse_models(g,cal$metadata,'linear',paste0('count_s',job$scenario),1L,FALSE)
    },error=function(e){message<<-conditionMessage(e);NULL})
    if(!is.null(message))errors[[length(errors)+1]]<-data.frame(scenario=job$scenario,replicate=rep,model='all',message=message)
    if(!is.null(ans)&&nrow(ans$messages)) {
      z<-ans$messages;z$scenario<-job$scenario;z$replicate<-rep;errors[[length(errors)+1]]<-z
    }
    tag<-function(z) {
      if(is.null(z)||!nrow(z))return(NULL)
      z$scenario<-job$scenario;z$replicate<-rep;z
    }
    if(!is.null(ans)) {
      checks[[length(checks)+1]]<-tag(ans$diagnostics)
      retries[[length(retries)+1]]<-tag(ans$retries)
    }
    trial<-list()
    for(model in models) {
      tab<-unavailable(genes,model,paste0('count_s',job$scenario),'filtered_out')
      tab$status[keep]<-'fit_failed'
      if(!is.null(ans)) {
        part<-ans$table[ans$table$model==model,];ix<-match(part$gene_symbol,genes)
        need(!anyNA(ix)&&!anyDuplicated(ix),'Unexpected simulation fit gene IDs')
        common<-intersect(names(tab),names(part));tab[ix,common]<-part[,common]
      }
      row<-performance_row(tab,generated$beta,generated$gamma,scenario,rep,model,keep)
      missing<-tab[tab$status!='available',c('gene_symbol','model','spec','status'),drop=FALSE]
      unavailable_rows[[length(unavailable_rows)+1]]<-tag(missing)
      inventory[[length(inventory)+1]]<-data.frame(scenario=job$scenario,replicate=rep,model=model,
        planned=length(genes),filtered=sum(!keep),available=row$available,
        fitted_unavailable=sum(keep & tab$status!='available'),diagnostic_scope='full_gene_diagnostics')
      row$seconds_total_so_far<-proc.time()[3]-start
      results[[length(results)+1]]<-row;trial[[model]]<-row
      # Keep one complete generated-data example, not all transcriptomes for all repetitions.
      if(rep==0L) {
        tab$true_beta<-generated$beta;tab$true_gamma<-generated$gamma;tab$q_planned<-bh(tab$p)
        wr(tab,file.path(job$output,paste0('example_',model,'.csv.gz')))
      }
    }
    reference<-trial$dream_native
    for(model in setdiff(models,'dream_native')) {
      z<-trial[[model]]
      pair_rows[[length(pair_rows)+1]]<-data.frame(scenario=job$scenario,replicate=rep,
        model=model,reference='dream_native',FDP_difference=z$FDP-reference$FDP,
        power_difference=z$power-reference$power,coverage_difference=z$coverage-reference$coverage,
        availability_difference=z$availability-reference$availability,
        both_usable=z$usable_run&&reference$usable_run)
    }
    save_json(list(id=job$id,done=rep-job$start+1L,total=job$end-job$start,replicate=rep),
              file.path(job$output,'progress.json'))
    cat('Count dataset completed',job$scenario,rep,'seconds=',proc.time()[3]-start,'\n')
  }
  wr(bind_rows(results),file.path(job$output,'metrics.csv'))
  wr(bind_rows(pair_rows),file.path(job$output,'paired_method_differences.csv'))
  wr(bind_rows(errors),file.path(job$output,'messages.csv'))
  wr(bind_rows(checks),file.path(job$output,'diagnostics.csv'))
  wr(bind_rows(retries),file.path(job$output,'retry_attempts.csv'))
  wr(bind_rows(unavailable_rows),file.path(job$output,'unavailable_genes.csv'))
  wr(bind_rows(inventory),file.path(job$output,'diagnostic_inventory.csv'))
  examples<-list.files(job$output,pattern='^example_.*csv.gz$')
  raw_files<-list.files(file.path(job$output,'failed_raw_fits'),recursive=TRUE,full.names=FALSE)
  finish_job(job,c('metrics.csv','paired_method_differences.csv','messages.csv',examples,
    'diagnostics.csv','retry_attempts.csv','unavailable_genes.csv','diagnostic_inventory.csv',
    if(length(raw_files))file.path('failed_raw_fits',raw_files)))
}
