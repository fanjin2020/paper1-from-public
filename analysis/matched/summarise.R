read_result <- function(path) {
  con<-if(endsWith(path,'.gz'))gzfile(path,'rt') else file(path,'rt')
  header<-tryCatch(readLines(con,n=2,warn=FALSE),finally=close(con))
  if(identical(header,'""'))return(data.frame())
  rd(path)
}

read_tasks <- function(out,prefix,name) {
  paths<-list.dirs(file.path(out,'tasks'),recursive=FALSE,full.names=TRUE)
  paths<-paths[startsWith(basename(paths),prefix)]
  bind_rows(lapply(paths,function(p){f<-file.path(p,name);if(file.exists(f))read_result(f) else NULL}))
}

mean_mc <- function(values,binary=FALSE) {
  v<-values[is.finite(values)];n<-length(v)
  if(!n)return(c(mean=NA_real_,MCSE=NA_real_,low=NA_real_,high=NA_real_,n=0))
  mu<-mean(v);se<-if(n>1)sd(v)/sqrt(n) else NA_real_
  low<-mu-1.96*se;high<-mu+1.96*se
  if(binary) {
    need(all(v%in%c(0,1)),'Wilson interval requires Bernoulli dataset outcomes')
    z<-qnorm(.975);den<-1+z^2/n
    center<-(mu+z^2/(2*n))/den
    half<-z*sqrt(mu*(1-mu)/n+z^2/(4*n^2))/den
    low<-max(0,center-half);high<-min(1,center+half)
  }
  c(mean=mu,MCSE=se,low=low,high=high,n=n)
}

summarise_all <- function(cfg) {
  out<-cfg$out;dest<-file.path(out,'results')
  reference_variant<-'shared__ordinal_nuisance__ordinary'
  h<-read_tasks(out,'human_','fits.csv.gz')
  h$q_joint<-ave(h$p,interaction(h$variant,h$model,drop=TRUE),FUN=bh)
  hp<-read_tasks(out,'human_','paired.csv.gz')
  hp$q_joint<-ave(hp$p,hp$paired_variant,FUN=bh)
  hp$q_HC3_joint<-ave(hp$p_HC3,hp$paired_variant,FUN=bh)
  hd<-read_tasks(out,'human_','designs.csv')
  wr(h,file.path(dest,'human_sensitivity_fits.csv.gz'))
  wr(hp,file.path(dest,'human_sensitivity_paired_changes.csv.gz'))
  wr(hd,file.path(dest,'human_sensitivity_designs.csv'))
  base<-h[h$variant==reference_variant & h$model=='base' & h$contrast=='clinical_status' &
    is.finite(h$q_joint) & h$q_joint<.05,]
  record_key<-function(z)paste(z$spec,z$gene_symbol,sep='::')
  selected<-record_key(base)
  summaries<-lapply(split(h,interaction(h$variant,h$model,drop=TRUE)),function(z) {
    data.frame(variant=z$variant[1],model=z$model[1],planned=nrow(z),available=sum(is.finite(z$p)),
      significant=sum(z$q_joint<.05,na.rm=TRUE),clinical_significant=sum(z$q_joint<.05 & z$contrast=='clinical_status',na.rm=TRUE),
      original_subset_size=length(selected),original_subset_significant=sum(z$q_joint<.05 & record_key(z)%in%selected,na.rm=TRUE))
  })
  wr(bind_rows(summaries),file.path(dest,'human_sensitivity_discoveries.csv'))
  psummary<-lapply(split(hp,hp$paired_variant),function(z) {
    sub<-z[record_key(z)%in%selected,]
    data.frame(paired_variant=z$paired_variant[1],planned=nrow(z),available=sum(is.finite(z$p)),
      full_model_based=sum(z$q_joint<.05,na.rm=TRUE),full_HC3=sum(z$q_HC3_joint<.05,na.rm=TRUE),
      full_overlap=sum(z$q_joint<.05&z$q_HC3_joint<.05,na.rm=TRUE),subset_size=nrow(sub),
      subset_model_based=sum(sub$q_joint<.05,na.rm=TRUE),subset_HC3=sum(sub$q_HC3_joint<.05,na.rm=TRUE))
  })
  wr(bind_rows(psummary),file.path(dest,'human_paired_family_and_selected_subset.csv'))
  current<-hp[hp$paired_variant=='shared__ordinal_nuisance',]
  selected_table<-current[record_key(current)%in%selected,]
  for(arm in c('base','quality')) {
    x<-h[h$variant==reference_variant & h$model==arm,]
    j<-match(record_key(selected_table),record_key(x))
    for(k in c('beta','SE','p','q_joint','low','high'))selected_table[[paste0(arm,'_',k)]]<-x[[k]][j]
  }
  selected_table$SE_ratio<-selected_table$quality_SE/selected_table$base_SE
  selected_table$interpretation<-rep('A non-significant paired difference does not establish coefficient equality',nrow(selected_table))
  wr(selected_table,file.path(dest,'human_original_discovery_complete_diagnostics.csv'))
  m<-read_tasks(out,'mouse_','background_fits.csv.gz')
  mp<-read_tasks(out,'mouse_','panel_fits.csv.gz')
  mp$q_panel_6spec_planned<-ave(mp$p,mp$model,FUN=bh)
  wr(m,file.path(dest,'mouse_seven_arm_background_fits.csv.gz'))
  wr(mp,file.path(dest,'mouse_seven_arm_panel_fits.csv.gz'))
  wr(read_tasks(out,'mouse_','designs.csv'),file.path(dest,'mouse_augmented_designs.csv'))
  wr(read_tasks(out,'mouse_','diagnostics.csv'),file.path(dest,'mouse_mixed_model_diagnostics.csv'))
  wr(read_tasks(out,'mouse_','retry_attempts.csv'),file.path(dest,'mouse_optimizer_retries.csv'))
  wr(read_tasks(out,'mouse_','messages.csv'),file.path(dest,'mouse_model_messages.csv'))
  contrasts<-list()
  for(spec in unique(m$spec))for(pair in list(c('mundlak_shared','dream_shared'),c('mundlak_native','dream_native'),
      c('mundlak_shared','cage_fixed'),c('mundlak_native','cage_fixed'),c('cage_fixed','baseline'))) {
    a<-m[m$spec==spec&m$model==pair[1],];b<-m[m$spec==spec&m$model==pair[2],]
    i<-match(a$gene_symbol,b$gene_symbol)
    good<-a$status=='available' & b$status[i]=='available'
    contrasts[[length(contrasts)+1]]<-data.frame(spec=spec,model=pair[1],reference=pair[2],
      gene_symbol=a$gene_symbol,paired_available=good,beta=a$beta,beta_reference=b$beta[i],
      beta_difference=a$beta-b$beta[i],SE_ratio=a$SE/b$SE[i],
      direction_changed=good & sign(a$beta)!=sign(b$beta[i]),
      q=a$q_background_spec,q_reference=b$q_background_spec[i],
      discovery_changed=good & ((a$q_background_spec<.05)!=(b$q_background_spec[i]<.05)))
  }
  mc<-bind_rows(contrasts);wr(mc,file.path(dest,'mouse_matched_method_comparisons.csv.gz'))
  ms<-lapply(split(m,interaction(m$spec,m$model,drop=TRUE)),function(z) {
    data.frame(spec=z$spec[1],model=z$model[1],planned=nrow(z),available=sum(z$status=='available'),
      significant_background=sum(z$q_background_spec<.05,na.rm=TRUE),median_SE=median(z$SE[z$status=='available'],na.rm=TRUE))
  })
  wr(bind_rows(ms),file.path(dest,'mouse_method_summary.csv'))
  mc_summary<-lapply(split(mc,interaction(mc$spec,mc$model,mc$reference,drop=TRUE)),function(z) {
    keep<-z$paired_available & is.finite(z$SE_ratio)
    data.frame(spec=z$spec[1],model=z$model[1],reference=z$reference[1],paired_available=sum(keep),
      median_beta_difference=median(z$beta_difference[keep],na.rm=TRUE),
      median_absolute_beta_difference=median(abs(z$beta_difference[keep]),na.rm=TRUE),
      median_SE_ratio=median(z$SE_ratio[keep],na.rm=TRUE),direction_changes=sum(z$direction_changed,na.rm=TRUE),
      discovery_changes=sum(z$discovery_changed,na.rm=TRUE))
  })
  wr(bind_rows(mc_summary),file.path(dest,'mouse_matched_method_summary.csv'))
  # Task ID provides the spec when joining repeated gene IDs across six designs.
  checks<-lapply(list.dirs(file.path(out,'tasks'),recursive=FALSE),function(path) {
    file<-file.path(path,'weighted_centering_identity.csv.gz')
    if(!file.exists(file))return(NULL)
    z<-read_result(file)
    data.frame(spec=sub('^mouse_','',basename(path)),genes=nrow(z),
      max_abs_weighted_identity_error=if(nrow(z))max(abs(z$difference_weighted)) else NA_real_,
      median_abs_arithmetic_difference=if(nrow(z))median(abs(z$difference_arithmetic)) else NA_real_,
      max_abs_arithmetic_difference=if(nrow(z))max(abs(z$difference_arithmetic)) else NA_real_)
  })
  wr(bind_rows(checks),file.path(dest,'weighted_centering_identity_summary.csv'))
  metrics<-read_tasks(out,'counts_','metrics.csv')
  need(nrow(metrics)==cfg$reps*8*5,'Count simulation is incomplete')
  need(!anyDuplicated(metrics[,c('scenario','replicate','model')]),'Duplicate count simulation dataset')
  wr(metrics,file.path(dest,'count_pipeline_per_dataset.csv'))
  wr(read_tasks(out,'counts_','messages.csv'),file.path(dest,'count_pipeline_messages.csv'))
  wr(read_tasks(out,'counts_','diagnostics.csv'),file.path(dest,'count_pipeline_diagnostics.csv'))
  wr(read_tasks(out,'counts_','retry_attempts.csv'),file.path(dest,'count_pipeline_optimizer_retries.csv'))
  wr(read_tasks(out,'counts_','unavailable_genes.csv'),file.path(dest,'count_pipeline_unavailable_genes.csv'))
  wr(read_tasks(out,'counts_','diagnostic_inventory.csv'),file.path(dest,'count_pipeline_diagnostic_inventory.csv'))
  paired<-read_tasks(out,'counts_','paired_method_differences.csv')
  wr(paired,file.path(dest,'count_pipeline_paired_per_dataset.csv'))
  sums<-lapply(split(metrics,interaction(metrics$scenario,metrics$model,drop=TRUE)),function(z) {
    ans<-list(scenario=z$scenario[1],model=z$model[1],repetitions=nrow(z),usable_runs=sum(z$usable_run),
      cage_log2_sd=z$cage_log2_sd[1],cage_mean_amplitude=z$cage_mean_amplitude[1],de_fraction=z$de_fraction[1],
      genes=z$genes[1],whole_run_failures=sum(!z$usable_run),
      FDR_lower_bound_all_runs=sum(z$FDP,na.rm=TRUE)/nrow(z),
      FDR_upper_bound_all_runs=(sum(z$FDP,na.rm=TRUE)+sum(!z$usable_run))/nrow(z))
    for(k in c('FDP','FWER','power','delivered_power','coverage','delivered_coverage','null_coverage',
               'alternative_coverage','availability','fit_availability_among_filtered','bias','RMSE','mean_width',
               'FDP_available_family')) {
      binary<-k=='FWER'||(k%in%c('FDP','FDP_available_family')&&z$de_fraction[1]==0)
      result<-mean_mc(z[[k]],binary=binary)
      label<-if(k=='FDP')'FDR' else k
      for(n in names(result))ans[[paste(label,n,sep='_')]]<-unname(result[n])
      ans[[paste0(label,'_interval_method')]]<-if(binary)'Wilson_95_dataset_Bernoulli' else 'normal_95_dataset_MCSE'
    }
    as.data.frame(ans)
  })
  sums<-bind_rows(sums);wr(sums,file.path(dest,'count_pipeline_summary.csv'))
  paired_sums<-lapply(split(paired,interaction(paired$scenario,paired$model,drop=TRUE)),function(z) {
    ans<-list(scenario=z$scenario[1],model=z$model[1],reference=z$reference[1],paired_repetitions=sum(z$both_usable))
    for(k in c('FDP_difference','power_difference','coverage_difference','availability_difference')) {
      res<-mean_mc(z[[k]][z$both_usable])
      for(n in names(res))ans[[paste(k,n,sep='_')]]<-unname(res[n])
    }
    as.data.frame(ans)
  })
  wr(bind_rows(paired_sums),file.path(dest,'count_pipeline_paired_summary.csv'))
  examples<-file.path(dest,'examples');dir.create(examples,showWarnings=FALSE)
  for(path in list.dirs(file.path(out,'tasks'),recursive=FALSE)) {
    for(f in list.files(path,pattern='^example_.*csv.gz$',full.names=TRUE))
      file.copy(f,file.path(examples,paste(basename(path),basename(f),sep='__')),overwrite=TRUE)
  }
  palette<-c(baseline='#666666',cage_fixed='#D55E00',cage_correlation='#0072B2',dream_shared='#56B4E9',
             dream_native='#009E73',mundlak_shared='#CC79A7',mundlak_native='#882255')
  save_json(as.list(palette),file.path(dest,'model_palette.json'))
  plot_count_summary(sums,file.path(out,'figures'),palette)
  note<-c('Matched human sensitivity, seven mouse analysis arms, and count-to-BH validation use a common prepared input set.',
    'The count generator uses the whole-section common gene background and observed allocation.',
    'Latent log2 cage SD is not the observed count ICC. Known truth is the conditional library-offset log2 mean coefficient.',
    'Voom uses a log-CPM approximation and estimated normalization; coverage therefore evaluates the complete pipeline against that generating target.',
    'FDR is mean per-dataset FDP, including zero discoveries. Entirely failed analyses have missing FDP, explicit availability, and all-run FDR bounds.',
    'Both planned-gene BH and available-gene BH are reported separately. Primary BH retains all prespecified simulated genes.',
    'MCSE uses variation across independent generated datasets, not across genes.',
    'All-null FDP and FWER use Wilson intervals for dataset-level Bernoulli outcomes. Non-null FDP and paired differences retain dataset-level normal MC intervals.',
    'Mundlak uses arithmetic cage means in shared and native-weight empirical arms; weighted-WLS identity checks are diagnostics, not additional mixed-model claims.',
    'Categorical Braak changes only the nuisance adjustment for the conditional clinical-status contrast; the Braak-slope estimand remains ordinal.',
    'Pairwise human tests condition on estimated weights and do not validate an eBayes-moderated paired test.',
    'This count study assesses genotype-coefficient inference, not a complete count-to-paired-difference FDR guarantee.',
    'Do not describe non-significant coefficient differences as evidence of equality.')
  writeLines(note,file.path(dest,'INTERPRETATION.txt'))
  flags<-character()
  diag<-bind_rows(list(read_tasks(out,'mouse_','diagnostics.csv'),read_tasks(out,'counts_','diagnostics.csv')))
  if(nrow(diag)&&any(diag$status=='retry_merge_failed'))flags<-c(flags,'RETRY_MERGE_FAILURES_REQUIRE_REVIEW')
  if(any(!metrics$usable_run) || any(vapply(ms,function(z)z$available==0,logical(1))))
    flags<-c(flags,'WHOLE_MODEL_FAILURES_REQUIRE_REVIEW')
  status<-if(length(flags))paste(flags,collapse=';') else 'ANALYSES_COMPLETE'
  save_json(list(status=status,validation_only=cfg$smoke,original_discovery_subset_size=length(selected),
    human_variants=length(unique(h$variant)),mouse_models=unique(m$model),
    count_scenarios=8,repetitions_per_scenario=cfg$reps,simulation_genes=unique(metrics$genes),
    interpretation=note),file.path(dest,'analysis_summary.json'))
}

plot_count_summary <- function(s,folder,palette) {
  dir.create(folder,recursive=TRUE,showWarnings=FALSE)
  methods<-c('baseline','cage_fixed','cage_correlation','dream_native','mundlak_native')
  labels<-c('No cage','Cage FE','Common corr.','dream','Mundlak')
  for(metric in c('FDR','power','coverage','availability')) {
    grDevices::cairo_pdf(file.path(folder,paste0('Count_pipeline_',metric,'.pdf')),width=12,height=7.5,family='sans')
    par(mfrow=c(2,4),mar=c(5.5,4,3.5,1),mgp=c(2.5,.7,0))
    for(scenario in sort(unique(s$scenario))) {
      z<-s[s$scenario==scenario,];z<-z[match(methods,z$model),]
      y<-z[[paste0(metric,'_mean')]];lo<-z[[paste0(metric,'_low')]];hi<-z[[paste0(metric,'_high')]]
      upper<-if(metric=='FDR')min(1,max(.1,hi[is.finite(hi)])) else 1
      plot(seq_along(methods),y,ylim=c(0,upper),xlim=c(.5,5.5),xaxt='n',xlab='',ylab=metric,
        pch=19,col=palette[methods],main=sprintf('SD %.1f; mean %.1f\nDE %.0f%%',z$cage_log2_sd[1],
          z$cage_mean_amplitude[1],100*z$de_fraction[1]))
      axis(1,at=seq_along(methods),labels=labels,las=2,cex.axis=.8)
      ok<-is.finite(y)&is.finite(lo)&is.finite(hi)
      if(any(ok))arrows(which(ok),pmax(0,lo[ok]),which(ok),pmin(1,hi[ok]),angle=90,code=3,length=.04,col=palette[methods][ok])
      if(metric=='FDR')abline(h=.05,lty=2,col='black')
      if(metric=='coverage')abline(h=.95,lty=2,col='black')
    }
    dev.off()
  }
}
