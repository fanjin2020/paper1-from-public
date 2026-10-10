# 逐基因加权协方差、dream区间和失败记录的核查函数。

# 用加权QR分解独立计算目标系数及其未缩放SE，检查lmFit的直接系数。
exact_wls_target <- function(y,w,x,target) {
  j<-match(target,colnames(x));need(!is.na(j),'设计缺少目标列')
  need(ncol(y)==nrow(x)&&identical(dim(y),dim(w)),'加权核查维度不一致')
  ans<-lapply(seq_len(nrow(y)),function(i) {
    need(all(is.finite(w[i,]))&&all(w[i,]>0),'权重不是正有限数')
    q<-qr(x*sqrt(w[i,]));need(q$rank==ncol(x),'加权设计秩不足')
    b<-qr.coef(q,y[i,]*sqrt(w[i,]))
    inv<-chol2inv(qr.R(q));position<-match(j,q$pivot)
    c(beta=b[j],unscaled=sqrt(inv[position,position]))
  })
  z<-as.data.frame(do.call(rbind,ans));names(z)<-c('beta','unscaled')
  z$gene_symbol<-rownames(y);z
}

# 在完全相同的响应和权重上交叉比较参数化与收缩设置，不按结果选择设置。
audit_se <- function(v,x,spec,branch) {
  need(inherits(v,'EList'),'SE核查需要含表达矩阵和逐观测权重的limma EList对象')
  truth<-exact_wls_target(v$E,v$weights,x,'genotype')
  raw<-limma::lmFit(v,x);j<-match('genotype',colnames(x))
  need(max(abs(truth$beta-raw$coefficients[,j]),na.rm=TRUE)<1e-7,'直接系数与QR解不一致')
  need(max(abs(truth$unscaled/raw$stdev.unscaled[,j]-1),na.rm=TRUE)<1e-7,
       '直接系数SE与逐基因加权协方差不一致')
  # 两个组均值的差与0/1系数相同；只对非缩放SE保留原对比变换作诊断。
  z<-cbind(groupWT=1-x[,j],group5X=x[,j],x[,setdiff(seq_len(ncol(x)),c(1,j)),drop=FALSE])
  need(qr(cbind(x,z))$rank==ncol(x),'两种参数化的模型空间不一致')
  ct<-setNames(rep(0,ncol(z)),colnames(z));ct[1:2]<-c(-1,1)
  legacy<-limma::contrasts.fit(limma::lmFit(v,z),ct)
  # 在原组均值参数化上逐基因直接计算 c'(X'WX)^-1c，不使用公共相关矩阵近似。
  exact_contrast<-vapply(seq_len(nrow(v$E)),function(i) {
    q<-qr(z*sqrt(v$weights[i,]));inv<-chol2inv(qr.R(q));c0<-ct[q$pivot]
    sqrt(drop(crossprod(c0,inv%*%c0)))
  },numeric(1))
  recoded<-limma::contrastAsCoef(z,matrix(ct,ncol=1))
  exact_refit<-limma::lmFit(v,recoded$design)
  exact_refit_se<-exact_refit$stdev.unscaled[,recoded$coef[1]]
  need(max(abs(exact_contrast/truth$unscaled-1))<1e-7,
       '精确线性对比方差与直接目标系数不一致')
  need(max(abs(exact_refit_se/truth$unscaled-1))<1e-7,
       'contrastAsCoef精确重拟合未通过QR核查')
  need(max(abs(raw$coefficients[,j]-legacy$coefficients[,1]),na.rm=TRUE)<1e-7,
       '两种参数化的目标系数不一致')
  out<-data.frame(gene_symbol=truth$gene_symbol,spec=spec,branch=branch,
    beta_direct=raw$coefficients[,j],beta_contrast=legacy$coefficients[,1],
    unscaled_exact=truth$unscaled,unscaled_direct=raw$stdev.unscaled[,j],
    unscaled_contrast=legacy$stdev.unscaled[,1],
    unscaled_contrast_exact=exact_contrast,unscaled_contrast_refit=exact_refit_se)
  for(robust in c(FALSE,TRUE)) {
    d<-moderate(raw,robust);l<-moderate(legacy,robust)
    id<-match(out$gene_symbol,rownames(d$coefficients));il<-match(out$gene_symbol,rownames(l$coefficients))
    suffix<-if(robust)'robust' else 'ordinary'
    out[[paste0('SE_direct_',suffix)]]<-d$stdev.unscaled[id,j]*sqrt(d$s2.post[id])
    out[[paste0('SE_contrast_',suffix)]]<-l$stdev.unscaled[il,1]*sqrt(l$s2.post[il])
    out[[paste0('p_direct_',suffix)]]<-d$p.value[id,j]
    out[[paste0('p_contrast_',suffix)]]<-l$p.value[il,1]
  }
  # 原稿线性年龄用中心化列。单独记录未中心化近似，证明参数化影响不来自计数。
  if('age'%in%colnames(z)&&!is.null(attr(x,'raw_age'))) {
    uncentered<-z;uncentered[,'age']<-attr(x,'raw_age')
    ul<-limma::contrasts.fit(limma::lmFit(v,uncentered),ct)
    ul<-moderate(ul,FALSE)
    out$SE_contrast_uncentered_ordinary<-ul$stdev.unscaled[,1]*sqrt(ul$s2.post)
  } else out$SE_contrast_uncentered_ordinary<-NA_real_
  out
}

# 以官方topTable区间反推SE并核对p值，避免混用dream的内部t转换量。
dream_rows <- function(fit,genes,spec,model,robust=FALSE,excluded=character()) {
  messages<-c(capture.output(print(attr(fit,'errors'))),capture.output(print(attr(fit,'error.initial'))))
  messages<-messages[!messages%in%c('NULL','character(0)')]
  out<-data.frame(gene_symbol=genes,spec=spec,model=model,beta=NA_real_,SE=NA_real_,
    df=NA_real_,p=NA_real_,low=NA_real_,high=NA_real_,status='fit_failed')
  keep<-!rownames(fit$coefficients)%in%excluded
  need(any(keep),'所有dream基因都存在未解决的数值问题')
  # 不收缩已定位为数值不稳定的拟合；计划基因仍保留在输出和BH分母中。
  f<-variancePartition::eBayes(fit[keep,],robust=robust)
  f<-f[,'genotype']
  tab<-variancePartition::topTable(f,coef=1,number=Inf,sort.by='none',confint=.95)
  ix<-match(rownames(tab),genes);need(!anyNA(ix)&&!anyDuplicated(ix),'dream返回未知或重复基因')
  df<-f$df.total;if(is.matrix(df))df<-df[,1]
  df<-rep_len(df,nrow(f$coefficients))[match(rownames(tab),rownames(f$coefficients))]
  se<-(tab$CI.R-tab$CI.L)/(2*qt(.975,df))
  ok<-is.finite(tab$logFC)&is.finite(se)&se>0&is.finite(df)&df>0&is.finite(tab$P.Value)
  if(any(ok))need(max(abs(2*pt(-abs(tab$logFC[ok]/se[ok]),df[ok])-tab$P.Value[ok]))<1e-7,
                 'dream官方区间与p值不一致')
  out[ix,c('beta','SE','df','p','low','high')]<-cbind(tab$logFC,se,df,tab$P.Value,tab$CI.L,tab$CI.R)
  out$status[ix]<-ifelse(ok,'available','fit_failed')
  out[out$status!='available',c('p','low','high')]<-NA_real_
  out$status[out$gene_symbol%in%excluded]<-'numerical_warning'
  list(table=out,messages=messages)
}

# 方差为零的边界解另行报告；负特征值或未收敛警告需要定位到基因。
numerical_dream_messages <- function(messages) {
  messages[grepl('converg|eigen|Hessian|positive definite|收敛|特征',messages,ignore.case=TRUE) &
             !grepl('boundary.*singular',messages,ignore.case=TRUE)]
}

# 先按128基因复查，再逐个检查有警告的块。没有复现的警告单独保留，不能直接称已解决。
locate_dream_warnings <- function(genes,fit_one,chunk_size=128L) {
  probe<-function(ids) {
    messages<-character();error<-FALSE
    tryCatch(withCallingHandlers(fit_one(ids),warning=function(w) {
      messages<<-c(messages,conditionMessage(w));invokeRestart('muffleWarning')
    }),error=function(e){error<<-TRUE;messages<<-c(messages,conditionMessage(e))})
    list(messages=numerical_dream_messages(messages),error=error,all=messages)
  }
  rows<-list()
  for(start in seq.int(1L,length(genes),by=chunk_size)) {
    ids<-genes[start:min(length(genes),start+chunk_size-1L)];block<-probe(ids)
    if(length(block$messages)||block$error) {
      n_before<-length(rows)
      for(gene in ids) {
        z<-probe(gene)
        if(length(z$messages)||z$error)rows[[length(rows)+1]]<-data.frame(gene_symbol=gene,
          status=if(z$error)'fit_error' else 'numerical_warning',message=paste(z$all,collapse=' | '))
      }
      if(length(rows)==n_before)rows[[length(rows)+1]]<-data.frame(gene_symbol=NA_character_,
        status='block_warning_unresolved',message=paste(c(paste(ids,collapse=';'),block$all),collapse=' | '))
    }
    if(start%%1024L==1L)cat('dream警告定位',min(length(genes),start+chunk_size-1L),'/',length(genes),'\n')
  }
  if(!length(rows))return(data.frame(gene_symbol=NA_character_,status='warning_not_reproduced',
      message='原运行有数值警告，分块复查未定位；仍需人工检查，未认证收敛'))
  do.call(rbind,rows)
}

# 收集官方接口的警告和错误；写入结果表后才拦截终端中的重复提示。
capture_dream_call <- function(fun) {
  messages<-character()
  value<-tryCatch(withCallingHandlers(fun(),warning=function(w) {
    messages<<-c(messages,conditionMessage(w));invokeRestart('muffleWarning')
  }),error=function(e)e)
  failed<-inherits(value,'error')
  if(failed)messages<-c(messages,conditionMessage(value))
  list(fit=if(failed)NULL else value,messages=unique(messages),failed=failed)
}

# 两个优化器和停止精度在查看重拟合结果前确定，不比较p值后再选算法。
dream_retry_controls <- function() {
  list(bobyqa=lme4::lmerControl(optimizer='bobyqa',calc.derivs=TRUE,
        optCtrl=list(maxfun=200000L,rhoend=1e-9)),
       nloptwrap=lme4::lmerControl(optimizer='nloptwrap',calc.derivs=TRUE,
        optCtrl=list(maxeval=200000L,xtol_abs=1e-10,ftol_abs=1e-10)))
}

# 比较未经经验贝叶斯收缩的目标量；包括方差和似然，不能只看系数相近。
dream_raw_metrics <- function(fit,gene) {
  i<-match(gene,rownames(fit$coefficients));j<-match('genotype',colnames(fit$coefficients))
  need(!is.na(i)&&!is.na(j),'重拟合没有返回目标基因或目标系数')
  v<-attr(fit,'varComp');need(!is.null(v)&&'cage'%in%colnames(v),'dream缺少笼位方差')
  k<-match(gene,rownames(v));need(!is.na(k),'方差组分没有对应基因')
  df<-if(is.matrix(fit$df.residual))fit$df.residual[i,j] else fit$df.residual[i]
  c(beta=unname(fit$coefficients[i,j]),
    SE_raw=unname(fit$stdev.unscaled[i,j]*fit$sigma[i]),df_raw=unname(df),
    logLik=unname(fit$logLik[i]),sigma=unname(fit$sigma[i]),
    cage_variance=unname(v[k,'cage']))
}

# 这是数值一致性容差，不是生物学等效界值；任何非有限结果均不能通过。
dream_metrics_agree <- function(a,b) {
  if(!identical(names(a),names(b))||any(!is.finite(c(a,b))))return(FALSE)
  tol<-c(beta=1e-5,SE_raw=1e-4,df_raw=1e-3,logLik=1e-6,sigma=1e-4,cage_variance=1e-4)
  if(!identical(names(a),names(tol)))return(FALSE)
  if(any(c(a['SE_raw'],b['SE_raw'],a['df_raw'],b['df_raw'],a['sigma'],b['sigma'])<=0))return(FALSE)
  if(any(c(a['cage_variance'],b['cage_variance'])<0))return(FALSE)
  all(abs(a-b)<=tol*(1+pmax(abs(a),abs(b))))
}

# 已定位的基因分别调用两个官方dream拟合；零方差本身不判失败，数值警告仍判失败。
probe_dream_retry <- function(gene,call_fit,controls) {
  fits<-list()
  probes<-lapply(names(controls),function(name) {
    z<-capture_dream_call(function()call_fit(gene,controls[[name]]))
    fits[[name]]<<-z$fit
    metrics<-setNames(rep(NA_real_,6),c('beta','SE_raw','df_raw','logLik','sigma','cage_variance'))
    if(!z$failed) {
      m<-tryCatch(dream_raw_metrics(z$fit,gene),error=function(e)e)
      if(inherits(m,'error'))z$messages<-c(z$messages,conditionMessage(m)) else metrics<-m
    }
    errors<-if(is.null(z$fit))NULL else attr(z$fit,'errors')
    initial<-if(is.null(z$fit))NULL else attr(z$fit,'error.initial')
    clean<-!z$failed&&!length(z$messages)&&!length(errors)&&!length(initial)&&all(is.finite(metrics))
    data.frame(gene_symbol=gene,optimizer=name,clean=clean,
      as.list(metrics),message=paste(c(z$messages,unlist(errors),unlist(initial)),collapse=' | '))
  })
  tab<-do.call(rbind,probes);cols<-c('beta','SE_raw','df_raw','logLik','sigma','cage_variance')
  one<-setNames(as.numeric(unlist(tab[1,cols],use.names=FALSE)),cols)
  two<-setNames(as.numeric(unlist(tab[2,cols],use.names=FALSE)),cols)
  accepted<-all(tab$clean)&&dream_metrics_agree(one,two)
  tab$accepted<-accepted
  list(table=tab,accepted=accepted,metrics=one,fit=if(accepted)fits[[1]] else NULL)
}

# 只替换已通过数值复核的原始拟合行，尚未进行经验贝叶斯收缩。
# dream没有公开的MArrayLM2行合并接口，因此逐项校验字段、维度和行名；
# 遇到未知结构就停止，不能把缺失的协方差或自由度用默认值补上。
replace_dream_raw_row <- function(background,retry_fit,gene) {
  need(inherits(background,'MArrayLM2')&&inherits(retry_fit,'MArrayLM2'),
       '原始dream合并需要两个MArrayLM2对象')
  forbidden<-c('s2.post','s2.prior','df.prior','var.prior','lods')
  need(!any(vapply(forbidden,function(k)!is.null(background[[k]])||!is.null(retry_fit[[k]]),logical(1))),
       '不能合并已经经验贝叶斯收缩的dream结果')
  ids<-rownames(background$coefficients);i<-match(gene,ids)
  need(!is.na(i)&&!anyDuplicated(ids)&&identical(rownames(retry_fit$coefficients),gene),
       '重试基因必须唯一且已存在于背景；缺失行仍保留为不可用')
  need(identical(names(background),names(retry_fit)),'dream原始拟合字段不一致')
  matrices<-c('coefficients','stdev.unscaled','df.residual','df.total','t','p.value',
              'z.std','residuals','hatvalues')
  vectors<-c('rdf','Amean','sigma','logLik','edf','F','F.p.value')
  shared<-c('design','contrasts','cov.coefficients','method','formula','data')
  known<-c(matrices,vectors,shared,'genes','cov.coefficients.list')
  need(all(names(background)%in%known),paste('未知dream字段：',paste(setdiff(names(background),known),collapse=', ')))
  out<-background
  for(k in names(background)) {
    a<-background[[k]];b<-retry_fit[[k]]
    if(is.null(a)&&is.null(b))next
    if(k%in%matrices||k=='genes') {
      need(length(dim(a))==2L&&length(dim(b))==2L&&nrow(a)==length(ids)&&nrow(b)==1L&&
           ncol(a)==ncol(b)&&identical(colnames(a),colnames(b)),paste('dream矩阵维度不一致：',k))
      need(identical(rownames(a),ids)&&identical(rownames(b),gene),paste('dream矩阵行名不一致：',k))
      out[[k]][i,]<-b[1,]
    } else if(k%in%vectors) {
      need(is.atomic(a)&&is.atomic(b)&&length(a)==length(ids)&&length(b)==1L,
           paste('dream逐基因向量长度不一致：',k))
      out[[k]][i]<-b[1]
    } else if(k=='cov.coefficients.list') {
      need(is.list(a)&&is.list(b)&&identical(names(a),ids)&&identical(names(b),gene),
           'dream逐基因协方差未按基因对应')
      need(identical(dimnames(a[[i]]),dimnames(b[[1]])),'dream协方差系数顺序不一致')
      out[[k]][[i]]<-b[[1]]
    } else need(identical(a,b),paste('dream共同设计字段发生改变：',k))
  }
  for(k in c('varComp','edf')) {
    a<-attr(background,k);b<-attr(retry_fit,k)
    need(!is.null(a)&&!is.null(b),paste('dream缺少原始属性：',k))
    if(k=='varComp') {
      need(identical(rownames(a),ids)&&identical(rownames(b),gene)&&identical(colnames(a),colnames(b)),
           'dream方差组分标识不一致')
      a[i,]<-b[1,]
    } else {
      need(length(a)==length(ids)&&length(b)==1L,'dream有效自由度长度不一致')
      a[i]<-b[1]
    }
    attr(out,k)<-a
  }
  errors<-attr(out,'errors')
  if(length(errors)) {
    need(!is.null(names(errors)),'dream错误属性缺少基因标识')
    attr(out,'errors')<-errors[names(errors)!=gene]
  }
  need(!length(attr(out,'error.initial')),'dream存在未解决的整体初始拟合错误')
  need(dream_metrics_agree(dream_raw_metrics(out,gene),dream_raw_metrics(retry_fit,gene)),
       '合并后原始估计与重试结果不一致')
  out
}

# 定位捕获的数值警告，并补入dream通过属性返回的基因级错误。
dream_problem_rows <- function(raw,genes,call_fit,control=NULL) {
  rows<-data.frame(gene_symbol=character(),status=character(),message=character())
  if(length(numerical_dream_messages(raw$messages)))
    rows<-locate_dream_warnings(genes,function(ids)call_fit(ids,control))
  errors<-attr(raw$fit,'errors')
  if(length(errors)) {
    for(k in seq_along(errors)) {
      gene<-names(errors)[k]
      if(is.null(gene)||is.na(gene)||!gene%in%genes)gene<-NA_character_
      rows<-rbind(rows,data.frame(gene_symbol=gene,status='fit_error',message=as.character(errors[[k]])))
    }
  }
  unique(rows)
}

# 保留正常基因的原始拟合，只替换通过复核的问题行，再对可用背景共同收缩。
checked_dream <- function(v,formula,d,genes,spec,model,param,robust=FALSE,
                          known_excluded=character(),retry=TRUE) {
  diagnostics<-data.frame(gene_symbol=character(),status=character(),message=character())
  attempts<-data.frame(gene_symbol=character(),optimizer=character(),clean=logical(),
    beta=numeric(),SE_raw=numeric(),df_raw=numeric(),logLik=numeric(),sigma=numeric(),
    cage_variance=numeric(),message=character(),accepted=logical())
  impact<-data.frame(gene_symbol=character(),spec=character(),model=character(),
    beta_with=numeric(),beta_without=numeric(),SE_with=numeric(),SE_without=numeric(),
    p_with=numeric(),p_without=numeric())
  all_messages<-character();recovered<-character();policy<-'default';default_table<-NULL
  call_fit<-function(ids,control=NULL) {
    args<-list(exprObj=v[ids,,drop=FALSE],formula=formula,data=d,REML=TRUE,
      ddf='Satterthwaite',BPPARAM=if(length(ids)==length(genes))param else BiocParallel::SerialParam(),
      hideErrorsInBackend=FALSE)
    if(!is.null(control))args$control<-control
    do.call(variancePartition::dream,args)
  }
  ans<-tryCatch({
    raw<-capture_dream_call(function()call_fit(genes))
    all_messages<-raw$messages;need(!raw$failed,paste(raw$messages,collapse=' | '))
    diagnostics<-dream_problem_rows(raw,genes,call_fit)
    extra<-setdiff(known_excluded,diagnostics$gene_symbol)
    if(length(extra))diagnostics<-rbind(diagnostics,data.frame(gene_symbol=extra,
      status='prior_numerical_warning',message='来源运行已记录数值警告，本轮必须复核'))
    bad<-unique(na.omit(diagnostics$gene_symbol));f<-raw$fit;excluded<-bad
    default_table<-dream_rows(f,genes,spec,model,robust,excluded)$table
    if(retry&&length(bad)&&!anyNA(diagnostics$gene_symbol)) {
      controls<-dream_retry_controls();probes<-list()
      for(gene in bad) {
        cat('dream定点复核',spec,model,gene,'\n')
        probes[[gene]]<-probe_dream_retry(gene,call_fit,controls)
        attempts<-rbind(attempts,probes[[gene]]$table)
      }
      accepted<-names(probes)[vapply(probes,function(x)x$accepted,logical(1))]
      if(length(accepted)) {
        for(gene in accepted) {
          merged<-tryCatch(replace_dream_raw_row(f,probes[[gene]]$fit,gene),error=function(e)e)
          if(inherits(merged,'error')) {
            diagnostics<-rbind(diagnostics,data.frame(gene_symbol=gene,status='retry_merge_failed',
              message=conditionMessage(merged)))
          } else {
            f<-merged;recovered<-c(recovered,gene)
            diagnostics$status[diagnostics$gene_symbol==gene]<-'recovered_optimizer_agreement'
          }
        }
        excluded<-setdiff(bad,recovered)
        if(length(recovered))policy<-'targeted_raw_retry_before_joint_moderation'
        cat('dream定点复核恢复',length(recovered),'/',length(bad),'；正常背景不重新拟合\n')
      }
    }
    moderation<-capture_dream_call(function()dream_rows(f,genes,spec,model,robust,excluded))
    all_messages<-unique(c(all_messages,moderation$messages))
    need(!moderation$failed,paste(moderation$messages,collapse=' | '))
    if(length(numerical_dream_messages(moderation$messages)))
      diagnostics<-rbind(diagnostics,data.frame(gene_symbol=NA_character_,
        status='postfit_warning_unresolved',message=paste(moderation$messages,collapse=' | ')))
    result<-moderation$fit
    if(length(recovered)) {
      without<-dream_rows(f,genes,spec,model,robust,union(excluded,recovered))$table
      with<-result$table
      impact<-data.frame(gene_symbol=genes,spec=spec,model=model,
        beta_with=with$beta,beta_without=without$beta,SE_with=with$SE,SE_without=without$SE,
        p_with=with$p,p_without=without$p)
    }
    result
  },error=function(e)list(table=data.frame(gene_symbol=genes,spec=spec,model=model,
    beta=NA_real_,SE=NA_real_,df=NA_real_,p=NA_real_,low=NA_real_,high=NA_real_,status='fit_failed'),
    messages=conditionMessage(e)))
  ans$messages<-unique(c(ans$messages,all_messages))
  diagnostics$spec<-rep(spec,nrow(diagnostics));diagnostics$model<-rep(model,nrow(diagnostics))
  attempts$spec<-rep(spec,nrow(attempts));attempts$model<-rep(model,nrow(attempts))
  ans$diagnostics<-diagnostics;ans$retry_attempts<-attempts;ans$background_impact<-impact
  ans$optimizer_policy<-policy;ans$default_table<-default_table;ans
}

# 两套权重分别提取方差。零方差的收敛解可用，非收敛解另记为失败。
checked_variance <- function(v,formula,d,genes,spec,policy,param) {
  warnings<-character()
  vc<-tryCatch(withCallingHandlers(variancePartition::fitVarPartModel(v,formula,d,
    REML=TRUE,fxn=variance_extract,BPPARAM=param),warning=function(w){
      warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')}),error=function(e)e)
  rows<-lapply(genes,function(gene) {
    value<-if(inherits(vc,'error'))vc else vc[[gene]]
    if(is.data.frame(value)) {
      cbind(gene_symbol=gene,spec=spec,weight_policy=policy,
        status=if(isTRUE(value$converged))'available' else 'failed',value)
    } else data.frame(gene_symbol=gene,spec=spec,weight_policy=policy,status='failed',
      cage_variance=NA_real_,residual_scale=NA_real_,mean_residual_variance=NA_real_,ICC=NA_real_,
      beta=NA_real_,mean_within_cage_correlation=NA_real_,
      at_zero=NA,singular=NA,converged=FALSE,
      message=if(inherits(value,'error'))conditionMessage(value) else paste(c('未返回有效方差记录',value),collapse=' | '))
  })
  list(table=do.call(rbind,rows),messages=unique(warnings))
}
