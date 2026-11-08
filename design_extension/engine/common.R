# 所有分析共用的输入检查、输出和统计函数；不包含任何历史结果。
options(stringsAsFactors=FALSE, digits=17, warn=1)
Sys.setenv(R_LIBS=paste(.libPaths(),collapse=.Platform$path.sep),
           R_LIBS_USER=paste(.libPaths(),collapse=.Platform$path.sep))

# 条件不成立时保留明确错误，不改变筛选门槛以继续计算。
need <- function(x, message) if (!isTRUE(x)) stop(message, call.=FALSE)
# 数值检查允许整数/双精度存储差异，但不允许改变维度、长度、数值或引入缺失。
# 标识符由调用处另外严格比较，不在这里忽略样本或基因的对应关系。
same_numeric_values <- function(x,y) {
  is.numeric(x) && is.numeric(y) && identical(dim(x),dim(y)) &&
    length(x)==length(y) && all(is.finite(x)) && all(is.finite(y)) &&
    all(as.numeric(x)==as.numeric(y))
}
# 保留原始标识符，避免R将连字符、数字样本名自动改写。
rd <- function(path, ...) read.csv(path, check.names=FALSE, stringsAsFactors=FALSE, ...)
# 使用临时文件发布完整表格，失败时不会把半张表当作结果。
wr <- function(x, path) {
  dir.create(dirname(path), recursive=TRUE, showWarnings=FALSE)
  temp <- paste0(path, '.tmp')
  con <- if (endsWith(path, '.gz')) gzfile(temp, 'wt') else file(temp, 'wt')
  tryCatch(write.csv(x, con, row.names=FALSE, na=''), finally=close(con))
  need(file.rename(temp,path), paste('不能保存',path))
}
# 缺失检验占据计划位置但仍输出NA，分母不会因失败而缩小。
bh <- function(p) {
  q <- p.adjust(ifelse(is.finite(p), p, 1), 'BH')
  q[!is.finite(p)] <- NA_real_; q
}
# 按计数矩阵读取整数输入；表达缺失不能被填成零。
counts <- function(path) {
  tab <- rd(path); x <- as.matrix(tab[,-1,drop=FALSE]); rownames(x) <- tab[[1]]
  need(!anyDuplicated(rownames(x)) && !anyDuplicated(colnames(x)), '计数行列编号重复')
  need(is.numeric(x) && all(is.finite(x)) && all(x>=0) && all(x==round(x)), '输入不是完整的非负整数计数')
  x
}
# 单进程为默认选择；指定多进程时显式使用SOCK，并在预检查中验证snow。
backend <- function(workers) {
  if (workers==1) BiocParallel::SerialParam() else
    BiocParallel::SnowParam(workers=workers,type='SOCK',progressbar=TRUE,stop.on.error=TRUE)
}
# 一个模型必须有足够残差自由度，不自动删除混杂协变量。
estimable <- function(x) qr(x)$rank==ncol(x) && nrow(x)-ncol(x)>=5
# 构造已有协变量空间；只有全组常数的性别或年龄不会产生冗余列。
human_design <- function(d, contrast, quality=FALSE) {
  target <- if (contrast=='braak_same_clinical_population') d$pathology else as.numeric(d$clinical_group=='Dementia')
  x <- cbind(Intercept=1, target=target)
  if (length(unique(d$age))>1) x <- cbind(x, age=d$age-mean(d$age))
  if (length(unique(d$sex))>1) x <- cbind(x, male=as.numeric(d$sex=='male'))
  if (contrast=='clinical_status_adjusted_for_braak') x <- cbind(x, Braak=d$pathology)
  if (quality) x <- cbind(x, PMI=d$PMI, RIN=d$RIN, pH=d$pH)
  x
}
# 同一响应的WLS差值使用两个影响向量的差，保留估计间协方差。
paired_wls <- function(y,w,x0,x1,j=2L) {
  need(estimable(x0)&&estimable(x1), '配对设计不可估计')
  need(qr(cbind(x1,x0))$rank==ncol(x1), '模型并非嵌套')
  a0 <- solve(crossprod(x0,w*x0),t(x0)*rep(w,each=ncol(x0)))
  a1 <- solve(crossprod(x1,w*x1),t(x1)*rep(w,each=ncol(x1)))
  d <- a1[j,]-a0[j,]; delta <- sum(d*y); df <- length(y)-ncol(x1)
  e <- as.numeric(y-x1%*%(a1%*%y)); h <- rowSums(x1*t(a1))
  unit <- sum(d^2/w); se <- sqrt(sum(w*e^2)/df*unit)
  hc3 <- if(max(h)<1-1e-10) sqrt(sum(d^2*(e/(1-h))^2)) else NA_real_
  valid <- is.finite(se)&&se>1e-12&&unit>1e-14
  data.frame(delta=delta, SE_delta=if(valid)se else NA_real_, df=df,
    p=if(valid)2*pt(-abs(delta/se),df) else NA_real_,
    low=if(valid)delta-qt(.975,df)*se else NA_real_,
    high=if(valid)delta+qt(.975,df)*se else NA_real_,
    SE_HC3=hc3, p_HC3=if(is.finite(hc3)&&hc3>1e-12)2*pt(-abs(delta/hc3),df) else NA_real_,
    max_leverage=max(h), status=if(valid)'available' else 'degenerate_difference')
}
# 提取limma或dream经验贝叶斯统计量，直接使用其系数对应的自由度。
fit_rows <- function(fit, genes, target, spec, model) {
  j <- match(target,colnames(fit$coefficients)); need(!is.na(j),'目标系数缺失')
  ix <- match(genes,rownames(fit$coefficients)); good <- !is.na(ix)
  answer <- data.frame(gene_symbol=genes,spec=spec,model=model,beta=NA_real_,SE=NA_real_,
    df=NA_real_,p=NA_real_,low=NA_real_,high=NA_real_,status='not_tested')
  if(any(good)) {
    k <- ix[good]; b <- fit$coefficients[k,j]
    se <- fit$stdev.unscaled[k,j]*sqrt(fit$s2.post[k])
    df <- if(is.matrix(fit$df.total)) fit$df.total[k,j] else rep_len(fit$df.total,nrow(fit$coefficients))[k]
    p <- fit$p.value[k,j]; valid <- is.finite(b)&is.finite(se)&se>0&is.finite(df)&df>0&is.finite(p)
    answer[good,c('beta','SE','df','p','low','high')] <- cbind(b,se,df,p,b-qt(.975,df)*se,b+qt(.975,df)*se)
    answer$status[good] <- ifelse(valid,'available','fit_failed')
    answer[answer$status!='available',c('p','low','high')] <- NA_real_
  }
  answer
}
# 收缩方式由参数文件明确指定；没有参数时使用普通经验贝叶斯。
moderate <- function(fit, robust=FALSE) {
  good <- is.finite(fit$sigma)&fit$sigma>0
  need(sum(good)>=100,'有效背景不足100基因')
  limma::eBayes(fit[good,],trend=FALSE,robust=robust)
}
# 提取随机截距方差；加权残差的平均相关与未加权ICC分别报告。
variance_extract <- function(fit) {
  v <- as.data.frame(lme4::VarCorr(fit)); tau <- v$vcov[v$grp=='cage'][1]
  s2 <- sigma(fit)^2; w <- weights(fit); if(is.null(w))w <- rep(1,nobs(fit))
  msg <- fit@optinfo$conv$lme4$messages; opt <- fit@optinfo$conv$opt
  nonboundary <- msg[!grepl('boundary.*singular',msg)]
  cg<-as.character(lme4::getME(fit,'flist')$cage)
  ij<-which(outer(cg,cg,'==')&upper.tri(matrix(FALSE,length(cg),length(cg))),arr.ind=TRUE)
  pair<-if(nrow(ij))tau/sqrt((tau+s2/w[ij[,1]])*(tau+s2/w[ij[,2]])) else numeric()
  data.frame(cage_variance=tau,residual_scale=s2,mean_residual_variance=mean(s2/w),
    beta=unname(lme4::fixef(fit)['genotype']),
    mean_within_cage_correlation=if(length(pair))mean(pair) else NA_real_,
    ICC=tau/(tau+mean(s2/w)),at_zero=tau<=1e-12,singular=lme4::isSingular(fit),
    converged=(is.null(opt)||all(opt==0))&&length(nonboundary)==0,
    message=paste(msg,collapse=' | '))
}
