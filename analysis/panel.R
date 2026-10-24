# 芯片面板从官方表达矩阵重算，阈值不依赖SEA-AD或小鼠结果。
a<-commandArgs(TRUE);work<-a[1];.libPaths(c(a[2],.libPaths()));root<-a[3]
source(file.path(root,'analysis','common.R'))
suppressPackageStartupMessages({library(limma);library(jsonlite)})
cfg<-fromJSON(file.path(work,'run_config.json'))
tab<-rd(file.path(work,'prepared','discovery_expression.csv.gz'))
x<-as.matrix(tab[,-1]);rownames(x)<-tab[[1]]
md<-rd(file.path(work,'prepared','discovery_metadata.csv'))
need(identical(colnames(x),md$unit_id),'发现数据表达与单位不匹配')
all<-list()
for(region in c('frontal cortex','entorhinal cortex')) {
  keep<-md$region==region & md$group%in%c('ad','control') & md$sex%in%c('male','female') & is.finite(md$age)
  d<-md[keep,,drop=FALSE];d$group<-factor(d$group,levels=c('control','ad'))
  need(all(table(d$group)>=10),'发现数据完整病例每组不足10')
  d$sex<-factor(d$sex,levels=c('female','male'));d$age<-d$age-mean(d$age)
  design<-model.matrix(~group+age+sex,d);need(estimable(design),'发现模型不可估计')
  fit<-eBayes(lmFit(x[,keep,drop=FALSE],design))
  z<-topTable(fit,coef='groupad',number=Inf,sort.by='none',adjust.method='BH')
  all[[region]]<-data.frame(gene_symbol=rownames(z),region=region,beta=z$logFC,p=z$P.Value,q=z$adj.P.Val,n=nrow(d))
}
answer<-do.call(rbind,all)
wr(answer,file.path(work,'results','panel_discovery.csv.gz'))
panel<-subset(answer,region=='frontal cortex' & q<cfg$panel_fdr)
need(nrow(panel)>0,'额叶面板为空；不使用旧面板替代')
wr(panel,file.path(work,'results','panel.csv'))
cat('本次额叶面板',nrow(panel),'基因；未强制为237\n')
