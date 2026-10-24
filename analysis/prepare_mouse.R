# 从公开Seurat对象的原始Spatial counts按作者SAMPLE编号聚合。
a<-commandArgs(TRUE); work<-a[1]; .libPaths(c(a[2],.libPaths())); root<-a[3]
source(file.path(root,'analysis','common.R'))
suppressPackageStartupMessages({library(SeuratObject);library(Matrix);library(jsonlite)})
cfg<-fromJSON(file.path(work,'run_config.json'))
md<-rd(file.path(work,'prepared','mouse_manifest.csv'),colClasses=c(SAMPLE='character',animal_id='character',cage='character'))
object_path<-readLines(file.path(work,'prepared','mouse_object_path.txt'),warn=FALSE)[1]
o<-readRDS(object_path)
need(inherits(o,'Seurat'),'公开对象不是Seurat')
meta<-o[[]]
need(all(c('SAMPLE','annotation','Condition','Age','Sex')%in%names(meta)),'Seurat身份字段缺失')
need(!anyDuplicated(rownames(meta)),'spot条码重复')
layers<-SeuratObject::Layers(o[['Spatial']])
need('counts'%in%layers,'Spatial没有唯一counts层；不能改用data层')
x<-SeuratObject::GetAssayData(o,assay='Spatial',layer='counts')
need(identical(colnames(x),rownames(meta)),'counts列与元数据条码不一致')
need(all(is.finite(x@x))&&all(x@x>=0)&&all(x@x==round(x@x)),'Spatial不是原始整数计数')
sample<-as.character(meta$SAMPLE); ix<-match(sample,md$SAMPLE)
eligible<-!is.na(ix); needed<-md$SAMPLE
need(all(needed%in%sample),'发表动物中有对象未覆盖样本')
for(i in seq_len(nrow(md))) {
  z<-meta[sample==md$SAMPLE[i],,drop=FALSE]
  need(all(toupper(z$Sex)==md$sex[i]),'Seurat与GEO性别冲突')
  need(all(as.character(z$Condition)==md$genotype[i]),'Seurat与GEO基因型冲突')
  age<-as.numeric(sub('mo$','',as.character(z$Age)))
  need(all(age==md$age[i]),'Seurat与GEO年龄冲突')
}
cortex<-as.character(meta$annotation)%in%cfg$cortex_labels
whole_n<-tabulate(ix[eligible],nbins=nrow(md))
cortex_n<-tabulate(ix[eligible&cortex],nbins=nrow(md))
md$whole_spots<-whole_n;md$cortex_spots<-cortex_n
md$included<-whole_n>0 & cortex_n>=cfg$minimum_cortex_spots
wr(md,file.path(work,'prepared','mouse_inclusion.csv'))
md<-md[md$included,,drop=FALSE];need(nrow(md)>10,'共同动物不足')
for(region in c('whole','cortex')) {
  keep<-sample%in%md$SAMPLE & (region=='whole'|cortex)
  group<-match(sample[keep],md$SAMPLE)
  map<-sparseMatrix(i=seq_along(group),j=group,x=1,dims=c(length(group),nrow(md)))
  result<-x[,keep,drop=FALSE]%*%map
  # 若对象有重复符号，按基因汇总counts，而不是保留第一条。
  genes<-rownames(result); unique_genes<-unique(genes)
  collapse<-sparseMatrix(i=match(genes,unique_genes),j=seq_along(genes),x=1)
  result<-collapse%*%result;rownames(result)<-unique_genes;colnames(result)<-md$animal_id
  wr(data.frame(gene_symbol=unique_genes,as.matrix(result),check.names=FALSE),
     file.path(work,'prepared',paste0('mouse_',region,'_counts.csv.gz')))
}
wr(md,file.path(work,'prepared','mouse_units.csv'))
cat('小鼠共同动物',nrow(md),'笼位',length(unique(md$cage)),'\n')
