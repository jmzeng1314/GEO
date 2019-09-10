### ---------------
###
### Create: Jianming Zeng
### Date: 2018-07-15 17:07:49
### Email: jmzeng1314@163.com
### Blog: http://www.bio-info-trainee.com/
### Forum:  http://www.biotrainee.com/thread-1376-1-1.html
### CAFS/SUSTC/Eli Lilly/University of Macau
### Update Log: 2018-07-09  First version
###
### ---------------

## 假如没有安装包，就运行下面被注释的代码
if(F){
  options()$repos  ## 查看使用install.packages安装时的默认镜像
  options()$BioC_mirror ##查看使用bioconductor的默认镜像
  options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/") 
  ##指定镜像，这个是中国科技大学镜像
  options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")) 
  ##指定install.packages安装镜像，这个是清华镜像
  options()$repos 
  options()$BioC_mirror
  
  if (!requireNamespace("BiocManager", quietly = TRUE)) 
    install.packages("BiocManager") ##判断是否存在BiocManager包，不存在的话安装
  
  library(BiocManager) 
  BiocManager::install(c('airway','DESeq2','edgeR','limma'),ask = F,update = F)
 
  
}

library(DESeq2)
library(edgeR)
library(limma)
library(airway)

## 表达矩阵来自于R包：  airway
if(F){
  library(airway)
  data(airway)
  exprSet=assay(airway)
  group_list=colData(airway)[,3]
  group_list=relevel(group_list,ref = 'trt')
  save(exprSet,group_list,file = 'airway_exprSet.Rdata')
}

rm(list = ls())
options(stringsAsFactors = F)
load(file = 'airway_exprSet.Rdata')

if(T){
  colnames(exprSet)
  pheatmap::pheatmap(cor(exprSet))
  group_list
  tmp=data.frame(g=group_list)
  rownames(tmp)=colnames(exprSet)
  # 组内的样本的相似性应该是要高于组间的！
  pheatmap::pheatmap(cor(exprSet),annotation_col = tmp)
  dim(exprSet)
  exprSet=exprSet[apply(exprSet,1, function(x) sum(x>1) > 5),]
  dim(exprSet)
  
  exprSet=log(edgeR::cpm(exprSet)+1)
  dim(exprSet)
  exprSet=exprSet[names(sort(apply(exprSet, 1,mad),decreasing = T)[1:500]),]
  dim(exprSet)
  M=cor(log2(exprSet+1))
  tmp=data.frame(g=group_list)
  rownames(tmp)=colnames(M)
  pheatmap::pheatmap(M,annotation_col = tmp)
  pheatmap::pheatmap(M,annotation_col = tmp,filename = 'cor.png')
  
  
  
  library(pheatmap)
  pheatmap(scale(cor(log2(exprSet+1))))
  
}

load(file = 'airway_exprSet.Rdata')

source('functions.R')
### ---------------
###
### Firstly run DEseq2 
###
### ---------------
suppressMessages(library(DESeq2)) 
group_list=c('untrt','trt' ,'untrt','trt' ,'untrt','trt','untrt','trt')
(colData <- data.frame(row.names=colnames(exprSet), 
                       group_list=group_list) )
dds <- DESeqDataSetFromMatrix(countData = exprSet,
                              colData = colData,
                              design = ~ group_list)
dds <- DESeq(dds)
res <- results(dds, 
               contrast=c("group_list","trt","untrt"))
resOrdered <- res[order(res$padj),]
head(resOrdered)
DEG =as.data.frame(resOrdered)
DEG = na.omit(DEG)
if(F){
  
  png("DESeq2_qc_dispersions.png", 1000, 1000, pointsize=20)
  plotDispEsts(dds, main="Dispersion plot")
  dev.off()
  
  rld <- rlogTransformation(dds)
  exprMatrix_rlog=assay(rld)
  
  x=apply(exprMatrix_rlog,1,mean)
  y=apply(exprMatrix_rlog,1,mad) 
  plot(x,y) 
  
  png("DESeq2_RAWvsNORM.png",height = 800,width = 800)
  par(cex = 0.7)
  n.sample=ncol(exprSet)
  if(n.sample>40) par(cex = 0.5)
  cols <- rainbow(n.sample*1.2)
  par(mfrow=c(2,2))
  boxplot(exprSet, col = cols,main="expression value",las=2)
  boxplot(exprMatrix_rlog, col = cols,main="expression value",las=2)
  hist(as.matrix(exprSet))
  hist(exprMatrix_rlog)
  dev.off()
  
  
  
}
nrDEG=DEG
DEseq_DEG=nrDEG
nrDEG=DEseq_DEG[,c(2,6)]
colnames(nrDEG)=c('log2FoldChange','pvalue') 
draw_h_v(exprSet,nrDEG,'DEseq2')
source('functions.R')
### ---------------
###
### Then run edgeR 
###
### ---------------
library(edgeR)
d <- DGEList(counts=exprSet,group=factor(group_list))
keep <- rowSums(cpm(d)>1) >= 2
table(keep)
d <- d[keep, , keep.lib.sizes=FALSE]
d$samples$lib.size <- colSums(d$counts)
d <- calcNormFactors(d)
d$samples

design <- model.matrix(~0+factor(group_list))
rownames(design)<-colnames(dge)
colnames(design)<-levels(factor(group_list))
dge=d
dge <- estimateGLMCommonDisp(dge,design)
dge <- estimateGLMTrendedDisp(dge, design)
dge <- estimateGLMTagwiseDisp(dge, design)

fit <- glmFit(dge, design)
# https://www.biostars.org/p/110861/
lrt <- glmLRT(fit,  contrast=c(-1,1)) 
nrDEG=topTags(lrt, n=nrow(dge))
nrDEG=as.data.frame(nrDEG)
head(nrDEG)
edgeR_DEG =nrDEG 
nrDEG=edgeR_DEG[,c(1,5)]
colnames(nrDEG)=c('log2FoldChange','pvalue')
draw_h_v(exprSet,nrDEG,'edgeR')

source('functions.R')
### ---------------
###
### Then run edgeR 
###
### --------------- 
suppressMessages(library(limma))
design <- model.matrix(~0+factor(group_list))
colnames(design)=levels(factor(group_list))
rownames(design)=colnames(exprSet)
design

dge <- DGEList(counts=exprSet)
dge <- calcNormFactors(dge)
logCPM <- cpm(dge, log=TRUE, prior.count=3)

v <- voom(dge,design,plot=TRUE, normalize="quantile")
fit <- lmFit(v, design)

group_list
cont.matrix=makeContrasts(contrasts=c('untrt-trt'),levels = design)
fit2=contrasts.fit(fit,cont.matrix)
fit2=eBayes(fit2)

tempOutput = topTable(fit2, coef='untrt-trt', n=Inf)
DEG_limma_voom = na.omit(tempOutput)
head(DEG_limma_voom)
nrDEG=DEG_limma_voom[,c(1,4)]
colnames(nrDEG)=c('log2FoldChange','pvalue')
draw_h_v(exprSet,nrDEG,'limma')


save(DEG_limma_voom,DEseq_DEG,edgeR_DEG,
     dds,exprSet,group_list,
     file = 'DEG_results.Rdata')


rm(list = ls())
load(file = 'DEG_results.Rdata')
source('functions.R')

nrDEG=DEG_limma_voom[,c(1,4)]
colnames(nrDEG)=c('log2FoldChange','pvalue')
draw_h_v(exprSet,nrDEG,'limma')

nrDEG=edgeR_DEG[,c(1,5)]
colnames(nrDEG)=c('log2FoldChange','pvalue')
draw_h_v(exprSet,nrDEG,'edgeR')


nrDEG=DEseq_DEG[,c(2,6)]
colnames(nrDEG)=c('log2FoldChange','pvalue') 
draw_h_v(exprSet,nrDEG,'DEseq2')




