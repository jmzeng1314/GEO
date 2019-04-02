rm(list = ls())  ## 魔幻操作，一键清空~
options(stringsAsFactors = F)
load(file = 'step1-output.Rdata')
if(F){
  library(utils)
  rforge <- "http://r-forge.r-project.org"
  install.packages("estimate", repos=rforge, dependencies=TRUE)
}
library(estimate)
help(package="estimate")

library(estimate)
dim(dat)
dat[1:4,1:4]
datf='tmp.txt'
write.table(dat,file=datf,sep = '\t',quote = F)
filterCommonGenes(input.f=datf, 
                  output.f="GSE39582.gct", 
                  id="GeneSymbol")
estimateScore(input.ds = "GSE39582.gct",
              output.ds="GSE39582_estimate_score.gct", 
              platform="affymetrix")
plotPurity(scores="GSE39582_estimate_score.gct", samples="GSM971957", 
           platform="affymetrix")
scores=read.table("GSE39582_estimate_score.gct",skip = 2,header = T)
rownames(scores)=scores[,1]
scores=t(scores[,3:ncol(scores)])
scores
dim(pd)
identical(rownames(pd),rownames(scores))
boxplot(scores[,4]~pd$source_name_ch1,las=2)
boxplot(scores[,4]~pd$ch,las=2)
pd=pd[,grepl('tnm',colnames(pd))]

save(pd,scores,file = 'output-by-estimate.Rdata')
library(ggpubr)
box <- lapply(colnames(pd),function(i) {
  dat <-  pd[,i,drop=F]  
  dat$group=dat[,1]
  dat$purity=scores[,4]
  ## 画boxplot 
  p <- ggboxplot(dat, x = "group", y = 'purity', 
                 add = "jitter" )
  p
})
library(cowplot)
plot_grid(plotlist=box, ncol=2 )





