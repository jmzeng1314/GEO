rm(list=ls())
### ---------------
###
### Create: Jianming Zeng
### Date: 
### Email: jmzeng1314@163.com
### Blog: http://www.bio-info-trainee.com/
### Forum:  http://www.biotrainee.com/thread-1376-1-1.html
### CAFS/SUSTC/Eli Lilly/University of Macau
### Update Log: 2018-07-09  First version
###
### ---------------
load(file='GSE3325_raw_exprSet.Rdata')
library(hgu133plus2.db)
eg2probe=toTable(hgu133plus2SYMBOL)
eg2probe[eg2probe$symbol=='TRAF4',]
raw_exprSet[1:4,1:4]
exprSet=log2(raw_exprSet)
dat=data.frame(gene= exprSet['211899_s_at',] ,
               mut= group_list)
head(dat)
if(require('ggpubr')){
  library(ggpubr)
  # google search : ggpubr boxplot add p-value
  # http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/76-add-p-values-and-significance-levels-to-ggplots/
  p <- ggboxplot(dat, x = "mut", y = "gene",
                 color = "mut", palette = "jco",
                 add = "jitter")
  #  Add p-value
  p + stat_compare_means(method = "t.test")
}

if(require('ggstatsplot')){
  library(ggstatsplot)
  ggbetweenstats(data = dat, x = mut,  y = gene)
}


if(require('ggplot2')){
  library(ggplot2)
  ggplot(dat,aes(x=mut,y=gene))+
    geom_boxplot()+
    theme_bw()
}



