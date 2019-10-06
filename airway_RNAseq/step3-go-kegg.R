rm(list = ls())
options(stringsAsFactors = F)
load('airway_DEG_results.Rdata')
source('run_DEG_RNA-seq.R')
deg=getDEGs(DEG_DEseq2,DEG_edgeR,DEG_limma_voom,thre_logFC=1,thre_p=0.05)

 
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
gene_up= bitr(unique(deg$up), fromType = "ENSEMBL",
           toType = c( "ENTREZID"),
           OrgDb = org.Hs.eg.db)[,2] 

gene_down= bitr(unique(deg$down), fromType = "ENSEMBL",
              toType = c( "ENTREZID"),
              OrgDb = org.Hs.eg.db)[,2] 
gene_diff=c(gene_up,gene_down) 

source('kegg_and_go_up_and_down.R')
# 同样的，里面包装了一些代码，比如setReadable函数
# 很有可能你使用的时候就发现过期了里面有一些参数
# 要学会调试代码，不要畏手畏脚。
run_kegg(gene_up,gene_down,pro='airway_test')
# 需要多go数据库的3个条目进行3次富集分析，非常耗时。
run_go(gene_up,gene_down,pro='airway_test')
# 很多绘图代码，都是依据数据本身特性需要调整的，而且高阶情况下需要AI等等。
go <- enrichGO(gene_up, OrgDb = "org.Hs.eg.db", ont="all") 
barplot(go, split="ONTOLOGY")+ facet_grid(ONTOLOGY~., scale="free") 
ggsave('gene_up_GO_all_barplot.png')
go <- enrichGO(gene_down, OrgDb = "org.Hs.eg.db", ont="all") 
barplot(go, split="ONTOLOGY",font.size =10)+ 
  facet_grid(ONTOLOGY~., scale="free") + 
  scale_x_discrete(labels=function(x) str_wrap(x, width=50))+
  ggsave('gene_down_GO_all_barplot.png')




