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
run_kegg(gene_up,gene_down,pro='npc_VS_normal')
# 需要多go数据库的3个条目进行3次富集分析，非常耗时。
run_go(gene_up,gene_down,pro='npc_VS_normal')

go <- enrichGO(gene_up, OrgDb = "org.Hs.eg.db", ont="all") 
barplot(go, split="ONTOLOGY")+ facet_grid(ONTOLOGY~., scale="free") 
ggsave('gene_up_GO_all_barplot.png')
go <- enrichGO(gene_down, OrgDb = "org.Hs.eg.db", ont="all") 
barplot(go, split="ONTOLOGY",font.size =10)+ 
  facet_grid(ONTOLOGY~., scale="free") + 
  scale_x_discrete(labels=function(x) str_wrap(x, width=50))+
  ggsave('gene_down_GO_all_barplot.png')




