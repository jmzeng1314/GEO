rm(list = ls())
options(stringsAsFactors = F)
load(file = 'airway_exprSet.Rdata')
group_list
group_list=relevel(group_list,ref = 'untrt')
source('run_DEG_RNA-seq.R')
run_DEG_RNAseq(exprSet,group_list,
               g1="untrt",g2="trt",
               pro='airway')
