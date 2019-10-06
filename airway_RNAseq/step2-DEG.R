rm(list = ls())
options(stringsAsFactors = F)
load(file = 'airway_exprSet.Rdata')
group_list
group_list=relevel(group_list,ref = 'untrt')
source('run_DEG_RNA-seq.R')
# 这个 run_DEG_RNAseq 函数，是我自定义的
# 主要是包装了3个RNA-seq数据分析的R包
# 以及部分可视化函数
run_DEG_RNAseq(exprSet,group_list,
               g1="untrt",g2="trt",
               pro='airway')
