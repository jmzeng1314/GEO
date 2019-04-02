rm(list = ls())  ## 魔幻操作，一键清空~
options(stringsAsFactors = F)
load(file = 'step1-output.Rdata')
# 每次都要检测数据
dat[1:4,1:4] 
table(group_list) #table函数，查看group_list中的分组个数
#通过为每个数据集绘制箱形图，比较数据集中的数据分布
boxplot(dat[1,]~group_list) #按照group_list分组画箱线图

bp=function(g){         #定义一个函数g，函数为{}里的内容
  library(ggpubr)
  df=data.frame(gene=g,stage=group_list)
  p <- ggboxplot(df, x = "stage", y = "gene",
                 color = "stage", palette = "jco",
                 add = "jitter")
  #  Add p-value
  p + stat_compare_means()
}
bp(dat[1,]) ## 调用上面定义好的函数，避免同样的绘图代码重复多次敲。
bp(dat[2,])
bp(dat[3,])
bp(dat[4,])
bp(dat['CDK9',])
dim(dat)
if(F){
  library(limma)
  design=model.matrix(~factor( group_list ))
  fit=lmFit(dat,design)
  fit=eBayes(fit)
  ## 上面是limma包用法的一种方式 
  options(digits = 4) #设置全局的数字有效位数为4
  #topTable(fit,coef=2,adjust='BH') 
  topTable(fit,coef=2,adjust='BH') 
  ## 但是上面的用法做不到随心所欲的指定任意两组进行比较
  
}
library(limma)
design <- model.matrix(~0+factor(group_list))
colnames(design)=levels(factor(group_list))
head(design)
exprSet=dat
rownames(design)=colnames(exprSet)
design
contrast.matrix<-makeContrasts("Tumor-Normal",
                               levels = design)
contrast.matrix ##这个矩阵声明，我们要把 Tumor 组跟 Normal 进行差异分析比较

deg = function(exprSet,design,contrast.matrix){
  ##step1
  fit <- lmFit(exprSet,design)
  ##step2
  fit2 <- contrasts.fit(fit, contrast.matrix) 
  ##这一步很重要，大家可以自行看看效果
  
  fit2 <- eBayes(fit2)  ## default no trend !!!
  ##eBayes() with trend=TRUE
  ##step3
  tempOutput = topTable(fit2, coef=1, n=Inf)
  nrDEG = na.omit(tempOutput) 
  #write.csv(nrDEG2,"limma_notrend.results.csv",quote = F)
  head(nrDEG)
  return(nrDEG)
}

deg = deg(exprSet,design,contrast.matrix)

head(deg)

save(deg,file = 'deg.Rdata')

load(file = 'deg.Rdata')

## for volcano 
if(T){
  nrDEG=deg
  head(nrDEG)
  attach(nrDEG)
  plot(logFC,-log10(P.Value))
  library(ggpubr)
  df=nrDEG
  df$v= -log10(P.Value) #df新增加一列'v',值为-log10(P.Value)
  ggscatter(df, x = "logFC", y = "v",size=0.5)
  
  df$g=ifelse(df$P.Value>0.01,'stable', #if 判断：如果这一基因的P.Value>0.01，则为stable基因
              ifelse( df$logFC >2,'up', #接上句else 否则：接下来开始判断那些P.Value<0.01的基因，再if 判断：如果logFC >1.5,则为up（上调）基因
                      ifelse( df$logFC < -2,'down','stable') )#接上句else 否则：接下来开始判断那些logFC <1.5 的基因，再if 判断：如果logFC <1.5，则为down（下调）基因，否则为stable基因
  )
  table(df$g)
  df$name=rownames(df)
  head(df)
  write.table(rownames(df[df$g == 'up',]),
       file = 'up_gene.txt',
       quote=F,row.names=F,col.names=F)
  write.table(rownames(df[df$g == 'down',]),
              file = 'down_gene.txt',
              quote=F,row.names=F,col.names=F)
  
  ggscatter(df, x = "logFC", y = "v",size=0.5,color = 'g')
  ggscatter(df, x = "logFC", y = "v", color = "g",size = 0.5,
            label = "name", repel = T,
            #label.select = rownames(df)[df$g != 'stable'] ,
            label.select = c('TTC9', 'AQP3', 'CXCL11','PTGS2'), #挑选一些基因在图中显示出来
            palette = c("#00AFBB", "#E7B800", "#FC4E07") )
  ggsave('volcano.png')
  
  ggscatter(df, x = "AveExpr", y = "logFC",size = 0.2)
  df$p_c = ifelse(df$P.Value<0.001,'p<0.001',
                  ifelse(df$P.Value<0.01,'0.001<p<0.01','p>0.01'))
  table(df$p_c )
  ggscatter(df,x = "AveExpr", y = "logFC", color = "p_c",size=0.2, 
            palette = c("green", "red", "black") )
  ggsave('MA.png')
  
  
}

## for heatmap 
if(T){ 
  load(file = 'step1-output.Rdata')
  # 每次都要检测数据
  dat[1:4,1:4]
  table(group_list)
  x=deg$logFC #deg取logFC这列并将其重新赋值给x
  names(x)=rownames(deg) #deg取probe_id这列，并将其作为名字给x
  cg=c(names(head(sort(x),100)),#对x进行从小到大排列，取前100及后100，并取其对应的探针名，作为向量赋值给cg
       names(tail(sort(x),100)))
  library(pheatmap)
  pheatmap(dat[cg,],show_colnames =F,show_rownames = F) #对dat按照cg取行，所得到的矩阵来画热图
  n=t(scale(t(dat[cg,])))#通过“scale”对log-ratio数值进行归一化，现在的dat是行名为探针，列名为样本名，由于scale这个函数应用在不同组数据间存在差异时，需要行名为样本，因此需要用t(dat[cg,])来转换，最后再转换回来
  
  n[n>2]=2
  n[n< -2]= -2
  n[1:4,1:4]
  pheatmap(n,show_colnames =F,show_rownames = F)
  ac=data.frame(g=group_list)
  rownames(ac)=colnames(n) #将ac的行名也就分组信息（是‘no TNBC’还是‘TNBC’）给到n的列名，即热图中位于上方的分组信息
  pheatmap(n,show_colnames =F,
           show_rownames = F,
          cluster_cols = F, 
           annotation_col=ac,filename = 'heatmap_top200_DEG.png') #列名注释信息为ac即分组信息
  
  
}

write.csv(deg,file = 'deg.csv')




