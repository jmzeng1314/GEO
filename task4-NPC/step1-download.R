## 
### ---------------
###
### Create: Jianming Zeng
### Date: 2018-12-20 15:43:52
### Email: jmzeng1314@163.com
### Blog: http://www.bio-info-trainee.com/
### Forum:  http://www.biotrainee.com/thread-1376-1-1.html
### CAFS/SUSTC/Eli Lilly/University of Macau
### Update Log: 2018-12-20  First version
###
### ---------------


rm(list = ls())  ## 魔幻操作，一键清空~
options(stringsAsFactors = F)#在调用as.data.frame的时，将stringsAsFactors设置为FALSE可以避免character类型自动转化为factor类型
# 注意查看下载文件的大小，检查数据 
f='GSE12452_eSet.Rdata'
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE12452
library(GEOquery)
# 这个包需要注意两个配置，一般来说自动化的配置是足够的。
#Setting options('download.file.method.GEOquery'='auto')
#Setting options('GEOquery.inmemory.gpl'=FALSE)
if(!file.exists(f)){
  gset <- getGEO('GSE12452', destdir=".",
                 AnnotGPL = F,     ## 注释文件
                 getGPL = F)       ## 平台文件
  save(gset,file=f)   ## 保存到本地
}
load('GSE12452_eSet.Rdata')  ## 载入数据
class(gset)  #查看数据类型
length(gset)  #
class(gset[[1]])
gset
# assayData: 44149 features, 41 samples

# 因为这个GEO数据集只有一个GPL平台，所以下载到的是一个含有一个元素的list
a=gset[[1]] #
dat=exprs(a) #a现在是一个对象，取a这个对象通过看说明书知道要用exprs这个函数
dim(dat)#看一下dat这个矩阵的维度
# GPL13667
dat[1:4,1:4] #查看dat这个矩阵的1至4行和1至4列，逗号前为行，逗号后为列
boxplot(dat,las=2)
pd=pData(a) #通过查看说明书知道取对象a里的临床信息用pData
## 挑选一些感兴趣的临床表型。

group_list= ifelse(grepl('normal',pd$title),'normal','npc')
table(group_list)

dat[1:4,1:4] 

# GPL570

if(F){
  library(GEOquery)
  #Download GPL file, put it in the current directory, and load it:
  gpl <- getGEO('GPL570', destdir=".")
  colnames(Table(gpl))  
  head(Table(gpl)[,c(1,15)]) ## you need to check this , which column do you need
  probe2gene=Table(gpl)[,c(1,15)]
  head(probe2gene)
  library(stringr)  
  save(probe2gene,file='probe2gene.Rdata')
}
# 
# load(file='probe2gene.Rdata')
# ids=probe2gene 

library(hgu133plus2.db)
ids=toTable(hgu133plus2SYMBOL) #toTable这个函数：通过看hgu133plus2.db这个包的说明书知道提取probe_id（探针名）和symbol（基因名）的对应关系的表达矩阵的函数为toTable
head(ids) #head为查看前六行

head(ids)
colnames(ids)=c('probe_id','symbol')  
ids=ids[ids$symbol != '',]
ids=ids[ids$probe_id %in%  rownames(dat),]

dat[1:4,1:4]   
dat=dat[ids$probe_id,] 

ids$median=apply(dat,1,median) #ids新建median这一列，列名为median，同时对dat这个矩阵按行操作，取每一行的中位数，将结果给到median这一列的每一行
ids=ids[order(ids$symbol,ids$median,decreasing = T),]#对ids$symbol按照ids$median中位数从大到小排列的顺序排序，将对应的行赋值为一个新的ids
ids=ids[!duplicated(ids$symbol),]#将symbol这一列取取出重复项，'!'为否，即取出不重复的项，去除重复的gene ，保留每个基因最大表达量结果s
dat=dat[ids$probe_id,] #新的ids取出probe_id这一列，将dat按照取出的这一列中的每一行组成一个新的dat
rownames(dat)=ids$symbol#把ids的symbol这一列中的每一行给dat作为dat的行名
dat[1:4,1:4]  #保留每个基因ID第一次出现的信息


save(dat,group_list,file = 'step1-output.Rdata')



 