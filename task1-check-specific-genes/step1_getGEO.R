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
rm(list=ls())
if(F){
  library(GEOquery)
  eSet <- getGEO('GSE3325', destdir=".",
                 AnnotGPL = F,
                 getGPL = F)
  save(eSet,file='GSE3325_eSet.Rdata')
}
load('GSE3325_eSet.Rdata')
b = eSet[[1]]
raw_exprSet=exprs(b) 
phe=pData(b)
library(stringr)
group_list= str_split(as.character(phe$title),' ',simplify = T)[,1]
save(raw_exprSet,group_list,
     file='GSE3325_raw_exprSet.Rdata')



rm(list=ls())
if(F){
  library(GEOquery)
  eSet <- getGEO('GSE17708', destdir=".",
                 AnnotGPL = F,
                 getGPL = F)
  save(eSet,file='GSE17708_eSet.Rdata')
}
load('GSE17708_eSet.Rdata')
b = eSet[[1]]
raw_exprSet=exprs(b) 
raw_exprSet[1:4,1:4]
phe=pData(b)
phe$title
library(stringr)
group_list= str_split(as.character(phe$title),' ',simplify = T)[,11]
group_list=paste0('group',group_list,'h')
save(raw_exprSet,group_list,
     file='GSE17708_raw_exprSet.Rdata')


