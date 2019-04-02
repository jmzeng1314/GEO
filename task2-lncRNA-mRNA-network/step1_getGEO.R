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
  eSet <- getGEO('GSE54238', destdir=".",
                 AnnotGPL = F,
                 getGPL = F)
  save(eSet,file='GSE54238_eSet.Rdata')
}
load('GSE54238_eSet.Rdata')
b = eSet[[1]]
raw_exprSet=exprs(b) 
phe=pData(b)
library(stringr)
group_list= str_split(as.character(phe$source_name_ch1),' ',simplify = T)[,1]
save(raw_exprSet,group_list,
     file='GSE54238_raw_exprSet.Rdata')



rm(list=ls())
if(F){
  library(GEOquery)
  eSet <- getGEO('GSE14520', destdir=".",
                 AnnotGPL = F,
                 getGPL = F)
  save(eSet,file='GSE14520_eSet.Rdata')
}
load('GSE14520_eSet.Rdata')
b = eSet[[1]]
raw_exprSet=exprs(b) 
raw_exprSet[1:4,1:4]
phe=pData(b)
phe$title
library(stringr)


