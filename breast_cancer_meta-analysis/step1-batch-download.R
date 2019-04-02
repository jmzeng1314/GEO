rm(list=ls())
options(stringsAsFactors = F)

### ---------------
###
### Create: Jianming Zeng
### Date: 2018-07-15 17:07:49
### Email: jmzeng1314@163.com
### Blog: http://www.bio-info-trainee.com/
### Forum:  http://www.biotrainee.com/thread-1376-1-1.html
### CAFS/SUSTC/Eli Lilly/University of Macau
### Update Log: 2018-07-31   First version
###
### ---------------

gse_list=read.table('gse_list.txt' )[,1]
gse_list
setwd('data/')
library(GEOquery)
if(F){
  
  for (gse in gse_list) {
    if(!file.exists(paste0(gse,'_eSet.Rdata'))){
      gset <- getGEO(gse, destdir=".",
                     AnnotGPL = F,
                     getGPL = F)
      save(gset,file=paste0(gse,'_eSet.Rdata'))
    }
   
  }
}
setwd('../')




