rm(list=ls())
options(stringsAsFactors = F)
gse_list=read.table('gse_list.txt' )[,1]
gse_list
library(GEOquery)  
list.files('data','*.Rdata')
a=list.files('data','*.Rdata')

load(file.path('data',a[1]))
gset[[1]] 
exprs(gset[[1]])['202048_s_at',]
data.frame(exprs(gset[[1]])['202048_s_at',])
if(F){
  cbx6=lapply(1:length(a),function(x){
    load(file.path('data',a[x]))
    data.frame(exprs(gset[[1]])['202048_s_at',])
  })
  cbx6
  save(cbx6,file = 'cbx6.rdata')
}
load(file = 'cbx6.rdata')


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

do.call(rbind,cbx6)

cxb6_g <- lapply(cbx6,function(x){
v=as.numeric(x[,1])
x$g=ifelse(v>median(v),'high','low')
return(x)
})
do.call(rbind,cxb6_g)

cxb6_expr <- do.call(rbind,cxb6_g)
View(cxb6_expr)
cxb6_expr$gsm=rownames(cxb6_expr)
head(cxb6_expr)
colnames(cxb6_expr)=c('value','group','gsm')
save(cxb6_expr,file = 'cxb6_expr.rdata')
write.csv(cxb6_expr,file = 'cxb6_expr.csv')



