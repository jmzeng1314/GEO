rm(list=ls())


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

options(stringsAsFactors = F)
library(stringr)

a=read.table('meta/meta1.txt',skip = 1,fill = T)
head(a)
a=a[,c(1,3,5)]
colnames(a)=c('gsm','month','event')
a$study='GSE12093'
meta=a

a=read.table('meta/meta2.txt',skip = 1,fill = T)
head(a)
a[,2]=as.numeric(str_split(a[,2],":",simplify = T)[,2])/30
a[,3]=as.numeric(str_split(a[,3],":",simplify = T)[,2]) 
colnames(a)=c('gsm','month','event')
a$study='GSE9195'
meta=rbind(meta,a)

a=read.table('meta/meta3.txt',skip = 1,fill = T)
head(a)
a=a[,c(1,5,3)] 
colnames(a)=c('gsm','month','event')
a$study='GSE1456'
head(a)
a$month=12*a$month
meta=rbind(meta,a)

a=read.table('meta/meta4.txt',skip = 1,fill = T)
head(a)
a=a[,c(1,3,5)]
a[,2]=a[,2]/30
colnames(a)=c('gsm','month','event')
head(a)
a$study='GSE6532'
meta=rbind(meta,a)

a=read.table('meta/meta5.txt',skip = 1,fill = T)
head(a)
a=a[,c(1,7,4)] 
colnames(a)=c('gsm','month','event')
head(a)
a$study='GSE2603'
meta=rbind(meta,a)

a=read.table('meta/meta6.txt',skip = 1,fill = T)
head(a)
a=na.omit(a)
a=a[,c(1,5,3)] 
a[,2]=a[,2]/30
colnames(a)=c('gsm','month','event')
head(a)
a$study='GSE16446'
meta=rbind(meta,a)

a=read.table('meta/meta7.txt',skip = 1,fill = T)
head(a)
a=na.omit(a)
a=a[,c(1,8,4)]  
colnames(a)=c('gsm','month','event')
head(a)
a$study='GSE21653'
meta=rbind(meta,a)

a=read.table('meta/meta8.txt',skip = 1,fill = T)
head(a)
a=na.omit(a)
a=a[,c(1,4,7)]  
a[,2]=a[,2]/30
colnames(a)=c('gsm','month','event')
head(a)
a$study='GSE2990'
meta=rbind(meta,a)

# a=read.table('meta/meta9.txt',skip = 1,fill = T)
# head(a)
# a=na.omit(a)
# a=a[,c(1,4,7)]  
# a[,2]=a[,2]/30
# colnames(a)=c('gsm','month','event')
# head(a)
# a$study='GSE11121'
# meta=rbind(meta,a)

a=read.table('meta/meta10.txt',skip = 1,fill = T)
head(a)
a=na.omit(a)
a=a[,c(1,3,5)]  
a[,2]=a[,2]/30
colnames(a)=c('gsm','month','event')
head(a)
a$study='GSE7390'
meta=rbind(meta,a)


a=read.table('meta/meta11.txt',skip = 1,fill = T)
head(a)
a=na.omit(a)
a=a[,c(1,3,2)]   
colnames(a)=c('gsm','month','event')
head(a)
a$study='GSE4922'
meta=rbind(meta,a)


a=read.table('meta/meta12.txt',skip = 1,fill = T)
head(a)
a=na.omit(a)
a=a[,c(1,3,2)]   
colnames(a)=c('gsm','month','event')
head(a)
a$study='GSE17705'
meta=rbind(meta,a)

a=read.table('meta/meta13.txt',skip = 1,fill = T)
head(a)
a=na.omit(a)
a=a[,c(1,3,2)]   
colnames(a)=c('gsm','month','event')
head(a)
a$study='GSE20685'
meta=rbind(meta,a)

meta[,2]=as.numeric(meta[,2])
meta[,3]=as.numeric(meta[,3])
meta=na.omit(meta)

boxplot(as.numeric(meta[,2]))
table(meta[,3])

save(meta,file = 'meta.rdata')













