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
a=read.table('~/Downloads/MAS5_1809_kerekitett_vegleges.txt',header = T)
a[1:4,1:4]
rownames(a)=a[,1]
a=a[,-1]
a[1:4,1:4]
b=read.table('~/Downloads/expdesc_1809.txt',header = T,fill = T,sep = '\t')[,c(1,3:6)]
head(b)
colnames(b)=c('gsm','rfs_event','rfs_time','os_event','os_time')
cxb6_expr <- data.frame(gsm=colnames(a),v=as.numeric(a['202048_s_at',]))
cxb6_expr$g=ifelse(cxb6_expr$v>median(cxb6_expr$v),'high','low')
head(cxb6_expr)
dat=merge(cxb6_expr,b,by='gsm')
dat$rfs_time=as.numeric(gsub(',','.',dat$rfs_time))
dat$os_time=as.numeric(gsub(',','.',dat$os_time))
head(dat)

library(survival)
library(survminer)
# 利用ggsurvplot快速绘制漂亮的生存曲线图
#dat=dat[dat$study==1,]
sfit <- survfit(Surv(rfs_time, rfs_event)~g, data=dat)
sfit <- survfit(Surv(os_time, os_event)~g, data=dat)
sfit
summary(sfit)
ggsurvplot(
  sfit, risk.table = TRUE, ggtheme = theme_bw(),
  pval = TRUE, pval.coord = c(0, 0.03)
)

ggsurvplot(sfit, conf.int=F, pval=TRUE)

ggsurvplot(sfit,palette = c("#E7B800", "#2E9FDF"),
           risk.table =TRUE,pval =TRUE,
           conf.int =TRUE,xlab ="Time in months", 
           ggtheme =theme_light(), 
           ncensor.plot = TRUE)
ggsurvplot(
  sfit, risk.table = TRUE, ggtheme = theme_bw(),
  pval = TRUE, pval.coord = c(0, 0.03)
) 



