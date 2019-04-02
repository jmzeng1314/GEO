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

load(file = 'cxb6_expr.rdata')
head(cxb6_expr)
load(file = 'meta.rdata')
# 0=right censored, 1=event at time, 2=left censored, 3=interval censored.
meta=meta[meta$month<120,]
head(meta)

dat=merge(cxb6_expr,meta,by='gsm')
head(dat)
table(dat$study)

library(survival)
library(survminer)
# 利用ggsurvplot快速绘制漂亮的生存曲线图
#dat=dat[dat$study==1,]
sfit <- survfit(Surv(month, event)~group, data=dat)
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

gse_list=unique(dat$study)

i=gse_list[13]
sub_dat=dat[dat$study==i,] 
sfit <- survfit(Surv(month, event)~group, data=sub_dat)
sfit
summary(sfit)
file=paste0('survival_study_',sub_dat[1,6],'.pdf')
p=ggsurvplot(sfit, conf.int=F, pval=TRUE)

pdf(file)
print(p$plot, newpage = FALSE)
dev.off()



