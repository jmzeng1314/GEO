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
options(stringsAsFactors = F)
load(file = 'step1-output.Rdata')
table(group_list)
# 每次都要检测数据
dat[1:4,1:4]
head(phe)
# R里面实现生存分析非常简单！

# 用my.surv <- surv(OS_MONTHS,OS_STATUS=='DECEASED')构建生存曲线。
# 用kmfit2 <- survfit(my.surv~TUMOR_STAGE_2009)来做某一个因子的KM生存曲线。
# 用 survdiff(my.surv~type, data=dat)来看看这个因子的不同水平是否有显著差异，其中默认用是的logrank test 方法。
# 用coxph(Surv(time, status) ~ ph.ecog + tt(age), data=lung) 来检测自己感兴趣的因子是否受其它因子(age,gender等等)的影响。


exprSet=dat
dim(exprSet)
colnames(phe)
colnames(phe)=c('event','grade','node','size','time')
phe = as.data.frame(apply(phe,2,as.numeric))
boxplot(phe$size)
phe$size=ifelse(phe$size>median(phe$size),'big','small')

library(survival)
library(survminer)
# 利用ggsurvplot快速绘制漂亮的生存曲线图
sfit <- survfit(Surv(time, event)~size, data=phe)
sfit
summary(sfit)
ggsurvplot(sfit, conf.int=F, pval=TRUE)
ggsurvplot(sfit,palette = c("#E7B800", "#2E9FDF"),
           risk.table =TRUE,pval =TRUE,
           conf.int =TRUE,xlab ="Time in months", 
           ggtheme =theme_light(), 
           ncensor.plot = TRUE)
## 多个 ggsurvplots作图生存曲线代码合并代码公布
sfit1=survfit(Surv(time, event)~size, data=phe)
sfit2=survfit(Surv(time, event)~grade, data=phe)
splots <- list()
splots[[1]] <- ggsurvplot(sfit1,pval =TRUE, data = phe, risk.table = TRUE)
splots[[2]] <- ggsurvplot(sfit2,pval =TRUE, data = phe, risk.table = TRUE)
# Arrange multiple ggsurvplots and print the output
arrange_ggsurvplots(splots, print = TRUE,  ncol = 2, nrow = 1, risk.table.height = 0.4)
# 可以看到grade跟生存显著相关，而size跟病人生存的关系并不显著。

## 挑选感兴趣的基因做差异分析
phe$CBX4=ifelse(exprSet['CBX4',]>median(exprSet['CBX4',]),'high','low')
table(phe$CBX4)
ggsurvplot(survfit(Surv(time, event)~CBX4, data=phe), conf.int=F, pval=TRUE)
phe$CBX6=ifelse(exprSet['CBX6',]>median(exprSet['CBX6',]),'high','low')
table(phe$CBX6)
ggsurvplot(survfit(Surv(time, event)~CBX6, data=phe), conf.int=F, pval=TRUE)
## 批量生存分析 使用  logrank test 方法
mySurv=with(phe,Surv(time, event))
log_rank_p <- apply(exprSet , 1 , function(gene){
  #gene=exprSet[1,]
  phe$group=ifelse(gene>median(gene),'high','low')  
  data.survdiff=survdiff(mySurv~group,data=phe)
  p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
  return(p.val)
})
log_rank_p=sort(log_rank_p)
head(log_rank_p)
boxplot(log_rank_p) 
phe$H6PD=ifelse(exprSet['H6PD',]>median(exprSet['H6PD',]),'high','low')
table(phe$H6PD)
ggsurvplot(survfit(Surv(time, event)~H6PD, data=phe), conf.int=F, pval=TRUE)

# 前面如果我们使用了WGCNA找到了跟grade相关的基因模块，然后在这里就可以跟生存分析的显著性基因做交集
# 这样就可以得到既能跟grade相关，又有临床预后意义的基因啦。

## 批量生存分析 使用 coxph 回归方法
colnames(phe)
mySurv=with(phe,Surv(time, event))
cox_results <-apply(exprSet , 1 , function(gene){
  group=ifelse(gene>median(gene),'high','low')
  survival_dat <- data.frame(group=group,grade=phe$grade,size=phe$size,stringsAsFactors = F)
  m=coxph(mySurv ~ grade + size + group, data =  survival_dat)
  
  beta <- coef(m)
  se <- sqrt(diag(vcov(m)))
  HR <- exp(beta)
  HRse <- HR * se
  
  #summary(m)
  tmp <- round(cbind(coef = beta, se = se, z = beta/se, p = 1 - pchisq((beta/se)^2, 1),
                     HR = HR, HRse = HRse,
                     HRz = (HR - 1) / HRse, HRp = 1 - pchisq(((HR - 1)/HRse)^2, 1),
                     HRCILL = exp(beta - qnorm(.975, 0, 1) * se),
                     HRCIUL = exp(beta + qnorm(.975, 0, 1) * se)), 3)
  return(tmp['grouplow',])
  
})
cox_results=t(cox_results)
table(cox_results[,4]<0.05)
cox_results[cox_results[,4]<0.05,]

length(setdiff(rownames(cox_results[cox_results[,4]<0.05,]),
               names(log_rank_p[log_rank_p<0.05])
))
length(setdiff( names(log_rank_p[log_rank_p<0.05]),
                rownames(cox_results[cox_results[,4]<0.05,])
))
length(unique( names(log_rank_p[log_rank_p<0.05]),
                rownames(cox_results[cox_results[,4]<0.05,])
))
save(log_rank_p,cox_results,exprSet,phe,file = 'surviva.Rdata')





