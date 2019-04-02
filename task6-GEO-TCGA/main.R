### ---------------
###
### Create: Jianming Zeng
### Date: 2018-10-14 16:57:30
### Email: jmzeng1314@163.com
### Blog: http://www.bio-info-trainee.com/
### Forum:  http://www.biotrainee.com/thread-1376-1-1.html
### CAFS/SUSTC/Eli Lilly/University of Macau
### Update Log: 2018-10-14 16:57:21  First version
###
### ---------------

### ---------------
###
### install the packages 
###
### ---------------
rm(list = ls())
if(!require('GEOquery')){
  source("https://bioconductor.org/biocLite.R") 
  options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/")  
  BiocInstaller::biocLite('GEOquery')
  install.packages('pheatmap')
}

library(clusterProfiler)
library(org.Hs.eg.db)
library(GEOquery)
# GSE64392
### ---------------
###
### step1:download data 
###
### ---------------
library(GEOquery)
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42872
eSet <- getGEO('GSE64392', destdir=".",
               AnnotGPL = F,
               getGPL = F)
save(eSet,file='GSE64392_eSet.Rdata')



### ---------------
###
### step2: visualization 
###
### ---------------
## hclust
rm(list = ls())
load(file='GSE64392_eSet.Rdata')
ob=eSet[[1]]
mat=exprs(ob)
meta=pData(ob)
exprSet=mat
group_list=ifelse(grepl('t',as.character(meta$title)),'tumor','normal')

colnames(exprSet)=paste(group_list,1:ncol(exprSet),sep='_')
# Define nodePar
nodePar <- list(lab.cex = 0.6, pch = c(NA, 19),
                cex = 0.7, col = "blue")
hc=hclust(dist(t(exprSet)))
par(mar=c(5,5,5,10))
png('hclust.png',width = 890,height = 1111,res=120)
plot(as.dendrogram(hc), nodePar = nodePar, horiz = TRUE)
dev.off()

## PCA

library(ggfortify)
df=as.data.frame(t(exprSet))
df$group=group_list
png('pca.png',res=120)
autoplot(prcomp( df[,1:(ncol(df)-1)] ), data=df,colour = 'group')+theme_bw()
dev.off()

M=cor(exprSet)
library(pheatmap)
pheatmap(M )


### ---------------
###
### step3: DEG 
###
### ---------------

rm(list = ls())
library(GEOquery)
load(file='GSE64392_eSet.Rdata')
ob=eSet[[1]]
mat=exprs(ob)
meta=pData(ob)
exprSet=mat
group_list=ifelse(grepl('t',as.character(meta$title)),'tumor','normal')

# DEG by limma 
suppressMessages(library(limma)) 
design <- model.matrix(~0+factor(group_list))
colnames(design)=levels(factor(group_list))
rownames(design)=colnames(exprSet)
design
contrast.matrix<-makeContrasts(paste0(unique(group_list),collapse = "-"),levels = design)
contrast.matrix<-makeContrasts("tumor-normal",levels = design)

contrast.matrix 
##这个矩阵声明，我们要把progres.组跟stable进行差异分析比较
##step1
fit <- lmFit(exprSet,design)
##step2
fit2 <- contrasts.fit(fit, contrast.matrix) ##这一步很重要，大家可以自行看看效果
fit2 <- eBayes(fit2)  ## default no trend !!!
##eBayes() with trend=TRUE
##step3
tempOutput = topTable(fit2, coef=1, n=Inf)
nrDEG = na.omit(tempOutput) 
#write.csv(nrDEG2,"limma_notrend.results.csv",quote = F)
head(nrDEG)
nrDEG$probe=rownames(nrDEG)
ids=read.table('ht20probe.txt',sep = '\t',header = T,fill=T)[,c(1,6)]
head(ids)
colnames(ids)=c('probe','accession')
ids=ids[ids$accession!='',]
library(org.Hs.eg.db)
id2acc=toTable(org.Hs.egACCNUM)
id2symbol=toTable(org.Hs.egSYMBOL)
tmp=merge(nrDEG,ids,by='probe')
tmp=merge(tmp,id2acc,by='accession')
nrDEG=merge(tmp,id2symbol,by='gene_id') 
save(nrDEG,file='GSE64392-nrDEG.Rdata')

attach(nrDEG)
plot(logFC,-log10(P.Value))
library(ggpubr)
df=nrDEG
df$v= -log10(P.Value)
ggscatter(df, x = "logFC", y = "v",size=0.5)

df$g=ifelse(df$P.Value>0.05,'stable',
            ifelse( df$logFC >1,'up',
                    ifelse( df$logFC < -1,'down','stable') )
)
table(df$g)
df$name=rownames(df)
ggscatter(df, x = "logFC", y = "v", color = "g",size = 0.5,
          label = "symbol", repel = T,
          #label.select = rownames(df)[df$g != 'stable'] ,
          label.select = c('CACHD1','SERTAD4'),
          palette = c("#00AFBB", "#E7B800", "#FC4E07") )

ggscatter(df, x = "AveExpr", y = "logFC",size = 0.2)
df$p_c = ifelse(df$P.Value<0.001,'p<0.001',
                ifelse(df$P.Value<0.01,'0.001<p<0.01','p>0.01'))
table(df$p_c )
ggscatter(df,x = "AveExpr", y = "logFC", color = "p_c",size=0.2, 
          palette = c("green", "red", "black") )



dat=exprSet[ df$probe[df$g != 'stable'] ,]
library(pheatmap)
pheatmap(dat,scale='row',fontsize_row =3)

exprSet[1:4,1:4]
exprSet=exprSet[df$probe,]
df$symbol

ids=df[,c(3,10)]
head(ids)
colnames(ids)=c('probe_id','symbol')
length(unique(ids$symbol))
tail(sort(table(ids$symbol)))
table(sort(table(ids$symbol)))
plot(table(sort(table(ids$symbol))))

table(rownames(exprSet) %in% ids$probe_id)
dim(exprSet)
exprSet=exprSet[rownames(exprSet) %in% ids$probe_id,]
dim(exprSet)

ids=ids[match(rownames(exprSet),ids$probe_id),]
head(ids)
exprSet[1:5,1:5]

jimmy <- function(exprSet,ids){
  tmp = by(exprSet,
           ids$symbol,
           function(x) rownames(x)[which.max(rowMeans(x))] )
  probes = as.character(tmp)
  print(dim(exprSet))
  exprSet=exprSet[rownames(exprSet) %in% probes ,]
  
  print(dim(exprSet))
  rownames(exprSet)=ids[match(rownames(exprSet),ids$probe_id),2]
  return(exprSet)
}

new_exprSet <- jimmy(exprSet,ids)
save(new_exprSet,file = 'new_exprSet.Rdata')


### ---------------
###
### step4: annotation
###
### ---------------

### ---------------
###
### step5: annotation
###
### ---------------

### ---------------
###
### step6: with TCGA 
###
### ---------------

load('new_exprSet.Rdata')
new_exprSet=as.data.frame(new_exprSet)
colnames(new_exprSet)=paste(group_list,1:ncol(exprSet),sep='_')
new_exprSet$SYMBOL=rownames(new_exprSet)

a=read.table('TCGA-COAD.htseq_counts.tsv.gz',header = T,stringsAsFactors = F)
a[1:4,1:4]
a$Ensembl_ID=as.character(a$Ensembl_ID)
a$Ensembl_ID=unlist(lapply(strsplit(a$Ensembl_ID,'[.]'),function(x) x[1]))

library(clusterProfiler)
library(org.Hs.eg.db)
df <- bitr(a$Ensembl_ID, fromType ='ENSEMBL' ,
           toType = c( "ENTREZID","SYMBOL"),
           OrgDb = org.Hs.eg.db)
head(df)
a=a[a$Ensembl_ID %in%  df$ENSEMBL ,]
a[1:4,1:4]
head(df)
tmp1=merge(a,df,by.x='Ensembl_ID',by.y='ENSEMBL')
tmp2=merge(tmp1,new_exprSet,by='SYMBOL')
expr=tmp2
#rownames(expr)=expr$ENTREZID
expr=expr[,grepl('TCGA',colnames(expr)) | grepl('normal',colnames(expr))  | grepl('tumor',colnames(expr)) ]

test=expr[,500:549]
save(test,file = 'test.Rdata')

load(file = 'test.Rdata')
boxplot(test,las=2)
source("http://bioconductor.org/biocLite.R")
options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/")
biocLite("sva")

library("sva")
dataMat=test
mod <- data.frame("(Intercept)"=rep(1, ncol(dataMat)))
rownames(mod) <- colnames(dataMat)
whichbatch <- as.factor(
  c(rep("tcga", 13), 
    rep("arrar", 37)))
combatout <- ComBat(as.matrix(dataMat), whichbatch, mod=mod)
boxplot(combatout,las=2)









