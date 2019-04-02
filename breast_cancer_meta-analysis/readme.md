# A brief example for meta-analysis 

Please first read the paper carefully: [Breast Cancer Res Treat.](https://www.ncbi.nlm.nih.gov/pubmed/20020197#) 2010  ， An online survival analysis tool to rapidly assess the effect of `22,277 genes` on breast cancer prognosis using microarray data of `1,809 patients.`

### GSE studies 

I find a list of GSE studies in that article, then find annother GSE studies in their website: http://kmplot.com/analysis/index.php?p=updates

```
GSE16446
GSE17907
GSE19615
GSE20685
GSE21653
GSE17705
GSE2603
GSE12276
GSE16391
GSE12093
GSE11121
GSE9195
GSE7390
GSE6532
GSE5327
GSE4922
GSE3494
GSE2990
GSE2034
GSE1456
```

### step1: download all of them by  GEOquery

```R
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
```

### step2: get the expression matrix by specific gene

```r
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
```

### step3: Try to understand the metadata of all samples.



### step4: Do survival analysis based on data from GEO 



### step5:Do survival analysis based on data from the supplementry of that paper 

we can download from :http://kmplot.com/analysis/index.php?p=download 

[![img](http://kmplot.com/analysis/pic/zip-logo.png)](http://kmplot.com/analysis/studies/@MAS5_1000_1809_rounded_final.zip) Supplemental database for our **2010 Breast Cancer Res Treatment **paper: a TXT file of the normalized gene expression data for 1,809 breast cancer microarrays.

[![img](http://kmplot.com/analysis/pic/txt-logo.png)](http://kmplot.com/analysis/studies/@expdesc_1809.txt) Supplemental table for our **2010 Breast Cancer Res Treatment **paper: a TXT file of the experimental descriptors for 1,809 breast cancer samples.

By using these information, we can do survival analysis.



