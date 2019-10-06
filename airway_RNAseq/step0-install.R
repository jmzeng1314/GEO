### ---------------
###
### Create: Jianming Zeng
### Date: 2018-07-15 17:07:49
### Email: jmzeng1314@163.com
### Blog: http://www.bio-info-trainee.com/
### Forum:  http://www.biotrainee.com/thread-1376-1-1.html
### CAFS/SUSTC/Eli Lilly/University of Macau
### Update Log: 2018-07-09  First version
### 2019-10-06 second version 
### ---------------

## 假如没有安装包，就运行下面被注释的代码
if(F){
  options()$repos  ## 查看使用install.packages安装时的默认镜像
  options()$BioC_mirror ##查看使用bioconductor的默认镜像
  options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/") 
  ##指定镜像，这个是中国科技大学镜像
  options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")) 
  ##指定install.packages安装镜像，这个是清华镜像
  options()$repos 
  options()$BioC_mirror
  
  if (!requireNamespace("BiocManager", quietly = TRUE)) 
    install.packages("BiocManager") ##判断是否存在BiocManager包，不存在的话安装
  
  library(BiocManager) 
  BiocManager::install(c('airway','DESeq2','edgeR','limma'),ask = F,update = F)
  
  
}

library(DESeq2)
library(edgeR)
library(limma)
library(airway)

