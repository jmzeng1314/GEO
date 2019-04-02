# R语言学徒作业

首先需要看完GEO视频并且尝试理解代码在：https://github.com/jmzeng1314/GEO 

视频在： <https://www.bilibili.com/video/av26731585/>

### 作业1 

看懂文章：https://www.jci.org/articles/view/96060/figure/1  看其C子图里面的TRAF4基因在4个数据集的表达量，画出更漂亮的boxplot。

提示：需要看完文章，了解作者所引用的数据并且下载对应的数据集，提取TRAF4基因对应的探针的表达量，根据对应的分组信息画boxplot。

- 2010-cancer cell MSKCC [GSE21032](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE21032)    [ProstateCancer](https://www.sciencedirect.com/topics/medicine-and-dentistry/prostate-cancer) Genomics Data Portalat [http://cbio.mskcc.org/prostate-portal](http://cbio.mskcc.org/prostate-portal/)  
- 2005-cancer cell GSE3325 Affymetrix U133 2.0 Plus arrays 
- 2007-BMCCancer. [GSE6919](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE6919) 
- 2012- Nature.  GEO([GSE35988](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE35988)).

### 作业2

了解数据集 ：[GSE17708](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE17708) 对应的文章： PMID: [20007254](https://www.ncbi.nlm.nih.gov/pubmed/20007254) 并且搞清楚该文章涉及的样本，实验设计。

找到最后一个时间点处理(72 h) 的 3个样本和3个untreated的A549 lung adenocarcinoma cell line的**差异表达基因集**，以及其**GO/KEGG富集分析**结果。

然后看看 BMC Systems Biology 2014 的文章是如何重新利用这个数据集的。  <https://doi.org/10.1186/1752-0509-8-55> 列出其分析点。

#### 进阶

还可以看看GSEA,GSVA是如何作用于整个表达矩阵，不局限于72小时的。

还可以看看这个R包和教程。https://blog.csdn.net/msw521sg/article/details/75452019  如何根据药物处理时间来找模块。

或者学习下面的几个R包：

```
Mfuzz
MaSigPro
ImpulseDE2
EBSeq-HMM
```

还可以使用WGCNA来分析这个数据集。https://github.com/jmzeng1314/my_WGCNA  

### 其它作业

下面这些芯片数据所依赖的文章看懂，查询到，然后下载数据集自己分析一波。

- GSE11072  2009-gastric cancer SBC Human 16K cDNA Microarray
- GSE42872 2015-melanoma-vemurafenib HuGene-1_0-st 
- GSE24673 2015-hub-gene-mcode-retinoblastoma  HuGene-1_0-st 
- GSE22863  2011-NSCLC HuGene-1_0-st 
- GSE622221, GSE4180414, GSE5140122 A total of 117 samples (54 cases and 63 controls) Affymetrix Human Genome U133 Plus 2.0 Array  2015-HCC
- GSE21815  2016-CRC Agilent-014850 Whole Human Genome Microarray 4x44K

 



