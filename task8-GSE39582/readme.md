# P值0.05界限为什么不好

感觉好像每隔几年学术界就掀起了一阵批判P值的讨论，正好我这里有个例子，分享给大家。

文献是：[Cancer Manag Res](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6149864/#). 2018 标题经过了更改：

- 最开始时  Tumour purity as a prognostic factor in colon cancer
- 然后是  Low tumor purity is associated with poor prognosis, heavy mutation burden, and intense immune phenotype in colon cancer

数据来源： [GSE17536](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE17536)/17537 and [GSE39582](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE39582) (n=794),  还有TCGA的 (n=454).Colorectal cancer (CRC) 

作者收集了 GSE17536/17537,GSE39582, and TCGA 共 1248病人，分析表明 Stage III and MMR-deficient
(dMMR) 病人具有比较低的肿瘤纯度，而较低的较低的肿瘤纯度是一个独立的病人预后风险因子。(大家去文章里面看P值)

以前的pan-cancer研究表明，在 idney renal clear cell carcinoma and lower-grade glioma 里面，高纯度的肿瘤预示着好的生存。

### 使用ESTIMATE算法看肿瘤纯度

首先下载表达矩阵，然后使用 算法处理它，得到每个样本的肿瘤纯度值，进行统计，如下：

![image-20190401210230509](http://www.bio-info-trainee.com/wp-content/uploads/2019/04/image-20190401210230509.png)

把肿瘤纯度来工具不同的分类指标进行量化，如下：

![image-20190401210328937](http://www.bio-info-trainee.com/wp-content/uploads/2019/04/image-20190401210328937.png)

对值得探索的因子绘森林图：

![image-20190401210403183](http://www.bio-info-trainee.com/wp-content/uploads/2019/04/image-20190401210403183.png)

看肿瘤纯度对生存的影响：

![image-20190401210422725](http://www.bio-info-trainee.com/wp-content/uploads/2019/04/image-20190401210422725.png)

### 使用CIBERSORT算法看细胞比例

作者 We applied CIBERSORT algorithm on `RNA-Seq and microarray data` to estimate the relative proportion of 22 immune cells of leukocytes for each CC sample.  

![image-20190401210633415](http://www.bio-info-trainee.com/wp-content/uploads/2019/04/image-20190401210633415.png)

发现：immunotherapy-associated markers (PD-1, PD-L1, CTLA-4, LAG-3, and TIM-3) were high expressed in low purity CC. 散点图展示相关性，如下：

![image-20190401210851460](http://www.bio-info-trainee.com/wp-content/uploads/2019/04/image-20190401210851460.png)

然后这些细胞比例和肿瘤纯度的相关性如下：

![image-20190401210757741](http://www.bio-info-trainee.com/wp-content/uploads/2019/04/image-20190401210757741.png)

