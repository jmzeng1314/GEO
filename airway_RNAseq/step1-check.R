library(airway)

## 表达矩阵来自于R包：  airway
# 如果当前工作目录不存在文件：'airway_exprSet.Rdata' 
# 就运行下面 if 代码块内容，载入R包airway及其数据集airway后拿到表达矩阵和表型信息
if(!file.exists('airway_exprSet.Rdata')){
  library(airway)
  data(airway)
  exprSet=assay(airway)
  group_list=colData(airway)[,3]
  save(exprSet,group_list,file = 'airway_exprSet.Rdata')
}
# 大家务必注意这样的代码技巧，多次存储Rdata文件。
# 存储后的Rdata文件，很容易就加载进入R语言工作环境。
load(file = 'airway_exprSet.Rdata')
table(group_list)
# 下面代码是为了检查
if(T){
  colnames(exprSet)
  pheatmap::pheatmap(cor(exprSet))
  group_list
  tmp=data.frame(g=group_list)
  rownames(tmp)=colnames(exprSet)
  # 组内的样本的相似性理论上应该是要高于组间的
  # 但是如果使用全部的基因的表达矩阵来计算样本之间的相关性
  # 是不能看到组内样本很好的聚集在一起。
  pheatmap::pheatmap(cor(exprSet),annotation_col = tmp)
  dim(exprSet)
  # 所以我这里初步过滤低表达量基因。
  exprSet=exprSet[apply(exprSet,1, function(x) sum(x>1) > 5),]
  dim(exprSet)
  
  exprSet=log(edgeR::cpm(exprSet)+1)
  dim(exprSet)
  # 再挑选top500的MAD值基因
  exprSet=exprSet[names(sort(apply(exprSet, 1,mad),decreasing = T)[1:500]),]
  dim(exprSet)
  M=cor(log2(exprSet+1))
  tmp=data.frame(g=group_list)
  rownames(tmp)=colnames(M)
  pheatmap::pheatmap(M,annotation_col = tmp)
  # 现在就可以看到，组内样本很好的聚集在一起
  # 组内的样本的相似性是要高于组间
  pheatmap::pheatmap(M,annotation_col = tmp,filename = 'cor.png')
  
  
  
  library(pheatmap)
  pheatmap(scale(cor(log2(exprSet+1))))
  
}

# 以上代码，就是为了检查R包airway及其数据集airway里面的表达矩阵和表型信息
