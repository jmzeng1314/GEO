draw_h_v <- function(exprSet,need_DEG,n='DEseq2'){
  ## we only need two columns of DEG, which are log2FoldChange and pvalue
  ## heatmap
  library(pheatmap)
  choose_gene=head(rownames(need_DEG),50) ## 50 maybe better
  choose_matrix=exprSet[choose_gene,]
  choose_matrix=t(scale(t(choose_matrix)))
  pheatmap(choose_matrix,filename = paste0(n,'_need_DEG_top50_heatmap.png'))
  
  logFC_cutoff <- with(need_DEG,mean(abs( log2FoldChange)) + 2*sd(abs( log2FoldChange)) )
  # logFC_cutoff=1
  
  need_DEG$change = as.factor(ifelse(need_DEG$pvalue < 0.05 & abs(need_DEG$log2FoldChange) > logFC_cutoff,
                                     ifelse(need_DEG$log2FoldChange > logFC_cutoff ,'UP','DOWN'),'NOT')
  )
  this_tile <- paste0('Cutoff for logFC is ',round(logFC_cutoff,3),
                      '\nThe number of up gene is ',nrow(need_DEG[need_DEG$change =='UP',]) ,
                      '\nThe number of down gene is ',nrow(need_DEG[need_DEG$change =='DOWN',])
  )
  library(ggplot2)
  g = ggplot(data=need_DEG, 
             aes(x=log2FoldChange, y=-log10(pvalue), 
                 color=change)) +
    geom_point(alpha=0.4, size=1.75) +
    theme_set(theme_set(theme_bw(base_size=20)))+
    xlab("log2 fold change") + ylab("-log10 p-value") +
    ggtitle( this_tile ) + theme(plot.title = element_text(size=15,hjust = 0.5))+
    scale_colour_manual(values = c('blue','black','red')) ## corresponding to the levels(res$change)
  print(g)
  ggsave(g,filename = paste0(n,'_volcano.png'))
}


run_DEG_RNAseq <- function(exprSet,group_list,
                           g1="control",g2="case",
                           pro='test'){
  print(table(group_list))
  colnames(exprSet)
  cat(paste0('Now process the project : ',pro))
  ### ---------------
  ###
  ### Firstly run DEseq2 
  ###
  ### ---------------
  suppressMessages(library(DESeq2))  
  (colData <- data.frame(row.names=colnames(exprSet), 
                         group_list=group_list) )
  dds <- DESeqDataSetFromMatrix(countData = exprSet,
                                colData = colData,
                                design = ~ group_list)
  dds <- DESeq(dds)
  res <- results(dds, 
                 contrast=c("group_list",g2,g1))
  resOrdered <- res[order(res$padj),]
  head(resOrdered)
  DEG =as.data.frame(resOrdered)
  DEG = na.omit(DEG)
  if(F){
    
    png("DESeq2_qc_dispersions.png", 1000, 1000, pointsize=20)
    plotDispEsts(dds, main="Dispersion plot")
    dev.off()
    
    rld <- rlogTransformation(dds)
    exprMatrix_rlog=assay(rld)
    
    x=apply(exprMatrix_rlog,1,mean)
    y=apply(exprMatrix_rlog,1,mad) 
    plot(x,y) 
    
    png("DESeq2_RAWvsNORM.png",height = 800,width = 800)
    par(cex = 0.7)
    n.sample=ncol(exprSet)
    if(n.sample>40) par(cex = 0.5)
    cols <- rainbow(n.sample*1.2)
    par(mfrow=c(2,2))
    boxplot(exprSet, col = cols,main="expression value",las=2)
    boxplot(exprMatrix_rlog, col = cols,main="expression value",las=2)
    hist(as.matrix(exprSet))
    hist(exprMatrix_rlog)
    dev.off()
    
    
    
  }
  nrDEG=DEG
  DEG_DEseq2=nrDEG
  nrDEG=DEG_DEseq2[,c(2,6)]
  colnames(nrDEG)=c('log2FoldChange','pvalue') 
  draw_h_v(exprSet,nrDEG,paste0(pro,'_DEseq2'))
   
  ### ---------------
  ###
  ### Then run edgeR 
  ###
  ### ---------------
  library(edgeR)
  g=factor(group_list)
  g=relevel(g,g1)
  d <- DGEList(counts=exprSet,group=g)
  keep <- rowSums(cpm(d)>1) >= 2
  table(keep)
  d <- d[keep, , keep.lib.sizes=FALSE]
  d$samples$lib.size <- colSums(d$counts)
  d <- calcNormFactors(d)
  d$samples
  
  dge=d
  design <- model.matrix(~0+factor(group_list))
  rownames(design)<-colnames(dge)
  colnames(design)<-levels(factor(group_list))
 
  dge <- estimateGLMCommonDisp(dge,design)
  dge <- estimateGLMTrendedDisp(dge, design)
  dge <- estimateGLMTagwiseDisp(dge, design)
  
  fit <- glmFit(dge, design)
  # https://www.biostars.org/p/110861/
  lrt <- glmLRT(fit,  contrast=c(-1,1)) 
  nrDEG=topTags(lrt, n=nrow(dge))
  nrDEG=as.data.frame(nrDEG)
  head(nrDEG)
  DEG_edgeR =nrDEG 
  nrDEG=DEG_edgeR[,c(1,5)]
  colnames(nrDEG)=c('log2FoldChange','pvalue')
  draw_h_v(exprSet,nrDEG,paste0(pro,'_edgeR')) 
   
  ### ---------------
  ###
  ### Then run limma 
  ###
  ### --------------- 
  suppressMessages(library(limma))
  design <- model.matrix(~0+factor(group_list))
  colnames(design)=levels(factor(group_list))
  rownames(design)=colnames(exprSet)
  design
  
  dge <- DGEList(counts=exprSet)
  dge <- calcNormFactors(dge)
  logCPM <- cpm(dge, log=TRUE, prior.count=3)
  
  v <- voom(dge,design,plot=TRUE, normalize="quantile")
  fit <- lmFit(v, design)
  
  group_list
  con=paste0(g2,'-',g1)
  cat(con)
  cont.matrix=makeContrasts(contrasts=c(con),levels = design)
  fit2=contrasts.fit(fit,cont.matrix)
  fit2=eBayes(fit2)
  
  tempOutput = topTable(fit2, coef=con, n=Inf)
  DEG_limma_voom = na.omit(tempOutput)
  head(DEG_limma_voom)
  nrDEG=DEG_limma_voom[,c(1,4)]
  colnames(nrDEG)=c('log2FoldChange','pvalue')
  draw_h_v(exprSet,nrDEG,paste0(pro,'_limma'))
  
  
  save(DEG_limma_voom,DEG_DEseq2,DEG_edgeR,
       dds,exprSet,group_list,
       file = paste0(pro,'_DEG_results.Rdata')) 
  
  allg=intersect(rownames(DEG_limma_voom),rownames(DEG_edgeR))
  allg=intersect(allg,rownames(DEG_DEseq2))
  
  nrDEG=cbind(DEG_limma_voom[allg,c(1,4)],
              DEG_edgeR[allg,c(1,5)],
              DEG_DEseq2[allg,c(2,6)]) 
  colnames(nrDEG)
  print(cor(nrDEG[,c(1,3,5)]))
  write.csv(nrDEG,file =  paste0(pro,'_DEG_results.csv'))
  return(nrDEG)
}

getDEGs <- function(DEG_DEseq2,DEG_edgeR,DEG_limma_voom,thre_logFC=1,thre_p=0.05){
  head(DEG_DEseq2)
  head(DEG_edgeR)
  head(DEG_limma_voom)
  # thre_logFC=1;thre_p=0.05
  u1=rownames(DEG_DEseq2[with(DEG_DEseq2,log2FoldChange>thre_logFC & padj<thre_p),])
  u2=rownames(DEG_edgeR[with(DEG_edgeR,logFC>thre_logFC & FDR<thre_p),])
  u3=rownames(DEG_limma_voom[with(DEG_limma_voom,logFC>thre_logFC & adj.P.Val<thre_p),])
  
  d1=rownames(DEG_DEseq2[with(DEG_DEseq2,log2FoldChange < -thre_logFC & padj<thre_p),])
  d2=rownames(DEG_edgeR[with(DEG_edgeR,logFC < -thre_logFC & FDR<thre_p),])
  d3=rownames(DEG_limma_voom[with(DEG_limma_voom,logFC < -thre_logFC & adj.P.Val<thre_p),])
  
  u=intersect(u1,u2);u=intersect(u3,u)
  d=intersect(d1,d2);d=intersect(d3,d)
  
  return(list(up=u,down=d))
}


run_GSEA_MSIGDB <- function(geneList,pro){
  library(ggplot2)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(org.Mm.eg.db)
  ### 对 MigDB中的全部基因集 做GSEA分析。
  # http://www.bio-info-trainee.com/2105.html
  # http://www.bio-info-trainee.com/2102.html 
  { 
    #选择gmt文件（MigDB中的全部基因集）
    d='~/biosoft/MSigDB/symbols/'
    gmts <- list.files(d,pattern = 'all')
    print(gmts)
    #GSEA分析
    library(GSEABase) # BiocManager::install('GSEABase')
    ## 下面使用lapply循环读取每个gmt文件，并且进行GSEA分析
    ## 如果存在之前分析后保存的结果文件，就不需要重复进行GSEA分析。
    f=paste0(pro,'gsea_results.Rdata')
    if(!file.exists(f)){
      gsea_results <- lapply(gmts, function(gmtfile){
        # gmtfile=gmts[2]
        geneset <- read.gmt(file.path(d,gmtfile)) 
        print(paste0('Now process the ',gmtfile))
        # minGSSize = 10, maxGSSize = 500, 
        # pvalueCutoff = 0.05, pAdjustMethod = "BH"
        egmt <- GSEA(geneList, TERM2GENE=geneset, verbose=FALSE)
        head(egmt)
        # gseaplot(egmt, geneSetID = rownames(egmt[1,]))
        
        return(egmt)
      })
      # 上面的代码耗时，所以保存结果到本地文件
      save(gsea_results,file = f)
    }
    load(file = f)
    #提取gsea结果，熟悉这个对象
    gsea_results_list<- lapply(gsea_results, function(x){
      cat(paste(dim(x@result)),'\n')
      x@result
    })
    gsea_results_df <- do.call(rbind, gsea_results_list)
    gseaplot(gsea_results[[2]],'KEGG_CELL_CYCLE') 
    gseaplot(gsea_results[[2]],'FARMER_BREAST_CANCER_CLUSTER_6') 
    gseaplot(gsea_results[[5]],'GO_CONDENSED_CHROMOSOME_OUTER_KINETOCHORE') 
    return(gsea_results_df)
  }
}

