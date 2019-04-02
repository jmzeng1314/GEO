
source("JTK_CYCLEv3.1.R")

project <- "Example3"

options(stringsAsFactors=FALSE)
annot <- read.delim("Example3_annot.txt")
data <- read.delim("Example3_data.txt")

rownames(data) <- data[,1]
data <- data[,-1]
jtkdist(10,3)

periods <- 6:6
jtk.init(periods,4)

cat("JTK analysis started on",date(),"\n")
flush.console()

st <- system.time({
  res <- apply(data,1,function(z) {
    jtkx(z, ampci=TRUE, conf=0.95)
    c(JTK.ADJP,JTK.PERIOD,JTK.LAG,JTK.AMP,JTK.AMP.CI,JTK.AMP.PVAL)
  })
  res <- as.data.frame(t(res))
  bhq <- p.adjust(unlist(res[,1]),"BH")
  res <- cbind(bhq,res)
  colnames(res) <- c("BH.Q", "ADJ.P","PER","LAG","AMP","AMP.LO","AMP.HI","AMP.PVAL")
  results <- cbind(annot,res,data)
  results <- results[order(res$ADJ.P,-res$AMP),]
})
print(st)

save(results,file=paste("JTK",project,"rda",sep="."))
write.table(results,file=paste("JTK",project,"txt",sep="."),row.names=F,col.names=T,quote=F,sep="\t")

