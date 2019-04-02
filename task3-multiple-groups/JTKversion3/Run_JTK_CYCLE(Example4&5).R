## Example4
## key features of dataset: 2h-resolution, no replicate at each time point but contain two missing time points
## main analysis request: periodic signals with period length between 20h and 28h

source("JTK_CYCLEv3.1.R")

project <- "Example4"

options(stringsAsFactors=FALSE)
data <- read.delim("Example4_data.txt")

Names <- data[,1]
data <- data[,-1]

jtkdist(ncol(data), reps=1)
periods <- 10:14
jtk.init(periods,2);

cat("JTK analysis started on",date(),"\n")
flush.console()

st <- system.time({
  res <- apply(data,1,function(z) {
	jtkx(z)
  c(JTK.ADJP,JTK.PERIOD,JTK.LAG,JTK.AMP)
  })
  res <- as.data.frame(t(res))
  bhq <- p.adjust(unlist(res[,1]),"BH")
  res <- cbind(bhq,res)
  colnames(res) <- c("BH.Q","ADJ.P","PER","LAG","AMP")
  results <- cbind(Names,res,data)
  results <- results[order(res$ADJ.P,-res$AMP),]
})
print(st)

write.table(results,file=paste("JTK",project,"txt",sep="."),row.names=F,col.names=T,quote=F,sep="\t")
	
## Example5
## key features of dataset: 4h-resolution, three replicates at each time point, and contain random missing values in each row
## main analysis request: 24h periodic signals 

source("JTK_CYCLEv3.1.R")

project <- "Example5"

options(stringsAsFactors=FALSE)
data <- read.delim("Example5_data.txt")

Names <- data[,1]
data <- data[,-1]

jtkdist(timepoints=10, reps=3)
periods <- 6:6
jtk.init(periods,4);

cat("JTK analysis started on",date(),"\n")
flush.console()

st <- system.time({
  res <- apply(data,1,function(z) {
    jtkx(z)
    c(JTK.ADJP,JTK.PERIOD,JTK.LAG,JTK.AMP)
  })
  res <- as.data.frame(t(res))
  bhq <- p.adjust(unlist(res[,1]),"BH")
  res <- cbind(bhq,res)
  colnames(res) <- c("BH.Q","ADJ.P","PER","LAG","AMP")
  results <- cbind(Names,res,data)
  results <- results[order(res$ADJ.P,-res$AMP),]
})
print(st)

write.table(results,file=paste("JTK",project,"txt",sep="."),row.names=F,col.names=T,quote=F,sep="\t")
