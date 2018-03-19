library(matrixStats)
library(NMF)
setwd("D:\\minÊı¾İ\\tcga\\prad_tcga\\analyse\\")
a<-read.csv(file="tcga_signature.csv",h=T,as.is=T,row.names=1)
a<-log2(a+1)
b<-a-rowMedians(data.matrix(a))


library(lattice)
d1<-b
win.graph();densityplot(unlist(d1))
quantile(unlist(d1),0.99)
win.graph();aheatmap(d1,color=colorRampPalette(c("blue", "white", "red"))(800), distfun="euclidean",
          # breaks=c(-1,1),
           hclustfun="complete",scale="none",labCol=NA,legend=T)
d2<-t(scale(t(b)))
aheatmap(d2,color=colorRampPalette(c("blue", "white", "red"))(802), distfun="euclidean",breaks=c(min(d2)-0.05,seq(-2,2,0.005),max(d2)+0.05),hclustfun="complete",scale="none",labCol=NA,legend=T)
