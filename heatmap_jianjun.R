library(matrixStats)
library(NMF)

a<-read.csv(file="tcga_signature.csv",h=T,as.is=T,row.names=1)
a<-log2(a+1)
b<-a-rowMedians(data.matrix(a))

d1<-b
aheatmap(d1,color=colorRampPalette(c("blue", "white", "red"))(802), distfun="euclidean",breaks=c(min(d1)-0.05,seq(-4,4,0.05),max(d1)+0.05),hclustfun="complete",scale="none",labCol=NA,legend=T)
d2<-t(scale(t(b)))
aheatmap(d2,color=colorRampPalette(c("blue", "white", "red"))(802), distfun="euclidean",breaks=c(min(d2)-0.05,seq(-2,2,0.005),max(d2)+0.05),hclustfun="complete",scale="none",labCol=NA,legend=T)