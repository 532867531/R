library(matrixStats)
library(NMF)


library(factoextra)
myclust<-function(x){
  hclust(x)
}
myDist<-function(x){
  get_dist(x,method = "euclidean")
}



a<-read.csv(file="tcga_signature.csv",h=T,as.is=T,row.names=1)
a<-log2(a+1)
b<-a-rowMedians(data.matrix(a))

d1<-b
library(lattice)
library(gplots)
win.graph();densityplot(unlist(d1))
fix(data_signatire)
quantile(data_signatire_1,0.005);quantile(data_signatire_1,0.955)
win.graph();heatmap.2(as.matrix(d1),col=colorRampPalette(c("blue", "white", "red"))(802), distfun=myDist,hclustfun=myclust,scale="none",
                      trace="none",
                      dendrogram = "both",
                      breaks = c(min(d1)-0.05,seq(-2,2,0.005),max(d1)+0.05))
d2<-t(scale(t(b)))
win.graph();aheatmap(d2,color=colorRampPalette(c("blue", "white", "red"))(802), distfun="euclidean",breaks=c(min(d2)-0.05,seq(-2,2,0.005),max(d2)+0.05),hclustfun="complete",scale="none",labCol=NA,legend=T)
