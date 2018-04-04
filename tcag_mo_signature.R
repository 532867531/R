##读取数据ran_seq_median
path="D:\\min数据\\tcga\\prad_tcga\\prad_tcga\\"
file=paste(path,"data_RNA_Seq_v2_expression_median.RDS",sep="")
useSymbol = TRUE
#file=paste(path,"data_RNA_Seq_v2_mRNA_median_Zscores.RDS",sep="")
if(file.exists(file)){
  rna_seq_data=readRDS(file=file)
}else{
  rna_data=read.csv(file=sub(pattern =".RDS",replacement = ".txt",x=file),sep="")
  saveRDS(rna_data,file=file)
}

library(dplyr)
library(sqldf)

##数据初步矫正去重，认为筛选的数据不在重复的基因中
ifelse(useSymbol,{
  d=duplicated(rna_seq_data$Hugo_Symbol);
  rna_seq_data=rna_seq_data[!d,];
  ("完成去重")
},
"跳过去重")
##选取基因
#geneSignatures=c("CXCR4","CXCR2","ITGAM","ITGAX","ANPEP","CD14","FUT4","CD33","CD34","CD38","ENTPD1","PTPRC","CEACAM8","CD80","CSF1R","IL4R","CSF3","CSF2","CXCL8","TNF","CXCL12","CSF1R","S100A8","S100A9","STAT1","STAT3","STAT5A","ARG1","NOS2","CD274","TLR3","TLR4","TGFB1","IL10","IDO1","PDCD1")
#geneSignatures=c("CEACAM8","S100A8","S100A9","CXCR2","ID01","CSF2","CSF3","PDCD1","TNF","PTPRC","CD38","CD80","CD274","STAT1","CXCR4","CD34","ENTPD1","TGFB1","IL4R","STAT3","STAT5A","FUT4","ITGAX","ITGAM","CD14","CSF1R","CD33","TLR4","CXCL12","TLR3","ANPEP")
##mo_signature

#reg=paste(geneSignatures,collapse = "|^")
#data_signatire=rna_seq_data[which(!is.na(stringi::stri_match(rna_seq_data$Hugo_Symbol,regex = reg))),]
#rna_seq_data=na.omit(rna_seq_data)
ifelse(useSymbol,{
  rna_seq_data=rna_seq_data[-which(is.na(rna_seq_data$Hugo_Symbol)),]
  rownames(rna_seq_data)<-rna_seq_data$Hugo_Symbol;
  "基因symbol应用完毕"
},{
  rownames(rna_seq_data)=rna_seq_data$Entrez_Gene_Id;
  "基因ID应用完毕"
})


##数据格式矫正
#rownames(data_signatire)=data_signatire$Hugo_Symbol

data_signatire=rna_seq_data[geneSignatures,]
data_signatire=na.omit((data_signatire))
tryCatch({
  data_signatire=data_signatire[,-which(colnames(data_signatire)=="Entrez_Gene_Id")]
  data_signatire=data_signatire[,-which(colnames(data_signatire)=="Hugo_Symbol")]
},finally = {
  "全部去除完毕!"
})
##查看未找到的名称
rownames(data_signatire)=geneSignatures
fix(data_signatire)
print(summary(is.na(rownames(data_signatire))))
data_signatire=log2(data_signatire+1)
# # 数据集群性评估，使用get_clust_tendency()计算Hopkins统计量
# res = get_clust_tendency(data_signatire, 40, graph = TRUE)
# res$hopkins_stat
# win.graph();res$plot
# ##由于k均值聚类需要指定要生成的聚类数量，因此我们将使用函数clusGap()来计算用于估计最优聚类数。函数fviz_gap_stat()用于可视化。
# set.seed(123)
# ## Compute the gap statistic
# gap_stat = clusGap(data_signatire, FUN = kmeans, nstart = 25, K.max = 10, B = 500)
# # Plot the result
# win.graph();fviz_gap_stat(gap_stat)
# ##对数据进行聚类分析
# res.hc=eclust(data_signatire,FUNcluster = "hclust"
#               ,hc_metric = "euclidean")
# 
# win.graph();fviz_dend(res.hc, rect = TRUE) # dendrogam
library(factoextra)
library(cluster)
# Compute k-means
res.km = eclust(t(data_signatire), "hclust")
# Gap statistic plot
win.graph();fviz_gap_stat(res.km$gap_stat)
# #下面的R代码生成Silhouette plot和分层聚类散点图。
# win.graph();fviz_silhouette(res.hc) # silhouette plot
# win.graph();fviz_cluster(res.hc) # scatter plot

fix(data_signatire)
#data_signatire=scale(data_signatire)
fix(data_signatire)
max(data_signatire)
min(data_signatire)
median(data_signatire)
mean(data_signatire)
##提取聚类轮廓图
sil = silhouette(res.km$cluster, dist(data_signatire))
rownames(sil) = rownames(USArrests)
head(sil[, 1:3])


#1、首先用dist()函数计算变量间距离
# dist.r = dist(data_signatire, method="euclidean")
# #其中method包括6种方法，表示不同的距离测度：”euclidean”, “maximum”, “manhattan”, “canberra”, “binary” or “minkowski”。相应的意义自行查找。
# #2、再用hclust()进行聚类
# hc.r = hclust(dist.r, method = "euclidean")
# #其中method包括7种方法，表示聚类的方法：”ward”, “single”, “complete”,”average”, “mcquitty”, “median” or “centroid”。相应的意义自行查找。
# #3、画图
# plot(hc.r, hang = -1,labels=NULL) 或者plot(hc.r, hang = 0.1,labels=F)

Cluster_Method<-c( "ward.D","ward.D2","single","complete","average" ,"mcquitty","median","centroid")
Dist_Methods<-c("euclidean", "maximum", "manhattan", "canberra", "binary" ,"minkowski")
for(amethods in  Dist_Methods){
  for(onemethod in Cluster_Method){
    for(k in c(3:5)){
      hc=hclust(d=dist(((data_signatire)),method = amethods),
                method = onemethod)
      hv=cutree(hc,k=k)
      out_put_dir="D:\\min数据\\tcga\\prad_tcga\\analyse\\"
      file=paste(out_put_dir,amethods,onemethod,k,".tiff",sep="_")
      tiff(filename = file,width = 1000, height = 800, units = "px")
      plot(hc,main = paste(amethods,onemethod,k,sep="_"))
      dev.off()
      #win.graph()
      #plot(hv,main = paste(onemethod,k,sep="_"))
      d=as.data.frame(hv[hc$order])
      colnames(d)="count"
      ##进行分类汇总
      library(sqldf)
      r=sqldf('select count(count) from d group by count')
      papaer=sort(c(118,206,174),decreasing = TRUE)
      ours=sort(c(r[,1]),decreasing = TRUE)
      
      matr=matrix(c(papaer,ours[c(1:3)]),byrow = TRUE,ncol = 3,nrow = 2)
      #matr=matrix(c(papaer,ours[c(1:3)]),byrow = TRUE,ncol = 3,nrow = 2)
      print(paste(cor(x=c(1:32),y=hc$order,method = "spearman"),amethods,onemethod,k,sep="_"))
      #print(paste("距离","_",k,amethods,onemethod,dist(matr)))
    }
  }
}


##删除离群数据
fix(data_signatire)

out_put_dir="D:\\min数据\\tcga\\prad_tcga\\analyse\\"
goHeatMap=function(cluster_method){
  myclust<-function(x){
    hclust(x,method=cluster_method)
  }
  library("RColorBrewer")
  library(gplots)
  ##win.graph();
  file=paste(out_put_dir,cluster_method,".tiff",sep="")
  devices=dev.list()
  tiff(filename = file,width = 8000, height = 2000, units = "px")
  hv=heatmap.2(as.matrix(data_signatire), col=rev( colorRampPalette(brewer.pal(11, "RdBu"))(256)[c((1+0):(128-20),(128+20):(256-0))]
  ), scale="row",
  distfun = dist,
  hclustfun = myclust,
  key=TRUE, keysize=1, symkey=FALSE, density.info="histogram", trace="none", cexRow=1,
  Rowv=FALSE, Colv=TRUE,dendrogram = "column",lmat=rbind(c(0,3,0),c(2,1,0),c(0,4,0)) ,lhei=c(3,10,3),lwid=c(2,9,2)
  ,key.title ="",cexCol = 1#,densadj = 0.25
  ,sepcolor="white"
  ,main = cluster_method
  )
  dev.off()
  hv
}


for(onedev in devices){
  dev.off()
}
for(one in Cluster_Method){
  hv=goHeatMap(one)
}