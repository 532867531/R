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
##geneSignatures=c("ANPEP","CD38","CSF3","CXCR2","CD80","IDO","SLEB2","CD274","PTPRC","TLR3","CEACAM8","STAT1","CXCR4","LOC116196","CSF1R","ITGAM","ITGAX","CD14","TGFB1","CXCL12","LOC199828","ENTPD1","FUT4","STAT3","IL4R","STAT5A","TLR4","CSF2","CXCL8","S100A8","S100A9","TNF")
geneSignatures=c("CD14","IL4R","PTPRC","ITGAM","LOC116196","ARG1","IL10","CD40","CD32","CD163","FCER2","CD200R1","PDCD1LG2","CD68","CSF1R","HLA-DRA","LY75","LOC90262","CCL2","JM2")  ##high=171,mid=260,low=67
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

fix(rna_seq_data)
##数据格式矫正
#rownames(data_signatire)=data_signatire$Hugo_Symbol

data_signatire=rna_seq_data[geneSignatures,]

fix(data_signatire)
tryCatch({
  data_signatire=data_signatire[,-which(colnames(data_signatire)=="Entrez_Gene_Id")]
  data_signatire=data_signatire[,-which(colnames(data_signatire)=="Hugo_Symbol")]
},finally = {
  "全部去除完毕!"
})

##查看未找到的名称
rownames(data_signatire)=geneSignatures
#write.csv(x=t(data_signatire),file = paste("D:\\min数据\\tcga\\prad_tcga\\analyse\\spss\\log2.csv",sep=""))
#shell.exec("D:\\min数据\\tcga\\prad_tcga\\analyse\\spss\\")
if(FALSE){
  {"是否手动执行行scale()，基因表达在不同的样品之间"}
  for(rowindex in c(1:nrow(data_signatire))){
    print(rowindex)
    data_signatire[rowindex,]=scale(data_signatire)
  }  
}
fix(data_signatire)
print(summary(is.na(rownames(data_signatire))))
if(TRUE){
  {"是否进行log2处理"}
  data_signatire=log2(data_signatire+1)  
  # write.csv(x=t(data_signatire),file = paste("D:\\min数据\\tcga\\prad_tcga\\analyse\\spss\\log2_true.csv",sep=""))
  #shell.exec("D:\\min数据\\tcga\\prad_tcga\\analyse\\spss\\")
}
if(FALSE){
  data_signatire=scale(data_signatire)
  fix(data_signatire)
}
if(FALSE){
  {"进行minMax归一化进行处理"}
  mapminmax=function(x,yMIN=-1,yMax=1){
    rowMin=min(x)
    rowMax=max(x)
    rowInterVal=rowMax-rowMin
    yMin=yMin;yMax=yMax;yInterVal=yMax-yMin
    yInterVal*(x-rowMin)/rowInterVal+yMin
  }
  fix(data_signatire)
  for(rowindex in c(1:nrow(data_signatire))){
    data_signatire[rowindex,]=mapminmax(data_signatire[rowindex,])
  }
  fix(data_signatire)
}

if(FALSE){
  {"是否进行sigmoid处理"}
  # 定义归一化方程
  sigmoid = function(x,a=1,centerX=0) {
    1 / (1 + exp(-(x+centerX)*a))-0.5
  }
  x <- seq(-10,10, 0.01)
  win.graph();plot(x,sigmoid(x),col='orange')
}

if(FALSE){
  {"查看每个基因表达数据的离群值"}
  library("outliers")
  win.graph()
  fix(data_signatire)
  boxplot(as.matrix(t(data_signatire)))
  ##使用5%，95% percentile代替离群值
  
}
if(FALSE){
  {"使用Matlab进行处理"}
  library(R.matlab)  
  
}
##
if(FALSE){
  {"确定最佳的聚类数目"}
  library(mclust)
  m_clust=mclust()
}
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
#res.km = eclust(t(data_signatire), "hclust")
# Gap statistic plot
#win.graph();fviz_gap_stat(res.km$gap_stat)
# #下面的R代码生成Silhouette plot和分层聚类散点图。
# win.graph();fviz_silhouette(res.hc) # silhouette plot
# win.graph();fviz_cluster(res.hc) # scatter plot

fix(data_signatire)
#data_signatire=scale(data_signatire)
# fix(data_signatire)
# max(data_signatire)
# min(data_signatire)
# median(data_signatire)
# mean(data_signatire)
# ##提取聚类轮廓图
# sil = silhouette(res.km$cluster, dist(data_signatire))
# rownames(sil) = rownames(USArrests)
# head(sil[, 1:3])

##删除离群数据
data_signatire=na.omit(data_signatire)
fix(data_signatire)
#data_signatire=data_signatire[,-which(colnames(data_signatire)=="TCGA.HC.7738.01")]

#http://blog.csdn.net/qazplm12_3/article/details/74516312  绘制热图前数据调整
output_comparision=data.frame()
out_put_dir="D:\\min数据\\tcga\\prad_tcga\\analyse\\"
k=5
goHeatMap=function(cluster_method="ward.D",dist_method="euclidean",k=5){
  
  myclust<-function(x){
    hclust(x,method=cluster_method)
  }
  myDist<-function(x){
    dist(x,method = dist_method)
  }
  library("RColorBrewer")
  library(gplots)
  
  ###初步聚类计算cophenetic距离
  hc=hclust(d=dist((t(data_signatire)),method = dist_method),
            method = cluster_method)
  
  summary(hc)
  ##进行cutree以及分组成标准格式列表
  hccut=as.data.frame(cutree(tree = hc,k=k))
  colnames(hccut)="group_index"
  hccut[,"SampleIndex"]=rownames(hccut)
  library(sqldf)
  r=sqldf("select count(`group_index`) from hccut group by `group_index`")
  ##调取某个基因的数据
  getGene="CXCL17"
  SAMPLES=rownames(hccut)
  SAMPLES_EXP=as.data.frame(t(rna_seq_data[getGene,SAMPLES]))
  SAMPLES_EXP[,"SampleIndex"]=rownames(SAMPLES_EXP)
  
  
  hccut=merge(hccut,SAMPLES_EXP,by.x="SampleIndex",by.y = "SampleIndex")
  #apply(X=hccut,1,FUN=function(line){print(line)})
  
  
  ##绘制统计信息
  library(ggplot2)
  library(reshape2)
  data_plot=hccut
  colnames(data_plot)
  ##画图前数据转换
  data_plot[,getGene]=log2(data_plot$CXCL17+1)
  data_plot$group_index=as.character(data_plot$group_index)
  ##开始绘图
  
  
  
  file=paste(out_put_dir,amethods,onemethod,k,"box_plot",".tiff",sep="_")
  tiff(filename = file,width = 1000, height = 800, units = "px")
  # win.graph()
  # p<-ggplot(data=data_plot, aes(x=group_index,y=getGene))
  # p=p+ geom_point()
  # p=p+geom_boxplot();
  # p + geom_jitter(aes(colour = group_index))
  ##data_plot数据的排序median，替换为low2high_rank;生成分组信息
  groupindexs=unique(data_plot$group_index)
  average=c()
  for(onekindindex in groupindexs){
    cat(onekindindex)
    average[onekindindex]=mean(data_plot[which(data_plot$group_index==onekindindex),getGene])
    #print(data_plot[which(data_plot$group_index==onekindindex),getGene])
  }
  average=as.data.frame(sort(average))
  average[,"group_index"]=rownames(average)
  index=1
  for(one in c(1:nrow(average))){
    average[one,"rowname_rank"]=paste("low2high",index,sep="_");index=index+1;
    average[one,"group_count"]=length(data_plot[which(data_plot$group_index==average[one,"group_index"]),getGene])
  }
  rownames(average)=average$rowname
  ##分组信息填充完毕，merge到data_plot
  data_plot=merge(x=data_plot,y=average,by="group_index")
  ##按照由小到大排序
  data_plot=data_plot[order(data_plot$rowname_rank),]
  
  my_comparision_matrix=as.data.frame(t(combn(x=unique(data_plot$rowname_rank),2)))
  colnames(my_comparision_matrix)=c("first","second")
  for(onerow in c(1:nrow(my_comparision_matrix))){
    my_comparision_matrix[onerow,"differ"]=as.numeric(stringi::stri_sub(str=my_comparision_matrix[onerow,"second"],from = stringi::stri_length(my_comparision_matrix[onerow,"second"])))-
      as.numeric(stringi::stri_sub(str=my_comparision_matrix[onerow,"first"],from = stringi::stri_length(my_comparision_matrix[onerow,"second"])))
  }
  ##按照差值进行排序
  my_comparision_matrix=my_comparision_matrix[order(my_comparision_matrix$first),]
  my_comparision_matrix=my_comparision_matrix[order(my_comparision_matrix$differ),]
  my_comparisons <- list()
  for(one in c(1:nrow(my_comparision_matrix))){
    oneterm=as.vector(c(as.character(my_comparision_matrix[one,"first"]),as.character(my_comparision_matrix[one,"second"])))
    my_comparisons[[one]]=oneterm
  }
  
  
  library(ggpubr)
  #win.graph()
  q=qplot(x=rowname_rank,y=CXCL17,data=data_plot,geom="boxplot",main = paste("cluster_method",cluster_method,"dist_method",dist_method,sep="_"))
  q=q+ geom_jitter(aes(colour = rowname_rank))
  q=q+geom_text(data = data_plot,aes(label=paste("n=",group_count,sep=""),y=-2))
  q=q+ggpubr::stat_compare_means(comparisons = my_comparisons,
                                 label.y = c(max(data_plot$CXCL17), max(data_plot$CXCL17)+1.5, max(data_plot$CXCL17)+3),
                                 stat_compare_means(label.y = 22),label = c("p.format"))
  print(q)
  dev.off()
  ##绘图结束
  ###输出随机组合的每组的个数，如果大于3
  n=k
  a=NULL
  a=combn(unique(data_plot$group_count),n-2,FUN = function(onea){
    return(c(unique(data_plot$group_count)[which(!unique(data_plot$group_count) %in% onea)],sum(onea)))
  })
  #结果转换为list
  a=as.data.frame(t(a))
  colnames(a)=c("first","second","third")
  for(rowIndex in c(1:nrow(a))){
    cat(unlist(sort(a[rowIndex,])))
    a[rowIndex,]=sort(a[rowIndex,])
  }
  a[,"amethods"]=amethods
  a[,"onemethod"]=onemethod
  a[,"ori_low"]=67;a[,"ori_middle"]=260;a[,"ori_high"]=117
  a[,"first_distance"]=a[,"first"]-a[,"ori_low"]
  a[,"second_distance"]=a[,"second"]-a[,"ori_middle"]
  a[,"third_distance"]=a[,"third"]-a[,"ori_high"]
  output_comparision=rbind(output_comparision,a)
  
  
  dc=cophenetic(hc)#计算系统聚类的cophenetic距离,h是hclust()函数生成的对象
  cor=cor(dist((t(data_signatire)),method = dist_method),dc)#d是dist()函数的距离，dc是cophenetic距离
  
  
  
  #win.graph()
  file=paste(out_put_dir,amethods,onemethod,k,"heatmap",".tiff",sep="_")
  tiff(filename = file,width = 1000, height = 1000, units = "px")
  #win.graph()
  hv=heatmap.2(as.matrix(data_signatire), col=rev( colorRampPalette(brewer.pal(11, "RdBu"))(256)[c((1+0):(128-20),(128+20):(256-0))]
  ), scale="row",
  distfun = myDist,
  hclustfun = myclust,
  key=TRUE, keysize=1, symkey=TRUE, density.info="histogram", trace="none", cexRow=1,
  Rowv=FALSE, Colv=TRUE,dendrogram = "column",lmat=rbind(c(0,3,0),c(2,1,0),c(0,4,0)) ,lhei=c(3,10,3),lwid=c(2,9,2)
  ,key.title ="",cexCol = 1#,densadj = 0.25
  ,sepcolor="white"
  ,main = paste("cluster_method",cluster_method,"dist_method",dist_method,"cophenetic:",cor,sep="_")
  )
  dev.off()
  hv
}
tryCatch({
  devices=dev.list()
  for(onedev in devices){
    dev.off()
  }  
},finally = {
  
})



##http://blog.sciencenet.cn/blog-651374-988817.html注意参数的配用
Cluster_Method<-c( "ward.D","ward.D2","single","complete","average" ,"mcquitty","median","centroid")
Dist_Methods<-  c("euclidean"
                  #, "maximum", "manhattan", "canberra", "binary" ,"minkowski"
)
for(amethods in  Cluster_Method){
  for(onemethod in Dist_Methods){
    hv=goHeatMap(amethods,onemethod)
  }
}
write.csv(x=output_comparision,file = paste(out_put_dir,"comparision.csv",sep=""))
shell.exec(out_put_dir)







###开始无图化的层次聚类带有
Cluster_Method<-c( "ward.D","ward.D2","single","complete","average" ,"mcquitty","median","centroid")
Dist_Methods<-c("euclidean", "maximum", "manhattan", "canberra", "binary" ,"minkowski")
for(amethods in  Dist_Methods){
  for(onemethod in Cluster_Method){
    for(k in c(3:3)){
      hc=hclust(d=dist((t(data_signatire)),method = amethods),
                method = onemethod)
      hv=cutree(hc,k=k)
      ##
      dc=cophenetic(hc)#计算系统聚类的cophenetic距离,h是hclust()函数生成的对象
      cor=cor(dist((t(data_signatire)),method = amethods),dc)#d是dist()函数的距离，dc是cophenetic距离
      print(cor)
      #通常认为该相关系数越接近1，说明聚类方法就越好
      ##
      
      out_put_dir="D:\\min数据\\tcga\\prad_tcga\\analyse\\"
      file=paste(out_put_dir,amethods,onemethod,k,".tiff",sep="_")
      tiff(filename = file,width = 1000, height = 800, units = "px")
      #win.graph()
      plot(hc,main = paste(amethods,onemethod,k,"cophenetic=",cor,sep="_"))
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
      print(paste(cor(x=c(1:nrow(t(data_signatire))),y=hc$order,method = "spearman"),amethods,onemethod,k
                  ,"cophenetic距离是：",cor,
                  sep="_"))
      #print(paste("距离","_",k,amethods,onemethod,dist(matr)))
      
      ###对聚类树进行剪枝
      hv=cutree(hc,k=3)
      hv
    }
  }
}
#####使用hclust函数结合不同的参数和dist函数和不同的参数进行聚类，和cutree,3
##以便于查看不同的样本编号
hc=hclust(d=dist(((data_signatire)),method = "euclidean"),
          method ="ward.D2")
hv=cutree(hc,k=k)
plot(hc)


rna_seq_data["CXCL17",]

shell.exec(out_put_dir)





######评价聚类效果
require(cluster)
library(factoextra)
par(mfrow=c(1,3))
for(onemethod in c("silhouette","wss","gap_stat")){
  r1=fviz_nbclust(t(data_signatire),
                  hcut,
                  method = onemethod
  )
  win.graph();plot(r1)
}

