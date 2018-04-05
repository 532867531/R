source("function_heatmap_cluster_combn_after_jianjun.R")
EDA=function(x){
  windows()
  par(mfrow=c(3,2))
  hist(x)
  dotchart(x)
  boxplot(x,horizontal = T)
  qqnorm(x);qqline(x)
  mtext("title",outer = TRUE)
  par(mfrow=c(1,1))
}


##定义参数列表
para_List=c()
#title0="t_cell_signature"
title0="PTEN_signature"
getGene0="ALL";
k_n0=3

##读取数据ran_seq_median
path="D:\\min数据\\tcga\\lihc_tcga\\lihc_tcga\\"
out_put_dir=paste(path,getGene0,"\\",sep="")
label_Func_readInRnaSeqFile=function(){print("读入RNAseq表达数据！")}
file=paste(path,list.files(path = path,pattern = "data.*?RNA_Seq.*?expression.*?median.*?txt"),sep="")
#file=paste(path,list.files(path = path,pattern = "data.*?RNA_Seq.*?Zscores.*?txt"),sep="")

RDSfile=sub(pattern = "txt",replacement = "RDS",perl=TRUE,x=file)
useSymbol = TRUE
library(data.table)
#file=paste(path,"data_RNA_Seq_v2_mRNA_median_Zscores.RDS",sep="")
if(file.exists(RDSfile)){
  rna_seq_data=readRDS(file=RDSfile)
}else{
  rna_data=read.csv(file = file,sep = "\t")
  saveRDS(rna_data,file=RDSfile)
  rna_seq_data=readRDS(file=RDSfile)
}
if(data.class(rna_seq_data)!="data.table"){
  rna_seq_data=data.table(rna_seq_data)
}

##rnaSEQ文件探针合并.max,average??
##这里我们使用data.table的.SD方法
receive=winDialogString(as.character(nrow(rna_seq_data)!=length(unique(rna_seq_data$Hugo_Symbol))),default = "发现重复基因名称")
rna_seq_data=rna_seq_data[,lapply(.SD,max),by=.(Hugo_Symbol),.SDcols=colnames(rna_seq_data)[c(2:ncol(rna_seq_data))]]
rna_seq_data=na.omit(as.data.frame(rna_seq_data))
rownames(rna_seq_data)=rna_seq_data$Hugo_Symbol
##选取基因
PTEN_geneSignatures=c("BZS","SLC30A1")##BZS=PTEN,
geneSignatures=PTEN_geneSignatures
#fix(rna_seq_data)
head(rna_seq_data,1)


ifelse(useSymbol,{
  rna_seq_data=na.omit(rna_seq_data)
  rownames(rna_seq_data)<-rna_seq_data$Hugo_Symbol;
  "基因symbol应用完毕"
},{
  rownames(rna_seq_data)=rna_seq_data$Entrez_Gene_Id;
  "基因ID应用完毕"
})

##开始获取signature基因
data_signatire=rna_seq_data[geneSignatures,]

tryCatch({
  if("Entrez_Gene_Id" %in% colnames(data_signatire)){
    data_signatire=data_signatire[,-which(colnames(data_signatire)=="Entrez_Gene_Id")]  
  }
  if("Hugo_Symbol" %in% colnames(data_signatire)){
    data_signatire=data_signatire[,-which(colnames(data_signatire)=="Hugo_Symbol")]
  }
},finally = {
  "全部去除完毕!"
  #  write.csv(x=data_signatire,file = paste(out_put_dir,"tcga_signature.csv",sep=""))
})


tryCatch({
  if("Entrez_Gene_Id" %in% colnames(rna_seq_data)){
    rna_seq_data=rna_seq_data[,-which(colnames(rna_seq_data)=="Entrez_Gene_Id")]  
  }
  if("Hugo_Symbol" %in% colnames(rna_seq_data)){
    rna_seq_data=rna_seq_data[,-which(colnames(rna_seq_data)=="Hugo_Symbol")]
  }
},finally = {
  "全部去除完毕!"
  #  write.csv(x=data_signatire,file = paste(out_put_dir,"tcga_signature.csv",sep=""))
})


##查看未找到的名称
rownames(data_signatire)=geneSignatures
#write.csv(x=t(data_signatire),file = paste("D:\\min数据\\tcga\\prad_tcga\\analyse\\spss\\log2.csv",sep=""))
#shell.exec("D:\\min数据\\tcga\\prad_tcga\\analyse\\spss\\")
fix(data_signatire)
print(summary(is.na(rownames(data_signatire))))

##非常必要加载高性能的包
library(matrixStats)
library(NMF)
if(TRUE){
  {"是否进行log2处理"}
  secondMin=min(as.matrix(data_signatire)[which(data_signatire!=0)])
  #print(secondMin)
  data_signatire=log2(data_signatire+1)  
  # write.csv(x=t(data_signatire),file = paste("D:\\min数据\\tcga\\prad_tcga\\analyse\\spss\\log2_true.csv",sep=""))
  #shell.exec("D:\\min数据\\tcga\\prad_tcga\\analyse\\spss\\")
  para_List[length(para_List)+1]="log2"
}
#!!!查看每个基因的数值分布
library(ggplot2)
library(lattice)
par(mfrow=c(2,1));
apply(data_signatire,MARGIN = 1, function(x){
win.graph();a=densityplot(unlist(x));print(a);
print(mean(unlist(x)))
})
par(mfrow=c(1,1))
if(TRUE){
  data_signatire=as.data.frame(t(scale(t(data_signatire))))
  para_List[length(para_List)+1]="rowWiseScale"
}
apply(data_signatire,MARGIN = 1, function(x){
  win.graph();a=densityplot(unlist(x));print(a);
  print(mean(unlist(x)))
})
if(FALSE){
  ###rowMeans进行中心化
  data_signatire=data_signatire-rowMedians(as.matrix(data_signatire))
  para_List[length(para_List)+1]="rowMeansCenteralized"
}

library(lattice)
curve=densityplot(unlist(data_signatire));win.graph();print(curve);
fix(data_signatire)
for(rowIndex in c(1:nrow(data_signatire))){
  data_signatire[rowIndex,]
  #win.graph();print(densityplot(unlist(data_signatire[rowIndex,])));
}

if(FALSE){
  {"确定最佳的聚类数目"}
  library(mclust)
  m_clust=mclust()
}

library(factoextra)
library(cluster)

label_Func_OutputDir=function(){print("打开输出文件文件夹！")}
if(!dir.exists(paste(out_put_dir,title0,paste(para_List,collapse = "_"),"r_heriachical\\",k_n0,"\\",sep=""))){
  dir.create(recursive=TRUE,paste(out_put_dir,title0,paste(para_List,collapse = "_"),"r_heriachical\\",k_n0,"\\",sep=""))
}
out_put_dir=paste(out_put_dir,title0,paste(para_List,collapse = "_"),"r_heriachical\\",k_n0,"\\",sep="")
shell.exec(path)
shell.exec(paste(out_put_dir,sep=""))
fix(data_signatire)
#######相关系数计算
data2cor=t(rna_seq_data)
##去除标准差为0的变量
data2cor=as.data.frame(data2cor[,which(colSds(data2cor)!=0)])
colSds(data2cor[,geneSignatures])==0
r_cor=cor(data2cor[,geneSignatures],data2cor)
fix(r_cor);win.graph();plot(t(r_cor));
power2_r_cor=data.frame(t(r_cor^2))
win.graph();plot(t(log2_r_cor));
##划分前百分25大的%
t_r_cor=data.frame(t(r_cor))
quantile(power2_r_cor$BZS,0.90)
quantile(power2_r_cor$SLC30A1,0.90)
filtered_gene=rownames(power2_r_cor[intersect(which(power2_r_cor$BZS>quantile(power2_r_cor$BZS,0.90))
                           ,which(power2_r_cor$SLC30A1>quantile(power2_r_cor$SLC30A1,0.90))),])

filtered_r_cor=r_cor[,filtered_gene]
library(ggplotgui)
library(dplyr)
ggplot_shiny(tbl_df(t(filtered_r_cor)))

##计算带P的相关系数

t_r_cor[,"rownames"]=rownames(t_r_cor)
t_r_cor1=apply(t_r_cor,MARGIN = 1,FUN = function(x){
  t_r_cor[x[[3]],"p.BZS"]=cor.test(data2cor[,x[[3]]],data2cor[,"BZS"])$p.value
  t_r_cor[x[[3]],"p.SLC30A1"]=cor.test(data2cor[,x[[3]]],data2cor[,"SLC30A1"])$p.value
  t(t_r_cor[x[[3]],])
})
t_r_cor1=as.data.frame(t(t_r_cor1))
colnames(t_r_cor1)=c("PTEN","SLC30A1","ROWNAMES","P.WITH.PTEN","P.WITH.SLC30A1")
fix(t_r_cor1)
write.csv(x=t_r_cor1,file=paste(out_put_dir,"correlation.csv",sep=""))
library("ggm")
pcor




##删除离群数据
data_signatire=na.omit(data_signatire)
#http://blog.csdn.net/qazplm12_3/article/details/74516312  绘制热图前数据调整
output_comparision=data.frame()

goHeatMap=function(d1=data_signatire,cluster_method,dist_method,k_n,getGene0,rna_seq_data_=rna_seq_data){
  myclust<-function(x){
    hclust(x,method=cluster_method)
  }
  myDist<-function(x){
    get_dist(x,method = dist_method)
  }
  
  
  hc=hclust(d=dist((t(data_signatire)),method = dist_method),
            method = cluster_method)
  
  summary(hc)
  ##进行cutree以及分组成标准格式列表
  hccut=as.data.frame(cutree(tree = hc,k=k_n0))
  colnames(hccut)="group_index"
  hccut[,"SampleIndex"]=rownames(hccut)
  ##cut数目的统计
  library(sqldf)
  r=sqldf("select count(`group_index`),group_index from hccut group by `group_index`")
  rownames(r)=r$group_index
  ##调取某个基因的数据
  getGene="CXCL17"
  SAMPLES=rownames(hccut)
  SAMPLES_EXP=as.data.frame(t(rna_seq_data[getGene,SAMPLES]))
  SAMPLES_EXP[,"SampleIndex"]=rownames(SAMPLES_EXP)
  
  
  hccut=merge(hccut,SAMPLES_EXP,by.x="SampleIndex",by.y = "SampleIndex")
  hccut=merge(hccut,r,by.x="group_index",by.y="group_index")
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
  
  
  
  file=paste(out_put_dir,amethods,onemethod,k_n0,"box_plot",".tiff",sep="_")
  tiff(filename = file,width =1000, height = 800, units = "px")
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
    print(data_plot[which(data_plot$group_index==onekindindex),getGene])
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
  
  
  library("RColorBrewer")
  library(gplots)
  
  ##绘制统计信息
  
  
  ##绘图结束
  ##
  rownames(hccut)=hccut$SampleIndex
  colnames(d1)=paste(colnames(d1),hccut[colnames(d1),4],hccut[colnames(d1),1],sep = "_")
  if(TRUE){
    file=paste(out_put_dir,cluster_method,dist_method,k_n,paste(para_List,collapse = "_"),"heatmap",".svg",sep="_");print(file)
    #tiff(filename = file,width = 1000, height = 5000, units = "px")
    svg(filename = file,width = 60, height = 7, pointsize = 12)
    #win.graph(width=10, height=7,pointsize=8)
    ##计算break
    PERCENTILE=0.005;lowQ=as.numeric(quantile(unlist(d1),PERCENTILE));highQ=as.numeric(quantile(unlist(d1),1-PERCENTILE))
    BREAKS=c(min(d1)-0.01,seq(lowQ,highQ,0.005),max(d1)+0.01)
    heatmap_2=TRUE
    if(heatmap_2){
      ##主要的绘图函数
      hm=heatmap.2(x=as.matrix(d1),col=colorRampPalette(c("blue", "white", "red"))(length(BREAKS)-1), distfun=myDist,hclustfun=myclust,
                   #hm=heatmap.2(x=as.matrix(d1),col=colorRampPalette(c("blue", "white", "red"))(256), distfun=myDist,hclustfun=myclust,
                   scale="none",
                   trace="none",
                   dendrogram = "column",
                   Rowv = FALSE,
                   Colv = TRUE,
                   sepcolor = "white",
                   symkey = TRUE,
                   #breaks = c(min(d1)-0.05,seq(-3,3,0.005),max(d1)+0.05)
                   # breaks = c(min(d1)-0.5,seq(-3,3,0.005),max(d1)+0.5)
                   breaks = BREAKS,
                   main = paste("hclust_method",cluster_method,"dist_method",dist_method,"ak_n",k_n0,sep="_"),
                   lmat=rbind(c(0,3,0),c(0,1,2),c(0,4,0)) ,lhei=c(3,10,2),lwid=c(1,9,1),
                   
      )
    }
    print(hm);dev.off();
  }
  hm
}

##http://blog.sciencenet.cn/blog-651374-988817.html注意参数的配用
Cluster_Method<-c( "ward.D","ward.D2","single","complete","average" ,"mcquitty","median","centroid")
Dist_Methods<-  c("euclidean"
                  #, "maximum", "manhattan", 
                  #"canberra", 
                  #"binary", "minkowski", "pearson", "spearman","kendall"
)


for(ak_n in c(k_n0:k_n0)){
  for(amethods in  Cluster_Method){
    for(onemethod in Dist_Methods){
      print(paste("onemethod",onemethod,"amethods",amethods,"ak_n",ak_n,sep="_"))
      hm=goHeatMap(data_signatire,amethods,onemethod
                   ,k_n=ak_n,rna_seq_data_=rna_seq_data
      )
    }
  }
}

