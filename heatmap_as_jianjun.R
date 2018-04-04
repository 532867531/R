library(matrixStats)
library(NMF)
##读取数据ran_seq_median
path="D:\\min数据\\tcga\\prad_tcga\\prad_tcga\\"
out_put_dir="D:\\min数据\\tcga\\prad_tcga\\analyse\\"
file=paste(path,"data_RNA_Seq_v2_expression_median.RDS",sep="")
useSymbol = TRUE
#file=paste(path,"data_RNA_Seq_v2_mRNA_median_Zscores.RDS",sep="")
if(file.exists(file)){
  rna_seq_data=readRDS(file=file)
}else{
  rna_data=read.csv(file=sub(pattern =".RDS",replacement = ".txt",x=file),sep="")
  saveRDS(rna_data,file=file)
  rna_seq_data=readRDS(file=file)
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

#数据格式矫正
#rownames(data_signatire)=data_signatire$Hugo_Symbol

data_signatire=rna_seq_data[geneSignatures,]

fix(data_signatire)
tryCatch({
  data_signatire=data_signatire[,-which(colnames(data_signatire)=="Entrez_Gene_Id")]
  data_signatire=data_signatire[,-which(colnames(data_signatire)=="Hugo_Symbol")]
},finally = {
  "全部去除完毕!"
  #write.csv(x=data_signatire,file = paste(out_put_dir,"tcga_signature.csv",sep=""))
})

##查看未找到的名称
rownames(data_signatire)=geneSignatures
#write.csv(x=t(data_signatire),file = paste("D:\\min数据\\tcga\\prad_tcga\\analyse\\spss\\log2.csv",sep=""))
#shell.exec("D:\\min数据\\tcga\\prad_tcga\\analyse\\spss\\")
isLog2=TRUE
if(isLog2){
  {"是否进行log2处理"}
 data_signatire=log2(data_signatire+1) 
}

isrowScale=TRUE
if(isrowScale){
  {"是否手动执行行scale()，基因表达在不同的样品之间"}
  data_signatire<-data_signatire-rowMedians(as.matrix(data_signatire))
}
fix(data_signatire)
print(summary(is.na(rownames(data_signatire))))

library(factoextra)
library(cluster)

data_signatire=na.omit(data_signatire)
fix(data_signatire)
#data_signatire=data_signatire[,-which(colnames(data_signatire)=="TCGA.HC.7738.01")]

#http://blog.csdn.net/qazplm12_3/article/details/74516312  绘制热图前数据调整
output_comparision=data.frame()
out_put_dir="D:\\min数据\\tcga\\prad_tcga\\analyse\\"

goHeatMap=function(cluster_method="ward.D",dist_method="euclidean",k_n=3){
  
  myclust<-function(x){
    hclust(x,method=cluster_method)
  }
  myDist<-function(x){
    get_dist(x,method = dist_method)
  }
  library("RColorBrewer")
  library(gplots)
  
  ###初步聚类计算cophenetic距离
  if(TRUE){
    hc=hclust(d=get_dist((t(data_signatire)),method = dist_method),
              method = cluster_method)
    
    summary(hc)
    ##进行cutree以及分组成标准格式列表
    hccut=as.data.frame(cutree(tree = hc,k=k_n));#win.graph();plot(cutree(hc,k=3));
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
  }
  
  ##绘制统计信息
  library(ggplot2)
  library(reshape2)
  data_plot=hccut
  colnames(data_plot)
  ##画图前数据转换
  data_plot[,getGene]=log2(data_plot$CXCL17+1)
  data_plot$group_index=as.character(data_plot$group_index)
  ##开始绘BOX图
  
  
  if(TRUE){
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
    if(TRUE){
     # print(q)
    }
    #dev.off()
  }
  ##绘图结束
  ###输出随机组合的每组的个数，如果大于3
  
  a=NULL
  input_1=na.omit(unique(data_plot$group_index))
  a=combn(input_1,k_n-2,FUN = function(onea){
    selectIndex=onea;print(paste("selectIndex:",length(selectIndex),selectIndex))
    unselectIndex=unique(data_plot$group_index[which(!data_plot$group_index %in% onea)]);print(paste("unselectedIndex",length(unselectIndex),unselectIndex))
    
    
    # print(c(unique(data_plot[which(data_plot$group_index %in% unselectIndex),"group_count"])
    #         ,sum(unique(data_plot[which(data_plot$group_index %in% selectIndex),"group_count"]))))
    
    ##未选中的求和
    unselected_Group_count=vector()
    for(one in unselectIndex){
      unselected_Group_count[length(unselected_Group_count)+1]=unique(data_plot[which(data_plot$group_index ==one),"group_count"])
    }
    selected_Group_count=vector()
    for(one in selectIndex){
      selected_Group_count[length(selected_Group_count)+1]=unique(data_plot[which(data_plot$group_index ==one),"group_count"])
    }
    
    return(c(sort(unselected_Group_count)
             ,sum(selected_Group_count)
    ))},simplify = TRUE)
  #结果转换为list
  a=as.data.frame(t(a))
  a=unique(a)
  colnames(a)=c("first","second","third")
  for(rowIndex in c(1:nrow(a))){
    print(unlist(sort(a[rowIndex,])))
    a[rowIndex,]=sort(a[rowIndex,])
    a[rowIndex,"k_n"]=k_n
  }
  a[,"amethods"]=amethods
  a[,"onemethod"]=onemethod
  a[,"ori_low"]=67;a[,"ori_middle"]=260;a[,"ori_high"]=171
  a[,"first_distance"]=a[,"first"]-a[,"ori_low"]
  a[,"second_distance"]=a[,"second"]-a[,"ori_high"]
  a[,"third_distance"]=a[,"third"]-a[,"ori_middle"]
  output_comparision=rbind(output_comparision,a)
  library(readr)
  if(TRUE){
    readr::write_csv(x=output_comparision,path= paste(out_put_dir,"comparision.csv",sep=""),append = TRUE,col_names = TRUE)  
  }
  
 ####开始绘制热图
  win.graph()
  aheatmap(x=data_signatire,color = colorRampPalette(c("blue", "white", "red"))(802),
           distfun=myDist,
           breaks=c(min(data_signatire)-0.05,seq(-5,5,0.05),max(data_signatire)+0.05),
           hclustfun=myclust,
           scale="none",
           labCol=NA,legend=T)
  }
  
tryCatch({
  # devices=dev.list()
  # for(onedev in devices){
  #   dev.off()
  # }  
},finally = {
  
})



##http://blog.sciencenet.cn/blog-651374-988817.html注意参数的配用
Cluster_Method<-c( "ward.D","ward.D2","single","complete","average" ,"mcquitty","median","centroid")
Dist_Methods<-  c("euclidean"
                  #, "maximum", "manhattan", 
                  #"canberra", 
                 # "binary", "minkowski", "pearson", "spearman","kendall"
)



  for(amethods in  Cluster_Method){
    for(onemethod in Dist_Methods){
      print(paste("onemethod",onemethod,"amethods",amethods,sep="_"))
      hv=goHeatMap(amethods,onemethod)
    }
  }


shell.exec(out_put_dir)

