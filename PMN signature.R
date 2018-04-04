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
ori_low=67;ori_middle=260;ori_high=171##CXCL17
ori_low=118;ori_middle=206;ori_high=174##cxcl15
getGene0="CXCL17"
getGene0="CXCL5"
k_n0=3
sigNature="PMN"

##读取数据ran_seq_median
path="D:\\min数据\\tcga\\prad_tcga\\prad_tcga\\"
out_put_dir=paste("D:\\min数据\\tcga\\prad_tcga\\analyse\\",getGene0,"\\",sigNature,"\\",sep="")
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
if(!dir.exists(paste(out_put_dir,"r_heriachical\\",k_n0,"\\",sep=""))){
  dir.create(recursive=TRUE,paste(out_put_dir,"r_heriachical\\",k_n0,"\\",sep=""))
}
out_put_dir=paste(out_put_dir,"r_heriachical\\",k_n0,"\\",sep="")
shell.exec(paste(out_put_dir,sep=""))




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
#Mo##geneSignatures=c("CD14","IL4R","PTPRC","ITGAM","LOC116196","ARG1","IL10","CD40","CD32","CD163","FCER2","CD200R1","PDCD1LG2","CD68","CSF1R","HLA-DRA","LY75","LOC90262","CCL2","JM2")  ##high=171,mid=260,low=67
#reg=paste(geneSignatures,collapse = "|^")
#data_signatire=rna_seq_data[which(!is.na(stringi::stri_match(rna_seq_data$Hugo_Symbol,regex = reg))),]
#rna_seq_data=na.omit(rna_seq_data)
geneSignatures=c("CEACAM8","S100A8","S100A9","CXCR2","ID01","CSF2","CSF3","PDCD1","TNF","PTPRC","CD38","CD80","CD274","STAT1","CXCR4","CD34","ENTPD1","TGFB1","IL4R","STAT3","STAT5A","FUT4","ITGAX","ITGAM","CD14","CSF1R","CD33","TLR4","CXCL12","TLR3","ANPEP")##PMNSignature
ifelse(useSymbol,{
  rna_seq_data=rna_seq_data[-which(is.na(rna_seq_data$Hugo_Symbol)),]
  rownames(rna_seq_data)<-rna_seq_data$Hugo_Symbol;
  "基因symbol应用完毕"
},{
  rownames(rna_seq_data)=rna_seq_data$Entrez_Gene_Id;
  "基因ID应用完毕"
})

#fix(rna_seq_data)
##数据格式矫正
#rownames(data_signatire)=data_signatire$Hugo_Symbol

data_signatire=rna_seq_data[geneSignatures,]

fix(data_signatire)
tryCatch({
  data_signatire=data_signatire[,-which(colnames(data_signatire)=="Entrez_Gene_Id")]
  data_signatire=data_signatire[,-which(colnames(data_signatire)=="Hugo_Symbol")]
},finally = {
  "全部去除完毕!"
  write.csv(x=data_signatire,file = paste(out_put_dir,"tcga_signature.csv",sep=""))
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
  secondMin=min(matrixStats::colMins(as.matrix(data_signatire))[which(matrixStats::colMins(as.matrix(data_signatire))!=0)])
  print(secondMin)
  data_signatire=log2(data_signatire+secondMin)  
  # write.csv(x=t(data_signatire),file = paste("D:\\min数据\\tcga\\prad_tcga\\analyse\\spss\\log2_true.csv",sep=""))
  #shell.exec("D:\\min数据\\tcga\\prad_tcga\\analyse\\spss\\")
  para_List[length(para_List)+1]="log2"
}
#!!!查看每个基因的数值分布

if(TRUE){
  data_signatire=as.data.frame(t(scale(t(data_signatire))))
  para_List[length(para_List)+1]="rowWiseScale"
}

if(FALSE){
  ###rowMeans进行中心化
  data_signatire=data_signatire-rowMedians(as.matrix(data_signatire))
  para_List[length(para_List)+1]="rowMeansCenteralized"
}

library(lattice)
curve=densityplot(unlist(data_signatire));win.graph(width = 15,height =10);plot(curve);
fix(data_signatire)
for(rowIndex in c(1:nrow(data_signatire))){
  data_signatire[rowIndex,]
  print(densityplot(unlist(data_signatire[rowIndex,])));win.graph()
}



if(FALSE){
  {"确定最佳的聚类数目"}
  library(mclust)
  m_clust=mclust()
}

library(factoextra)
library(cluster)


fix(data_signatire)

##删除离群数据
data_signatire=na.omit(data_signatire)
#data_signatire=data_signatire[,-which(colnames(data_signatire)=="TCGA.HC.7738.01")]

#http://blog.csdn.net/qazplm12_3/article/details/74516312  绘制热图前数据调整
output_comparision=data.frame()

goHeatMap=function(d1,cluster_method,dist_method,k_n=k_n0,getGene=getGene0,rna_seq_data_){
  
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
    hc=hclust(d=get_dist((t(d1)),method = dist_method),
              method = cluster_method)
    
    summary(hc)
    ##进行cutree以及分组成标准格式列表
    hccut=as.data.frame(cutree(tree = hc,k=k_n));#win.graph();plot(cutree(hc,k=3));
    colnames(hccut)="group_index"
    hccut[,"SampleIndex"]=rownames(hccut)
    library(sqldf)
    r=sqldf("select count(`group_index`) from hccut group by `group_index`")
    ##调取某个基因的数据
    SAMPLES=rownames(hccut)
    SAMPLES_EXP=as.data.frame(t(rna_seq_data_[getGene,SAMPLES]))
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
  print(paste("提取数据基因",getGene,sep=""))
  EDA(log2(data_plot[,getGene]+secondMin));dev.off()
  data_plot[,getGene]=log2(data_plot[,getGene]+secondMin)
  data_plot$group_index=as.character(data_plot$group_index)
  ##开始绘图
  
  
  if(TRUE){
    file=paste(out_put_dir,amethods,onemethod,k_n,paste(para_List,collapse = "_"),"box_plot",".tiff",sep="_")
    tiff(filename = file,width = 1000, height = 800, units = "px")
    #win.graph()
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
    q=qplot(x=rowname_rank,y=data_plot[,getGene],data=data_plot,geom="boxplot",main = paste("cluster_method",cluster_method,"dist_method",dist_method,sep="_")
            ,outlier.colour = "black")
    q=q+ geom_jitter(aes(colour = rowname_rank))
    q=q+geom_text(data = data_plot,aes(label=paste("n=",group_count,sep=""),y=min(data_plot[,getGene])-0.5))
    q=q+stat_compare_means(comparisons = my_comparisons,paired = FALSE,label = "pb.format",
                           hide.ns = FALSE,symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")),
                           label.y = c(seq(max(data_plot[,getGene]),max(data_plot[,getGene])+0.75*k_n,0.75)[c(1:k_n)]),
                           method = "t.test")
    q=q+stat_compare_means(label.y= max(c(seq(max(data_plot[,getGene]),max(data_plot[,getGene])+0.75*k_n,0.75)[c(1:k_n)]))+log10(max(c(seq(max(data_plot[,getGene]),max(data_plot[,getGene])+0.75*k_n,0.75)[c(1:k_n)]))),
                           method = "anova")
    
    
    plot(q)
    dev.off()
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
  a[,"a_clust_methods"]=amethods
  a[,"one_dist_method"]=onemethod
  a[,"ori_low"]=ori_low;a[,"ori_middle"]=ori_middle;a[,"ori_high"]=ori_high
  a[,"first_distance"]=a[,"first"]-a[,"ori_low"]
  a[,"second_distance"]=a[,"second"]-a[,"ori_high"]
  a[,"third_distance"]=a[,"third"]-a[,"ori_middle"]
  output_comparision=rbind(output_comparision,a)
  library(readr)
  if(TRUE){
    readr::write_csv(x=output_comparision,path= paste(out_put_dir,"comparision.csv",sep=""),append = TRUE,col_names = TRUE)  
  }
  
  cor="NULL"
  if(FALSE){
    dc=cophenetic(hc)#计算系统聚类的cophenetic距离,h是hclust()函数生成的对象
    cor=cor(get_dist((t(d1)),method = dist_method),dc)#d是dist()函数的距离，dc是cophenetic距离
  }
  
  
  
  ##绘制热图
  if(TRUE){
    file=paste(out_put_dir,amethods,onemethod,k_n,paste(para_List,collapse = "_"),"heatmap",".tiff",sep="_");print(file)
    tiff(filename = file,width = 1000, height = 1000, units = "px")
    ##计算break
    PERCENTILE=0.005;lowQ=as.numeric(quantile(unlist(d1),PERCENTILE));highQ=as.numeric(quantile(unlist(d1),1-PERCENTILE))
    BREAKS=c(min(d1)-0.01,seq(lowQ,highQ,0.005),max(d1)+0.01)
    ##主要的绘图函数
    hm=heatmap.2(x=as.matrix(d1),col=colorRampPalette(c("blue", "white", "red"))(length(BREAKS)-1), distfun=myDist,hclustfun=myclust,
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
                 main = paste("onemethod",onemethod,"amethods",amethods,"ak_n",ak_n,sep="_"),
                 lmat=rbind(c(0,3,0),c(0,1,2),c(0,4,0)) ,lhei=c(3,10,2),lwid=c(1,9,1)
    )
    
    dev.off()
  }
  #hv
}


for(onedev in dev.list()){
  dev.off()
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
      hv=goHeatMap(data_signatire,amethods,onemethod
                   ,k_n=ak_n,rna_seq_data_=rna_seq_data
      )
    }
  }
}


for(one in dev.list()){
  dev.off()
}

shell.exec(out_put_dir)






if(FALSE){
  ###开始无图化的层次聚类带有
  Cluster_Method<-c( "ward.D","ward.D2","single","complete","average" ,"mcquitty","median","centroid")
  Dist_Methods<-c("euclidean", "maximum", "manhattan", "canberra", "binary" ,"minkowski")
  for(amethods in  Dist_Methods){
    for(onemethod in Cluster_Method){
      for(k in c(3:3)){
        hc=hclust(d=dist((t(data_signatire)),method = amethods),
                  method = onemethod)
        hv=cutree(hc,k=k_n)
        
        ##
        dc=cophenetic(hc)#计算系统聚类的cophenetic距离,h是hclust()函数生成的对象
        cor=cor(get_dist((t(data_signatire)),method = amethods),dc)#d是dist()函数的距离，dc是cophenetic距离
        print(cor)
        #通常认为该相关系数越接近1，说明聚类方法就越好
        ##
        
        out_put_dir="D:\\min数据\\tcga\\prad_tcga\\analyse\\"
        file=paste(out_put_dir,amethods,onemethod,k,paste(para_List,collapse = "_"),".tiff",sep="_")
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
  hv=cutree(hc,k=k_n)
  plot(hc)
  
  
  rna_seq_data[getGene,]
  
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
    win.graph(width = 12,height =10);plot(r1)
  }
}