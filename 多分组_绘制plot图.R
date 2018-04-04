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
ori_low=37;ori_middle=46;ori_high=34##cxcl15
getGene0="CXCL17"
getGene0="CXCL5"
#getGene0="PML"
ref_Gene="ZBTB7A"
ref_Gene1="PTEN"
k_n0=3

##读取数据ran_seq_median
##path="D:\\min数据\\tcga\\prad_tcga\\prad_tcga\\"
path="D:\\min数据\\tcga\\prad_su2c_2015\\prad_su2c_2015\\"
out_put_dir=paste(path,getGene0,"\\",sep="")
label_Func_readInRnaSeqFile=function(){print("读入RNAseq表达数据！")}
file=paste(path,list.files(path = path,pattern = "data.*?RNA_Seq.*?expression.*?median.*?txt"),sep="")
RDSfile=sub(pattern = "txt",replacement = "RDS",perl=TRUE,x=file)
useSymbol = TRUE
#file=paste(path,"data_RNA_Seq_v2_mRNA_median_Zscores.RDS",sep="")
if(file.exists(RDSfile)){
  rna_seq_data=readRDS(file=RDSfile)
}else{
  rna_data=read.csv(file = file,sep = "\t")
  saveRDS(rna_data,file=RDSfile)
  rna_seq_data=readRDS(file=RDSfile)
}




label_Func_FileCNA=function(){print("读取CNA文件！")}
CNAfile=paste(path,list.files(path = path,pattern = "data_CNA.*txt"),sep="")
CNARDSfile=sub(pattern = "txt",replacement = "RDS",perl=TRUE,x=CNAfile)
if(file.exists(CNARDSfile)){
  data_CNA=readRDS(file=CNARDSfile)
}else{
  data_CNA=read.csv(file = CNAfile,sep = "\t")
  saveRDS(data_CNA,file=CNARDSfile)
  data_CNA=readRDS(file=CNARDSfile)
}
fix(data_CNA)



label_Func_OutputDir=function(){print("打开输出文件文件夹！")}
if(!dir.exists(paste(out_put_dir,"r_heriachical\\",k_n0,"\\",sep=""))){
  dir.create(recursive=TRUE,paste(out_put_dir,"r_heriachical\\",k_n0,"\\",sep=""))
}
out_put_dir=paste(out_put_dir,"r_heriachical\\",k_n0,"\\",sep="")
shell.exec(path)
shell.exec(paste(out_put_dir,sep=""))


library(dplyr)
library(sqldf)
##数据初步矫正去重，认为筛选的数据不在重复的基因中;或者进行median、average处理
ifelse(useSymbol,{
  library(matrixStats)
  library(NMF)
  Samples=colnames(rna_seq_data)[which(colnames(rna_seq_data)!="Hugo_Symbol")]
  sql=paste("select",samples,"from rna_seq_data")
  sqldf("select ")
  d=duplicated(rna_seq_data$Hugo_Symbol);
  dup_Genes=as.character(rna_seq_data[d,"Hugo_Symbol"])
  Samples=colnames(rna_seq_data)[which(colnames(rna_seq_data)!="Hugo_Symbol")]
  for (oneGene in dup_Genes) {
    for(oneSample in Samples){
      rna_seq_data[which(rna_seq_data$Hugo_Symbol==oneGene),oneSample]=max(na.omit(rna_seq_data[rna_seq_data$Hugo_Symbol==oneGene,oneSample]))
    }
  }
  
  
  rna_seq_data=rna_seq_data[!d,];
  ("完成去重")
},
"跳过去重")
fix(rna_seq_data)
ifelse(useSymbol,{
  rna_seq_data=na.omit(rna_seq_data)
  rownames(rna_seq_data)<-rna_seq_data$Hugo_Symbol;
  "基因symbol应用完毕"
},{
  rownames(rna_seq_data)=rna_seq_data$Entrez_Gene_Id;
  "基因ID应用完毕"
})



##提取zbtb7a
data_refGene=rna_seq_data[ref_Gene,]
if("Hugo_Symbol" %in% colnames(rna_seq_data)){
  data_refGene=data_refGene[,-which(colnames(data_refGene)=="Hugo_Symbol")] 
}
fix(data_refGene)
##进行参考基因的分组
criteria_refGene.high=as.numeric(quantile(data_refGene,1-0.205))
criteria_refGene.low=as.numeric(quantile(data_refGene,1-0.795))

Data_refGene.high=data_refGene[,which(data_refGene>=criteria_refGene.high)]
Data_refGene.low=data_refGene[,which(data_refGene<=criteria_refGene.low)]

Sample_refGene.high=colnames(Data_refGene.high)
Sample_refGene.high=Sample_refGene.high[which(Sample_refGene.high %in% colnames(rna_seq_data))]
Sample_refGene.low=colnames(Data_refGene.low)
Sample_refGene.low=Sample_refGene.low[which(Sample_refGene.low %in% colnames(rna_seq_data))]


##开始获取signature基因
data_signatire=rna_seq_data[getGene0,]

tryCatch({
  if("Entrez_Gene_Id" %in% colnames(rna_seq_data)){
    data_signatire=data_signatire[,-which(colnames(data_signatire)=="Entrez_Gene_Id")]  
  }
  if("Hugo_Symbol" %in% colnames(rna_seq_data)){
    data_signatire=data_signatire[,-which(colnames(data_signatire)=="Hugo_Symbol")]
  }
},finally = {
  "全部去除完毕!"
})

##查看未找到的名称
rownames(data_signatire)=getGene0
fix(data_signatire)
print(summary(is.na(rownames(data_signatire))))
##寻找criteria基因的样品进行分组,绘图时加标签进行分组
lable_Func_preProcess=function(){print("对signature基因进行分组处理！")}
Sample_Alt_Normal=data_CNA[which(data_CNA$Hugo_Symbol==ref_Gene1),]
rownames(Sample_Alt_Normal)=Sample_Alt_Normal$Hugo_Symbol
if( "Hugo_Symbol"%in%colnames(Sample_Alt_Normal)){
  Sample_Alt_Normal=Sample_Alt_Normal[,-which(colnames(Sample_Alt_Normal)=="Hugo_Symbol")]
}
if("Entrez_Gene_Id"%in%colnames(Sample_Alt_Normal)){
  Sample_Alt_Normal=Sample_Alt_Normal[,-which(colnames(Sample_Alt_Normal)=="Entrez_Gene_Id")]
}

Sample_Alt=rownames(as.data.frame(t(Sample_Alt_Normal)[which(Sample_Alt_Normal==-2),]))
Sample_Alt=Sample_Alt[which(Sample_Alt %in% colnames(rna_seq_data))]
Sample_Normal=rownames(as.data.frame(t(Sample_Alt_Normal)[which(abs(Sample_Alt_Normal)<2),]))
Sample_Normal=Sample_Normal[which(Sample_Normal %in% colnames(rna_seq_data))]
###画图
data_plot=as.data.frame(t(data_signatire))
data_plot[Sample_Alt,"ALT"]=paste(ref_Gene1,"ALT",sep = "_")
data_plot[Sample_Normal,"ALT"]=paste(ref_Gene1,"Normal",sep = "_")
data_plot[Sample_refGene.low,paste(ref_Gene,"Level",sep="_")]=paste(ref_Gene,"Low",sep = "_")
data_plot[Sample_refGene.high,paste(ref_Gene,"Level",sep="_")]=paste(ref_Gene,"High",sep = "_")
library(ggpubr);library(ggrepel)
secondMin=min(data_plot[,getGene0][which(data_plot[,getGene0]!=0)])
data_plot[,getGene0]=log2(data_plot[,getGene0]+secondMin)
write.csv(x=data_plot,file = paste(out_put_dir,"tcga_signature.csv",sep=""))


#可视化DEPDC1基因表达谱
ggboxplot(data=data_plot,x="ALT",y="CXCL5", 
     color = "ALT",add = "jitter", legend="none")+ 
  
  geom_hline(yintercept = mean(data_plot[,getGene0]), linetype=2)+# Add horizontal line at base mean 
  stat_compare_means(method = "t.test") # Add global annova p-value 
win.graph();plot(q)
