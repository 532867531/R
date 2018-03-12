library(DEGseq)
##先使用杨翔师兄的seq数据进行测试
path="G:\\YX-两次RNAseq\\T1_VS_C1\\Differentially_Expressed_Gene\\"
path="G:\\YX-两次RNAseq\\T2_VS_C2\\Differentially_Expressed_Gene\\"
t1c1=paste(path,"T1_C1.anno.csv",sep="")
t1c1=paste(path,"T2_C2.anno.csv",sep="")
t1c1_data=read.csv(file = t1c1)
colnames(t1c1_data)[1]="Gene"
colnames(t1c1_data)


t1c1_no=as.data.frame(t1c1_data[which(t1c1_data$Significant=="no"),])
##摘出官方结果
outputDir_=paste(path,"test\\",sep="")
DEGseq::DEGexp(
  #geneExpMatrix1 = as.matrix(t1c1_data[,c("Gene","T1_count")]),
  #geneExpMatrix2 = as.matrix(t1c1_data[,c("Gene","C1_count")]),
  geneExpMatrix1 = as.matrix(t1c1_data[,c("Gene","T2_count")]),
  geneExpMatrix2 = as.matrix(t1c1_data[,c("Gene","C2_count")]),
               #geneCol1 = 1,expCol1 = 2,
               groupLabel1 = "T1_LABLE",
               #geneCol2=1,expCol2 = 2,
               groupLabel2 = "C1_LABLE",
              rawCount = TRUE,##似乎影响不大
              normalMethod = "none",
              method="MARS",
              qValue = 0.05,
              foldChange = 2,
              thresholdKind = 5,
              outputDir = outputDir_
               )
##读取我们的结果
our_r=read.csv(file=paste(outputDir_,"output_score.txt",sep=""),sep="\t")
##生成对比列表
compare_r=merge(t1c1_data,our_r,by.x="Gene",by.y="GeneNames")
compare_r=compare_r[,c("Gene","Log2FoldChange","log2.Fold_change.","log2.Fold_change..normalized")]
##看吧
fix(compare_r)

