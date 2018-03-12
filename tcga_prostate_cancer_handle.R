##读取数据ran_seq_median
path="D:\\min数据\\tcga\\prad_tcga\\prad_tcga\\"
file=paste(path,"data_RNA_Seq_v2_expression_median.RDS",sep="")
if(file.exists(file)){
  rna_seq_data=readRDS(file=file)
}else{
  rna_data=read.csv(file=file,sep="")
}

index = duplicated(rna_seq_data$Hugo_Symbol)
rna_seq_data_single = rna_seq_data[!index,]
library(dplyr)
dplyr::group_by(rna_seq_data,)