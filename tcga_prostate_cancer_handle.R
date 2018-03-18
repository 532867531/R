##读取数据ran_seq_median
path="D:\\min数据\\tcga\\prad_tcga\\prad_tcga\\"
file=paste(path,"data_RNA_Seq_v2_expression_median.RDS",sep="")
if(file.exists(file)){
  rna_seq_data=readRDS(file=file)
}else{
  rna_data=read.csv(file=file,sep="")
}


rna_seq_data=rna_seq_data[,-which(colnames(rna_seq_data)=="Entrez_Gene_Id")]
library(dplyr)
rna_seq_data=tbl_df(rna_seq_data)
##选取基因
geneSignatures=c("CXCR4","CXCR2","ITGAM","ITGAX","ANPEP","CD14","FUT4","CD33","CD34","CD38","ENTPD1","PTPRC","CEACAM8","CD80","CSF1R","IL4R","CSF3","CSF2","CXCL8","TNF","CXCL12","CSF1R","S100A8","S100A9","STAT1","STAT3","STAT5A","ARG1","NOS2","CD274","TLR3","TLR4","TGFB1","IL10","IDO1","PDCD1")
reg=paste(geneSignatures,collapse = "|")
data_signatire=rna_seq_data[which(!is.na(stringi::stri_match(rna_seq_data$Hugo_Symbol,regex = reg))),]
##数据格式矫正

