library(factoextra)
library(cluster)
# 数据准备
# 使用内置的R数据集USArrests
data("USArrests")
# remove any missing value (i.e, NA values for not available)
USArrests = na.omit(USArrests) #view the first 6 rows of the data
head(USArrests, n=6) 
# 显示测试数据示例如下
# 在聚类之前我们可以先进行一些必要的数据检查即数据描述性统计，如平均值、标准差等
desc_stats = data.frame( Min=apply(USArrests, 2, min),#minimum
                         Med=apply(USArrests, 2, median),#median
                         Mean=apply(USArrests, 2, mean),#mean
                         SD=apply(USArrests, 2, sd),#Standard deviation
                         Max=apply(USArrests, 2, max)#maximum
)
desc_stats = round(desc_stats, 1)#保留小数点后一位head(desc_stats)
desc_stats
# 变量有很大的方差及均值时需进行标准化
df = scale(USArrests)




# 数据集群性评估，使用get_clust_tendency()计算Hopkins统计量
res = get_clust_tendency(df, 40, graph = TRUE)
res$hopkins_stat
win.graph();res$plot
##由于k均值聚类需要指定要生成的聚类数量，因此我们将使用函数clusGap()来计算用于估计最优聚类数。函数fviz_gap_stat()用于可视化。
set.seed(123)
## Compute the gap statistic
gap_stat = clusGap(df, FUN = kmeans, nstart = 25, K.max = 10, B = 500)
# Plot the result
win.graph();fviz_gap_stat(gap_stat)
##图中显示最佳为聚成四类（k=4）kmeans进行聚类  kmeans按四组进行聚类，选择25个随机集
km.res = kmeans(df, 4, nstart = 25)
# Visualize clusters using factoextra
win.graph();fviz_cluster(km.res, USArrests)
##提取聚类轮廓图
sil = silhouette(km.res$cluster, dist(df))
rownames(sil) = rownames(USArrests)
head(sil[, 1:3])
##四个cluster的基本信息
# Visualize
fviz_silhouette(sil)
##图片尺寸宽900 dpi较适合微信手机端阅读  图中可以看出有负值，可以通过函数silhouette()确定是哪个观测值
neg_sil_index = which(sil[, "sil_width"] < 0)
sil[neg_sil_index, , drop = FALSE]
##显示为负的观测值


#eclust():增强的聚类分析
#与其他聚类分析包相比，eclust()有以下优点： 简化了聚类分析的工作流程，可以用于计算层次聚类和分区聚类，eclust()自动计算最佳聚类簇数。 自动提供Silhouette plot，可以结合ggplot2绘制优美的图形，使用eclust()的K均值聚类
# Compute k-means
res.km = eclust(df, "kmeans")
# Gap statistic plot
fviz_gap_stat(res.km$gap_stat)
##使用eclust()的层次聚类
# Enhanced hierarchical clustering
res.hc = eclust(df, "hclust") # compute hclust
fviz_dend(res.hc, rect = TRUE) # dendrogam
#下面的R代码生成Silhouette plot和分层聚类散点图。
fviz_silhouette(res.hc) # silhouette plot
fviz_cluster(res.hc) # scatter plot


############################################
############################################
############################################
############################################
############################################
##################下面开始分析我们自己的数据###############
############################################
############################################
############################################
##读取数据ran_seq_median
path="D:\\min数据\\tcga\\prad_tcga\\prad_tcga\\"
file=paste(path,"data_RNA_Seq_v2_expression_median.RDS",sep="")
data=readRDS(file=file)


