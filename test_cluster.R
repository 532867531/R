library(factoextra)
library(cluster)
# ����׼��
# ʹ�����õ�R���ݼ�USArrests
data("USArrests")
# remove any missing value (i.e, NA values for not available)
USArrests = na.omit(USArrests) #view the first 6 rows of the data
head(USArrests, n=6) 
# ��ʾ��������ʾ������
# �ھ���֮ǰ���ǿ����Ƚ���һЩ��Ҫ�����ݼ�鼴����������ͳ�ƣ���ƽ��ֵ����׼���
desc_stats = data.frame( Min=apply(USArrests, 2, min),#minimum
                         Med=apply(USArrests, 2, median),#median
                         Mean=apply(USArrests, 2, mean),#mean
                         SD=apply(USArrests, 2, sd),#Standard deviation
                         Max=apply(USArrests, 2, max)#maximum
)
desc_stats = round(desc_stats, 1)#����С�����һλhead(desc_stats)
desc_stats
# �����кܴ�ķ����ֵʱ����б�׼��
df = scale(USArrests)




# ���ݼ�Ⱥ��������ʹ��get_clust_tendency()����Hopkinsͳ����
res = get_clust_tendency(df, 40, graph = TRUE)
res$hopkins_stat
win.graph();res$plot
##����k��ֵ������Ҫָ��Ҫ���ɵľ���������������ǽ�ʹ�ú���clusGap()���������ڹ������ž�����������fviz_gap_stat()���ڿ��ӻ���
set.seed(123)
## Compute the gap statistic
gap_stat = clusGap(df, FUN = kmeans, nstart = 25, K.max = 10, B = 500)
# Plot the result
win.graph();fviz_gap_stat(gap_stat)
##ͼ����ʾ���Ϊ�۳����ࣨk=4��kmeans���о���  kmeans��������о��࣬ѡ��25�������
km.res = kmeans(df, 4, nstart = 25)
# Visualize clusters using factoextra
win.graph();fviz_cluster(km.res, USArrests)
##��ȡ��������ͼ
sil = silhouette(km.res$cluster, dist(df))
rownames(sil) = rownames(USArrests)
head(sil[, 1:3])
##�ĸ�cluster�Ļ�����Ϣ
# Visualize
fviz_silhouette(sil)
##ͼƬ�ߴ��900 dpi���ʺ�΢���ֻ����Ķ�  ͼ�п��Կ����и�ֵ������ͨ������silhouette()ȷ�����ĸ��۲�ֵ
neg_sil_index = which(sil[, "sil_width"] < 0)
sil[neg_sil_index, , drop = FALSE]
##��ʾΪ���Ĺ۲�ֵ


#eclust():��ǿ�ľ������
#�����������������ȣ�eclust()�������ŵ㣺 ���˾�������Ĺ������̣��������ڼ����ξ���ͷ������࣬eclust()�Զ�������Ѿ�������� �Զ��ṩSilhouette plot�����Խ��ggplot2����������ͼ�Σ�ʹ��eclust()��K��ֵ����
# Compute k-means
res.km = eclust(df, "kmeans")
# Gap statistic plot
fviz_gap_stat(res.km$gap_stat)
##ʹ��eclust()�Ĳ�ξ���
# Enhanced hierarchical clustering
res.hc = eclust(df, "hclust") # compute hclust
fviz_dend(res.hc, rect = TRUE) # dendrogam
#�����R��������Silhouette plot�ͷֲ����ɢ��ͼ��
fviz_silhouette(res.hc) # silhouette plot
fviz_cluster(res.hc) # scatter plot


############################################
############################################
############################################
############################################
############################################
##################���濪ʼ���������Լ�������###############
############################################
############################################
############################################
##��ȡ����ran_seq_median
path="D:\\min����\\tcga\\prad_tcga\\prad_tcga\\"
file=paste(path,"data_RNA_Seq_v2_expression_median.RDS",sep="")
data=readRDS(file=file)

