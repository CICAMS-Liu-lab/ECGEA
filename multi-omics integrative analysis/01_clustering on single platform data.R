#####Signal Omic Cluster heatmap
library(ConsensusClusterPlus)
library(dplyr)
library(tidyverse)
library(pheatmap)
df <- read.table("Tumor_Signal_Omic_exp.txt",header = T,row.names = 1)
title = "./mRNA/"

conf <- read_tsv("All_155_key_phenotype_PNI.conf")
conf$ID <- paste("X",conf$ID,sep = "")
df_new = sweep(log2(df),1, apply(log2(df),1,median,na.rm=T))
df_new <- as.matrix(df_new)
results = ConsensusClusterPlus(df_new,maxK=10,reps=1000,pItem=0.80,pFeature=1,
          title="ward.D_ward.D_hc_pearson_sweep_log_new",clusterAlg="hc", innerLinkage = "ward.D",finalLinkage = "ward.D",
          distance="pearson",seed=123.456,writeTable = T,plot = "pdf")
res <- calcICL(results,title="ward.D_ward.D_pam_pearson_sweep_log_consensus_cluster",plot="pdf",writeTable=T)
write.table(res[["itemConsensus"]],file = "Signal_Omic_sweep_log_itemConsensus.txt",col.names = T,row.names = F,sep = "\t",quote = F)

annon <- data.frame(results[[4]]$consensusClass)
colnames(annon) <- "ConsenseCluster"
write.table(annon,file = "Signal_Omic_cluster.txt",col.names = T,row.names = T,sep = "\t",quote = F)
annon$ConsenseCluster <- as.character(annon$ConsenseCluster)
cols = results[[4]]$consensusClass
data=df_new[,order(cols)]
data <- apply(data, 2, function(x) ifelse(x > 3, 3, x))
data <- apply(data, 2, function(x) ifelse(x < -3, -3, x))
bk = unique(c(seq(-3,3, length=100)))
annon$id <- rownames(annon)
annon_col <- annon %>% left_join(conf,by=c("id"="ID")) %>% select(ConsenseCluster,Gender,N_info,TNM,Smoking_history,Drinking_history)
rownames(annon_col) <- annon$id
annon_col$N_info <- as.character(annon_col$N_info)
ann_colors = list(ConsenseCluster = c('1' ="#CC0C00FF", '2'="#5C88DAFF", '3'="#008B8B", '4'="#84BD00FF", '5'="#B0E2FF"),
                  TNM = c("I"="grey70","II"="grey50","III"="grey30","IV"="black"),Gender = c("Female"="#AA6B98","Male"="#496283"),
                  Smoking_history = c("never"= "#C6E3EE", "light"="#87C3D9", "moderate"="#5097B3","heavy"="#496283"),
                  Drinking_history = c("never"= "#C6E3EE", "light"="#87C3D9", "moderate"="#5097B3","heavy"="#496283"),
                  Location = c("Upper"="#8B6914","Middle"="#8B4C39","Lower"="#8B2252"),grade = c("G1"="#7570B3","G2"="#E7298A","G3"="#66A61E"),
                  ICE_PNI = c("0"="grey75","1"="#A020F0","2"="#CDCD00","NA"="grey75"),
                  N_info = c("0"="#F8CCCC","1"="#E88B8A", "2"="#C35454","3"="#80343E"))
pdf("Signal_Omic_cluster.pdf",width=10,height=8)
p1<-pheatmap(data, scale="row",breaks = bk,cluster_rows=T,cluster_cols=F,color = colorRampPalette(c("blue", "white", "red"))(100), fontsize_col = 4,clustering_distance_rows="euclidean",clustering_distance_cols="correlation", 
             clustering_method="ward.D2",show_rownames=F,show_colnames=T,annotation_col = annon_col,annotation_colors = ann_colors)
p1
dev.off()

