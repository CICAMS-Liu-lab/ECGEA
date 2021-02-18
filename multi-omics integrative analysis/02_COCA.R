library(ConsensusClusterPlus)
library(dplyr)
library(tidyverse)
library(pheatmap)
library(philentropy)
library(coca)
library(ggpubr)
library(RColorBrewer)
df <- read.table("Coca_155_mRNA_promoter_miRNA_cnv_subtype_T.txt",header = T,row.names = 1)
title = "./COCA/Consensus/"
conf <- read_tsv("All_155_key_phenotype_PNI.conf")
conf$ID <- paste("X",conf$ID,sep = "")
df_value <- read_tsv("ESCC_155_Apopec_TMB_PD1_PDL1.txt")
df_value$Sample_id <- paste("X",df_value$Sample_id,sep = "")
conf <- conf %>% left_join(df_value,by=c("ID"="Sample_id"))
conf$ICE_PNI <- gsub("2","1",conf$ICE_PNI)
df_new <- as.matrix(cor(df),method = "pearson")
results = ConsensusClusterPlus(df_new,maxK=10,reps=1000,pItem=0.85,pFeature=1,
          title="ward.D_ward.D_pam_cor_pearson",clusterAlg="pam", innerLinkage = "ward.D",finalLinkage = "ward.D",
          distance="pearson",seed=123.456,writeTable = T,plot = "pdf")

annon <- data.frame(results[[4]]$consensusClass)
colnames(annon) <- "ConsenseCluster"
write.table(annon,file = "Coca_155_pam_cor_pearson_4_re.txt",col.names = T,row.names = T,sep = "\t",quote = F)
annon$ConsenseCluster <- as.character(annon$ConsenseCluster)
cols = results[[4]]$consensusClass
data=df[,order(cols)]
annon$id <- rownames(annon)
annon_col <- annon %>% left_join(conf,by=c("id"="ID")) %>% select(ConsenseCluster,N_info,TNM)
rownames(annon_col) <- annon$id
annon_col$N_info <- as.character(annon_col$N_info)

ann_colors = list(Con_pam_dist_5 = c('1' ="#CC0C00FF", '2'="#5C88DAFF", '3'="#008B8B", '4'="#84BD00FF", '5'="#B0E2FF"),
                ConsenseCluster = c('1' ="#CC0C00FF", '2'="#5C88DAFF", '3'="#008B8B", '4'="#84BD00FF", '5'="#B0E2FF"),
                  TNM = c("I"="grey70","II"="grey50","III"="grey30","IV"="black"),Gender = c("Female"="#008B8B","Male"="#8B1C62"),
                  Smoking_history = c("never"= "#FFAEB9", "light"="#EEA2AD", "moderate"="#CD8C95","heavy"="#8B5F65"),
                  Drinking_history = c("never"= "#FFAEB9", "light"="#EEA2AD", "moderate"="#CD8C95","heavy"="#8B5F65"),
                  Location = c("Upper"="#8B6914","Middle"="#8B4C39","Lower"="#8B2252"),grade = c("G1"="#C6E3EE","G2"="#87C3D9", "G3"="#5097B3"),
                  N_info = c("0"="#F8CCCC","1"="#E88B8A", "2"="#C35454","3"="#80343E"))

write.table(data,file = "ESCC_COCA_matrix.xls",col.names = T,row.names = T,sep = "\t",quote = F)

data_new <- read.table("ESCC_COCA_matrix.xls",header = T,row.names = 1)

#####COCA#####
datasetNames <- rownames(df)
M = length(table(datasetNames))
sumK = dim(data)[1]
Indicator <- as.integer(c(rep(1, 4), rep(2, 4), rep(3, 4), rep(4,4)))
for(i in 1:sumK){
  data[i,] <- data[i,]*Indicator[i]
}
mycols <- c("#EDEDED","#EE0000","#0000FF","#080808","#00CD00")
mycols1 <- c("#EDEDED",(RColorBrewer::brewer.pal(n = 4, name = "Set3")))
mycols_set <- c("#EDEDED","#FB8072","#80B1D3","#FDB462","#BEBADA")
mycols_set1 <- c("#EDEDED","#D75B5A","#60B9B5","#4A6486","#AA6B98")
#data_test <- data[c("CNV3","Promoter3","miRNA3","mRNA1","CNV1","Promoter2","miRNA4","miRNA1","CNV2","Promoter4","mRNA2","miRNA2","mRNA4","Promoter1","CNV4","mRNA3"),]
#data_test <- data[c("Promoter1","miRNA3","mRNA1","CNV1","Promoter2","miRNA4","miRNA1","CNV2","CNV3","Promoter3","Promoter4","mRNA2","miRNA2","mRNA4","CNV4","mRNA3"),]
data_test <- data[c("miRNA3","Promoter3","CNV3","Promoter2","mRNA1","miRNA4","CNV1","CNV2","miRNA1","mRNA3","Promoter4","mRNA2","miRNA2","Promoter1","mRNA4","CNV4"),]
pdf("COCA_Omic_cluster.pdf",width=10,height=8)
p1 <- pheatmap(data_test, scale="none",color = mycols_set1, cluster_rows=F,cluster_cols=F, fontsize_col = 6,clustering_distance_rows="binary",
             clustering_method="complete",gaps_row = c(3, 7,10,14),gaps_col = c(39, 77,107),show_rownames=T,show_colnames=F,annotation_col = annon_col,annotation_colors = ann_colors, drop_levels = FALSE, na_col = "seashell2")
dev.off()


