##########################################################
library(caret)
setwd("/Volumes/LYH/Project/Baidu")
count = read.delim("expression/gene_count_155pair.xls",row.names = 1)
colnames(count) = gsub("X","",colnames(count))
tumour.count = count[,grepl("T",colnames(count))]
colnames(tumour.count) = gsub("^00|^0","",colnames(tumour.count))
tumour.count = tumour.count[,!(colnames(tumour.count)%in%"8T")]
annotation.col = read.delim("expression/used2.annotation.test.xls")
tumour.count = tumour.count[,row.names(annotation.col)]
library(edgeR)
dge <- DGEList(tumour.count,group= factor(annotation.col$coca4))
keep <- rowSums(cpm(dge)>1) >= 20
dge <- dge[keep, , keep.lib.sizes=FALSE]
dge <- calcNormFactors(dge, method='TMM')
design <- model.matrix(~0+factor(annotation.col$coca4))
colnames(design) <- factor(unique(annotation.col$coca4))
dge <- estimateGLMCommonDisp(dge, design)
dge <- estimateGLMTrendedDisp(dge, design)
dge <- estimateGLMTagwiseDisp(dge, design)
cpm <- cpm(dge, normalized.lib.sizes=TRUE)
# test
fit <- glmFit(dge, design)
lrt.4vs3 <- glmLRT(fit, contrast = c(1,-1,0,0))
lrt.4vs2 <- glmLRT(fit, contrast = c(1,0,-1,0))
lrt.4vs1 <- glmLRT(fit, contrast = c(1,0,0,-1))
lrt.3vs2 <- glmLRT(fit, contrast = c(0,1,-1,0))
lrt.3vs1 <- glmLRT(fit, contrast = c(0,1,0,-1))
lrt.2vs1 <- glmLRT(fit, contrast = c(0,0,1,-1))

degs.res.4vs3 <- topTags(lrt.4vs3, n = Inf, adjust.method = 'BH', sort.by = 'PValue')     
degs.res.4vs2 <- topTags(lrt.4vs2, n = Inf, adjust.method = 'BH', sort.by = 'PValue')     
degs.res.4vs1 <- topTags(lrt.4vs1, n = Inf, adjust.method = 'BH', sort.by = 'PValue')     
degs.res.3vs2 <- topTags(lrt.3vs2, n = Inf, adjust.method = 'BH', sort.by = 'PValue')     
degs.res.3vs1 <- topTags(lrt.3vs1, n = Inf, adjust.method = 'BH', sort.by = 'PValue')     
degs.res.2vs1 <- topTags(lrt.2vs1, n = Inf, adjust.method = 'BH', sort.by = 'PValue')     

degs_sig_4vs3 = row.names(degs.res.4vs3)[degs.res.4vs3$table$FDR<=0.05 & abs(degs.res.4vs3$table$logFC)>=1]
degs_sig_4vs2 = row.names(degs.res.4vs2)[degs.res.4vs2$table$FDR<=0.05 & abs(degs.res.4vs2$table$logFC)>=1]
degs_sig_4vs1 = row.names(degs.res.4vs1)[degs.res.4vs1$table$FDR<=0.05 & abs(degs.res.4vs1$table$logFC)>=1]
degs_sig_3vs2 = row.names(degs.res.3vs2)[degs.res.3vs2$table$FDR<=0.05 & abs(degs.res.3vs2$table$logFC)>=1]
degs_sig_3vs1 = row.names(degs.res.3vs1)[degs.res.3vs1$table$FDR<=0.05 & abs(degs.res.3vs1$table$logFC)>=1]
degs_sig_2vs1 = row.names(degs.res.2vs1)[degs.res.2vs1$table$FDR<=0.05 & abs(degs.res.2vs1$table$logFC)>=1]


allgene = c(degs_sig_4vs3,degs_sig_4vs2,degs_sig_4vs1,degs_sig_3vs2,degs_sig_3vs1,degs_sig_2vs1)
allgene = allgene[!duplicated(allgene)]

expr2 = read.delim("expression/Used.baidu.155tumour.TPM.txt",row.names = 1)
mat = expr2[allgene,]
mat = data.frame(t(mat))
row.names(mat) = gsub("X","",row.names(mat))
mat$coca4 = annotation.col[row.names(mat),"coca4"]


up_sig_4vs3 = row.names(degs.res.4vs3)[degs.res.4vs3$table$FDR<=0.05 & degs.res.4vs3$table$logFC>=1]
up_sig_4vs2 = row.names(degs.res.4vs2)[degs.res.4vs2$table$FDR<=0.05 & degs.res.4vs2$table$logFC>=1]
up_sig_4vs1 = row.names(degs.res.4vs1)[degs.res.4vs1$table$FDR<=0.05 & degs.res.4vs1$table$logFC>=1]
up_sig_3vs2 = row.names(degs.res.3vs2)[degs.res.3vs2$table$FDR<=0.05 & degs.res.3vs2$table$logFC>=1]
up_sig_3vs1 = row.names(degs.res.3vs1)[degs.res.3vs1$table$FDR<=0.05 & degs.res.3vs1$table$logFC>=1]
up_sig_2vs1 = row.names(degs.res.2vs1)[degs.res.2vs1$table$FDR<=0.05 & degs.res.2vs1$table$logFC>=1]

down_sig_4vs3 = row.names(degs.res.4vs3)[degs.res.4vs3$table$FDR<=0.05 & degs.res.4vs3$table$logFC<=-1]   
down_sig_4vs2 = row.names(degs.res.4vs2)[degs.res.4vs2$table$FDR<=0.05 & degs.res.4vs2$table$logFC<=-1]   
down_sig_4vs1 = row.names(degs.res.4vs1)[degs.res.4vs1$table$FDR<=0.05 & degs.res.4vs1$table$logFC<=-1]   
down_sig_3vs2 = row.names(degs.res.3vs2)[degs.res.3vs2$table$FDR<=0.05 & degs.res.3vs2$table$logFC<=-1]   
down_sig_3vs1 = row.names(degs.res.3vs1)[degs.res.3vs1$table$FDR<=0.05 & degs.res.3vs1$table$logFC<=-1]   
down_sig_2vs1 = row.names(degs.res.2vs1)[degs.res.2vs1$table$FDR<=0.05 & degs.res.2vs1$table$logFC<=-1]   
##########################################################
# 1 check batch
expr = read.delim("expression/gene_TPM_mat.xls",row.names = 1)
expr   = expr[,grepl("EP",colnames(expr))]
row.names(expr) <- gsub("_.*", "", row.names(expr ))
colnames(expr) = gsub(".*\\.","",colnames(expr))
sampleID=read.delim("expression/immune therapy sample ID-ZYH.csv",sep=",")
expr = expr[,as.character(sampleID[,2][sampleID[,2] %in%   colnames(expr) ] )]
colnames(expr) = sampleID[,1][sampleID[,2] %in%   colnames(expr) ] 
expr2 = read.delim("expression/Used.baidu.155tumour.TPM.txt",row.names = 1)
# expr2 = read.delim("expression/baidu.ESCC.TPM.txt",row.names = 1)
colnames(expr2) = gsub("X","",colnames(expr2))
row.names(expr2 ) <- gsub("_.*", "", row.names(expr2))
mat = merge(expr,expr2,by="row.names")
row.names(mat)=mat[,1]
mat = mat[,-1]
mat = mat[,!(colnames(mat) %in% c("ITPD31","ITPD32","ITPD28"))]
mat <- subset(mat, apply(mat, 1, median) > 0)
batch = c(rep("ITPD",40),c(rep("before",155)))

mat_batch =  limma::removeBatchEffect(log2(mat+1),batch)
mat_batch.pca = prcomp(t(mat_batch), center = TRUE,scale = TRUE)
library(factoextra)
fviz_pca_ind(mat_batch.pca,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = F,     # Avoid text overlapping
             title = "median across samples > 0"
)
annotation.col = read.delim('expression//used2.annotation.test.xls')
drugs.info = read.delim("expression/drugs_info.xls",row.names = 1)
anno = data.frame("drugs"=annotation.col[,c("coca4")],row.names = row.names(annotation.col))
anno = rbind(drugs.info,anno)
sample.order = row.names(anno)[order(anno$drugs)]
anno = anno[order(anno$drugs),]

plot.temp = cor(mat_batch[row.names(mat_batch) %in% gsub("_.*","",allgene),])
plot.temp = plot.temp[sample.order,sample.order]
ComplexHeatmap::Heatmap(plot.temp,show_row_names = F,show_column_names = F,cluster_rows = F,
                        cluster_columns = F,
                        left_annotation = ComplexHeatmap::rowAnnotation(
                          "coca4" = anno))
mat_batch = mat_batch[row.names(mat_batch) %in% gsub("_.*","",allgene),]
coca1 = apply(mat_batch[,colnames(mat_batch) %in% row.names(annotation.col)[annotation.col$coca4 == 1]],1,median)
coca2 = apply(mat_batch[,colnames(mat_batch) %in% row.names(annotation.col)[annotation.col$coca4 == 2]],1,median)
coca3 = apply(mat_batch[,colnames(mat_batch) %in% row.names(annotation.col)[annotation.col$coca4 == 3]],1,median)
coca4 = apply(mat_batch[,colnames(mat_batch) %in% (row.names(annotation.col)[annotation.col$coca4 == 4])],1,median)
pr = apply(mat_batch[,colnames(mat_batch) %in% row.names(drugs.info)[drugs.info$drugs == "PR"]],1,median)
pd = apply(mat_batch[,colnames(mat_batch) %in% row.names(drugs.info)[drugs.info$drugs == "PD"]],1,median)
sd = apply(mat_batch[,colnames(mat_batch) %in% row.names(drugs.info)[drugs.info$drugs == "SD"]],1,median)

cormat = data.frame("coca1"= coca1,"coca2"=coca2,"coca3"=coca3,"coca4"=coca4,"PR"=pr,"PD"=pd,"SD"=sd)
ComplexHeatmap::Heatmap(cor(cormat),#method = "spearman")   ,
                        show_row_names = T,show_column_names = T,cluster_rows = T,
                        col  = circlize::colorRamp2(c(0.935, 0.945,0.95), c("cornflowerblue","white", "red")),
                        cluster_columns = T)
library(flexclust)
nutrient.scale<-cormat
d<-dist(t(nutrient.scale))
fit<-hclust(d,method='complete')

plot(fit,hang=-1,cex=0.8)#hang的作用在于让名字显示得更加规整

clus <-cutree(fit, 4)
d = data.frame(label=names(clus),member=factor(clus))
library(tidyverse)
library(ggtree)
pdf("cluster.coca.drugs.intergration.pdf",h=2,w=4)
ggtree(fit, linetype='dashed', color = "#487AA1") %<+% 
  d + layout_dendrogram() + 
  geom_tiplab(aes(color = member),angle=90,hjust=1,size =6) + 
  theme_dendrogram(plot.margin=margin(6,6,80,6))              
dev.off()


ComplexHeatmap::Heatmap(scale(1-dist(t(cormat))))






