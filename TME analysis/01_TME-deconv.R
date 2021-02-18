# 1.0 ssgsea
######################################################################################
library(GSVA)
library(tidyverse)
signature = "signature/Bindea.cell.2013_Charoentong.Cell_report.2017.MADS.TREG_ShiLM.CCR.2019_fibro_end.rds"
signature <- subset(signature, signature$Ensembl != "NA")
signature_list <- split(signature[, 3], signature[, 1])
mat<- subset(mat, apply(mat, 1, median) > 0)
scrore <- gsva(as.matrix(mat),signature_list, method = "ssgsea",kcdf = "Gaussian", tau = 0.25, ssgsea.norm = TRUE)
write.table(scrore, "ssgsea/all_ssgsea_gsva_NES.txt", sep = "\t")
######################################################################################
# 2 cluster
######################################################################################
# 2.0 load package and function#####
library(ComplexHeatmap)
library(ConsensusClusterPlus)
scale_rows = function(x){
  m = apply(x, 1, mean, na.rm = T)
  s = apply(x, 1, sd, na.rm = T)
  return((x - m) / s)
}
scale_mat = function(mat, scale){
  if(!(scale %in% c("none", "row", "column"))){
    stop("scale argument shoud take values: 'none', 'row' or 'column'")
  }
  mat = switch(scale, none = mat, row = scale_rows(mat), column = t(scale_rows(t(mat))))
  return(mat)
}

######################
ssgsea = read.delim("ssgsea/all_ssgsea_gsva_NES.txt")
colnames(ssgsea) = gsub("X","",colnames(ssgsea))
TILs.ssgsea = ssgsea
sHc <- hclust(ddist <- dist(t(TILs.ssgsea)), method ="ward.D2")
ans <- ConsensusClusterPlus(ddist, maxK = 10, pItem =0.9,
                            reps = 1000, title = "mc_consensus_k7_1000", clusterAlg = "hc",
                            innerLinkage = "ward.D2", finalLinkage = "ward.D2",
                            plot = "pdf", writeTable = TRUE)

res = ans
annotation.col=data.frame("TME2.subtype"=sapply(res[[2]]$consensusClass, function(x) paste("Cluster",x,sep = "")),
                          "TME3.subtype"=sapply(res[[3]]$consensusClass, function(x) paste("Cluster",x,sep = "")),
                          "TME4.subtype"=sapply(res[[4]]$consensusClass, function(x) paste("Cluster",x,sep = "")))
lymp.ratio = read.delim("annotation/lymphocyte_info.csv",sep = ',',row.names = 1)
row.names(lymp.ratio)=gsub("^0","",gsub("^00","",row.names(lymp.ratio)))
annotation.col =merge(annotation.col,lymp.ratio,by="row.names",all.x = T)
row.names(annotation.col) = annotation.col$Row.names
annotation.col = annotation.col[,-1]
drugs = read.delim("annotation/drugs_info.xls",row.names = 1)
annotation.col =merge(annotation.col,drugs,by="row.names",all.x = T)
row.names(annotation.col) = annotation.col$Row.names
annotation.col = annotation.col[,-1]
purity = read.delim("annotation/absolute-estimate-PAMES-purity.xls",row.names = 1)
row.names(purity) = gsub("^0","",gsub("^00","",gsub("X","",row.names(purity))))
annotation.col =merge(annotation.col,purity,by="row.names",all.x = T)
row.names(annotation.col) = annotation.col$Row.names
annotation.col = annotation.col[,-1]
annnnotation.clinical =  read.delim("annotation/final.clinical.csv",sep = ",",row.names = 1)
row.names(annnnotation.clinical) = gsub("^0","",gsub("^00","",gsub("X","",row.names(annnnotation.clinical))))
row.names(annnnotation.clinical) =paste(row.names(annnnotation.clinical),"T",sep = "")
colnames(annnnotation.clinical)
annotation.col =merge(annotation.col,annnnotation.clinical,by="row.names",all.x = T)
row.names(annotation.col) = annotation.col$Row.names
annotation.col = annotation.col[,-1]
used.anno.col = read.delim("annotation/used.annotation.xls")
used.anno.col = used.anno.col[,c(2,4:11,20, 23:29,31:32,34:38)]
annotation.col =merge(annotation.col,used.anno.col,by="row.names",all.x = T)
row.names(annotation.col) = annotation.col$Row.names
annotation.col = annotation.col[,-1]

TILs.ssgsea[TILs.ssgsea < 0] <- 0.0000000001
annotation.col$shannon.index = -sapply(1:ncol(TILs.ssgsea),function(y){
  sum(apply(as.matrix(TILs.ssgsea[,y]), 2, 
            function(x) {(x/sum(TILs.ssgsea[,y]))*log(x/sum(TILs.ssgsea[,y]),2)}))
}
)

annotation.col$TME3.subtype.used=sapply(1:nrow(annotation.col),function(x) {
  if(annotation.col$TME4.subtype[x]=="Cluster3"){
    "Cluster1"
  }else if(annotation.col$TME4.subtype[x]=="Cluster1"){
    "Cluster2"
  }else if(annotation.col$TME4.subtype[x]=="Cluster2"){
    "Cluster3"
  }else {
    "Cluster3"
  }
})
row.names(ssgsea) = c("aDC","B cells","CD8 T cells","Cytotoxic cells","DC","Endothelial",
                      "Eosinophils","Fibroblasts", "iDC","Macrophages","Mast cells", "MDSC",
                      "Neutrophils","NK CD56bright cells","NK CD56dim cells","NK cells","pDC",
                      "T helper cells ","naive CD4 T cells","CD4 Tcm cells","CD4 Tem cells",
                      "TFH","Tgd","Th1","Th17", "Th2","Tregs" )    

htmat= scale_mat(ssgsea,scale="row")
annotation.col = annotation.col[order(annotation.col$TME2.subtype),] 
annotation.row = read.delim("annotation/cellType.txt",row.names = 1,sep = ",",header=F)

col_ha = HeatmapAnnotation(df=annotation.col[,c(8,17,13,14,29,30,31,32,37:39,40,26,3,2,1,49)],
                           col=list(
                             Gender = c(Female = "skyblue","Male"="seashell3","NA"="gray"),
                             # Drinking.history = c("heavy"="black","light"="gray","never"="black"),
                             TME2.subtype=c("Cluster1"="#D65F4C","Cluster2"="#689e46"),
                             TME3.subtype=c("Cluster1"="red","Cluster2"="blue","Cluster3"="green"),
                             TME4.subtype=c("Cluster1"="red","Cluster2"="blue","Cluster3"="green","Cluster4"="lightgreen"),
                             TME3.subtype.used =c("Cluster1"="#D65F4C","Cluster2"="slateblue","Cluster3"="#689e46"),
                             "lymphocyte_ratio" =circlize::colorRamp2(c(-2, 0, 2), c("#5b5985", "white", "#0e0e0e")),
                             coca4=c('1' ="#CC0C00FF", '2'="#5C88DAFF", '3'="#008B8B", '4'="#84BD00FF"),
                             drugs=c("PR"="#00AFBB","SD"= "#E7B800", "PD"="#FC4E07","missing"= "lightgray","Without drugs"="lightgray"),
                             "os.state"=c("0" = "pink",'1'="#98A0A7","loss"="lightgray"),
                             "psf.state"=c("0" = "pink",'1'="#98A0A7","loss"="lightgray"),
                             "cytotoxic.score" = circlize::colorRamp2(c(-2, 0, 2), c("#5b5985", "white", "#2079c6")),
                             "exhaustion.score" = circlize::colorRamp2(c(-2, 0, 2), c("#5b5985", "white", "#2079c6")),
                             PD1=circlize::colorRamp2(c(-2, 0, 2), c("#598ac6", "white", "#ffbf53")),
                             PDL1=circlize::colorRamp2(c(-2, 0, 2), c("#598ac6", "white", "#ffbf53")),
                             "MSI_score"=circlize::colorRamp2(c(-2, 0, 2), c("#598ac6", "white", "#ffbf53")),
                             TMB.score=circlize::colorRamp2(c(-2, 0, 2), c("#598ac6", "white", "#ffbf53")),
                             APOBEC.score = circlize::colorRamp2(c(-2, 0, 2), c("#598ac6", "white", "#ffbf53"))
                             # "shannon.index"=circlize::colorRamp2(c(-2, 0, 2), c("#598ac6", "white", "#ffbf53"))
                             
                           ),
                           na_col = "lightgray",annotation_name_gp = gpar(fontsize = 6),
                           simple_anno_size = unit(3, "mm"), 
                           annotation_legend_param = list(title_gp = gpar(fontsize = 6, fontface = "bold"),
                                                          labels_gp = gpar(fontsize = 6)
                                                          ,direction = "horizontal")
)

pdf("fig/TME.Cluster-bindea-3.pdf",width=8,height=6)
Heatmap(htmat[,row.names(annotation.col)[order(annotation.col$TME2.subtype)]],top_annotation = col_ha,
        row_split = annotation.row[row.names(htmat),1],
        cluster_rows = T,cluster_columns = FALSE,height = unit(6, "cm"),
        left_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = 2:3),
                                                         labels = c("stromal","adaptive", "innate"), 
                                                         labels_gp = gpar(col = "white", fontsize = 10))),
        row_names_gp = gpar(fontsize = 6),show_column_names = F,show_parent_dend_line = FALSE,row_title = NULL,
        col = circlize::colorRamp2(c(-2, 0, 2), c("navyblue", "white", "firebrick")),name = "NES",
        heatmap_legend_param = list(title_gp = gpar(fontsize = 6, fontface = "bold"),labels_gp = gpar(fontsize = 6),direction = "horizontal")
        
)

dev.off()
write.table(annotation.col,"annotation/used2.annotation.xls",sep = "\t")
######################################################################################

#######################预后分析#################
library(survival)   
library(survminer)
library(stringr)
library(dplyr)
annotation.col =read.delim("annotation/used2.annotation.xls",stringsAsFactors =F)
annotation.col.surv = annotation.col[!grepl("^ITPD",row.names(annotation.col)),]
annotation.col.surv = annotation.col.surv[annotation.col.surv$os.state!="loss",]
annotation.col.surv = annotation.col.surv[annotation.col.surv$TME3.subtype.used!='Cluster2',]
annotation.col.surv$s <- grepl("1", annotation.col.surv$os.state, ignore.case = TRUE)
survival.data <- annotation.col.surv[, c("os.time", "s", "TME3.subtype.used")]
f.m <- formula(Surv(as.numeric(os.time), event = survival.data$s) ~
                 survival.data$TME3.subtype.used)
fit <- do.call(survfit, list(formula = f.m, data = survival.data))
label.add.n <- function(x) {
  na.idx <- is.na(survival.data [, "os.time"])
  negative.idx <- survival.data [, "os.time"] < 0
  idx <- !(na.idx | negative.idx)
  return(paste0(x, " (n = ", sum(survival.data [idx, "TME3.subtype.used"] == x),
                ")"))
}

d <- survminer::surv_summary(fit, data = survival.data)
order <- unname(sapply(levels(d$strata), function(x) unlist(str_split(x,
                                                                      "="))[2]))
labels <- sapply(order, label.add.n)
pdf("fig/prognosis.os.bindea.2.pdf",w= 6,h=6)
ggsurvplot(fit, risk.table = TRUE, pval = TRUE,legend.labs = labels, 
           conf.int = F, #xlim = xlim, main = main, xlab = xlab,
           #legend.title = legend, 
           palette = c("#D65F4C","#689e46"),
)

dev.off()
#################################

mat = read.delim("expression/Used.baidu.155tumour.TPM.txt",row.names = 1)
colnames(mat) = gsub("X","",colnames(mat))
row.names(mat) <- gsub("_.*", "", row.names(mat))

exhaustion = c("CTLA4", "HAVCR2", "LAG3", "PDCD1","TIGIT")
exhaustion = c("ENSG00000163599","ENSG00000135077","ENSG00000089692","ENSG00000188389","ENSG00000181847")
cytotoxic = c("CST7", "GZMA", "GZMB", "IFNG", "NKG7" , "PRF1","CD8A","CD8B","GNLY","GZMH")
cytotoxic = c("ENSG00000077984","ENSG00000145649","ENSG00000100453","ENSG00000111537",
              "ENSG00000105374","ENSG00000180644","ENSG00000153563","ENSG00000172116",
              "ENSG00000115523","ENSG00000100450")
cytotoxic.score = apply(mat[cytotoxic,],2,mean)
exhaustion.score = apply(mat[exhaustion,],2,mean)
annotation.col.surv$cytotoxic.score2 = cytotoxic.score[row.names(annotation.col.surv)]
annotation.col.surv$exhaustion.score2 = exhaustion.score[row.names(annotation.col.surv)]
annotation.col.surv$cytotoxic.score = scale(annotation.col.surv$cytotoxic.score2)
annotation.col.surv$exhaustion.score = scale(annotation.col.surv$exhaustion.score2)

p=ggboxplot(reshape2::melt(annotation.col.surv[,c("TME3.subtype.used","cytotoxic.score","exhaustion.score")]),
            x = "variable", y = "value", color = "TME3.subtype.used",
            palette = c("Cluster1"="#D65F4C","Cluster2"="slateblue","Cluster3"="#689e46"),
            add = "jitter")
my_comparisons <- list( c("Cluster1","Cluster2"),c("Cluster2","Cluster3"), c("Cluster1","Cluster3"))
p + stat_compare_means(comparisons = my_comparisons,method = "wilcox.test")
ggsave("fig/barplot_cytotoxic_exhaustion.2.score.tiff",width = 6,height = 3)
p=ggboxplot(annotation.col.surv,
            x = "TME3.subtype.used", y = "cytotoxic.score2", color = "TME3.subtype.used",
            palette = c("Cluster1"="#D65F4C","Cluster2"="slateblue","Cluster3"="#689e46"),
            add = "jitter")
my_comparisons <- list( c("Cluster1","Cluster2"),c("Cluster2","Cluster3"), c("Cluster1","Cluster3"))
p + stat_compare_means(comparisons = my_comparisons,method = "wilcox.test")
ggsave("fig/barplot_bindea_cytotoxic.2.score.tiff",width = 3,height = 4)

p=ggboxplot(annotation.col.surv,
            x = "TME3.subtype.used", y = "exhaustion.score", color = "TME3.subtype.used",
            palette = c("Cluster1"="#D65F4C","Cluster2"="slateblue","Cluster3"="#689e46"),
            add = "jitter")
my_comparisons <- list( c("Cluster1","Cluster2"),c("Cluster2","Cluster3"), c("Cluster1","Cluster3"))
p + stat_compare_means(comparisons = my_comparisons,method = "wilcox.test")
ggsave("fig/barplot_bindea_exhaustion.2.score.tiff",width = 3,height = 4)
