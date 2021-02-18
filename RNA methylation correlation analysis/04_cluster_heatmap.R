options(bitmapType='cairo')
args <- commandArgs(TRUE)
if(length(args)!=2){
  cat("Usage: Rscript cluster_heatmap.R gene_CPM_TMM.txt condition.conf\n")
  q()
}


fin=args[1]
condition=args[2]
prefix <- unlist(strsplit(basename(fin), split="\\."))[1]
data=read.table(fin,header=T,row.names=1,check.names=F)
con=read.table(condition,header=T,sep="\t",row.names=1)
if(ncol(data)<=1){
    print("No heatmap plot for less than 3 samples")
    quit(status=0)
}

Tumor=grepl("T",names(data))
data_new=data[,Tumor]
#data=data_new[,!names(data_new)=="218T"]
data=data_new

mostVar <- function(data, n, i_want_most_var = TRUE) {
  data.var <- apply(data, 1, stats::var)
  data[order(data.var, decreasing = i_want_most_var)[1:n],]
}

#keep <- rowSums(data>0)>10
#data <- data[keep,]
dim(data)
#data <- mostVar(data,n=1000)

library("pheatmap")

con$Coca[which(con$Con_Cor_4=="1")] <- "COCA1"
con$Coca[which(con$Con_Cor_4=="2")] <- "COCA2"
con$Coca[which(con$Con_Cor_4=="3")] <- "COCA3"
con$Coca[which(con$Con_Cor_4=="4")] <- "COCA4"
con$CIMP <- gsub("CIMP_1$","CIMP+",con$CIMP)
con$CIMP <- gsub("CIMP_2$","CIMP-",con$CIMP)

annotation_col =con[,c(21,17,5,9:10,3:4)]

ann_colors = list(
	GLM6 =c("CIMP+" ="black", "CIMP-"="white"),
        CIMP_biomarker=c("+" ="black", "-"="white"),
	TNM = c("I"="grey75","II"="grey55","III"="grey35","IV"="grey15"),
	Gender = c("Female"="#AA6B98","Male"="#496283"),
	Smoking_history = c("NA"="#B5B1B1", "never"= "#C6E3EE", "light"="#87C3D9", "moderate"="#5097B3","heavy"="#496283"),
	Drinking_history = c("NA"="#B5B1B1", "never"= "#C6E3EE", "light"="#87C3D9", "moderate"="#5097B3","heavy"="#496283"),
	grade = c("G1"="#C6E3EE","G2"="#87C3D9", "G3"="#5097B3"),
        N_info = c("0"="#F8CCCC","1"="#E88B8A", "2"="#C35454","3"="#80343E")

)
pdf(file=paste0(prefix,"_heatmap.pdf"),width=15,height=9, onefile=FALSE)
res=pheatmap(data,cluster_rows=T,cluster_cols=T,
         fontsize_col = 5,fontsize_row = 10,
	 clustering_distance_cols="euclidean",clustering_method="ward.D2",
	 legend=T,annotation_legend = T,show_rownames=1,show_colnames=T,cutree_cols=2,
	 annotation_col = annotation_col,annotation_colors = ann_colors)
order_row = res$tree_row$order
order_col = res$tree_col$order
datat = data.frame(data[order_row,order_col])
datat = data.frame(rownames(datat),datat,check.names =F)
colnames(datat)[1] = "id"
write.table(datat,file=paste0(prefix,"_reorder.txt"),row.names=FALSE,quote = FALSE,sep='\t')
dev.off()
