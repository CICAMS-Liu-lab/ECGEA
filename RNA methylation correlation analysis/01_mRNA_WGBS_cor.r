library(ggplot2)
library(RColorBrewer)
ComputeCorrForMiM<-function(MM_matrix,method,OutPutFileName){
	CCbox<-vector()
	Pbox<-vector()
	N1<-dim(MM_matrix)[1]
	N2<-dim(MM_matrix)[2] #col num
	N_sample=(N2-2)/2
	WGBS_col<-3:(2+N_sample)  #col of WGBS
	mRNA_col<-(3+N_sample):N2  #col of mRNA
	for(i in 1:N1){
		vector_WGBS<-as.numeric(MM_matrix[i,WGBS_col])
		vector_mRNA<-as.numeric(MM_matrix[i,mRNA_col])
		name=MM_matrix[i,1]
		Symbol=MM_matrix[i,2]
		pdf(file = paste("plot_",name,".pdf", sep = ""),width=10,height=10)
		if(method=="pearson"){
			if(length(vector_WGBS)<3){
				print("No enough values for caculating the pvalues,please input at least 3 pairs of values!")
				break
			}
			else{
				CorrCoeffcient<-cor.test(vector_WGBS/100,vector_mRNA,method = c("pearson"))
				plot(vector_WGBS/100,vector_mRNA,pch=16,col=brewer.pal(308,"Paired")[4:5])
				dev.off()
			}
		}
		if(method=="spearman"){
			CorrCoeffcient<-cor.test(vector_WGBS/100,vector_mRNA,method = c("spearman"))
			par(mar=c(1,1,3,1),pin=c(8,8))
			plot(vector_WGBS/100,vector_mRNA,pch=16,main=Symbol,cex.main=2.5,
			     font.main=1,xlab="Methylation level",ylab="Expression level",
			     cex.lab = 2,cex.axis=1.5,col=brewer.pal(308,"Paired")[4:5])
			legend("topright",c("Normal","Tumor"),pch=16,col=brewer.pal(12,"Paired")[4:5],cex=2)
			dev.off()
		}
		if(method=="kendall"){
			CorrCoeffcient<-cor.test(vector_WGBS/100,vector_mRNA,method = c("kendall"))
			plot(vector_WGBS/100,vector_mRNA,pch=16,col=brewer.pal(308,"Paired")[4:5])
			
			dev.off()
		}
		CCbox[i]<-CorrCoeffcient$estimate
		Pbox[i]<-CorrCoeffcient$p.value
		
	}
	FDRbox<-p.adjust(as.numeric(Pbox),"BH") ## "BH":Benjamini & Hochberg (1995) ("BH" or its alias "fdr")
	if (N_sample==2){
		ResultBox<-cbind(MM_matrix,CCbox)
		colnames(ResultBox)<-c(names(MM_matrix),"CorrC")
		ResultBox<-ResultBox[order(ResultBox[,1],decreasing=FALSE),]
		write.table(ResultBox,OutPutFileName,sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
		return(ResultBox)
	}else{
		ResultBox<-cbind(MM_matrix,CCbox,Pbox,FDRbox)
		colnames(ResultBox)<-c(names(MM_matrix),"CorrC","P_value","FDR")
		ResultBox<-ResultBox[order(ResultBox[["FDR"]],decreasing=FALSE),]
		write.table(ResultBox,OutPutFileName,sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
		return(ResultBox) 
	}
}

args<-commandArgs(trailingOnly=TRUE)
if(length(args) < 4){
	stop("
	Usage:
	Rscript cor.r --args <dif_WGBS-dif_mRNA.exp.xls> <method> <out_prefix>
	####################################
		 -  spearman
		|
	method - -  kendall
		|		
		 -  pearson
	
	Notice: when sample_num less than 3, \"pearson\" can't be used !!!
	#### dif_WGBS-dif_mRNA.exp.xls ####
	WGBS		Targetgene	Sample1_m	Sample2_m		... 	Sample1_g		Sample1_g 	...
	miR-1		gene1		19.45		4.34			...	3.538			22.259		...
	miR-2		gene2		350.61		68.33			...	0.521			4.799		...
	####################################
	")
}
MM_matrix<-unique(read.delim(args[2],header=TRUE,sep="\t"))
OutPutFileName<-paste(args[4],".xls",sep="")
Result_CC<-ComputeCorrForMiM(MM_matrix,method=args[3],OutPutFileName)
N_sample=(dim(MM_matrix)[2]-2)/2
if (N_sample==2){
	significant_result<-subset(Result_CC,as.numeric(as.character(CorrC))==-1)
	Sig_OutPutFileName<-paste(args[4],".sig.xls",sep="") ## select -1
	write.table(significant_result,Sig_OutPutFileName,sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
}else{
	significant_result<-Result_CC[which(Result_CC[["FDR"]]<=0.05),]
	Sig_OutPutFileName<-paste(args[4],".sig.xls",sep="") ## select sig
	write.table(significant_result,Sig_OutPutFileName,sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
}
