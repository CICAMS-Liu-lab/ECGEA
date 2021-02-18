args=commandArgs(TRUE)
if(length(args)!=2){
        cat("Useage: Rscript getMiRNACount.R work_evn sample_name")
        q()
}

cwd=args[1]
setwd(cwd)
sample=args[2]
file=dir(pattern="miRNAs_expressed_all_samples")
data=read.csv(file[1],sep="\t",header=T,comment.char = "",check.names = F)
colnames(data)=c("miRNA","read_count","precursor","total","seq","seq(norm)")
data$miRNA=as.factor(data$miRNA)
expression=as.data.frame(floor(as.matrix(tapply(data$read_count,data$miRNA,mean))))
colnames(expression)<-args[2]
write.table(expression,file=paste(args[2],"mean.count.txt",sep="."),sep="\t",quote=F,col.names=NA)

