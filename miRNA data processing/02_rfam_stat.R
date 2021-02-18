args<-commandArgs(TRUE)
if(length(args)!=1){
	cat("Usage: Rscript rfam_stat.R rfam.stat\n")
	q()
}

ab<-read.table(args[1],sep="\t")
sample<-unlist(strsplit(args[1],'.',fixed = TRUE))[1]

absum<-tapply(ab$V1,ab$V2,sum)
data<-as.data.frame(absum[order(absum,decreasing=T)])
colnames(data)<-"Abundance"
data$pct<-data$Abundance/sum(data$Abundance)
write.table(data,paste(sample,".proportion.stat.txt",sep=""),col.names=NA,sep="\t",quote=F)

main<-data[data$pct>=0.05,]
other<-data[data$pct<0.05,]
ot<-data.frame(row.names = "others",Abundance=sum(other$Abundance),pct=sum(other$pct))
draw<-rbind(main,ot)
draw$plab<-paste(format(round(100*draw$pct, 2),nsmall=2), "%")
labels<-rownames(draw)

png(paste(sample,".proportion.png",sep=""),width=6.25,height=6.25,units="in",res=1200)

pie(draw$pct,main=paste(basename(sample),".proportion",sep=""),label=draw$plab,col=rainbow(7))
legend("topright",labels,bty="n",fill=rainbow(7),cex = 0.8)
dev.off()

