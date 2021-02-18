library(scales)
library(RColorBrewer)
library(ggplot2)
library(data.table)

###########
#    Normal	Tumor
#       3.600      79
#       1.800      54
#       3.333      74
#       2.283      62
#       4.533      85
#       2.883      55
#       4.700      88

#data<-read.csv("2D_density_input.xls",header=T,sep="\t")
#data <- fread("2D_density_input.xls", stringsAsFactors = T, encoding = "UTF-8")
options(bitmapType='cairo')
args <- commandArgs(TRUE)
if(length(args)!=1){
  cat("Usage: Rscript 2D_density.r promoter_CG_ratio.xls\n")
  q()
}


fin=args[1]
prefix <- unlist(strsplit(basename(fin), split="\\."))[1]
data<-read.csv(fin,header=T,sep="\t")
#data <- fread(fin, stringsAsFactors = T, encoding = "UTF-8")

###以下适用1M test.xls ###
#keep1 <- rowSums(data>0)>1
#data <- data[keep1,]
#keep2 <- rowSums(data<1)>1
#data <- data[keep2,]
myPalette <- colorRampPalette(rev(brewer.pal(8, "Spectral")))
pdf(file=paste0(prefix,"_dimensional.pdf"),width=8,height=8)
ggplot(data,aes(x=Normal,y=Tumor))+
	stat_density2d(aes(fill=..density..),geom="raster",contour=FALSE)+
	scale_fill_gradientn(colours = myPalette(8)) +
	xlab("CpG methylation (normal)")+
	ylab("CpG methylation (tumor)")+
	theme_bw()+
	theme(legend.title = element_blank(),legend.text = element_text(size=15),
	axis.title.x = element_text(size=22),axis.title.y = element_text(size=20),
	axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
dev.off()
