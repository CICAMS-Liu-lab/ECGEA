library(ggplot2)
library(ggsci)
library(ggpubr)

data<-read.csv("DMRs_pie_cgi_input.xls",header = T,sep="\t")
hyper=data[which(data$DMR =='recurrent hyperDMRs'|data$DMR == 'hyper DMRs'|data$DMR =='Recurrent_hyperDMRs'),]
#hyper=data[which(data$DMR =='hyper_DMRs'|data$DMR =='recurrent_hyperDMRs'),]
hyper$per=paste0(round(hyper$Ratio*100,1),"%")

hypo=data[which(data$DMR =='hypo DMRs'|data$DMR =='Recurrent hypoDMRs'),]
hypo$per=paste0(round(hypo$Ratio*100,1),"%")

#hyper=hyper[order(hyper$DMR,hyper$Feature,decreasing=T),]
#order <- sort(hyper$Feature,index.return=TRUE,decreasing = TRUE)
#hyper$Feature <- factor(hyper$Feature, levels = hyper$Feature[order$ix])

#pdf(file="hyper_DMRs_pie_CGI.pdf",width=15,height=8, onefile=FALSE)
#ggplot(hyper,aes(x=DMR,y=Ratio,fill=factor(Feature)))+
p1=ggplot(hyper,aes(x=DMR,y=Ratio,fill=factor(Feature,levels=c("Promoter","Intergenic","Genebody","Enhancer","CGI","Promoter_CGI","Intergenic_CGI","Genebody_CGI","Enhancer_CGI"))))+
	geom_bar(stat = "identity",position='stack')+
	#scale_fill_npg()+
	scale_fill_manual(values=c("#CD3333","#6CA6CD","#388E8E","#00008B","#CD8500","#CD5C5C","#CAE1FF","#66CDAA","#7171C6"))+
	geom_text(aes(label = per), colour ="white",size=5,position = position_stack(vjust = 0.5))+
	xlab("")+
	ylab("Frequency")+
        theme_bw()+
        theme(legend.title=element_blank(),legend.text = element_text(size = 20),
        axis.text.x = element_text(size= 25),axis.text.y = element_text(size= 20),axis.title.y = element_text(size= 25))+
#	+guides(fill=FALSE)+
	theme(legend.title=element_blank(),legend.text = element_text(size = 15),
	axis.text.x = element_text(size= 15))+
	theme_classic()+
	coord_polar(theta = "y")+
	theme(legend.title=element_blank()) 
#dev.off()

p2=ggplot(hypo,aes(x=DMR,y=Ratio,fill=factor(Feature,levels=c("Promoter","Intergenic","Genebody","Enhancer","CGI","Promoter_CGI","Intergenic_CGI","Genebody_CGI","Enhancer_CGI"))))+
        geom_bar(stat = "identity",position='stack')+
	#scale_fill_npg()+
	scale_fill_manual(values=c("#CD3333","#6CA6CD","#388E8E","#00008B","#CD8500","#CD5C5C","#CAE1FF","#66CDAA","#7171C6"))+
	geom_text(aes(label = per),colour ="white",size=5,position = position_stack(vjust = 0.5))+
        xlab("")+
        ylab("Frequency")+
	#labs(y="")+
	theme_bw()+
        theme(legend.title=element_blank(),legend.text = element_text(size = 20),
        axis.text.x = element_text(size= 25),axis.text.y = element_text(size= 20),axis.title.y = element_text(size= 25))+
	#+guides(fill=FALSE)+
	theme(legend.title=element_blank(),legend.text = element_text(size = 15),
        axis.text.x = element_text(size= 15))+
        theme_classic()+
        coord_polar(theta = "y")+
	theme(legend.title=element_blank())
	
pdf(file="hypo_hyper_DMRs_pie_new.pdf",width=15,height=8, onefile=FALSE)
ggarrange(p2,p1,ncol = 2, nrow = 1,common.legend = FALSE,legend="bottom")
dev.off()
