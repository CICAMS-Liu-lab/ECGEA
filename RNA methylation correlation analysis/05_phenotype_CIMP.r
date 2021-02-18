library(ggplot2)
#annotation=read.table("All_phenotype_CIMP_split.conf",header=T,sep="\t",check.names=F)
annotation=read.table("GLM_0.7_6+.conf",header=T,sep="\t",check.names=F)
annotation=na.omit(annotation)
annotation$CIMP <- gsub("CIMP_1$","CIMP+",annotation$CIMP)
annotation$CIMP <- gsub("CIMP_2$","CIMP-",annotation$CIMP)
annotation$CIMP_biomarker <- gsub("^\\+$","CIMP+",annotation$CIMP_biomarker)
annotation$CIMP_biomarker <- gsub("^\\-$","CIMP-",annotation$CIMP_biomarker)
annotation$NN_info[which(annotation$N_info > 0)] <- "LNMets"
annotation$NN_info[which(annotation$N_info=="0")] <- "No LNMets"
annotation$TNM <- gsub("III$","III-IV",annotation$TNM)
annotation$TNM <- gsub("^IV","III-IV",annotation$TNM)
annotation$TNM <- gsub("^II$","I-II",annotation$TNM)
annotation$TNM <- gsub("^I$","I-II",annotation$TNM)

## CIMP vs TNM ###
freTable = as.data.frame(prop.table(table(annotation[,c("CIMP","TNM")]),1))
freTable$TNM=freTable[,2]
pvale=fisher.test(matrix(table(annotation[,c("CIMP","TNM")]),ncol=4,nrow=2))
pdf("CIMP_TNM.pdf",height=7,w=7)
ggplot(data=freTable, aes(x=CIMP, y=Freq, fill=TNM)) +
       geom_bar(stat="identity",width=0.4) +
       theme_bw()+
       theme(axis.text.x = element_text(angle = 0,color="black",size=25,hjust=0.5,lineheight=1))+
       theme(axis.text.y = element_text(color="black",size=17,hjust=1,lineheight=1),
       axis.title.x = element_blank(),axis.title.y = element_text(size=25),
       legend.position = "bottom",
       legend.text = element_text(size=15),legend.title = element_blank())+
       labs(title=paste("TNM", paste0("p_value = ",round(pvale$p.value,4)),sep="\n"),y="Frequence")+
       scale_fill_manual(values=c("LightSalmon1","SteelBlue"))
dev.off()

## CIMP vs N_info ###
freTable = as.data.frame(prop.table(table(annotation[,c("CIMP","NN_info")]),1))
freTable$NN_info=freTable[,2]
pvale=fisher.test(matrix(table(annotation[,c("CIMP","NN_info")]),ncol=4,nrow=2))
pdf("CIMP_N_info.pdf",height=7,w=7)
ggplot(data=freTable, aes(x=CIMP, y=Freq, fill=NN_info)) +
       geom_bar(stat="identity",width=0.4) +
       theme_bw()+
       theme(axis.text.x = element_text(angle = 0,color="black",size=25,hjust=0.5,lineheight=1))+
       theme(axis.text.y = element_text(color="black",size=17,hjust=1,lineheight=1),
       axis.title.x = element_blank(),axis.title.y = element_text(size=25),
       legend.position = "bottom",
       legend.text = element_text(size=15),legend.title = element_blank())+
       labs(title=paste("N_info", paste0("p_value = ",round(pvale$p.value,4)),sep="\n"),y="Frequence")+
       scale_fill_manual(values=c("IndianRed2","MediumSlateBlue"))
dev.off()

