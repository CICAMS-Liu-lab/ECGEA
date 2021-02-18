library(ggplot2)
library(gridExtra)
library(ggpubr)
library(grid)

qc1<-read.csv("expression-ESCC_COCA_cluster",header = T,  check.names = F,sep="\t")
#qc1<-read.csv("methylation-ESCC_COCA_cluster",header = T,  check.names = F,sep="\t")
names(qc1)[298] <- "iCluster"
path <- unlist(strsplit(basename(fin), split="\\."))[1]
for (s in 1:5) {
sym=c("NFE2L2_Nrf2 pathway", "KEAP1_Nrf2 pathway", "CUL3_Nrf2 pathway", "TP63_GI differentiation pathway", "SOX2_GI differentiation pathway")
prefix <- names(qc1[sym[s]])
gene <- unlist(strsplit(prefix, split="\\_"))[1]
max<-max(range(qc1[sym[s]],na.rm=T))
qc1$genename<-qc1[sym[s]][,1]
compare <- compare_means(genename~iCluster, data=qc1,method = "anova")
merge<-merge(prefix,compare$p)

qc1$iCluster <- gsub("cluster", "C", qc1$iCluster)
pdf(file=paste0(gene,"-COCA_155_4-cluster.pdf"),width=10,height=12)
p <- ggplot(qc1, aes(x=iCluster, y=qc1[sym[s]][,1],fill=iCluster)) + geom_boxplot(outlier.colour = NA)+
  xlab("")+ylab("CPM_log2FC")+ labs(title=prefix)+theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.title=element_blank(),legend.key.size=unit(2,'cm'))+geom_jitter(shape=16, position=position_jitter(0.2))+
  stat_compare_means(method = "anova",label.y =max+1, size=12)+
  ylim(-4, 6)+
  scale_fill_manual(values=c("cluster1"="#A020F0","cluster2"="#CDCD00","cluster3"="#DDA0DD","cluster4"="#66CDAA","cluster5"="blue",
                             "C1"="#CC0C00FF","C2"="#5C88DAFF","C3"="#008B8B","C4"="#84BD00FF","5"="blue","6"="white"))+
  theme(axis.line.y = element_line(colour="black"),
        axis.line.x = element_line(colour="black"),
        axis.text.x =element_text(size=30,colour = "black"),
        axis.text.y=element_text(size=24,colour = "black"),
        axis.title.y=element_text(size = 26, vjust=0),
        #axis.ticks.x = element_blank(),
        legend.text = element_text(size = 24),
        title = element_text(size = 24),
        panel.background = element_blank())

grid.draw(p)
dev.off()
}