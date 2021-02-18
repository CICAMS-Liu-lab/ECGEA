library(dplyr)
library(tidyverse)
library(stringr)
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(grid)

df_mut <- read_tsv("oncomatrix.xls")
df_cnv <- read_tsv("cnv_status.xls")
df_wgbs_promoter <- read_tsv("promoter_methylation_status.xls")
df_rna <- read_tsv("expression_status.xls")
coca_cluster <- read_tsv("ESCC_COCA_cluster.txt")

colnames(df_mut) <- gsub(" ","-",colnames(df_mut)) 
colnames(df_cnv) <- gsub(" ","-",colnames(df_cnv))
colnames(df_wgbs_promoter) <- gsub(" ","-",colnames(df_wgbs_promoter))
colnames(df_rna) <- gsub(" ","-",colnames(df_rna))

df_mut_long <- df_mut %>% gather(gene, mutation, -Sample_id)
df_cnv_long <- df_cnv %>% gather(gene, cnv, -Sample_id)
df_wgbs_promoter_long <- df_wgbs_promoter %>% gather(gene, wgbs_promoter, -Sample_id)
df_rna_long <- df_rna %>% gather(gene, rna, -Sample_id)

df_mut_rna <- df_rna_long %>% left_join(df_mut_long,by=c("Sample_id",'gene'))
df_mut_cnv_rna <- df_mut_rna %>% left_join(df_cnv_long,by=c("Sample_id",'gene'))
df_mut_cnv_wgbs_gene_promoter_rna <- df_mut_cnv_rna %>% left_join(df_wgbs_promoter_long,by=c("Sample_id",'gene'))
df_mut_cnv_wgbs_gene_promoter_rna_long <- df_mut_cnv_wgbs_gene_promoter_rna %>% gather(type,alter, rna:mutation:cnv:wgbs_promoter, -Sample_id) %>% separate(col = gene, into = c('gene_id','pathway'), sep = '_')

NRF2.pathway <- "NRF2_pathway"
Differentiation.pathway <- "Differentiation_pathway"
Cellcycle.pathway <- "Cell_cycle_pathway"
RTK.pathway <- "RTK-RAS-PI3K_pathway"
Proliferation.pathway <- "Proliferation_pathway"
Chromatin.pathway <- "Chromatin_pathway"

NRF2.gene <- c("NFE2L2","KEAP1","CUL3")
Differentiation.gene <- c("TP63","SOX2","NOTCH1","ZNF750")
Cellcycle.gene <- c("TP53","CDKN2A","CCND1","CDK6","CCNE1","RB1")
RTK.gene <- c("ERBB2","EGFR","FGFR1","FGFR2","KRAS","PTEN","PIK3CA","MET")
Proliferation.gene <- c("MYC","SMAD4","SMAD2","FBXW7","APC","CTNNB1","PTCH1")
Chromatin.gene <- c("SMARCA4","KDM6A","KMT2D","KMT2C","PBRM1","ARID1A")
gene <- c(NRF2.gene, Differentiation.gene, Cellcycle.gene, RTK.gene, Proliferation.gene, Chromatin.gene)

for(n in gene){
  df_mut_cnv_wgbs_gene_promoter_rna_long$gene_id <- gsub(paste0("^",n,"$"),paste0(n,".pathway",sep=""),df_mut_cnv_wgbs_gene_promoter_rna_long$gene_id)
}

df_path <- df_mut_cnv_wgbs_gene_promoter_rna_long %>% filter(str_detect(gene_id,pattern="pathway")) %>% mutate(hgt=case_when(type=="cnv" ~ 0.90, type=="mutation" ~ 0.60, type=="wgbs_promoter" ~ 0.30, type=="rna" ~ 0.90))

NRF2.term <- paste0(NRF2.gene,".pathway",sep="")
NRF2.df_path <- df_path %>% filter(gene_id %in% NRF2.term)
Differentiation.term <- paste0(Differentiation.gene,".pathway",sep="")
Differentiation.df_path <- df_path %>% filter(gene_id %in% Differentiation.term)
Cellcycle.term <- paste0(Cellcycle.gene,".pathway",sep="")
Cellcycle.df_path <- df_path %>% filter(gene_id %in% Cellcycle.term)
RTK.term <- paste0(RTK.gene,".pathway",sep="")
RTK.df_path <- df_path %>% filter(gene_id %in% RTK.term)
Proliferation.term <- paste0(Proliferation.gene,".pathway",sep="")
Proliferation.df_path <- df_path %>% filter(gene_id %in% Proliferation.term)
Chromatin.term <- paste0(Chromatin.gene,".pathway",sep="")
Chromatin.df_path <- df_path %>% filter(gene_id %in% Chromatin.term)

NRF2.gene_order_all <- NRF2.df_path %>% filter(type=="rna"|type=="wgbs_promoter"|type=="mutation"|type=="cnv") %>% filter(alter=="NA"|alter=="0") %>% group_by(gene_id)%>%summarise(n())%>%arrange(desc(`n()`)) %>%select(gene_id) %>%unlist()
Differentiation.gene_order_all <- Differentiation.df_path %>% filter(type=="rna"|type=="wgbs_promoter"|type=="mutation"|type=="cnv") %>% filter(alter=="NA"|alter=="0") %>% group_by(gene_id)%>%summarise(n())%>%arrange(desc(`n()`)) %>%select(gene_id) %>%unlist()
Cellcycle.gene_order_all <- Cellcycle.df_path %>% filter(type=="rna"|type=="wgbs_promoter"|type=="mutation"|type=="cnv") %>% filter(alter=="NA"|alter=="0") %>% group_by(gene_id)%>%summarise(n())%>%arrange(desc(`n()`)) %>%select(gene_id) %>%unlist()
RTK.gene_order_all <- RTK.df_path %>% filter(type=="rna"|type=="wgbs_promoter"|type=="mutation"|type=="cnv") %>% filter(alter=="NA"|alter=="0") %>% group_by(gene_id)%>%summarise(n())%>%arrange(desc(`n()`)) %>%select(gene_id) %>%unlist()
Proliferation.gene_order_all <- Proliferation.df_path %>% filter(type=="rna"|type=="wgbs_promoter"|type=="mutation"|type=="cnv") %>% filter(alter=="NA"|alter=="0") %>% group_by(gene_id)%>%summarise(n())%>%arrange(desc(`n()`)) %>%select(gene_id) %>%unlist()
Chromatin.gene_order_all <- Chromatin.df_path %>% filter(type=="rna"|type=="wgbs_promoter"|type=="mutation"|type=="cnv") %>% filter(alter=="NA"|alter=="0") %>% group_by(gene_id)%>%summarise(n())%>%arrange(desc(`n()`)) %>%select(gene_id) %>%unlist()

Cellcycle.gene_order_all[6] <- "TP53.pathway"
Cellcycle.gene_order_all[5] <- "CDKN2A.pathway"
Cellcycle.gene_order_all[4] <- "CCND1.pathway"
Cellcycle.gene_order_all[3] <- "CCNE1.pathway"
Cellcycle.gene_order_all[2] <- "CDK6.pathway"
Cellcycle.gene_order_all[1] <- "RB1.pathway"

NRF2.gene_order <- gsub(".pathway","",NRF2.gene_order_all)
Differentiation.gene_order <- gsub(".pathway","",Differentiation.gene_order_all)
Cellcycle.gene_order <- gsub(".pathway","",Cellcycle.gene_order_all)
RTK.gene_order <- gsub(".pathway","",RTK.gene_order_all)
Proliferation.gene_order <- gsub(".pathway","",Proliferation.gene_order_all)
Chromatin.gene_order <- gsub(".pathway","",Chromatin.gene_order_all)

myfunction <- function(gene,gene_order_all,df_mut,df_cnv,df_wgbs_promoter){
  sample_order_cluster <- c()
  for(n in gene){
    colnames(df_mut) <- gsub(n,paste0(n,".pathway",sep=""),colnames(df_mut))
    colnames(df_cnv) <- gsub(n,paste0(n,".pathway",sep=""),colnames(df_cnv))
    colnames(df_wgbs_promoter) <- gsub(n,paste0(n,".pathway",sep=""),colnames(df_wgbs_promoter))
  }
  df_mut_sam <- df_mut %>% select(Sample_id, contains("pathway"))
  
  df_mut_sam <- df_mut_sam %>% mutate(ATG7.pathway=NA)
  df_mut_sam <- df_mut_sam %>% mutate(SOX2.pathway=NA)
  df_mut_sam <- df_mut_sam %>% mutate(CCND1.pathway=NA)
  df_mut_sam <- df_mut_sam %>% mutate(MET.pathway=NA)
  df_mut_sam <- df_mut_sam %>% mutate(VEGFA.pathway=NA)
  df_mut_sam <- df_mut_sam %>% mutate(FGFR2.pathway=NA)
  df_mut_sam <- df_mut_sam %>% mutate(MYC.pathway=NA)
  df_mut_sam <- df_mut_sam %>% mutate(SMAD2.pathway=NA)
  df_mut_sam <- df_mut_sam %>% mutate(CTNNB1.pathway=NA)
  
  df_cnv_sam <- df_cnv %>% select(Sample_id, contains("pathway"))
  df_wgbs_promoter_sam <- df_wgbs_promoter %>% select(Sample_id, contains("pathway"))
  df_sam_tmp <- df_mut_sam %>% inner_join(df_cnv_sam, by=c("Sample_id"="Sample_id"))
  df_sam_tmp <- df_sam_tmp %>% inner_join(df_wgbs_promoter_sam, by=c("Sample_id"="Sample_id"))
  df_sam_tmp <- df_sam_tmp %>% mutate_at(vars(matches("pathway")), function(x) factor(x, levels=c(1,-1,0)))
  sample_order <- df_sam_tmp %>% arrange_(.dots=paste0("`",c(paste(rep(rev(gene_order_all),each=2),c(paste('y',sep = '.'), paste('x',sep = '.')), sep="."),rev(as.character(gene_order_all))), "`")) %>% select(Sample_id) %>%unlist() 
  sample_order <- as.character(sample_order)
  sample_order_cluster <- c(sample_order_cluster, sample_order)
  return(sample_order_cluster)
}

NRF2.sample_order_cluster <- myfunction(NRF2.gene,NRF2.gene_order_all,df_mut,df_cnv,df_wgbs_promoter)
Differentiation.sample_order_cluster <- myfunction(Differentiation.gene,Differentiation.gene_order_all,df_mut,df_cnv,df_wgbs_promoter)
Cellcycle.sample_order_cluster <- myfunction(Cellcycle.gene,Cellcycle.gene_order_all,df_mut,df_cnv,df_wgbs_promoter)
RTK.sample_order_cluster <- myfunction(RTK.gene,RTK.gene_order_all,df_mut,df_cnv,df_wgbs_promoter)
Proliferation.sample_order_cluster <- myfunction(Proliferation.gene,Proliferation.gene_order_all,df_mut,df_cnv,df_wgbs_promoter)
Chromatin.sample_order_cluster <- myfunction(Chromatin.gene,Chromatin.gene_order_all,df_mut,df_cnv,df_wgbs_promoter)

for(n in NRF2.gene){
  NRF2.df_path$gene_id <- gsub(paste0(n,".pathway",sep=""),n,NRF2.df_path$gene_id)
}
NRF2.df_sam_order <- NRF2.df_path %>% mutate(Sample_id=factor(Sample_id, levels=NRF2.sample_order_cluster))%>% mutate(gene_id=factor(gene_id, levels=NRF2.gene_order))
NRF2.cc <- NRF2.df_sam_order %>% filter(str_detect(type,pattern = "rna"))
NRF2.cc$alter <- gsub("-3","#2d00f7",NRF2.cc$alter)
NRF2.cc$alter <- gsub("3","#f20089",NRF2.cc$alter)
NRF2.cc$alter[is.na(NRF2.cc$alter)] <- "#e9ecef"

for(n in Differentiation.gene){
  Differentiation.df_path$gene_id <- gsub(paste0(n,".pathway",sep=""),n,Differentiation.df_path$gene_id)
}
Differentiation.df_sam_order <- Differentiation.df_path %>% mutate(Sample_id=factor(Sample_id, levels=Differentiation.sample_order_cluster))%>% mutate(gene_id=factor(gene_id, levels=Differentiation.gene_order))
Differentiation.cc <- Differentiation.df_sam_order %>% filter(str_detect(type,pattern = "rna"))
Differentiation.cc$alter <- gsub("-3","#2d00f7",Differentiation.cc$alter)
Differentiation.cc$alter <- gsub("3","#f20089",Differentiation.cc$alter)
Differentiation.cc$alter[is.na(Differentiation.cc$alter)] <- "#e9ecef"

for(n in Cellcycle.gene){
  Cellcycle.df_path$gene_id <- gsub(paste0(n,".pathway",sep=""),n,Cellcycle.df_path$gene_id)
}
Cellcycle.df_sam_order <- Cellcycle.df_path %>% mutate(Sample_id=factor(Sample_id, levels=Cellcycle.sample_order_cluster))%>% mutate(gene_id=factor(gene_id, levels=Cellcycle.gene_order))
Cellcycle.cc <- Cellcycle.df_sam_order %>% filter(str_detect(type,pattern = "rna"))
Cellcycle.cc$alter <- gsub("-3","#2d00f7",Cellcycle.cc$alter)
Cellcycle.cc$alter <- gsub("3","#f20089",Cellcycle.cc$alter)
Cellcycle.cc$alter[is.na(Cellcycle.cc$alter)] <- "#e9ecef"

for(n in RTK.gene){
  RTK.df_path$gene_id <- gsub(paste0(n,".pathway",sep=""),n,RTK.df_path$gene_id)
}
RTK.df_sam_order <- RTK.df_path %>% mutate(Sample_id=factor(Sample_id, levels=RTK.sample_order_cluster))%>% mutate(gene_id=factor(gene_id, levels=RTK.gene_order))
RTK.cc <- RTK.df_sam_order %>% filter(str_detect(type,pattern = "rna"))
RTK.cc$alter <- gsub("-3","#2d00f7",RTK.cc$alter)
RTK.cc$alter <- gsub("3","#f20089",RTK.cc$alter)
RTK.cc$alter[is.na(RTK.cc$alter)] <- "#e9ecef"

for(n in Proliferation.gene){
  Proliferation.df_path$gene_id <- gsub(paste0(n,".pathway",sep=""),n,Proliferation.df_path$gene_id)
}
Proliferation.df_sam_order <- Proliferation.df_path %>% mutate(Sample_id=factor(Sample_id, levels=Proliferation.sample_order_cluster))%>% mutate(gene_id=factor(gene_id, levels=Proliferation.gene_order))
Proliferation.cc <- Proliferation.df_sam_order %>% filter(str_detect(type,pattern = "rna"))
Proliferation.cc$alter <- gsub("-3","#2d00f7",Proliferation.cc$alter)
Proliferation.cc$alter <- gsub("3","#f20089",Proliferation.cc$alter)
Proliferation.cc$alter[is.na(Proliferation.cc$alter)] <- "#e9ecef"

for(n in Chromatin.gene){
  Chromatin.df_path$gene_id <- gsub(paste0(n,".pathway",sep=""),n,Chromatin.df_path$gene_id)
}
Chromatin.df_sam_order <- Chromatin.df_path %>% mutate(Sample_id=factor(Sample_id, levels=Chromatin.sample_order_cluster))%>% mutate(gene_id=factor(gene_id, levels=Chromatin.gene_order))
Chromatin.cc <- Chromatin.df_sam_order %>% filter(str_detect(type,pattern = "rna"))
Chromatin.cc$alter <- gsub("-3","#2d00f7",Chromatin.cc$alter)
Chromatin.cc$alter <- gsub("3","#f20089",Chromatin.cc$alter)
Chromatin.cc$alter[is.na(Chromatin.cc$alter)] <- "#e9ecef"

NRF2.df_sam_order_tmp <- NRF2.df_sam_order
NRF2.gene_order_tmp <- NRF2.gene_order
for (i in 1:length(unique(NRF2.df_sam_order$gene_id))) {
  gene <- unique(NRF2.df_sam_order$gene_id)[i]
  tmp <- NRF2.df_sam_order[,-6]
  tmp <- spread(tmp, type, alter)
  gene_frq <- tmp %>% filter(gene_id == gene)
  gene_frq_ratio <- (1-(dim(gene_frq[which(is.na(gene_frq$mutation) & gene_frq$cnv == 0 & is.na(gene_frq$wgbs_promoter)),])[1]/155))*100
  paste0(gene,"(",round(gene_frq_ratio, 0),"%)",sep="")
  NRF2.df_sam_order_tmp$gene_id <- gsub(gene, paste0(gene,"(",round(gene_frq_ratio, 0),"%)",sep=""), NRF2.df_sam_order_tmp$gene_id)
  NRF2.gene_order_tmp <- gsub(gene, paste0(gene,"(",round(gene_frq_ratio, 0),"%)",sep=""), NRF2.gene_order_tmp)
}
NRF2.df_sam_order_tmp$gene_id <-  factor(NRF2.df_sam_order_tmp$gene_id, levels=NRF2.gene_order_tmp)

NRF2.p = NRF2.df_sam_order_tmp %>%
  ggplot(aes(y=gene_id,x= Sample_id,fill=alter,height=hgt))+
  geom_tile(data=NRF2.df_sam_order_tmp%>%filter(type=="cnv"))+
  geom_tile(data=NRF2.df_sam_order_tmp%>%filter(type=="mutation"))+
  geom_tile(data=NRF2.df_sam_order_tmp%>%filter(type=="wgbs_promoter"))+
  xlab("")+ylab("")+
scale_fill_manual(values=c("0"="#e9ecef","NA"="none","Missense_Mutation"="#4B7AB1", 
                           "Nonsense_Mutation"="#ED9147",
                           "Frame_Shift_Ins"="#96B458",
                           "Frame_Shift_Del"="#824CA7",
                           "In_Frame_Ins"="#00A087FF",
                           "In_Frame_Del"="#B03060",
                           "Splice_Site"="#006400",
                           "Nonstop_Mutation"="red",
                           "Multi_Hit"="#7E6148FF",
                           "1"="#F17C67","-1"="#4cb4e7",
                           "2"="#FF0080","-2"="#0080FF"))+
  theme(
    axis.text.x=element_blank(),
    axis.text.y = element_text(color="black",size=18),
    axis.ticks = element_blank(),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    legend.position = "none",
    legend.title=element_blank())+
  theme(plot.title = element_text(hjust = 0.5,size=20))

NRF2.df_sam_order_onco <- NRF2.df_sam_order_tmp[,c(1,2,5)]
NRF2.df_sam_order_count = NRF2.df_sam_order_onco %>% 
  group_by(gene_id, alter) %>% 
  summarise(n= n())
NRF2.df_sam_order_count <- NRF2.df_sam_order_count %>% filter(alter != 0)
NRF2.gene_order_all <- gsub(".pathway", "", NRF2.gene_order_all)
NRF2.df_sam_order_count$gene_id <- factor(NRF2.df_sam_order_count$gene_id, levels=(NRF2.gene_order_all))
NRF2.df_sam_order_count$alter <- factor(NRF2.df_sam_order_count$alter, levels = rev(c(1, -1, 2, -2, "Missense_Mutation",
                                                                                      "Nonsense_Mutation", "Frame_Shift_Ins",
                                                                                      "Frame_Shift_Del", "In_Frame_Ins",
                                                                                      "In_Frame_Del", "Splice_Site",
                                                                                      "Nonstop_Mutation", "Multi_Hit", 3, -3)))
NRF2.df_sam_order_count_cnv <- NRF2.df_sam_order_count %>% filter(alter == -1 | alter == 1)
NRF2.df_sam_order_count_mut <- NRF2.df_sam_order_count %>% filter(alter != -1 & alter != 1 & alter != -2 & alter != 2 & alter != -3 & alter != 3)
NRF2.df_sam_order_count_rna <- NRF2.df_sam_order_count %>% filter(alter == -3 | alter == 3)
NRF2.df_sam_order_count_methyl <- NRF2.df_sam_order_count %>% filter(alter == -2 | alter == 2)

NRF2_p_count_cnv = ggplot(NRF2.df_sam_order_count_cnv, aes(x=gene_id, y=n, fill=alter))+
  geom_bar(stat="identity")+
  scale_fill_manual(values=c("Missense_Mutation"="#4B7AB1", "Nonsense_Mutation"="#ED9147", "Frame_Shift_Ins"="#96B458",
                             "Frame_Shift_Del"="#824CA7", "In_Frame_Ins"="#00A087FF", "In_Frame_Del"="#B03060",
                             "Splice_Site"="#006400", "Nonstop_Mutation"="red", "Multi_Hit"="#7E6148FF",
                             "-3"="#2d00f7", "3"="#f20089", "1"="#F17C67", "-1"="#4cb4e7", "2"="#FF0080", "-2"="#0080FF"))+
  coord_flip()+
  ylim(c(0,60))+
  xlab("")+
  ylab("")+
  theme(legend.position = "none",
        axis.line = element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y=element_blank(),
        panel.background = element_blank(),
        axis.ticks.y=element_blank(),
        axis.ticks.x=element_blank())

NRF2_p_count_mut = ggplot(NRF2.df_sam_order_count_mut, aes(x=gene_id, y=n, fill=alter))+
  geom_bar(stat="identity")+
  scale_fill_manual(values=c("Missense_Mutation"="#4B7AB1", "Nonsense_Mutation"="#ED9147", "Frame_Shift_Ins"="#96B458",
                             "Frame_Shift_Del"="#824CA7", "In_Frame_Ins"="#00A087FF", "In_Frame_Del"="#B03060",
                             "Splice_Site"="#006400", "Nonstop_Mutation"="red", "Multi_Hit"="#7E6148FF",
                             "-3"="#2d00f7", "3"="#f20089", "1"="#F17C67", "-1"="#4cb4e7", "2"="#FF0080", "-2"="#0080FF"))+
  coord_flip()+
  ylim(c(0,150))+
  xlab("")+
  ylab("")+
  theme(legend.position = "none",
        axis.line = element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y=element_blank(),
        panel.background = element_blank(),
        axis.ticks.y=element_blank(),
        axis.ticks.x=element_blank())

NRF2_p_count_rna = ggplot(NRF2.df_sam_order_count_rna, aes(x=gene_id, y=n, fill=alter))+
  geom_bar(stat="identity")+
  scale_fill_manual(values=c("Missense_Mutation"="#4B7AB1", "Nonsense_Mutation"="#ED9147", "Frame_Shift_Ins"="#96B458",
                             "Frame_Shift_Del"="#824CA7", "In_Frame_Ins"="#00A087FF", "In_Frame_Del"="#B03060",
                             "Splice_Site"="#006400", "Nonstop_Mutation"="red", "Multi_Hit"="#7E6148FF",
                             "-3"="#2d00f7", "3"="#f20089", "1"="#F17C67", "-1"="#4cb4e7", "2"="#FF0080", "-2"="#0080FF"))+
  coord_flip()+
  ylim(c(0,100))+
  xlab("")+
  ylab("")+
  theme(legend.position = "none",
        axis.line = element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y=element_blank(),
        panel.background = element_blank(),
        axis.ticks.y=element_blank(),
        axis.ticks.x=element_blank())

NRF2_p_count_methyl = ggplot(NRF2.df_sam_order_count_methyl, aes(x=gene_id, y=n, fill=alter))+
  geom_bar(stat="identity")+
  scale_fill_manual(values=c("Missense_Mutation"="#4B7AB1", "Nonsense_Mutation"="#ED9147", "Frame_Shift_Ins"="#96B458",
                             "Frame_Shift_Del"="#824CA7", "In_Frame_Ins"="#00A087FF", "In_Frame_Del"="#B03060",
                             "Splice_Site"="#006400", "Nonstop_Mutation"="red", "Multi_Hit"="#7E6148FF",
                             "-3"="#2d00f7", "3"="#f20089", "1"="#F17C67", "-1"="#4cb4e7", "2"="#FF0080", "-2"="#0080FF"))+
  coord_flip()+
  ylim(c(0,30))+
  xlab("")+
  ylab("")+
  theme(legend.position = "none",
        axis.line = element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y=element_blank(),
        panel.background = element_blank(),
        axis.ticks.y=element_blank(),
        axis.ticks.x=element_blank())

Differentiation.df_sam_order_tmp <- Differentiation.df_sam_order
Differentiation.gene_order_tmp <- Differentiation.gene_order
for (i in 1:length(unique(Differentiation.df_sam_order$gene_id))) {
  gene <- unique(Differentiation.df_sam_order$gene_id)[i]
  tmp <- Differentiation.df_sam_order[,-6]
  tmp <- spread(tmp, type, alter)
  gene_frq <- tmp %>% filter(gene_id == gene)
  gene_frq_ratio <- (1-(dim(gene_frq[which(is.na(gene_frq$mutation) & gene_frq$cnv == 0 & is.na(gene_frq$wgbs_promoter)),])[1]/155))*100
  paste0(gene,"(",round(gene_frq_ratio, 0),"%)",sep="")
  Differentiation.df_sam_order_tmp$gene_id <- gsub(gene, paste0(gene,"(",round(gene_frq_ratio, 0),"%)",sep=""), Differentiation.df_sam_order_tmp$gene_id)
  Differentiation.gene_order_tmp <- gsub(gene, paste0(gene,"(",round(gene_frq_ratio, 0),"%)",sep=""), Differentiation.gene_order_tmp)
}
Differentiation.df_sam_order_tmp$gene_id <-  factor(Differentiation.df_sam_order_tmp$gene_id, levels=Differentiation.gene_order_tmp)

Differentiation.p = Differentiation.df_sam_order_tmp %>%
  ggplot(aes(y=gene_id,x= Sample_id,fill=alter,height=hgt))+
  geom_tile(data=Differentiation.df_sam_order_tmp%>%filter(type=="cnv"))+
  geom_tile(data=Differentiation.df_sam_order_tmp%>%filter(type=="mutation"))+
  geom_tile(data=Differentiation.df_sam_order_tmp%>%filter(type=="wgbs_promoter"))+
  xlab("")+ylab("")+
  scale_fill_manual(values=c("0"="#e9ecef","NA"="none","Missense_Mutation"="#4B7AB1", 
                             "Nonsense_Mutation"="#ED9147",
                             "Frame_Shift_Ins"="#96B458",
                             "Frame_Shift_Del"="#824CA7",
                             "In_Frame_Ins"="#00A087FF",
                             "In_Frame_Del"="#B03060",
                             "Splice_Site"="#006400",
                             "Nonstop_Mutation"="red",
                             "Multi_Hit"="#7E6148FF",
                             "1"="#F17C67","-1"="#4cb4e7",
                             "2"="#FF0080","-2"="#0080FF"))+
  theme(
    axis.text.x=element_blank(),
    axis.text.y = element_text(color="black",size=18),
    axis.ticks = element_blank(),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    legend.position = "none",
    legend.title=element_blank())+
  theme(plot.title = element_text(hjust = 0.5,size=20))

Differentiation.df_sam_order_onco <- Differentiation.df_sam_order_tmp[,c(1,2,5)]
Differentiation.df_sam_order_count = Differentiation.df_sam_order_onco %>% 
  group_by(gene_id, alter) %>% 
  summarise(n= n())
Differentiation.df_sam_order_count <- Differentiation.df_sam_order_count %>% filter(alter != 0)
Differentiation.gene_order_all <- gsub(".pathway", "", Differentiation.gene_order_all)
Differentiation.df_sam_order_count$gene_id <- factor(Differentiation.df_sam_order_count$gene_id, levels=(Differentiation.gene_order_all))
Differentiation.df_sam_order_count$alter <- factor(Differentiation.df_sam_order_count$alter, levels = rev(c(1, -1, 2, -2, "Missense_Mutation",
                                                                                                            "Nonsense_Mutation", "Frame_Shift_Ins",
                                                                                                            "Frame_Shift_Del", "In_Frame_Ins",
                                                                                                            "In_Frame_Del", "Splice_Site",
                                                                                                            "Nonstop_Mutation", "Multi_Hit", 3, -3)))
Differentiation.df_sam_order_count_cnv <- Differentiation.df_sam_order_count %>% filter(alter == -1 | alter == 1)
Differentiation.df_sam_order_count_mut <- Differentiation.df_sam_order_count %>% filter(alter != -1 & alter != 1 & alter != -2 & alter != 2 & alter != -3 & alter != 3)
Differentiation.df_sam_order_count_rna <- Differentiation.df_sam_order_count %>% filter(alter == -3 | alter == 3)
Differentiation.df_sam_order_count_methyl <- Differentiation.df_sam_order_count %>% filter(alter == -2 | alter == 2)

Differentiation_p_count_cnv = ggplot(Differentiation.df_sam_order_count_cnv, aes(x=gene_id, y=n, fill=alter))+
  geom_bar(stat="identity")+
  scale_fill_manual(values=c("Missense_Mutation"="#4B7AB1", "Nonsense_Mutation"="#ED9147", "Frame_Shift_Ins"="#96B458",
                             "Frame_Shift_Del"="#824CA7", "In_Frame_Ins"="#00A087FF", "In_Frame_Del"="#B03060",
                             "Splice_Site"="#006400", "Nonstop_Mutation"="red", "Multi_Hit"="#7E6148FF",
                             "-3"="#2d00f7", "3"="#f20089", "1"="#F17C67", "-1"="#4cb4e7", "2"="#FF0080", "-2"="#0080FF"))+
  coord_flip()+
  ylim(c(0,60))+
  xlab("")+
  ylab("")+
  theme(legend.position = "none",
        axis.line = element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y=element_blank(),
        panel.background = element_blank(),
        axis.ticks.y=element_blank(),
        axis.ticks.x=element_blank())

Differentiation_p_count_mut = ggplot(Differentiation.df_sam_order_count_mut, aes(x=gene_id, y=n, fill=alter))+
  geom_bar(stat="identity")+
  scale_fill_manual(values=c("Missense_Mutation"="#4B7AB1", "Nonsense_Mutation"="#ED9147", "Frame_Shift_Ins"="#96B458",
                             "Frame_Shift_Del"="#824CA7", "In_Frame_Ins"="#00A087FF", "In_Frame_Del"="#B03060",
                             "Splice_Site"="#006400", "Nonstop_Mutation"="red", "Multi_Hit"="#7E6148FF",
                             "-3"="#2d00f7", "3"="#f20089", "1"="#F17C67", "-1"="#4cb4e7", "2"="#FF0080", "-2"="#0080FF"))+
  coord_flip()+
  ylim(c(0,150))+
  xlab("")+
  ylab("")+
  theme(legend.position = "none",
        axis.line = element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y=element_blank(),
        panel.background = element_blank(),
        axis.ticks.y=element_blank(),
        axis.ticks.x=element_blank())

Differentiation_p_count_rna = ggplot(Differentiation.df_sam_order_count_rna, aes(x=gene_id, y=n, fill=alter))+
  geom_bar(stat="identity")+
  scale_fill_manual(values=c("Missense_Mutation"="#4B7AB1", "Nonsense_Mutation"="#ED9147", "Frame_Shift_Ins"="#96B458",
                             "Frame_Shift_Del"="#824CA7", "In_Frame_Ins"="#00A087FF", "In_Frame_Del"="#B03060",
                             "Splice_Site"="#006400", "Nonstop_Mutation"="red", "Multi_Hit"="#7E6148FF",
                             "-3"="#2d00f7", "3"="#f20089", "1"="#F17C67", "-1"="#4cb4e7", "2"="#FF0080", "-2"="#0080FF"))+
  coord_flip()+
  ylim(c(0,100))+
  xlab("")+
  ylab("")+
  theme(legend.position = "none",
        axis.line = element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y=element_blank(),
        panel.background = element_blank(),
        axis.ticks.y=element_blank(),
        axis.ticks.x=element_blank())

Differentiation_p_count_methyl = ggplot(Differentiation.df_sam_order_count_methyl, aes(x=gene_id, y=n, fill=alter))+
  geom_bar(stat="identity")+
  scale_fill_manual(values=c("Missense_Mutation"="#4B7AB1", "Nonsense_Mutation"="#ED9147", "Frame_Shift_Ins"="#96B458",
                             "Frame_Shift_Del"="#824CA7", "In_Frame_Ins"="#00A087FF", "In_Frame_Del"="#B03060",
                             "Splice_Site"="#006400", "Nonstop_Mutation"="red", "Multi_Hit"="#7E6148FF",
                             "-3"="#2d00f7", "3"="#f20089", "1"="#F17C67", "-1"="#4cb4e7", "2"="#FF0080", "-2"="#0080FF"))+
  coord_flip()+
  ylim(c(0,30))+
  xlab("")+
  ylab("")+
  theme(legend.position = "none",
        axis.line = element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y=element_blank(),
        panel.background = element_blank(),
        axis.ticks.y=element_blank(),
        axis.ticks.x=element_blank())

Cellcycle.df_sam_order_tmp <- Cellcycle.df_sam_order
Cellcycle.gene_order_tmp <- Cellcycle.gene_order
for (i in 1:length(unique(Cellcycle.df_sam_order$gene_id))) {
  gene <- unique(Cellcycle.df_sam_order$gene_id)[i]
  tmp <- Cellcycle.df_sam_order[,-6]
  tmp <- spread(tmp, type, alter)
  gene_frq <- tmp %>% filter(gene_id == gene)
  gene_frq_ratio <- (1-(dim(gene_frq[which(is.na(gene_frq$mutation) & gene_frq$cnv == 0 & is.na(gene_frq$wgbs_promoter)),])[1]/155))*100
  paste0(gene,"(",round(gene_frq_ratio, 0),"%)",sep="")
  Cellcycle.df_sam_order_tmp$gene_id <- gsub(gene, paste0(gene,"(",round(gene_frq_ratio, 0),"%)",sep=""), Cellcycle.df_sam_order_tmp$gene_id)
  Cellcycle.gene_order_tmp <- gsub(gene, paste0(gene,"(",round(gene_frq_ratio, 0),"%)",sep=""), Cellcycle.gene_order_tmp)
}
Cellcycle.df_sam_order_tmp$gene_id <-  factor(Cellcycle.df_sam_order_tmp$gene_id, levels=Cellcycle.gene_order_tmp)

Cellcycle.p = Cellcycle.df_sam_order_tmp %>%
  ggplot(aes(y=gene_id,x= Sample_id,fill=alter,height=hgt))+
  geom_tile(data=Cellcycle.df_sam_order_tmp%>%filter(type=="cnv"))+
  geom_tile(data=Cellcycle.df_sam_order_tmp%>%filter(type=="mutation"))+
  geom_tile(data=Cellcycle.df_sam_order_tmp%>%filter(type=="wgbs_promoter"))+
  xlab("")+ylab("")+
  scale_fill_manual(values=c("0"="#e9ecef","NA"="none","Missense_Mutation"="#4B7AB1", 
                             "Nonsense_Mutation"="#ED9147",
                             "Frame_Shift_Ins"="#96B458",
                             "Frame_Shift_Del"="#824CA7",
                             "In_Frame_Ins"="#00A087FF",
                             "In_Frame_Del"="#B03060",
                             "Splice_Site"="#006400",
                             "Nonstop_Mutation"="red",
                             "Multi_Hit"="#7E6148FF",
                             "1"="#F17C67","-1"="#4cb4e7",
                             "2"="#FF0080","-2"="#0080FF"))+
  theme(
    axis.text.x=element_blank(),
    axis.text.y = element_text(color="black",size=18),
    axis.ticks = element_blank(),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    legend.position = "none",
    legend.title=element_blank())+
  theme(plot.title = element_text(hjust = 0.5,size=20))

Cellcycle.df_sam_order_onco <- Cellcycle.df_sam_order[,c(1,2,5)]
Cellcycle.df_sam_order_count = Cellcycle.df_sam_order_onco %>% 
  group_by(gene_id, alter) %>% 
  summarise(n= n())
Cellcycle.df_sam_order_count <- Cellcycle.df_sam_order_count %>% filter(alter != 0)
Cellcycle.gene_order_all <- gsub(".pathway", "", Cellcycle.gene_order_all)
Cellcycle.df_sam_order_count$gene_id <- factor(Cellcycle.df_sam_order_count$gene_id, levels=(Cellcycle.gene_order_all))
Cellcycle.df_sam_order_count$alter <- factor(Cellcycle.df_sam_order_count$alter, levels = rev(c(1, -1, 2, -2, "Missense_Mutation",
                                                                                                "Nonsense_Mutation", "Frame_Shift_Ins",
                                                                                                "Frame_Shift_Del", "In_Frame_Ins",
                                                                                                "In_Frame_Del", "Splice_Site",
                                                                                                "Nonstop_Mutation", "Multi_Hit", 3, -3)))
Cellcycle.df_sam_order_count_cnv <- Cellcycle.df_sam_order_count %>% filter(alter == -1 | alter == 1)
Cellcycle.df_sam_order_count_mut <- Cellcycle.df_sam_order_count %>% filter(alter != -1 & alter != 1 & alter != -2 & alter != 2 & alter != -3 & alter != 3)
Cellcycle.df_sam_order_count_rna <- Cellcycle.df_sam_order_count %>% filter(alter == -3 | alter == 3)
Cellcycle.df_sam_order_count_methyl <- Cellcycle.df_sam_order_count %>% filter(alter == -2 | alter == 2)

Cellcycle_p_count_cnv = ggplot(Cellcycle.df_sam_order_count_cnv, aes(x=gene_id, y=n, fill=alter))+
  geom_bar(stat="identity")+
  scale_fill_manual(values=c("Missense_Mutation"="#4B7AB1", "Nonsense_Mutation"="#ED9147", "Frame_Shift_Ins"="#96B458",
                             "Frame_Shift_Del"="#824CA7", "In_Frame_Ins"="#00A087FF", "In_Frame_Del"="#B03060",
                             "Splice_Site"="#006400", "Nonstop_Mutation"="red", "Multi_Hit"="#7E6148FF",
                             "-3"="#2d00f7", "3"="#f20089", "1"="#F17C67", "-1"="#4cb4e7", "2"="#FF0080", "-2"="#0080FF"))+
  coord_flip()+
  ylim(c(0,60))+
  xlab("")+
  ylab("")+
  theme(legend.position = "none",
        axis.line = element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y=element_blank(),
        panel.background = element_blank(),
        axis.ticks.y=element_blank(),
        axis.ticks.x=element_blank())

Cellcycle_p_count_mut = ggplot(Cellcycle.df_sam_order_count_mut, aes(x=gene_id, y=n, fill=alter))+
  geom_bar(stat="identity")+
  scale_fill_manual(values=c("Missense_Mutation"="#4B7AB1", "Nonsense_Mutation"="#ED9147", "Frame_Shift_Ins"="#96B458",
                             "Frame_Shift_Del"="#824CA7", "In_Frame_Ins"="#00A087FF", "In_Frame_Del"="#B03060",
                             "Splice_Site"="#006400", "Nonstop_Mutation"="red", "Multi_Hit"="#7E6148FF",
                             "-3"="#2d00f7", "3"="#f20089", "1"="#F17C67", "-1"="#4cb4e7", "2"="#FF0080", "-2"="#0080FF"))+
  coord_flip()+
  ylim(c(0,150))+
  xlab("")+
  ylab("")+
  theme(legend.position = "none",
        axis.line = element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y=element_blank(),
        panel.background = element_blank(),
        axis.ticks.y=element_blank(),
        axis.ticks.x=element_blank())

Cellcycle_p_count_rna = ggplot(Cellcycle.df_sam_order_count_rna, aes(x=gene_id, y=n, fill=alter))+
  geom_bar(stat="identity")+
  scale_fill_manual(values=c("Missense_Mutation"="#4B7AB1", "Nonsense_Mutation"="#ED9147", "Frame_Shift_Ins"="#96B458",
                             "Frame_Shift_Del"="#824CA7", "In_Frame_Ins"="#00A087FF", "In_Frame_Del"="#B03060",
                             "Splice_Site"="#006400", "Nonstop_Mutation"="red", "Multi_Hit"="#7E6148FF",
                             "-3"="#2d00f7", "3"="#f20089", "1"="#F17C67", "-1"="#4cb4e7", "2"="#FF0080", "-2"="#0080FF"))+
  coord_flip()+
  ylim(c(0,100))+
  xlab("")+
  ylab("")+
  theme(legend.position = "none",
        axis.line = element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y=element_blank(),
        panel.background = element_blank(),
        axis.ticks.y=element_blank(),
        axis.ticks.x=element_blank())

Cellcycle_p_count_methyl = ggplot(Cellcycle.df_sam_order_count_methyl, aes(x=gene_id, y=n, fill=alter))+
  geom_bar(stat="identity")+
  scale_fill_manual(values=c("Missense_Mutation"="#4B7AB1", "Nonsense_Mutation"="#ED9147", "Frame_Shift_Ins"="#96B458",
                             "Frame_Shift_Del"="#824CA7", "In_Frame_Ins"="#00A087FF", "In_Frame_Del"="#B03060",
                             "Splice_Site"="#006400", "Nonstop_Mutation"="red", "Multi_Hit"="#7E6148FF",
                             "-3"="#2d00f7", "3"="#f20089", "1"="#F17C67", "-1"="#4cb4e7", "2"="#FF0080", "-2"="#0080FF"))+
  coord_flip()+
  ylim(c(0,30))+
  xlab("")+
  ylab("")+
  theme(legend.position = "none",
        axis.line = element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y=element_blank(),
        panel.background = element_blank(),
        axis.ticks.y=element_blank(),
        axis.ticks.x=element_blank())

RTK.df_sam_order_tmp <- RTK.df_sam_order
RTK.gene_order_tmp <- RTK.gene_order
for (i in 1:length(unique(RTK.df_sam_order$gene_id))) {
  gene <- unique(RTK.df_sam_order$gene_id)[i]
  tmp <- RTK.df_sam_order[,-6]
  tmp <- spread(tmp, type, alter)
  gene_frq <- tmp %>% filter(gene_id == gene)
  gene_frq_ratio <- (1-(dim(gene_frq[which(is.na(gene_frq$mutation) & gene_frq$cnv == 0 & is.na(gene_frq$wgbs_promoter)),])[1]/155))*100
  paste0(gene,"(",round(gene_frq_ratio, 0),"%)",sep="")
  RTK.df_sam_order_tmp$gene_id <- gsub(gene, paste0(gene,"(",round(gene_frq_ratio, 0),"%)",sep=""), RTK.df_sam_order_tmp$gene_id)
  RTK.gene_order_tmp <- gsub(gene, paste0(gene,"(",round(gene_frq_ratio, 0),"%)",sep=""), RTK.gene_order_tmp)
}
RTK.df_sam_order_tmp$gene_id <-  factor(RTK.df_sam_order_tmp$gene_id, levels=RTK.gene_order_tmp)

RTK.p = RTK.df_sam_order_tmp %>%
  ggplot(aes(y=gene_id,x= Sample_id,fill=alter,height=hgt))+
  geom_tile(data=RTK.df_sam_order_tmp%>%filter(type=="cnv"))+
  geom_tile(data=RTK.df_sam_order_tmp%>%filter(type=="mutation"))+
  geom_tile(data=RTK.df_sam_order_tmp%>%filter(type=="wgbs_promoter"))+
  xlab("")+ylab("")+
  scale_fill_manual(values=c("0"="#e9ecef","NA"="none","Missense_Mutation"="#4B7AB1", 
                             "Nonsense_Mutation"="#ED9147",
                             "Frame_Shift_Ins"="#96B458",
                             "Frame_Shift_Del"="#824CA7",
                             "In_Frame_Ins"="#00A087FF",
                             "In_Frame_Del"="#B03060",
                             "Splice_Site"="#006400",
                             "Nonstop_Mutation"="red",
                             "Multi_Hit"="#7E6148FF",
                             "1"="#F17C67","-1"="#4cb4e7",
                             "2"="#FF0080","-2"="#0080FF"))+
  theme(
    axis.text.x=element_blank(),
    axis.text.y = element_text(color="black",size=18),
    axis.ticks = element_blank(),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    legend.position = "none",
    legend.title=element_blank())+
  theme(plot.title = element_text(hjust = 0.5,size=20))

RTK.df_sam_order_onco <- RTK.df_sam_order[,c(1,2,5)]
RTK.df_sam_order_count = RTK.df_sam_order_onco %>% 
  group_by(gene_id, alter) %>% 
  summarise(n= n())
RTK.df_sam_order_count <- RTK.df_sam_order_count %>% filter(alter != 0)
RTK.gene_order_all <- gsub(".pathway", "", RTK.gene_order_all)
RTK.df_sam_order_count$gene_id <- factor(RTK.df_sam_order_count$gene_id, levels=(RTK.gene_order_all))
RTK.df_sam_order_count$alter <- factor(RTK.df_sam_order_count$alter, levels = rev(c(1, -1, 2, -2, "Missense_Mutation",
                                                                                    "Nonsense_Mutation", "Frame_Shift_Ins",
                                                                                    "Frame_Shift_Del", "In_Frame_Ins",
                                                                                    "In_Frame_Del", "Splice_Site",
                                                                                    "Nonstop_Mutation", "Multi_Hit", 3, -3)))
RTK.df_sam_order_count_cnv <- RTK.df_sam_order_count %>% filter(alter == -1 | alter == 1)
RTK.df_sam_order_count_mut <- RTK.df_sam_order_count %>% filter(alter != -1 & alter != 1 & alter != -2 & alter != 2 & alter != -3 & alter != 3)
RTK.df_sam_order_count_rna <- RTK.df_sam_order_count %>% filter(alter == -3 | alter == 3)
RTK.df_sam_order_count_methyl <- RTK.df_sam_order_count %>% filter(alter == -2 | alter == 2)

RTK_p_count_cnv = ggplot(RTK.df_sam_order_count_cnv, aes(x=gene_id, y=n, fill=alter))+
  geom_bar(stat="identity")+
  scale_fill_manual(values=c("Missense_Mutation"="#4B7AB1", "Nonsense_Mutation"="#ED9147", "Frame_Shift_Ins"="#96B458",
                             "Frame_Shift_Del"="#824CA7", "In_Frame_Ins"="#00A087FF", "In_Frame_Del"="#B03060",
                             "Splice_Site"="#006400", "Nonstop_Mutation"="red", "Multi_Hit"="#7E6148FF",
                             "-3"="#2d00f7", "3"="#f20089", "1"="#F17C67", "-1"="#4cb4e7", "2"="#FF0080", "-2"="#0080FF"))+
  coord_flip()+
  ylim(c(0,60))+
  xlab("")+
  ylab("")+
  theme(legend.position = "none",
        axis.line = element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y=element_blank(),
        panel.background = element_blank(),
        axis.ticks.y=element_blank(),
        axis.ticks.x=element_blank())

RTK_p_count_mut = ggplot(RTK.df_sam_order_count_mut, aes(x=gene_id, y=n, fill=alter))+
  geom_bar(stat="identity")+
  scale_fill_manual(values=c("Missense_Mutation"="#4B7AB1", "Nonsense_Mutation"="#ED9147", "Frame_Shift_Ins"="#96B458",
                             "Frame_Shift_Del"="#824CA7", "In_Frame_Ins"="#00A087FF", "In_Frame_Del"="#B03060",
                             "Splice_Site"="#006400", "Nonstop_Mutation"="red", "Multi_Hit"="#7E6148FF",
                             "-3"="#2d00f7", "3"="#f20089", "1"="#F17C67", "-1"="#4cb4e7", "2"="#FF0080", "-2"="#0080FF"))+
  coord_flip()+
  ylim(c(0,150))+
  xlab("")+
  ylab("")+
  theme(legend.position = "none",
        axis.line = element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y=element_blank(),
        panel.background = element_blank(),
        axis.ticks.y=element_blank(),
        axis.ticks.x=element_blank())

RTK_p_count_rna = ggplot(RTK.df_sam_order_count_rna, aes(x=gene_id, y=n, fill=alter))+
  geom_bar(stat="identity")+
  scale_fill_manual(values=c("Missense_Mutation"="#4B7AB1", "Nonsense_Mutation"="#ED9147", "Frame_Shift_Ins"="#96B458",
                             "Frame_Shift_Del"="#824CA7", "In_Frame_Ins"="#00A087FF", "In_Frame_Del"="#B03060",
                             "Splice_Site"="#006400", "Nonstop_Mutation"="red", "Multi_Hit"="#7E6148FF",
                             "-3"="#2d00f7", "3"="#f20089", "1"="#F17C67", "-1"="#4cb4e7", "2"="#FF0080", "-2"="#0080FF"))+
  coord_flip()+
  ylim(c(0,100))+
  xlab("")+
  ylab("")+
  theme(legend.position = "none",
        axis.line = element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y=element_blank(),
        panel.background = element_blank(),
        axis.ticks.y=element_blank(),
        axis.ticks.x=element_blank())

RTK_p_count_methyl = ggplot(RTK.df_sam_order_count_methyl, aes(x=gene_id, y=n, fill=alter))+
  geom_bar(stat="identity")+
  scale_fill_manual(values=c("Missense_Mutation"="#4B7AB1", "Nonsense_Mutation"="#ED9147", "Frame_Shift_Ins"="#96B458",
                             "Frame_Shift_Del"="#824CA7", "In_Frame_Ins"="#00A087FF", "In_Frame_Del"="#B03060",
                             "Splice_Site"="#006400", "Nonstop_Mutation"="red", "Multi_Hit"="#7E6148FF",
                             "-3"="#2d00f7", "3"="#f20089", "1"="#F17C67", "-1"="#4cb4e7", "2"="#FF0080", "-2"="#0080FF"))+
  coord_flip()+
  ylim(c(0,30))+
  xlab("")+
  ylab("")+
  theme(legend.position = "none",
        axis.line = element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y=element_blank(),
        panel.background = element_blank(),
        axis.ticks.y=element_blank(),
        axis.ticks.x=element_blank())

Proliferation.df_sam_order_tmp <- Proliferation.df_sam_order
Proliferation.gene_order_tmp <- Proliferation.gene_order
for (i in 1:length(unique(Proliferation.df_sam_order$gene_id))) {
  gene <- unique(Proliferation.df_sam_order$gene_id)[i]
  tmp <- Proliferation.df_sam_order[,-6]
  tmp <- spread(tmp, type, alter)
  gene_frq <- tmp %>% filter(gene_id == gene)
  gene_frq_ratio <- (1-(dim(gene_frq[which(is.na(gene_frq$mutation) & gene_frq$cnv == 0 & is.na(gene_frq$wgbs_promoter)),])[1]/155))*100
  paste0(gene,"(",round(gene_frq_ratio, 0),"%)",sep="")
  Proliferation.df_sam_order_tmp$gene_id <- gsub(gene, paste0(gene,"(",round(gene_frq_ratio, 0),"%)",sep=""), Proliferation.df_sam_order_tmp$gene_id)
  Proliferation.gene_order_tmp <- gsub(gene, paste0(gene,"(",round(gene_frq_ratio, 0),"%)",sep=""), Proliferation.gene_order_tmp)
}
Proliferation.df_sam_order_tmp$gene_id <-  factor(Proliferation.df_sam_order_tmp$gene_id, levels=Proliferation.gene_order_tmp)

Proliferation.p = Proliferation.df_sam_order_tmp %>%
  ggplot(aes(y=gene_id,x= Sample_id,fill=alter,height=hgt))+
  geom_tile(data=Proliferation.df_sam_order_tmp%>%filter(type=="cnv"))+
  geom_tile(data=Proliferation.df_sam_order_tmp%>%filter(type=="mutation"))+
  geom_tile(data=Proliferation.df_sam_order_tmp%>%filter(type=="wgbs_promoter"))+
  xlab("")+ylab("")+
  scale_fill_manual(values=c("0"="#e9ecef","NA"="none","Missense_Mutation"="#4B7AB1", 
                             "Nonsense_Mutation"="#ED9147",
                             "Frame_Shift_Ins"="#96B458",
                             "Frame_Shift_Del"="#824CA7",
                             "In_Frame_Ins"="#00A087FF",
                             "In_Frame_Del"="#B03060",
                             "Splice_Site"="#006400",
                             "Nonstop_Mutation"="red",
                             "Multi_Hit"="#7E6148FF",
                             "1"="#F17C67","-1"="#4cb4e7",
                             "2"="#FF0080","-2"="#0080FF"))+
  theme(
    axis.text.x=element_blank(),
    axis.text.y = element_text(color="black",size=18),
    axis.ticks = element_blank(),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    legend.position = "none",
    legend.title=element_blank())+
  theme(plot.title = element_text(hjust = 0.5,size=20))

Proliferation.df_sam_order_onco <- Proliferation.df_sam_order_tmp[,c(1,2,5)]
Proliferation.df_sam_order_count = Proliferation.df_sam_order_onco %>% 
  group_by(gene_id, alter) %>% 
  summarise(n= n())
Proliferation.df_sam_order_count <- Proliferation.df_sam_order_count %>% filter(alter != 0)
Proliferation.gene_order_all <- gsub(".pathway", "", Proliferation.gene_order_all)
Proliferation.df_sam_order_count$gene_id <- factor(Proliferation.df_sam_order_count$gene_id, levels=(Proliferation.gene_order_all))
Proliferation.df_sam_order_count$alter <- factor(Proliferation.df_sam_order_count$alter, levels = rev(c(1, -1, 2, -2, "Missense_Mutation",
                                                                                                        "Nonsense_Mutation", "Frame_Shift_Ins",
                                                                                                        "Frame_Shift_Del", "In_Frame_Ins",
                                                                                                        "In_Frame_Del", "Splice_Site",
                                                                                                        "Nonstop_Mutation", "Multi_Hit", 3, -3)))
Proliferation.df_sam_order_count_cnv <- Proliferation.df_sam_order_count %>% filter(alter == -1 | alter == 1)
Proliferation.df_sam_order_count_mut <- Proliferation.df_sam_order_count %>% filter(alter != -1 & alter != 1 & alter != -2 & alter != 2 & alter != -3 & alter != 3)
Proliferation.df_sam_order_count_rna <- Proliferation.df_sam_order_count %>% filter(alter == -3 | alter == 3)
Proliferation.df_sam_order_count_methyl <- Proliferation.df_sam_order_count %>% filter(alter == -2 | alter == 2)

Proliferation_p_count_cnv = ggplot(Proliferation.df_sam_order_count_cnv, aes(x=gene_id, y=n, fill=alter))+
  geom_bar(stat="identity")+
  scale_fill_manual(values=c("Missense_Mutation"="#4B7AB1", "Nonsense_Mutation"="#ED9147", "Frame_Shift_Ins"="#96B458",
                             "Frame_Shift_Del"="#824CA7", "In_Frame_Ins"="#00A087FF", "In_Frame_Del"="#B03060",
                             "Splice_Site"="#006400", "Nonstop_Mutation"="red", "Multi_Hit"="#7E6148FF",
                             "-3"="#2d00f7", "3"="#f20089", "1"="#F17C67", "-1"="#4cb4e7", "2"="#FF0080", "-2"="#0080FF"))+
  coord_flip()+
  ylim(c(0,60))+
  xlab("")+
  ylab("")+
  theme(legend.position = "none",
        axis.line = element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y=element_blank(),
        panel.background = element_blank(),
        axis.ticks.y=element_blank(),
        axis.ticks.x=element_blank())

Proliferation_p_count_mut = ggplot(Proliferation.df_sam_order_count_mut, aes(x=gene_id, y=n, fill=alter))+
  geom_bar(stat="identity")+
  scale_fill_manual(values=c("Missense_Mutation"="#4B7AB1", "Nonsense_Mutation"="#ED9147", "Frame_Shift_Ins"="#96B458",
                             "Frame_Shift_Del"="#824CA7", "In_Frame_Ins"="#00A087FF", "In_Frame_Del"="#B03060",
                             "Splice_Site"="#006400", "Nonstop_Mutation"="red", "Multi_Hit"="#7E6148FF",
                             "-3"="#2d00f7", "3"="#f20089", "1"="#F17C67", "-1"="#4cb4e7", "2"="#FF0080", "-2"="#0080FF"))+
  coord_flip()+
  ylim(c(0,150))+
  xlab("")+
  ylab("")+
  theme(legend.position = "none",
        axis.line = element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y=element_blank(),
        panel.background = element_blank(),
        axis.ticks.y=element_blank(),
        axis.ticks.x=element_blank())

Proliferation_p_count_rna = ggplot(Proliferation.df_sam_order_count_rna, aes(x=gene_id, y=n, fill=alter))+
  geom_bar(stat="identity")+
  scale_fill_manual(values=c("Missense_Mutation"="#4B7AB1", "Nonsense_Mutation"="#ED9147", "Frame_Shift_Ins"="#96B458",
                             "Frame_Shift_Del"="#824CA7", "In_Frame_Ins"="#00A087FF", "In_Frame_Del"="#B03060",
                             "Splice_Site"="#006400", "Nonstop_Mutation"="red", "Multi_Hit"="#7E6148FF",
                             "-3"="#2d00f7", "3"="#f20089", "1"="#F17C67", "-1"="#4cb4e7", "2"="#FF0080", "-2"="#0080FF"))+
  coord_flip()+
  ylim(c(0,100))+
  xlab("")+
  ylab("")+
  theme(legend.position = "none",
        axis.line = element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y=element_blank(),
        panel.background = element_blank(),
        axis.ticks.y=element_blank(),
        axis.ticks.x=element_blank())

Proliferation_p_count_methyl = ggplot(Proliferation.df_sam_order_count_methyl, aes(x=gene_id, y=n, fill=alter))+
  geom_bar(stat="identity")+
  scale_fill_manual(values=c("Missense_Mutation"="#4B7AB1", "Nonsense_Mutation"="#ED9147", "Frame_Shift_Ins"="#96B458",
                             "Frame_Shift_Del"="#824CA7", "In_Frame_Ins"="#00A087FF", "In_Frame_Del"="#B03060",
                             "Splice_Site"="#006400", "Nonstop_Mutation"="red", "Multi_Hit"="#7E6148FF",
                             "-3"="#2d00f7", "3"="#f20089", "1"="#F17C67", "-1"="#4cb4e7", "2"="#FF0080", "-2"="#0080FF"))+
  coord_flip()+
  ylim(c(0,30))+
  xlab("")+
  ylab("")+
  theme(legend.position = "none",
        axis.line = element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y=element_blank(),
        panel.background = element_blank(),
        axis.ticks.y=element_blank(),
        axis.ticks.x=element_blank())

Chromatin.df_sam_order_tmp <- Chromatin.df_sam_order
Chromatin.gene_order_tmp <- Chromatin.gene_order
for (i in 1:length(unique(Chromatin.df_sam_order$gene_id))) {
  gene <- unique(Chromatin.df_sam_order$gene_id)[i]
  tmp <- Chromatin.df_sam_order[,-6]
  tmp <- spread(tmp, type, alter)
  gene_frq <- tmp %>% filter(gene_id == gene)
  gene_frq_ratio <- (1-(dim(gene_frq[which(is.na(gene_frq$mutation) & gene_frq$cnv == 0 & is.na(gene_frq$wgbs_promoter)),])[1]/155))*100
  paste0(gene,"(",round(gene_frq_ratio, 0),"%)",sep="")
  Chromatin.df_sam_order_tmp$gene_id <- gsub(gene, paste0(gene,"(",round(gene_frq_ratio, 0),"%)",sep=""), Chromatin.df_sam_order_tmp$gene_id)
  Chromatin.gene_order_tmp <- gsub(gene, paste0(gene,"(",round(gene_frq_ratio, 0),"%)",sep=""), Chromatin.gene_order_tmp)
}
Chromatin.df_sam_order_tmp$gene_id <-  factor(Chromatin.df_sam_order_tmp$gene_id, levels=Chromatin.gene_order_tmp)

Chromatin.p = Chromatin.df_sam_order_tmp %>%
  ggplot(aes(y=gene_id,x= Sample_id,fill=alter,height=hgt))+
  geom_tile(data=Chromatin.df_sam_order_tmp%>%filter(type=="cnv"))+
  geom_tile(data=Chromatin.df_sam_order_tmp%>%filter(type=="mutation"))+
  geom_tile(data=Chromatin.df_sam_order_tmp%>%filter(type=="wgbs_promoter"))+
  xlab("")+ylab("")+
  scale_fill_manual(values=c("0"="#e9ecef","NA"="none","Missense_Mutation"="#4B7AB1", 
                             "Nonsense_Mutation"="#ED9147",
                             "Frame_Shift_Ins"="#96B458",
                             "Frame_Shift_Del"="#824CA7",
                             "In_Frame_Ins"="#00A087FF",
                             "In_Frame_Del"="#B03060",
                             "Splice_Site"="#006400",
                             "Nonstop_Mutation"="red",
                             "Multi_Hit"="#7E6148FF",
                             "1"="#F17C67","-1"="#4cb4e7",
                             "2"="#FF0080","-2"="#0080FF"))+
  theme(
    axis.text.x=element_blank(),
    axis.text.y = element_text(color="black",size=18),
    axis.ticks = element_blank(),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    legend.position = "none",
    legend.title=element_blank())+
  theme(plot.title = element_text(hjust = 0.5,size=20))

Chromatin.df_sam_order_onco <- Chromatin.df_sam_order[,c(1,2,5)]
Chromatin.df_sam_order_count = Chromatin.df_sam_order_onco %>% 
  group_by(gene_id, alter) %>% 
  summarise(n= n())
Chromatin.df_sam_order_count <- Chromatin.df_sam_order_count %>% filter(alter != 0)
Chromatin.gene_order_all <- gsub(".pathway", "", Chromatin.gene_order_all)
Chromatin.df_sam_order_count$gene_id <- factor(Chromatin.df_sam_order_count$gene_id, levels=(Chromatin.gene_order_all))
Chromatin.df_sam_order_count$alter <- factor(Chromatin.df_sam_order_count$alter, levels = rev(c(1, -1, 2, -2, "Missense_Mutation",
                                                                                                "Nonsense_Mutation", "Frame_Shift_Ins",
                                                                                                "Frame_Shift_Del", "In_Frame_Ins",
                                                                                                "In_Frame_Del", "Splice_Site",
                                                                                                "Nonstop_Mutation", "Multi_Hit", 3, -3)))
Chromatin.df_sam_order_count_cnv <- Chromatin.df_sam_order_count %>% filter(alter == -1 | alter == 1)
Chromatin.df_sam_order_count_mut <- Chromatin.df_sam_order_count %>% filter(alter != -1 & alter != 1 & alter != -2 & alter != 2 & alter != -3 & alter != 3)
Chromatin.df_sam_order_count_rna <- Chromatin.df_sam_order_count %>% filter(alter == -3 | alter == 3)
Chromatin.df_sam_order_count_methyl <- Chromatin.df_sam_order_count %>% filter(alter == -2 | alter == 2)

Chromatin_p_count_cnv = ggplot(Chromatin.df_sam_order_count_cnv, aes(x=gene_id, y=n, fill=alter))+
  geom_bar(stat="identity")+
  scale_fill_manual(values=c("Missense_Mutation"="#4B7AB1", "Nonsense_Mutation"="#ED9147", "Frame_Shift_Ins"="#96B458",
                             "Frame_Shift_Del"="#824CA7", "In_Frame_Ins"="#00A087FF", "In_Frame_Del"="#B03060",
                             "Splice_Site"="#006400", "Nonstop_Mutation"="red", "Multi_Hit"="#7E6148FF",
                             "-3"="#2d00f7", "3"="#f20089", "1"="#F17C67", "-1"="#4cb4e7", "2"="#FF0080", "-2"="#0080FF"))+
  coord_flip()+
  ylim(c(0,60))+
  xlab("")+
  ylab("")+
  theme(legend.position = "none",
        axis.line = element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y=element_blank(),
        panel.background = element_blank(),
        axis.ticks.y=element_blank(),
        axis.ticks.x=element_blank())

Chromatin_p_count_mut = ggplot(Chromatin.df_sam_order_count_mut, aes(x=gene_id, y=n, fill=alter))+
  geom_bar(stat="identity")+
  scale_fill_manual(values=c("Missense_Mutation"="#4B7AB1", "Nonsense_Mutation"="#ED9147", "Frame_Shift_Ins"="#96B458",
                             "Frame_Shift_Del"="#824CA7", "In_Frame_Ins"="#00A087FF", "In_Frame_Del"="#B03060",
                             "Splice_Site"="#006400", "Nonstop_Mutation"="red", "Multi_Hit"="#7E6148FF",
                             "-3"="#2d00f7", "3"="#f20089", "1"="#F17C67", "-1"="#4cb4e7", "2"="#FF0080", "-2"="#0080FF"))+
  coord_flip()+
  ylim(c(0,150))+
  xlab("")+
  ylab("")+
  theme(legend.position = "none",
        axis.line = element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y=element_blank(),
        panel.background = element_blank(),
        axis.ticks.y=element_blank(),
        axis.ticks.x=element_blank())

Chromatin_p_count_rna = ggplot(Chromatin.df_sam_order_count_rna, aes(x=gene_id, y=n, fill=alter))+
  geom_bar(stat="identity")+
  scale_fill_manual(values=c("Missense_Mutation"="#4B7AB1", "Nonsense_Mutation"="#ED9147", "Frame_Shift_Ins"="#96B458",
                             "Frame_Shift_Del"="#824CA7", "In_Frame_Ins"="#00A087FF", "In_Frame_Del"="#B03060",
                             "Splice_Site"="#006400", "Nonstop_Mutation"="red", "Multi_Hit"="#7E6148FF",
                             "-3"="#2d00f7", "3"="#f20089", "1"="#F17C67", "-1"="#4cb4e7", "2"="#FF0080", "-2"="#0080FF"))+
  coord_flip()+
  ylim(c(0,100))+
  xlab("")+
  ylab("")+
  theme(legend.position = "none",
        axis.line = element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y=element_blank(),
        panel.background = element_blank(),
        axis.ticks.y=element_blank(),
        axis.ticks.x=element_blank())

Chromatin_p_count_methyl = ggplot(Chromatin.df_sam_order_count_methyl, aes(x=gene_id, y=n, fill=alter))+
  geom_bar(stat="identity")+
  scale_fill_manual(values=c("Missense_Mutation"="#4B7AB1", "Nonsense_Mutation"="#ED9147", "Frame_Shift_Ins"="#96B458",
                             "Frame_Shift_Del"="#824CA7", "In_Frame_Ins"="#00A087FF", "In_Frame_Del"="#B03060",
                             "Splice_Site"="#006400", "Nonstop_Mutation"="red", "Multi_Hit"="#7E6148FF",
                             "-3"="#2d00f7", "3"="#f20089", "1"="#F17C67", "-1"="#4cb4e7", "2"="#FF0080", "-2"="#0080FF"))+
  coord_flip()+
  ylim(c(0,30))+
  xlab("")+
  ylab("")+
  theme(legend.position = "none",
        axis.line = element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y=element_blank(),
        panel.background = element_blank(),
        axis.ticks.y=element_blank(),
        axis.ticks.x=element_blank())

p_NRF2.p = ggplot_gtable(ggplot_build(NRF2.p))
p_NRF2_p_count_cnv = ggplot_gtable(ggplot_build(NRF2_p_count_cnv))
p_NRF2_p_count_methyl = ggplot_gtable(ggplot_build(NRF2_p_count_methyl))
p_NRF2_p_count_mut = ggplot_gtable(ggplot_build(NRF2_p_count_mut))
p_NRF2_p_count_rna = ggplot_gtable(ggplot_build(NRF2_p_count_rna))

p_Differentiation.p = ggplot_gtable(ggplot_build(Differentiation.p))
p_Differentiation_p_count_cnv = ggplot_gtable(ggplot_build(Differentiation_p_count_cnv))
p_Differentiation_p_count_methyl = ggplot_gtable(ggplot_build(Differentiation_p_count_methyl))
p_Differentiation_p_count_mut = ggplot_gtable(ggplot_build(Differentiation_p_count_mut))
p_Differentiation_p_count_rna = ggplot_gtable(ggplot_build(Differentiation_p_count_rna))

p_Cellcycle.p = ggplot_gtable(ggplot_build(Cellcycle.p))
p_Cellcycle_p_count_cnv = ggplot_gtable(ggplot_build(Cellcycle_p_count_cnv))
p_Cellcycle_p_count_methyl = ggplot_gtable(ggplot_build(Cellcycle_p_count_methyl))
p_Cellcycle_p_count_mut = ggplot_gtable(ggplot_build(Cellcycle_p_count_mut))
p_Cellcycle_p_count_rna = ggplot_gtable(ggplot_build(Cellcycle_p_count_rna))

p_RTK.p = ggplot_gtable(ggplot_build(RTK.p))
p_RTK_p_count_cnv = ggplot_gtable(ggplot_build(RTK_p_count_cnv))
p_RTK_p_count_methyl = ggplot_gtable(ggplot_build(RTK_p_count_methyl))
p_RTK_p_count_mut = ggplot_gtable(ggplot_build(RTK_p_count_mut))
p_RTK_p_count_rna = ggplot_gtable(ggplot_build(RTK_p_count_rna))

p_Proliferation.p = ggplot_gtable(ggplot_build(Proliferation.p))
p_Proliferation_p_count_cnv = ggplot_gtable(ggplot_build(Proliferation_p_count_cnv))
p_Proliferation_p_count_methyl = ggplot_gtable(ggplot_build(Proliferation_p_count_methyl))
p_Proliferation_p_count_mut = ggplot_gtable(ggplot_build(Proliferation_p_count_mut))
p_Proliferation_p_count_rna = ggplot_gtable(ggplot_build(Proliferation_p_count_rna))

p_Chromatin.p = ggplot_gtable(ggplot_build(Chromatin.p))
p_Chromatin_p_count_cnv = ggplot_gtable(ggplot_build(Chromatin_p_count_cnv))
p_Chromatin_p_count_methyl = ggplot_gtable(ggplot_build(Chromatin_p_count_methyl))
p_Chromatin_p_count_mut = ggplot_gtable(ggplot_build(Chromatin_p_count_mut))
p_Chromatin_p_count_rna = ggplot_gtable(ggplot_build(Chromatin_p_count_rna))

p_Differentiation.p$widths = p_NRF2.p$widths
p_Cellcycle.p$widths = p_NRF2.p$widths
p_RTK.p$widths = p_NRF2.p$widths
p_Proliferation.p$widths = p_NRF2.p$widths
p_Chromatin.p$widths = p_NRF2.p$widths

p_Differentiation_p_count_cnv$widths = p_NRF2_p_count_cnv$widths
p_Cellcycle_p_count_cnv$widths = p_NRF2_p_count_cnv$widths
p_RTK_p_count_cnv$widths = p_NRF2_p_count_cnv$widths
p_Proliferation_p_count_cnv$widths = p_NRF2_p_count_cnv$widths
p_Chromatin_p_count_cnv$widths = p_NRF2_p_count_cnv$widths

p_Differentiation_p_count_methyl$widths = p_NRF2_p_count_methyl$widths
p_Cellcycle_p_count_methyl$widths = p_NRF2_p_count_methyl$widths
p_RTK_p_count_methyl$widths = p_NRF2_p_count_methyl$widths
p_Proliferation_p_count_methyl$widths = p_NRF2_p_count_methyl$widths
p_Chromatin_p_count_methyl$widths = p_NRF2_p_count_methyl$widths

p_Differentiation_p_count_mut$widths = p_NRF2_p_count_mut$widths
p_Cellcycle_p_count_mut$widths = p_NRF2_p_count_mut$widths
p_RTK_p_count_mut$widths = p_NRF2_p_count_mut$widths
p_Proliferation_p_count_mut$widths = p_NRF2_p_count_mut$widths
p_Chromatin_p_count_mut$widths = p_NRF2_p_count_mut$widths

p_Differentiation_p_count_rna$widths = p_NRF2_p_count_rna$widths
p_Cellcycle_p_count_rna$widths = p_NRF2_p_count_rna$widths
p_RTK_p_count_rna$widths = p_NRF2_p_count_rna$widths
p_Proliferation_p_count_rna$widths = p_NRF2_p_count_rna$widths
p_Chromatin_p_count_rna$widths = p_NRF2_p_count_rna$widths

p_NRF2.p$heights[c(1,8,9,10,12)] = unit(rep(0,5),"cm")
p_NRF2_p_count_rna$heights[c(1,8,9,10,12)] = unit(rep(0,5),"cm")
p_NRF2_p_count_mut$heights[c(1,8,9,10,12)] = unit(rep(0,5),"cm")
p_NRF2_p_count_methyl$heights[c(1,8,9,10,12)] = unit(rep(0,5),"cm")
p_NRF2_p_count_cnv$heights[c(1,8,9,10,12)] = unit(rep(0,5),"cm")

p_Differentiation.p$heights[c(1,8,9,10,12)] = unit(rep(0,5),"cm")
p_Differentiation_p_count_rna$heights[c(1,8,9,10,12)] = unit(rep(0,5),"cm")
p_Differentiation_p_count_mut$heights[c(1,8,9,10,12)] = unit(rep(0,5),"cm")
p_Differentiation_p_count_methyl$heights[c(1,8,9,10,12)] = unit(rep(0,5),"cm")
p_Differentiation_p_count_cnv$heights[c(1,8,9,10,12)] = unit(rep(0,5),"cm")

p_Cellcycle.p$heights[c(1,8,9,10,12)] = unit(rep(0,5),"cm")
p_Cellcycle_p_count_rna$heights[c(1,8,9,10,12)] = unit(rep(0,5),"cm")
p_Cellcycle_p_count_mut$heights[c(1,8,9,10,12)] = unit(rep(0,5),"cm")
p_Cellcycle_p_count_methyl$heights[c(1,8,9,10,12)] = unit(rep(0,5),"cm")
p_Cellcycle_p_count_cnv$heights[c(1,8,9,10,12)] = unit(rep(0,5),"cm")

p_RTK.p$heights[c(1,8,9,10,12)] = unit(rep(0,5),"cm")
p_RTK_p_count_rna$heights[c(1,8,9,10,12)] = unit(rep(0,5),"cm")
p_RTK_p_count_mut$heights[c(1,8,9,10,12)] = unit(rep(0,5),"cm")
p_RTK_p_count_methyl$heights[c(1,8,9,10,12)] = unit(rep(0,5),"cm")
p_RTK_p_count_cnv$heights[c(1,8,9,10,12)] = unit(rep(0,5),"cm")

p_Proliferation.p$heights[c(1,8,9,10,12)] = unit(rep(0,5),"cm")
p_Proliferation_p_count_rna$heights[c(1,8,9,10,12)] = unit(rep(0,5),"cm")
p_Proliferation_p_count_mut$heights[c(1,8,9,10,12)] = unit(rep(0,5),"cm")
p_Proliferation_p_count_methyl$heights[c(1,8,9,10,12)] = unit(rep(0,5),"cm")
p_Proliferation_p_count_cnv$heights[c(1,8,9,10,12)] = unit(rep(0,5),"cm")

p_Chromatin.p$heights[c(1,8,9,10,12)] = unit(rep(0,5),"cm")
p_Chromatin_p_count_rna$heights[c(1,8,9,10,12)] = unit(rep(0,5),"cm")
p_Chromatin_p_count_mut$heights[c(1,8,9,10,12)] = unit(rep(0,5),"cm")
p_Chromatin_p_count_methyl$heights[c(1,8,9,10,12)] = unit(rep(0,5),"cm")
p_Chromatin_p_count_cnv$heights[c(1,8,9,10,12)] = unit(rep(0,5),"cm")


m2 = grid.arrange(p_NRF2.p, p_Differentiation.p, p_Cellcycle.p, p_RTK.p, p_Proliferation.p, p_Chromatin.p,
                  p_NRF2_p_count_cnv, p_Differentiation_p_count_cnv, p_Cellcycle_p_count_cnv, p_RTK_p_count_cnv, p_Proliferation_p_count_cnv, p_Chromatin_p_count_cnv,
                  p_NRF2_p_count_mut, p_Differentiation_p_count_mut, p_Cellcycle_p_count_mut, p_RTK_p_count_mut, p_Proliferation_p_count_mut, p_Chromatin_p_count_mut,
                  p_NRF2_p_count_methyl, p_Differentiation_p_count_methyl, p_Cellcycle_p_count_methyl, p_RTK_p_count_methyl, p_Proliferation_p_count_methyl, p_Chromatin_p_count_methyl,
                  layout_matrix=rbind(c(rep(3,9),9,15,21),c(rep(3,9),9,15,21),c(rep(3,9),9,15,21),c(rep(3,9),9,15,21),c(rep(3,9),9,15,21),c(rep(3,9),9,15,21),
                                      c(rep(4,9),10,16,22),c(rep(4,9),10,16,22),c(rep(4,9),10,16,22),c(rep(4,9),10,16,22),c(rep(4,9),10,16,22),c(rep(4,9),10,16,22),c(rep(4,9),10,16,22),c(rep(4,9),10,16,22),
                                      c(rep(2,9),8,14,20),c(rep(2,9),8,14,20),c(rep(2,9),8,14,20),c(rep(2,9),8,14,20),
                                      c(rep(6,9),12,18,24),c(rep(6,9),12,18,24),c(rep(6,9),12,18,24),c(rep(6,9),12,18,24),c(rep(6,9),12,18,24),c(rep(6,9),12,18,24),
                                      c(rep(1,9),7,13,19),c(rep(1,9),7,13,19),c(rep(1,9),7,13,19),
                                      c(rep(5,9),11,17,23),c(rep(5,9),11,17,23),c(rep(5,9),11,17,23),c(rep(5,9),11,17,23),c(rep(5,9),11,17,23),c(rep(5,9),11,17,23),c(rep(5,9),11,17,23)
                  ))

remain <- Cellcycle.df_sam_order
remain$alter <- gsub(-1,"Deletion",remain$alter)
remain$alter <- gsub(1,"Amplification",remain$alter)
remain$alter <- gsub(-2,"zHypo-methylation",remain$alter)
remain$alter <- gsub(2,"zHyper-methylation",remain$alter)
remain$alter <- gsub(-3,"zLow-expression",remain$alter)
remain$alter <- gsub(3,"zOver-expression",remain$alter)
remain$alter <- gsub(0,NA,remain$alter)

p_legend = remain %>%
  ggplot(aes(y=gene_id,x= Sample_id,fill=alter,height=hgt))+
  geom_tile(data=remain%>%filter(type=="cnv"))+
  geom_tile(data=remain%>%filter(type=="mutation"))+
  geom_tile(data=remain%>%filter(type=="wgbs_promoter"))+
  geom_tile(data=remain%>%filter(type=="rna"),colour=Cellcycle.cc$alter,width=0.85,size=0.4)+
  xlab("")+ylab("")+
  scale_fill_manual(values=c("NA"="none","Missense_Mutation"="#4B7AB1", 
                             "Nonsense_Mutation"="#ED9147",
                             "Frame_Shift_Ins"="#96B458",
                             "Frame_Shift_Del"="#824CA7",
                             "In_Frame_Ins"="#00A087FF",
                             "In_Frame_Del"="#B03060",
                             "Splice_Site"="#006400",
                             "Nonstop_Mutation"="red",
                             "Multi_Hit"="#7E6148FF",
                             "Amplification"="#F17C67","Deletion"="#4cb4e7",
                             "zHyper-methylation"="#FF0080","zHypo-methylation"="#0080FF",
                             "zLow-expression"="#2d00f7","zOver-expression"="#f20089"))+
  theme(title= element_text(size=13, color="black"))+
  theme(legend.text= element_text(size=14, color="black"))

p_legend = ggplot_gtable(ggplot_build(p_legend))
p_legend = p_legend$grobs[[15]]
legend.m = grid.arrange(p_legend,nrow=1)