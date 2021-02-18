# author: Lixng, derive from Jiawei and zehua's code
# last modified date: 20180109
args <- commandArgs(TRUE)
if(length(args)<7){
  cat("Usage: Rscript DEG_edgeR.R mat_count.txt ctrl_col_idx case_col_idx ctrl_label case_label FC pvalue BCV\n
       The last argument BCV is optional. If your comparison design has replicates, ignore this argument. Otherwise, please provide the empirical BCV to estimate the replicate variance. Typical values for the common BCV (squareroot-dispersion) for datasets arising from well-controlled experiments are 0.4 for human
data, 0.1 for data on genetically identical model organisms or 0.01 for technical replicates.\n")
  q()
}

library(edgeR)
library(gridExtra)
library(grid)
library(scales)
library(ggplot2)


case <- args[5]
ctrl <- args[4]

fc.cut <- as.numeric(args[6])
padj.cut <- as.numeric(args[7])

count.data <- read.table(args[1], header=T, row.names=1,check.names=F)
ctrl.idx <- unlist(strsplit(args[2], ","))
case.idx <- unlist(strsplit(args[3], ","))
#count.data <- count.data[, c(ctrl.idx, case.idx)]
#nsample <- ncol(count.data)
nsample=length(c(ctrl.idx,case.idx))
condition <- data.frame(group=c(rep(args[4], length(ctrl.idx)), rep(args[5], length(case.idx))))
condition$group <- relevel(condition$group, ctrl)
#rownames(condition) <- colnames(count.data)
rownames(condition) <- c(ctrl.idx,case.idx)

#condition$group2 <- c(1,1,1,1,1,2)
#dge <- DGEList(counts=count.data, group=condition$group)
dge <- DGEList(counts=count.data)

#keep <- rowSums(cpm(dge)>1) >= 1
#dge <- dge[keep, keep.lib.sizes=FALSE]


design <- model.matrix(~group, data=condition)
colnames(design) <- levels(condition$group)
# estimate library size
dge <- calcNormFactors(dge, method='TMM')
dge = dge[,c(ctrl.idx,case.idx)]
dge$samples$group=condition$group

keep <- rowSums(cpm(dge)>1) >= 1
dge <- dge[keep, keep.lib.sizes=TRUE]

nct <- cpm(dge, normalized.lib.sizes=TRUE)
nct <- data.frame(round(nct, 3))

if(nsample>2){
  # esimate dispersion
  dge <- estimateGLMCommonDisp(dge, design)
  dge <- estimateGLMTrendedDisp(dge, design)
  dge <- estimateGLMTagwiseDisp(dge, design)
  # test
  fit <- glmFit(dge, design)
  lrt <- glmLRT(fit)
}else{
  print("No replicate, using empirical dispersion! BCV = 0.1")
  if (length(args) != 8){
    bcv <- 0.1
  }else{
    bcv <- as.numeric(args[8])
  }
  fit <- glmFit(dge, design, dispersion=bcv^2)
  lrt <- glmLRT(fit)
}
# summary results
result.table <- lrt$table
#pvals[rowSums(count.data) == 0] <- NA
result.table$padj <- p.adjust(result.table$PValue, method="BH")
result.table$padj [is.na(result.table$padj )] <- 1
result.table <- merge(nct, result.table, by='row.names')
edrout <- data.frame(EnsemblGene_GeneSymbol=result.table$Row.names,
                     mean_logCPM=result.table$logCPM,
                     log2FC=result.table$logFC,
                     pvalue=result.table$PValue,
                     padj=result.table$padj)

edrout <- cbind(edrout, result.table[,c(2:(1+nsample))])
edrout <- edrout[order(edrout$padj), ]
write.table(edrout, paste0('DEG_',case,'_vs_',ctrl,'_all.xls'), sep="\t", row.names=F, quote=F)

result.table$sig <- (result.table$logFC > log2(fc.cut) | result.table$logFC < log2(1/fc.cut)) & result.table$padj<padj.cut

edrout.sig <- subset(edrout, padj<padj.cut & (log2FC>log2(fc.cut) | log2FC<log2(1/fc.cut)))
write.table(edrout.sig, paste0('DEG_', case, '_vs_', ctrl, '_sig.xls'), sep="\t", row.names=F, quote=F)
edrout.sig.up <- data.frame(edrout.sig[edrout.sig$log2FC>0,]$EnsemblGene_GeneSymbol)
edrout.sig.down <- data.frame(edrout.sig[edrout.sig$log2FC<0,]$EnsemblGene_GeneSymbol)
write.table(edrout.sig.up, paste0('DEG_', case, '_vs_', ctrl, '_sig.up.txt'), sep="\t",col.names=F, row.names=F, quote=F)
write.table(edrout.sig.down, paste0('DEG_', case, '_vs_', ctrl, '_sig.down.txt'), sep="\t", col.names=F, row.names=F, quote=F)
df <- data.frame(up_regulated=nrow(edrout.sig.up), down_regulated=nrow(edrout.sig.down))
file <- paste0('DEG_', case, '_vs_', ctrl, '_diff_stats.txt')
cat(paste0("# Criteria for identifying DEGs: FC > ",fc.cut,' and FDR < ', padj.cut, '\n'),  file = file)
write.table(df,file,append=T, quote=F, row.names = F)

## MAplot
p <- ggplot(result.table, aes(result.table$logCPM, result.table$logFC)) +
     geom_point(data=result.table[result.table$sig=="FALSE",], aes(logCPM, logFC), colour="#BFBFBF", size=2) +
     geom_point(data=result.table[result.table$sig=="TRUE",], aes(logCPM, logFC), colour="#E25668", size=2) +
     geom_hline(aes(yintercept=0), colour="#4B7AB1",linetype="dashed", size=1) +
     geom_hline(aes(yintercept=log2(1/fc.cut)),colour="#96B458", linetype="dashed", size=1) +
     geom_hline(aes(yintercept=log2(fc.cut)), colour="#96B458", linetype="dashed", size=1) +
     xlab("Average logCPM") +
     ylab("log2FC") +
     ggtitle(paste0("MA Plot for ", case, "_vs_", ctrl))+
     theme(legend.position = "none",
           axis.text.x = element_text(size=14),
           axis.text.y = element_text(size=14),
           title = element_text(face="bold", size=14))
ggsave(paste0('DEG_', case, '_vs_', ctrl, '_MAplot.pdf'), p, width=8, height=6)
ggsave(paste0('DEG_', case, '_vs_', ctrl, '_MAplot.pdf'), p, width=8, height=6)

# Volcano plot
pdf(paste0('DEG_', case, '_vs_', ctrl, '_Volcano.pdf'), width=10, height=8, onefile=F)
result.table.up <- result.table[result.table$sig=="TRUE"&result.table$logFC>0,]
result.table.up10 <- head(result.table.up[order(result.table.up$logFC, decreasing=T),], 10)
n.up <- nrow(result.table.up10)
if(n.up>0){
    result.table.up10$num <- seq(1:n.up)
}

result.table.down <- result.table[result.table$sig=="TRUE"&result.table$logFC<0,]
result.table.down10 <- head(result.table.down[order(result.table.down$logFC, decreasing=F),], 10)
n.down <- nrow(result.table.down10)
if(n.down>0){
    result.table.down10$num <- seq(1:n.down)
}

x.lim <- max(c(abs(result.table.down$logFC),result.table.up$logFC))*1.1
p <- ggplot(result.table, aes(result.table$logCPM, result.table$logFC)) +
  geom_point(data=result.table[result.table$sig=="FALSE",], aes(logFC, -log10(padj)), colour="#BFBFBF", size=2) +
  geom_point(data=result.table.up, aes(logFC, -log10(padj)), colour="#E25668", size=4) +
  geom_point(data=result.table.down, aes(logFC, -log10(padj)), colour="#4B7AB1", size=4) +
  ylim(0, NA)+
  xlim(-x.lim,x.lim)+
  geom_hline(aes(yintercept=-log10(padj.cut)), colour="gray65",linetype="dashed", size=1) +
  geom_vline(aes(xintercept=log2(1/fc.cut)),colour="gray65", linetype="dashed", size=1) +
  geom_vline(aes(xintercept=log2(fc.cut)), colour="gray65", linetype="dashed", size=1) +
  xlab("log2FC") +
  ylab("-log10(padj)") +
  ggtitle(paste0("Volcano Plot for ", case, "_vs_", ctrl))+
  theme(legend.position = "none",
        axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        title = element_text(face="bold", size=14))

if(n.up>0){
 p = p+geom_text(data = result.table.up10, aes(logFC, -log10(padj), label = num),colour="black",size=3)
}
if(n.down>0){
 p = p+geom_text(data = result.table.down10, aes(logFC, -log10(padj), label = num),colour="black",size=3)
}

deg_table_up <- data.frame("Up_regulated"= gsub("ENS[a-zA-Z0-9]+_","",result.table.up10$Row.names))
deg_table_down <- data.frame("Down_regulated"=gsub("ENS[a-zA-Z0-9]+_","",result.table.down10$Row.names))

if(n.up>0){
  attach(result.table.up10)
  q.up <- ggplot(result.table.up10, aes(x=rep(1,n.up),y=seq(10,11-n.up)))+
  geom_point(colour="#E25668", size=4)+
  geom_text(data = result.table.up10, aes(x=rep(1,n.up), y=seq(10,11-n.up), label = num, hjust=0.5, vjust=0.5),colour="black",size=3)+
  geom_text(data = result.table.up10, aes(x=rep(1.2,n.up), y=seq(10,11-n.up), hjust=0,label = gsub("ENS[a-zA-Z0-9]+_","",result.table.up10$Row.names)) ,colour="black",size=3)+
  xlim(1,3)+
  ylim(1,10)+
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())
  detach(result.table.up10)
}else{
   q.up<- grid.rect(gp=gpar(col=NA))
}

if(n.down>0){
  attach(result.table.down10)
  q.down <- ggplot(result.table.down10, aes(x=rep(1,n.down),y=seq(10, 11-n.down)))+
  geom_point(colour="#4B7AB1", size=4)+
  geom_text(data = result.table.down10, aes(x=rep(1,n.down), y=seq(10, 11-n.down), label = num, hjust=0.5, vjust=0.5),colour="black",size=3)+
  geom_text(data = result.table.down10, aes(x=rep(1.2,n.down), y=seq(10, 11-n.down), hjust=0,label = gsub("ENS[a-zA-Z0-9]+_","",result.table.down10$Row.names)),colour="black",size=3)+
  xlim(1,3)+
  ylim(1,10)+
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())
  detach(result.table.down10)
}else{
   q.down<- grid.rect(gp=gpar(col=NA))
}
t.up <- textGrob("Top 10 Up-regulated Genes",gp=gpar(fontsize=10),hjust=0.5,vjust=1)
t.down <- textGrob("Top 10 Down-regulated Genes",gp=gpar(fontsize=10),hjust=0.5,vjust=1)
m <- grid.arrange(p,  grid.arrange(t.up, q.up, t.down, q.down, nrow=4, heights=c(0.5,4.5,0.5,4.5)), ncol=2, widths=c(8,2))
grid.draw(m)
dev.off()

png(paste0('DEG_', case, '_vs_', ctrl, '_Volcano.png'), width=10, height=8, units="in",res=300)
grid.draw(m)
dev.off()

