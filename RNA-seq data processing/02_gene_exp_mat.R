args <- commandArgs(TRUE)
if(length(args)<1){
	cat("Usage: Rscript gene_exp_mat.R [wgcid.genes.results...]\n")
	q()
}

library(edgeR)

df = read.csv(args[1], sep='\t', header=T)
df[,1] <- gsub("\\.\\d+\\_","\\_",df[,1])
count_mat = df[,c(1,5)]
TPM_mat = df[,c(1,6)]
FPKM_mat = df[,c(1,7)]

for (sample in args[-1]){
	df  = read.csv(sample, '\t', header=T)
	count_mat = cbind(count_mat, df[,5])
	TPM_mat = cbind(TPM_mat, df[,6])
	FPKM_mat = cbind(FPKM_mat, df[,7])
}

col = c('EnsemblGene_GeneSymbol',sapply(args, function(x){gsub("\\.genes\\.results","",basename(x))}))
colnames(count_mat) = col
colnames(TPM_mat) = col
colnames(FPKM_mat) = col
count_mat[,-1]=round(count_mat[,-1])

write.table(count_mat, "gene_count_mat.xls",sep="\t", quote=F, row.names=F)
write.table(TPM_mat, "gene_TPM_mat.xls",sep="\t", quote=F, row.names=F)
write.table(FPKM_mat, "gene_FPKM_mat.xls",sep="\t", quote=F, row.names=F)


# normalize count with TMM
sample.number = length(args)
if(sample.number > 1){
	artificial.condition <- data.frame(group=c(0, rep(1, (sample.number-1))))
	rownames(artificial.condition) <- col[-1]
	dge <- DGEList(counts=count_mat[,-1], group=artificial.condition$group)
	dge <- calcNormFactors(dge, method='TMM')
	nct <- cpm(dge, normalized.lib.sizes=TRUE)
	nct.out <- cbind(count_mat[,1], nct)
	colnames(nct.out)[1] <- 'EnsemblGene_GeneSymbol'
	write.table(nct.out, "gene_CPM_TMM.xls",sep="\t", quote=F, row.names=F)
}


