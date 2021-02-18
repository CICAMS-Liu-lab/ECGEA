#/nextcode/sge_software/anaconda2/bin/R −−max−mem− s i z e =20Gb
library(GenomicRanges)
library(ggplot2)
library(SomaticSignatures)
library(SomaticCancerAlterations)
library(BSgenome.Hsapiens.1000genomes.hs37d5)
library(BSgenome)
library(MutationalPatterns)
require('pheatmap')
require(maftools)

###Statistics were performed for each mutation type in each sample ###

#args=commandArgs(TRUE)
#if(length(args)!=2){
# cat("Usage: Rscript signature.primary.R maf output.pdf")
#}
#input=args[1]
#output=args[2]

#laml = read.maf(maf = input)
#laml.titv = titv(maf = laml, plot = FALSE, useSyn = TRUE)
#pdf(output)
#plotTiTv(res = laml.titv,showBarcodes = TRUE,textSize = 8)
#dev.off()


### Build Vranges data frame which be used next step ###
maf=read.table("ESCC_155_pairs.maf",header=T)
maf$ID=paste(maf$Tumor_Sample_Barcode,"_vs_",maf$Matched_Norm_Sample_Barcode,sep="")
df=as.data.frame(maf)

vinfo=makeGRangesFromDataFrame(df,
                         keep.extra.columns=TRUE,
                         ignore.strand=TRUE,
                         seqinfo=NULL,
                         seqnames.field=c("seqnames", "seqname",
                                          "Chromosome", "chrom",
                                          "chr", "chromosome_name",
                                          "seqid"),
                         start.field="Start_Position",
                         end.field=c("End_Position", "stop"),
                         strand.field="strand",
                         starts.in.df.are.0based=FALSE)
sca_data = unname(subset(vinfo, Variant_Type %in% "SNP"))
sca_vr = VRanges(
    seqnames = seqnames(sca_data),
    ranges = ranges(sca_data),
    ref = sca_data$Reference_Allele,
    alt = sca_data$Tumor_Seq_Allele2,
    sampleNames = sca_data$Tumor_Sample_Barcode,
    seqinfo = seqinfo(sca_data),
    study = sca_data$Tumor_Sample_Barcode)
###Extract triples and creat Matrix ###
sca_motifs = mutationContext(sca_vr, BSgenome.Hsapiens.1000genomes.hs37d5)
sca_mm = motifMatrix(sca_motifs, group = "study", normalize = TRUE)
#pdf("sample.pdf")
#plotMutationSpectrum(sca_motifs, "study",colorby = "alteration")
#dev.off()

###Evaluate the number of signature,each of number from sigs 2 to sigs 8  run 30 turn ###

#n_sigs = 2:30
#gof_nmf = assessNumberSignatures(sca_mm, n_sigs, nReplicates = 30)
#pdf("select.sig.pdf")
#plotNumberSignatures(gof_nmf)
#dev.off()

###decompositaion when  signature equal 5 
n_sigs = 10
sigs_nmf = identifySignatures(sca_mm, n_sigs, nmfDecomposition)
sig_matrics=signatures(sigs_nmf)
#nmf_res <- extract_signatures(mut_mat, rank = 5,nrun=30)
#sig_matrics=as.matrix(nmf_res$signatures)
write.table(samples(sigs_nmf),"sample_sigs.xls",row.names = TRUE, col.names = TRUE,sep="\t")
write.table(sig_matrics,"sig_matrics.xls",row.names = TRUE, col.names = TRUE,sep="\t")


pdf("motif.sig.heatmap.pdf")
plotSignatureMap(sigs_nmf) + ggtitle("Somatic Signatures: NMF - Heatmap")
dev.off()

pdf("sigs.pdf")
SomaticSignatures::plotSignatures(sigs_nmf) + ggtitle("Somatic Signatures: NMF - Barchart")
dev.off()


pdf("sample.sig.heatmap.pdf")
plotSampleMap(sigs_nmf)
dev.off()

pdf("sample.sig.ratio.pdf")
plotSamples(sigs_nmf)
dev.off()


###Calculate the cosine similarity of cosmic signatures ###
row.name=row.names(sig_matrics)
change.name=""
for (name in row.name){
        first.third=strsplit(name, ' ')[[1]][2]
        middle=strsplit(name, ' ')[[1]][1]
        first=strsplit(first.third, '[.]')[[1]][1]
        third=strsplit(first.third, '[.]')[[1]][2]
        tmp_middle=strsplit(middle, '')[[1]]
        change.name=paste(change.name,paste(first,"[",tmp_middle[1],">",tmp_middle[2],"]",third,sep=""),sep=" ")
}
new.name=strsplit(change.name, ' ')[[1]][-1]

row.names(sig_matrics)=new.name
ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
library(ref_genome, character.only = TRUE)
cancer_signatures = read.table("./signatures_probabilities.txt", sep = "\t", header = TRUE)
new_order = match(row.names(sig_matrics), cancer_signatures$Somatic.Mutation.Type)
cancer_signatures = cancer_signatures[as.vector(new_order),]
row.names(cancer_signatures) = cancer_signatures$Somatic.Mutation.Type
cancer_signatures = as.matrix(cancer_signatures[,4:33])
matrix=matrix(nrow=ncol(sig_matrics),ncol=ncol(cancer_signatures))
for(i in 1:ncol(sig_matrics)){
        result.list<-list()
        for(j in 1:ncol(cancer_signatures)){
                matrix[i,j]=cos_sim(sig_matrics[,i], cancer_signatures[,j])
                result.list<-c(result.list,cos_sim(sig_matrics[,i], cancer_signatures[,j]))
        }
}
rownames(matrix)=colnames(sig_matrics)
colnames(matrix)=colnames(cancer_signatures)
pdf("cosmic.sig.pdf")
pheatmap::pheatmap(mat = matrix, cluster_rows = FALSE, main = "cosine similarity against validated signatures")
dev.off()

factor=""
row.names(matrix)=colnames(sig_matrics)
colnames(matrix)=colnames(cancer_signatures)
similarity=apply(matrix,1,max)
cosmic=apply(matrix, 1, function(t) colnames(matrix)[which.max(t)])
for(i in similarity){if(i>0.9) factor=c(factor,"PASS") else factor=c(factor,"NOT PASS")}
result=cbind(similarity,cosmic,factor[-1])
result=as.data.frame(result)
colnames(result)[3]="similarity>0.8"
table=read.table("./database/cosmic.note.txt",sep="\t",header=T)
pos=match(result[,2],table[,1])
dat=cbind(result,table[pos,])

write.table(dat[,-c(4,9)],"result.cosmic.xls",row.names = TRUE, col.names = TRUE,sep="\t")
