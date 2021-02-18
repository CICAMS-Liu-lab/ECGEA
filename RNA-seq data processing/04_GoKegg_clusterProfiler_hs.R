######### encrich KEGG & GO 
####clusterprofiler  (entrzid prefered)
#R version > 3.4.0
#library("Matrix")
library("getopt")
command=matrix(c("input_dir","i",1,"character",
                 "output_dir","o",1,"character",
                 "help","h",0,"logical"),byrow=T,ncol=4)
args=getopt(command)
if (!is.null(args$help) || is.null(args$input_dir) || is.null(args$output_dir)) {
  cat(paste(getopt(command, usage = T), "\n"))
  q()
}

inputdir <- args$input_dir
outputdir <- args$output_dir

setwd(inputdir)
######################################################################################
library(clusterProfiler)
library(org.Hs.eg.db)

######################################################################################

######background genes for GO enrichment
entrez_object <- org.Hs.egGO
mapped_genes <- mappedkeys(entrez_object)
universe <- mapped_genes
####transfer geneid to genesymbol
etrz2sym <- org.Hs.egSYMBOL
mappedgenes <- mappedkeys(etrz2sym)
entrz_symb <- as.list(etrz2sym[mappedgenes])
####kegg database#####
keggdb<-read.table(file='/glusterfs/home/wang_beibei0103/database/kegg/kegg_hs.db',sep='\t',header=T)
kegg2gene<-keggdb[,c(1,2)]
kegg2name<-read.table(file='/glusterfs/home/wang_beibei0103/database/kegg/kegg2name.conf',sep='\t',header=T)


for(file in dir(pattern="gene_list.txt")){
  gene= readLines(file)
  sample <- unlist(strsplit(file,"\\."))[1]
  
  #####Go enrichment
  if(!any(gene %in% universe)){
    next
    }else{
      ego <- enrichGO(gene = gene, 
                      OrgDb = org.Hs.eg.db, 
                      universe = universe, 
                      keyType = "ENTREZID", 
                      ont = "ALL",
                      pAdjustMethod = "BH",
                      minGSSize = 1,
                      maxGSSize = 10000000000000,
                      pvalueCutoff  = 1,
                      qvalueCutoff  = 1,
                      readable = TRUE)
      ego <- ego@result
      ego_sigBP <- ego[(ego$ONTOLOGY == 'BP')&(ego$pvalue < 0.01),]
      ego_sigMF <- ego[(ego$ONTOLOGY == 'MF')&(ego$pvalue < 0.01),]
      ego_sigCC <- ego[(ego$ONTOLOGY == 'CC')&(ego$pvalue < 0.01),]
      write.table(ego_sigBP,file=paste(outputdir,"/",sample,"_BP.xls",sep=""),row.names=F,quote=F,sep="\t")
      write.table(ego_sigMF,file=paste(outputdir,"/",sample,"_MF.xls",sep=""),row.names=F,quote=F,sep="\t")
      write.table(ego_sigCC,file=paste(outputdir,"/",sample,"_CC.xls",sep=""),row.names=F,quote=F,sep="\t")
    }
  
  
  if(!any(gene %in% universe)){
    next
    }else{
      kk <- enricher(gene = gene,
                     TERM2GENE=kegg2gene,
                     TERM2NAME=kegg2name,
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 1,
                     qvalueCutoff  = 1,
                     minGSSize=1,
                     maxGSSize = 10000000000000)
      kk <- kk@result
      kk <- kk[kk$pvalue < 0.1,]##########or 0.1?
      #####convert geneid to gene symbol
      idconvert<-kk$geneID
      idconvert_list<-strsplit(idconvert,'/')
      for (i in 1:length(idconvert_list)){
        smyb_str<-c()
        for (j in 1:length(idconvert_list[[i]])){
          symb<-idconvert_list[[i]][j]
          smyb1<-entrz_symb[symb][[1]]
          smyb_str<-c(smyb_str,smyb1)
          symbid<-paste(smyb_str,collapse = "/")
        }
       kk$Symbolid[i]<- symbid
      }
      write.table(kk,file=paste(outputdir,"/",sample,"_KEGG.xls",sep=""),row.names=F,quote=F,sep="\t")
    }
}
