args <- commandArgs(TRUE)
if(length(args)!=4){
	cat("Usage: Rscript deg_GOstats_Hs_WGBS.R enrich_dir <type> <case_idx> <ctrl_idx>")
	q()
}
cwd=args[1]
setwd(cwd)
type=args[2]
case=args[3]
ctrl=args[4]

library("org.Hs.eg.db")
library("GO.db")
library("GOstats")
library("annotate")
library("DBI")

entrez_object <- org.Hs.egGO
mapped_genes <- mappedkeys(entrez_object)
universe <- mapped_genes
anno="org.Hs.eg.db"
dataset="org.Hs.eg"

go_categories=c("BP","MF","CC")
filepost="_MF.xls"

# for(file in dir(pattern="gene_list.txt"))
#{
file=paste(case,ctrl,type,"dmr.entrez.gene_list.txt",sep="_")
	gene <- readLines(file)
	#sample <- unlist(strsplit(file,"\\."))[1]
	
	for(go_category in go_categories){
		if(!any(gene %in% universe)){
		next
		}else{
		params <- new("GOHyperGParams",geneIds=gene,universeGeneIds=universe,ontology=go_category,pvalueCutoff=0.01,conditional=F,testDirection="over",annotation=anno)
		hgOver <- hyperGTest(params)
		
		origGeneIds <- geneIds(params)
		selected <- intersect(geneIds(params), universeGeneIds(params))
		cat2Entrez <- categoryToEntrezBuilder(params)
		## get the gene (Entrez ID) in the category
		geneInCat <- lapply(as.list(summary(hgOver)[,1]),
			function(goid) {
			selected[selected %in% cat2Entrez[[goid]]]
			} )
			
		new=NULL
		for (i in 1:length(geneInCat)){
			symbs=getSYMBOL(geneInCat[[i]], data=dataset)
			new=c(new,toString(as.list(symbs)))}
		
		GeneSymbol=as.matrix(new)
		
		gGhyp.pv <- pvalues(hgOver)
		gGhyp.fdr <- p.adjust(gGhyp.pv, "BH")
		FDR=as.matrix(gGhyp.fdr[1:length(geneInCat)])
		
		addInfo=cbind(FDR,GeneSymbol)
		colnames(addInfo)=c("FDR","GeneSymbol")
		
		write.table(cbind(summary(hgOver),addInfo),file=paste(case,"_",ctrl,"_",type,"_",go_category,".xls",sep=""),row.names=F,quote=F,sep="\t")
		}	
	}
	if(!any(gene %in% universe)){
	next
	}else{
	params <- new("KEGGHyperGParams",
				geneIds=gene,
				universeGeneIds=universe,
				pvalueCutoff=0.1,
				testDirection="over",
				annotation=anno
				)
	hgOver <- hyperGTest(params)
	
	origGeneIds <- geneIds(params)
	selected <- intersect(geneIds(params), universeGeneIds(params))
	cat2Entrez <- categoryToEntrezBuilder(params)
	## get the gene (Entrez ID) in the category
	geneInCat <- lapply(as.list(summary(hgOver)[,1]),
		function(goid) {
		selected[selected %in% cat2Entrez[[goid]]]
		} )
		
	new=NULL
	for (i in 1:length(geneInCat)){
		symbs=getSYMBOL(geneInCat[[i]], data=dataset)
		new=c(new,toString(as.list(symbs)))}
	
	GeneSymbol=as.matrix(new)
	
	gGhyp.pv <- pvalues(hgOver)
	gGhyp.fdr <- p.adjust(gGhyp.pv, "BH")
	FDR=as.matrix(gGhyp.fdr[1:length(geneInCat)])
	
	addInfo=cbind(FDR,GeneSymbol)
	colnames(addInfo)=c("FDR","GeneSymbol")
	
	write.table(cbind(summary(hgOver),addInfo),file=paste(case,ctrl,type,"KEGG.xls",sep="_"),row.names=F,quote=F,sep="\t")	
	}
#}
