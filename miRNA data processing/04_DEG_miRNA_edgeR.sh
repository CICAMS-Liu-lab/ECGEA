${Rscript} ${DEG_edgeR} \
../summary/miRNA.matrix.count.txt \
${ctrl_idx} \
${case_idx} \
${ctrl} \
${case} \
${fc_cut} \
${fdr_cut}
${perl} ${target_select} -i DEG_${prefix}_sig.up.txt -t ${miRGate_miRNA_target} -n 3 -head F |cut -f1 |sort -u > DEG_${prefix}_sig.up.tar.name.txt
${perl} ${target_select} -i DEG_${prefix}_sig.down.txt -t ${miRGate_miRNA_target} -n 3 -head F |cut -f1 |sort -u > DEG_${prefix}_sig.down.tar.name.txt

${perl} ${target_select} -i DEG_${prefix}_sig.up.txt -t ${miRGate_miRNA_target} -n 3 -head F|cut -f1-7 |awk 'BEGIN{OFS="\t";print "Gene_Name","Gene_ID","miRNA_ID","Source","Info","Method","Confidence"}{print $0}' > DEG_${prefix}_sig.up.txt.tar.xls
${perl} ${target_select} -i DEG_${prefix}_sig.down.txt -t ${miRGate_miRNA_target} -n 3 -head F|cut -f1-7 |awk 'BEGIN{OFS="\t";print "Gene_Name","Gene_ID","miRNA_ID","Source","Info","Method","Confidence"}{print $0}' > DEG_${prefix}_sig.down.txt.tar.xls



${Rscript} ${Gene_sybmol_to_entrez_id} ./

rename _sig.up.tar.name.txt.entrez.gene_list.txt _sig_up_tar_gene_list.txt *
rename _sig.down.tar.name.txt.entrez.gene_list.txt _sig_down_tar_gene_list.txt *
${Rscript} \
${GO_enrich} \
. \

${Rscript}  ../script/CategoryPlot.r  DEG_${prefix}_sig_up_tar_gene_list_BP.xls DEG_${prefix}_sig_up_tar_gene_list_MF.xls DEG_${prefix}_sig_up_tar_gene_list_CC.xls
${Rscript}  ../script/CategoryPlot.r  DEG_${prefix}_sig_down_tar_gene_list_BP.xls DEG_${prefix}_sig_down_tar_gene_list_MF.xls DEG_${prefix}_sig_down_tar_gene_list_CC.xls
${Rscript}  ../script/bubble.R DEG_${prefix}_sig_down_tar_gene_list_KEGG.xls
${Rscript}  ../script/bubble.R DEG_${prefix}_sig_up_tar_gene_list_KEGG.xls
