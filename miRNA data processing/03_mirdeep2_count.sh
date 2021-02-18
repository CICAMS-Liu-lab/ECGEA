PATH=${PATH}

${perl} \
${mirdeep2} \
${prefix}_collapsed.fa \
${ref_genome} \
${prefix}.arf \
${ref_mature_miRNA} \
none \
${ref_hairpin_miRNA} \
-t ${species}

${Rscript} \
../script/Mean_count.R \
. \
${prefix}

#cat ${prefix}.count.txt |awk 'BEGIN{OFS="\t"}{if($1!~/hsa-(.*)-(.*)-(.*)-/) print $0}' > ${prefix}.mean.count.txt
