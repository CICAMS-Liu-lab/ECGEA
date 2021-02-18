PATH=${PATH}
gunzip -c ${prefix}*R1.fastq.gz > ${prefix}_1.fastq 
${perl} \
${mirdeep_mapper} \
${prefix}_1.fastq  \
-e \
-h \
-j  \
-k TGGAATTCTCGGGTGCC  \
-l 18  \
-m  \
-p ${bowtie_index} \
-s ${prefix}_collapsed.fa  \
-t ${prefix}.arf \
-u \
-n \
-o 16
