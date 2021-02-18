${RSEM_calcuate_expression} \
--bam \
--no-bam-output \
--estimate-rspd \
--paired-end \
--append-names \
-p 8 \
-time \
${prefix}.transcriptome.bam \
${RSEM_ref} \
${prefix}
