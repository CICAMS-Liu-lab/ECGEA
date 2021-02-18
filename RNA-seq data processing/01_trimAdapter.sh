### skewer trimming adapter
${skewer} \
-m pe \
-f sanger \
-r 0.1 \
-d 0.03 \
-q 20 \
-Q 20 \
-l 75 \
-n \
-o ${prefix} \
--quiet \
--compress \
-X \
-t 4 \
${prefix}_R1.fastq.gz \
${prefix}_R2.fastq.gz

mv ${prefix}-trimmed-pair1.fastq.gz ${prefix}_R1.trimmed.fastq.gz
mv ${prefix}-trimmed-pair2.fastq.gz ${prefix}_R2.trimmed.fastq.gz
mv ${prefix}-untrimmed-excluded-pair1.fastq.gz ${prefix}_R1.dropped.fastq.gz
mv ${prefix}-untrimmed-excluded-pair2.fastq.gz ${prefix}_R2.dropped.fastq.gz
