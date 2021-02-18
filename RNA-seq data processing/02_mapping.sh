${STAR} \
--genomeDir ${STAR_index} \
--sjdbGTFfile ${gene_annotation} \
--limitBAMsortRAM 40000000000 \
--runThreadN 24 \
--limitIObufferSize 500000000 \
--outFilterType BySJout \
--outFilterMismatchNmax 999 \
--outFilterMismatchNoverLmax 0.04 \
--outFilterMultimapNmax 20 \
--outFilterMatchNminOverLread 0.66 \
--outFilterIntronMotifs None \
--outSJfilterReads All \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMstrandField intronMotif \
--outSAMattrRGline ID:${prefix} SM:${prefix} PL:${platform} \
--alignSJoverhangMin 8 \
--alignSJDBoverhangMin 1 \
--alignIntronMin 20 \
--alignIntronMax 1000000 \
--alignMatesGapMax 1000000 \
--chimSegmentMin 15 \
--chimJunctionOverhangMin 15 \
--chimScoreMin 0 \
--chimScoreDropMax 20 \
--chimScoreSeparation 10 \
--chimScoreJunctionNonGTAG -1 \
--quantMode TranscriptomeSAM \
--quantTranscriptomeBan IndelSoftclipSingleend \
--outReadsUnmapped Fastx \
--readFilesIn ${prefix}_R1.trimmed.fastq.gz ${prefix}_R2.trimmed.fastq.gz \
--readFilesCommand zcat \
--outFileNamePrefix  ${prefix}.

mv ${prefix}.Aligned.sortedByCoord.out.bam ${prefix}_sorted.bam

mv ${prefix}.Aligned.toTranscriptome.out.bam ${prefix}.transcriptome.bam

${sambamba} \
index \
-t 24 \
${prefix}_sorted.bam
