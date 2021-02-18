${bwa} \
mem \
-M \
-t 20 \
-R "@RG\\tID:${prefix}\\tPL:ILLUMINA\\tSM:${prefix}" \
-K 10000000 \
${bwa_index} \
${prefix}_R1.fastq.gz \
${prefix}_R2.fastq.gz \
| \
${sentieon} util sort \
-o ${prefix}_sorted.bam \
-t 20 \
--sam2bam \
-i -

${sentieon} driver \
-t 20 -i ${prefix}_sorted.bam \
--algo LocusCollector \
--fun score_info ${prefix}_score.txt

${sentieon} driver \
-t 20 -i ${prefix}_sorted.bam \
--algo Dedup  \
--score_info ${prefix}_score.txt \
--metrics ${prefix}_dedup.metrics \
${prefix}_sorted_dedup.bam

${sentieon} driver \
-t 20 \
-r ${genome} \
-i ${prefix}_sorted_dedup.bam \
--algo Realigner \
-k ${indel_mills} \
-k ${indel_1kg} \
${prefix}_sorted_dedup_realign.bam

${sentieon} driver \
-t 20 \
-r ${genome} \
-i ${prefix}_sorted_dedup_realign.bam \
--algo QualCal \
-k ${indel_1kg} \
-k ${indel_mills} \
-k ${dbsnp} \
${prefix}_recal_data.table
cd ../
wd=`pwd`
if [ ! -d "/native/$USER" ];then
   mkdir /native/$USER
fi

rsync -hltr  Sample_${case_idx}/*sorted_dedup_realign.bam* /native/$USER/Sample_${case_idx}/
rsync -hltr  Sample_${case_idx}/*recal_data.table /native/$USER/Sample_${case_idx}/
rsync -hltr  Sample_${ctrl_idx}/*sorted_dedup_realign.bam* /native/$USER/Sample_${ctrl_idx}/
rsync -hltr  Sample_${ctrl_idx}/*recal_data.table /native/$USER/Sample_${ctrl_idx}/
mkdir /native/$USER/Pair_${case}_vs_${ctrl}
cd /native/$USER/Pair_${case}_vs_${ctrl}

${sentieon} driver \
-t 16 \
-i ../Sample_${case_idx}/${case_idx}_sorted_dedup_realign.bam \
-q ../Sample_${case_idx}/${case_idx}_recal_data.table \
-i ../Sample_${ctrl_idx}/${ctrl_idx}_sorted_dedup_realign.bam \
-q ../Sample_${ctrl_idx}/${ctrl_idx}_recal_data.table \
-r ${genome} \
--algo TNhaplotyper \
--tumor_sample ${case_idx} \
--normal_sample ${ctrl_idx} \
--dbsnp ${dbsnp} \
--cosmic ${cosmic} \
${case_idx}_vs_${ctrl_idx}_TN.vcf

## exclude variants in black region, then extract PASS and clustered_events variants as final somatics variants
${bedtools} subtract -header -a ${case_idx}_vs_${ctrl_idx}_TN.vcf -b ${black_region} | awk -F '\t' '($7=="PASS" || $7=="clustered_events") || $1~"^#"' | uniq > ${case_idx}_vs_${ctrl_idx}_TN.PASS.vcf

mkdir annovar_tmp

${python} \
$wd/${vcf2anno} \
-a ${annovar} \
-s ${rtg} \
-c ${avdblist} \
-i ${case_idx}_vs_${ctrl_idx}_TN.PASS.vcf \
-o ./annovar_tmp \
-p ${case_idx}_vs_${ctrl_idx} \
-g refGene \
-f 1000g2015aug_all,1000g2015aug_afr,1000g2015aug_amr,1000g2015aug_eas,1000g2015aug_eur,1000g2015aug_sas,exac03,avsnp150,clinvar_20180603,cosmic86,dbnsfp33a \
-r cytoband,rmsk,simpleRepeat,targetScanS,tfbsConsSites \
-t 16 \
-v somatic

mv annovar_tmp/${case_idx}_vs_${ctrl_idx}_TN.PASS.anno.vcf .
mv annovar_tmp/${case_idx}_vs_${ctrl_idx}_TN.PASS.anno.vcf.stat .

rsync  -hltr --exclude *bam  /native/$USER/Pair_${case}_vs_${ctrl}/* $wd/Pair_${case}_vs_${ctrl}/
rm  -rf /native/$USER/Pair_${case}_vs_${ctrl}
cd $wd/Pair_${prefix}
