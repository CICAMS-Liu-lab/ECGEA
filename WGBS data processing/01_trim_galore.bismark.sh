cd ../
wd=`pwd`
if [ ! -d "/native/$USER" ];then
   mkdir /native/$USER
fi

rsync -hltr  Sample_${prefix} /native/$USER/
cd /native/$USER/Sample_${prefix}

export PATH=/nextcode/sge_software/anaconda2/bin:$PATH


${trim_galore}   --output_dir ./ --fastqc --fastqc_args "-t 16" --paired --clip_R1 6 --clip_R2 6 --trim1 --adapter AGATCGGAAGAGC --adapter2 AGATCGGAAGAGC --path_to_cutadapt ${cutadapt} ${prefix}_R1.fastq.gz ${prefix}_R2.fastq.gz


${bismark} ${methyl_database} -1 ${prefix}_R1_val_1.fq.gz -2 ${prefix}_R2_val_2.fq.gz --path_to_bowtie ${path_to_bowtie} --samtools_path ${samtools_path} -I 0 -X 1000 -p 2 --parallel 5  --score_min L,0,-0.6 --bowtie2 -N 1 -L 30 -o . --temp_dir ./temp_01/

${sentieon} util sort -o ${prefix}_sorted.bam -t 16 -i ${prefix}_R1_val_1_bismark_bt2_pe.bam
/nextcode/sge_software/sambamba/0.5.4/sambamba_v0.5.4 markdup -t 16  -r ${prefix}_sorted.bam ${prefix}_sorted_redup.bam
/nextcode/sge_software/sambamba/0.5.4/sambamba_v0.5.4 sort -n -t 16 -o ${prefix}_sorted_redup_sortn.bam ${prefix}_sorted_redup.bam


${bismark_methylation_extractor}  ${prefix}_sorted_redup_sortn.bam -p --no_overlap --bedGraph --cytosine_report --CX --genome_folder ${methyl_database} --buffer_size 30G --counts --parallel 5 --gzip --merge_non_CpG 

gunzip *CX_report.txt.gz

${perl} ${script_dir}/perl_cpvalue_v2.txt ${prefix}_R1_val_1_bismark_bt2_pe.CX_report.txt ${prefix}.CX_report.txt 0.01
${perl} ${script_dir}/split_cx.txt ${prefix}.CX_report.txt 2>split_cx.e
${perl} ${script_dir}/produce_density_work.txt ./ 
${perl} ${script_dir}/density_10k_figure.txt ./ 

rsync  -hltr --exclude *txt.gz /native/$USER/Sample_${prefix}/* $wd/Sample_${prefix}/
rm  -rf /native/$USER/Sample_${prefix}
cd ../
wd=`pwd`
if [ ! -d "/native/$USER" ];then
   mkdir /native/$USER
fi

mkdir -p /native/$USER/Sample_${prefix}_1
rsync  -hltr Sample_${prefix}/${prefix}.CX_report.txt  /native/$USER/Sample_${prefix}_1/
rsync  -hltr Sample_${prefix}/${prefix}_R1_val_1_bismark_bt2_PE_report.txt     /native/$USER/Sample_${prefix}_1/

cd /native/$USER/Sample_${prefix}_1/

export PATH=/nextcode/sge_software/anaconda2/bin:$PATH

${perl} ${script_dir}/C_profile.txt ${prefix}.CX_report.txt >${prefix}.C_profile.stat 

grep chr_lambda ${prefix}.CX_report.txt >chr_lambda.CX_report.txt
${perl} ${script_dir}/stat_conversion.txt chr_lambda.CX_report.txt >${prefix}.conversion.res 

${perl} ${script_dir}/sam_report_pe_stat.txt   ${prefix}_R1_val_1_bismark_bt2_PE_report.txt  > ${prefix}.mapped.res

${perl} ${script_dir}/genome_meth_levle_profile.txt  ${prefix}.CX_report.txt >${prefix}.genome_methy_level.res 
${perl} ${script_dir}/proportion_pie.txt ${prefix}.genome_methy_level.res ${prefix}

rsync  -hltr  /native/$USER/Sample_${prefix}_1/*res $wd/Sample_${prefix}/
rsync  -hltr  /native/$USER/Sample_${prefix}_1/*stat $wd/Sample_${prefix}/

rm  -rf /native/$USER/Sample_${prefix}_1 
