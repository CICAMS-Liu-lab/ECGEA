cd ../
wd=`pwd`
if [ ! -d "/native/$USER" ];then
   mkdir /native/$USER
fi

mkdir /native/$USER/Pair_${case_idx}_vs_${ctrl_idx}

rsync -hltr  Sample_${case_idx}/${case_idx}.CX_report.txt /native/$USER/Pair_${case_idx}_vs_${ctrl_idx}/
rsync -hltr  Sample_${ctrl_idx}/${ctrl_idx}.CX_report.txt /native/$USER/Pair_${case_idx}_vs_${ctrl_idx}/

cd /native/$USER/Pair_${case_idx}_vs_${ctrl_idx}

export PATH=/nextcode/sge_software/BEDTools/currentVersion/bin/:$PATH

cat ${case_idx}.CX_report.txt  | awk 'BEGIN{FS="\t";OFS="\t"}{if($1!~/_/){if($6=="CG"){if($4+$5!=0){per=100*$4/($4+$5)}else{per=0};end=$2+1;print $1,$2,end,per}}}' > ${case_idx}_bedgraph.bed
cat ${ctrl_idx}.CX_report.txt  | awk 'BEGIN{FS="\t";OFS="\t"}{if($1!~/_/){if($6=="CG"){if($4+$5!=0){per=100*$4/($4+$5)}else{per=0};end=$2+1;print $1,$2,end,per}}}' > ${ctrl_idx}_bedgraph.bed


bedtools sort -i ${case_idx}_bedgraph.bed >${case_idx}.sorted.bed
bedtools sort -i ${ctrl_idx}_bedgraph.bed >${ctrl_idx}.sorted.bed

perl ${metilene}/metilene_input.pl --in1 ${case_idx}.sorted.bed --in2 ${ctrl_idx}.sorted.bed --out ${case_idx}_${ctrl_idx}.metilene.input --h1 ${case_idx} --h2 ${ctrl_idx} 
${metilene}/metilene -t 8 -M 100 -m 5 -d 0.1  -a ${case_idx} -b ${ctrl_idx} ${case_idx}_${ctrl_idx}.metilene.input >${case_idx}_${ctrl_idx}.metilene.result

grep "^chr" ${case_idx}_${ctrl_idx}.metilene.result > ${case_idx}_${ctrl_idx}.metilene.result.filter

export PATH=/nextcode/sge_software/anaconda2/bin:$PATH

perl ${metilene}/metilene_output.pl -q ${case_idx}_${ctrl_idx}.metilene.result.filter -o ${case_idx}_${ctrl_idx}_metilene_CG  -p 0.05 -c 5 -d 0.1  -a ${case_idx} -b ${ctrl_idx}

sort -g -k4 ${case_idx}_${ctrl_idx}_metilene_CG_qval.0.05.out |awk 'BEGIN{FS="\t";OFS="\t";print "chr","start","stop","qvalue","mean_methylation_difference" ,"CpGs","mean_in_Case","mean_in_Control"}{print $0}' > ${case_idx}_${ctrl_idx}.metilene.result.xls


export PATH=/nextcode/sge_software/anaconda2/bin/:$PATH

sort -k1,1 -k2,2n ${case_idx}_${ctrl_idx}.metilene.result.xls |sed '1d'  > metilene.sorted.bed

bedtools intersect -a metilene.sorted.bed -b /nextcode/nfs_test/users/song_yunjie/WGBS/hg19_database/genome/hg19_refgene.sorted.gtf -wao > ${case_idx}_${ctrl_idx}_metilene_DMR.annotation.xls

sed -i '1i chr\tstart\tstop\tqvalue\tmean_methylation_difference\tCpGs\tmean_in_Case\tmean_in_Control\tchr\tstart\tend\tfeature\tstrand\ttranscript_id\toverlapped_base'   ${case_idx}_${ctrl_idx}_metilene_DMR.annotation.xls

cat  ${case_idx}_${ctrl_idx}_metilene_DMR.annotation.xls |awk '{print $(NF-1)}' |grep -v "transcript_id" |grep -v "\."|sort |uniq >${case_idx}_${ctrl_idx}_metilene_uniq.transcript.list



${Rscript} $wd/script/transcript_to_gene.R ${case_idx}_${ctrl_idx}_metilene_uniq.transcript.list

${Rscript} $wd/script/Gene_sybmol_to_entrez_id.R ./ metilene ${case_idx} ${ctrl_idx}
${Rscript} $wd/script/deg_GOstats_Hs_WGBS.R ./ metilene ${case_idx} ${ctrl_idx}
${Rscript} $wd/script/GO_bar_plot.R BP.xls MF.xls CC.xls metilene ${case_idx} ${ctrl_idx}
${Rscript} $wd/script/KEGG.bubble_plot.R KEGG.xls metilene ${case_idx} ${ctrl_idx}

rsync  -hltr  /native/$USER/Pair_${case_idx}_vs_${ctrl_idx}/*.pdf $wd/Pair_${case_idx}_vs_${ctrl_idx}/
rsync  -hltr  /native/$USER/Pair_${case_idx}_vs_${ctrl_idx}/*.xls $wd/Pair_${case_idx}_vs_${ctrl_idx}/
rsync  -hltr /native/$USER/Pair_${case_idx}_vs_${ctrl_idx}/*dmr.gene.list $wd/Pair_${case_idx}_vs_${ctrl_idx}/
rsync  -hltr /native/$USER/Pair_${case_idx}_vs_${ctrl_idx}/*tiff $wd/Pair_${case_idx}_vs_${ctrl_idx}/


rm  -rf /native/$USER/Pair_${case_idx}_vs_${ctrl_idx}


