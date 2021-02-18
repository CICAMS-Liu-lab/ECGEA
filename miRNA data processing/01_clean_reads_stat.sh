>$2
printf "sample\ttotal_reads\tclean_reads\tpct_clean_reads\tunique_reads\tpct_unique_reads\tmapped_reads\tpct_mapped_reads\n"  >>$2
for path in `find $1 -type d -name "dir_mapper*"`;do
	cd $path
        sample=`dirname $path`
	echo "***$path"
        sample=${sample##\../Sample_}
	unique=$(( `wc -l reads_nr.fa|cut -f1  -d" "`/2 )) 
	clean=$(( `wc -l reads_no_short.fa|cut -f1  -d" "`/2 )) 
	total=$(( `wc -l reads.fa|cut -f1  -d" "`/2 )) 
	mapped_count=$((`cat *fastq_mapped|grep ">"|awk -F "x" '{sum+=$2}END{print sum}'`))
#	mapped=$(( `wc -l *fastq_mapped|cut -f1  -d" "`/2 )) 
	cd -
        printf "%s\t%s\t%s\t%.4f\t%s\t%.4f\t%s\t%.4f\n" $sample $total $clean `echo $clean/$total|bc -l` $unique `echo $unique/$clean|bc -l` $mapped_count `echo $mapped_count/$clean|bc -l` >>$2
done

