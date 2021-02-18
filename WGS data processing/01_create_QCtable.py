#!/usr/bin/env python
import os
import argparse
import textwrap

CHECK_VERSION = "v0.1"
AUTHOR = "lierhan"

def parse_raw_data_metrics(raw_data_metrics):
    raw_data_metrics_list = []
    with open(raw_data_metrics, 'r') as myReaderFile:
	lines = myReaderFile.readlines()
	for line in lines:
	    line = line.rstrip()
	    if line.startswith("Yield"):
		continue
	    else:
		raw_data_metrics_list.append(line.split('\t'))
    raw_total_throughput = float(raw_data_metrics_list[0][0]) + float(raw_data_metrics_list[1][0])
    raw_total_reads = int(raw_data_metrics_list[0][1]) + int(raw_data_metrics_list[1][1])
    return raw_total_throughput, raw_total_reads

def parse_alignment_summary_metrics(alignment_summary_metrics):
    alignment_summary_metrics_dict = {"FIRST_OF_PAIR":{}, "SECOND_OF_PAIR":{}, "PAIR":{}}
    with open(alignment_summary_metrics, 'r') as myReaderFile:
	lines = myReaderFile.readlines()
	for line in lines:
	    line = line.rstrip()
	    if line.startswith("CATEGORY"):
		header_list = line.split('\t')
	    elif line.startswith("FIRST_OF_PAIR"):
		FIRST_OF_PAIR_list = line.split('\t')
		for i in range(1, len(FIRST_OF_PAIR_list)):
		    item = header_list[i]
		    alignment_summary_metrics_dict["FIRST_OF_PAIR"][item] = FIRST_OF_PAIR_list[i]
	    elif line.startswith("SECOND_OF_PAIR"):
		SECOND_OF_PAIR_list = line.split('\t')
		for i in range(1, len(SECOND_OF_PAIR_list)):
		    item = header_list[i]
		    alignment_summary_metrics_dict["SECOND_OF_PAIR"][item] = SECOND_OF_PAIR_list[i]
	    elif line.startswith("PAIR"):
		PAIR_list = line.split('\t')
		for i in range(1, len(PAIR_list)):
		    item = header_list[i]
		    alignment_summary_metrics_dict["PAIR"][item] = PAIR_list[i]
    return alignment_summary_metrics_dict	

def parse_common_metrics(common_metrics, flag):
    common_metrics_dict = {}
    with open(common_metrics, 'r') as myReaderFile:
        while True:
            line = myReaderFile.readline()
            line = line.rstrip()
            if line.startswith(flag):
                header_list = line.split('\t')
                value_line = myReaderFile.readline()
                value_list = value_line.rstrip().split('\t')
                for i in range(len(value_list)):
                    item = header_list[i]
                    common_metrics_dict[item] = value_list[i]
                break
    return common_metrics_dict

def parse_flagstat(flagstat):
    flagstat_dict = {}
    value_list = []
    with open(flagstat, 'r') as myReaderFile:
        lines = myReaderFile.readlines()
        for line in lines:
            line = line.rstrip()
            value_list.append(line.split(' ')[0])
    flagstat_dict["total_reads"] = value_list[0]
    flagstat_dict["duplicates"] = value_list[1]
    flagstat_dict["mapped"] = value_list[2]
    flagstat_dict["paired"] = value_list[3]
    flagstat_dict["read1"] = value_list[4]
    flagstat_dict["read2"] = value_list[5]
    flagstat_dict["properly_paired"] = value_list[6]
    flagstat_dict["with_itself_and_mate_mapped"] = value_list[7]
    flagstat_dict["singletons"] = value_list[8]
    flagstat_dict["with_mate_mapped_to_a_different_chr"] = value_list[9]
    flagstat_dict["with_mate_mapped_to_a_different_chr_maq5"] = value_list[10]
    return flagstat_dict

def parse_target_base_coverage(targetCoverage):
    coverage = []
    bases = []
    TARGET_TERRITORY = 0
    ratio = []
    with open(targetCoverage, 'r') as myreadfile:
	lines = myreadfile.readlines()
	for line in lines:
	    line = line.rstrip()
	    if line.startswith("all"):
		line_list = line.split('\t')
		coverage.append(int(line_list[1]))
		bases.append(int(line_list[2]))
		TARGET_TERRITORY = int(line_list[3])
		ratio.append(float(line_list[4]))
    return coverage, bases, TARGET_TERRITORY, ratio

def calculate_capture_efficiency(reads_TargetReadRegion, reads_TargetReadRegion150, reads_TargetReadRegion500):
    with open(reads_TargetReadRegion, 'r') as target, open(reads_TargetReadRegion150, 'r') as target150, open(reads_TargetReadRegion500, 'r') as target500:
	num_target = int(target.readline().rstrip())
	num_target150 = int(target150.readline().rstrip())
	num_target500 = int(target500.readline().rstrip())
    return num_target, num_target150, num_target500

def main(raw_data_metrics, alignment_summary_metrics, coverageMetrics, flagstat, insert_size_metrics, prefix, analysis_type, *args):
    ## parse raw data metrics
    fqstat_raw_total_throughput, fqstat_raw_total_reads = parse_raw_data_metrics(raw_data_metrics)

    ## load alignment_summary, insert_size, flagstat and dedup metrics, then prepare metrics dict
    alignment_summary_metrics_dict = parse_alignment_summary_metrics(alignment_summary_metrics)
    insert_size_metrics_dict = parse_common_metrics(insert_size_metrics, 'MEDIAN_INSERT_SIZE')
    flagstat_dict  = parse_flagstat(flagstat)

    ## extract value from alignment_summary, insert_size, flagstat and dedup metrics dict
    picard_total_reads = int(alignment_summary_metrics_dict["PAIR"]["TOTAL_READS"])
    picard_pct_reads_aligned = float(alignment_summary_metrics_dict["PAIR"]["PCT_PF_READS_ALIGNED"]) * 100
    picard_mean_read_length = float(alignment_summary_metrics_dict["PAIR"]["MEAN_READ_LENGTH"])
    picard_pct_reads_improper_pairs = float(alignment_summary_metrics_dict["PAIR"]["PCT_PF_READS_IMPROPER_PAIRS"]) * 100
    picard_mean_insert_size = float(insert_size_metrics_dict["MEAN_INSERT_SIZE"])
    flagstat_total_reads = int(flagstat_dict["total_reads"])
    flagstat_singletons = int(flagstat_dict["singletons"])
    flagstat_mate_mapped_diff_chr = int(flagstat_dict["with_mate_mapped_to_a_different_chr_maq5"])

    ## calculate final value to present in QC table
    picard_pct_reads_proper_pairs = (100.0 - picard_pct_reads_improper_pairs) * picard_pct_reads_aligned / 100
    picard_total_throughput = (picard_total_reads * picard_mean_read_length) / 1000000000
    flagstat_pct_singleton = flagstat_singletons * 100.0 / flagstat_total_reads
    flagstat_pct_mate_mapping_diff_chr = flagstat_mate_mapped_diff_chr * 100.0 / flagstat_total_reads

    if analysis_type == "WGS":
	# load raw wgs metrics
        raw_wgs_metrics_dict = parse_common_metrics(coverageMetrics, 'GENOME_TERRITORY')        
	picard_pct_duplication = float(raw_wgs_metrics_dict["PCT_EXC_DUPE"]) * 100
        picard_mean_coverage_depth = float(raw_wgs_metrics_dict["MEAN_COVERAGE"])
        picard_pct_coverage_1X = float(raw_wgs_metrics_dict["PCT_1X"]) * 100
        picard_pct_coverage_5X = float(raw_wgs_metrics_dict["PCT_5X"]) * 100
        picard_pct_coverage_10X = float(raw_wgs_metrics_dict["PCT_10X"]) * 100
        picard_pct_coverage_20X = float(raw_wgs_metrics_dict["PCT_20X"]) * 100

        ## write value to QC table
        output =  prefix + "_sorted_dedup_QCtable.txt"
        with open(output, 'w') as myWriterFile:
#            myWriterFile.write("QC Statistics Items\t%s\n" % prefix)
            myWriterFile.write("Total raw data yield(G)\t%.3f\n" % (fqstat_raw_total_throughput / 1000))
            myWriterFile.write("Total raw data reads(M)\t%.3f\n" % (fqstat_raw_total_reads * 1.0 /1000000))
            myWriterFile.write("Raw data read length\t150\n")
            myWriterFile.write("Mean insert size\t%.3f\n" % picard_mean_insert_size)
            myWriterFile.write("Percentage mapped reads\t%.3f%%\n" % picard_pct_reads_aligned)
            myWriterFile.write("Percentage properly paired reads\t%.3f%%\n" % picard_pct_reads_proper_pairs)
            myWriterFile.write("Percentage singletons\t%.3f%%\n" % flagstat_pct_singleton)
            myWriterFile.write("Percentage with mate mapped to a different chr(mapQ>=5)\t%.3f%%\n" % flagstat_pct_mate_mapping_diff_chr)
            myWriterFile.write("Mean coverage sequencing depth on genome region\t%.3f\n" % picard_mean_coverage_depth)
            myWriterFile.write("Percentage PCR duplicates\t%.3f%%\n" % picard_pct_duplication)
            myWriterFile.write("Fraction of genome region covered with at least 1X\t%.3f%%\n" % picard_pct_coverage_1X)
            myWriterFile.write("Fraction of genome region covered with at least 5X\t%.3f%%\n" % picard_pct_coverage_5X)
            myWriterFile.write("Fraction of genome region covered with at least 10X\t%.3f%%\n" % picard_pct_coverage_10X)
            myWriterFile.write("Fraction of genome region covered with at least 20X\t%.3f%%\n" % picard_pct_coverage_20X)

    if analysis_type == "WES":
        # calculate capture efficency
        num_target, num_target150, num_target500 = calculate_capture_efficiency(args[0], args[1], args[2])
        sambamba_capture_efficiency = num_target * 100.0 / picard_total_reads
        sambamba_capture_efficiency_150 = num_target150 * 100.0 / picard_total_reads
        sambamba_capture_efficiency_500 = num_target500 * 100.0 / picard_total_reads
 
        ## load captrue metrics and calculate value to present in the QC table
        capture_metrics_dict = parse_common_metrics(coverageMetrics, "BAIT_SET")
	picard_pct_duplication = float(capture_metrics_dict["PCT_EXC_DUPE"]) * 100
        picard_mean_target_coverage_depth = float(capture_metrics_dict["MEAN_BAIT_COVERAGE"])
        picard_pct_coverage_1X = float(capture_metrics_dict["PCT_TARGET_BASES_1X"]) * 100
        picard_pct_coverage_2X = float(capture_metrics_dict["PCT_TARGET_BASES_2X"]) * 100
        picard_pct_coverage_10X = float(capture_metrics_dict["PCT_TARGET_BASES_10X"]) * 100
        picard_pct_coverage_20X = float(capture_metrics_dict["PCT_TARGET_BASES_20X"]) * 100
 
        ## write value to QC table
        output =  prefix + "_sorted_dedup_QCtable.txt"
        with open(output, 'w') as myWriterFile:
#            myWriterFile.write("QC Statistics Items\t%s\n" % prefix)
            myWriterFile.write("Total raw data yield(G)\t%.3f\n" % (fqstat_raw_total_throughput / 1000))
            myWriterFile.write("Total raw data reads(M)\t%.3f\n" % (fqstat_raw_total_reads * 1.0 /1000000))
            myWriterFile.write("Raw data read length\t150\n")
            myWriterFile.write("Total clean data yield(G)\t%.3f\n" % picard_total_throughput)
            myWriterFile.write("Total clean data reads(M)\t%.3f\n" % (picard_total_reads * 1.0 / 1000000))
            myWriterFile.write("Mean clean data read length\t%.3f\n" % picard_mean_read_length)
            myWriterFile.write("Mean insert size\t%.3f\n" % picard_mean_insert_size)
            myWriterFile.write("Percentage mapped reads\t%.3f%%\n" % picard_pct_reads_aligned)
            myWriterFile.write("Percentage properly paired reads\t%.3f%%\n" % picard_pct_reads_proper_pairs)
            myWriterFile.write("Percentage singletons\t%.3f%%\n" % flagstat_pct_singleton)
            myWriterFile.write("Percentage with mate mapped to a different chr(mapQ>=5)\t%.3f%%\n" % flagstat_pct_mate_mapping_diff_chr)
            myWriterFile.write("Capture efficiency rate on target regions\t%.3f%%\n" % sambamba_capture_efficiency)
            myWriterFile.write("Capture efficiency rate on or near +-150 target regions\t%.3f%%\n" % sambamba_capture_efficiency_150)
            myWriterFile.write("Capture efficiency rate on or near +-500 target regions\t%.3f%%\n" % sambamba_capture_efficiency_500)
            myWriterFile.write("Mean coverage sequencing depth on official target region\t%.3f\n" % picard_mean_target_coverage_depth)
            myWriterFile.write("Percentage PCR duplicates\t%.3f%%\n" % picard_pct_duplication)
            myWriterFile.write("Fraction of official target region covered with at least 1X\t%.3f%%\n" % picard_pct_coverage_1X)
            myWriterFile.write("Fraction of official target region covered with at least 2X\t%.3f%%\n" % picard_pct_coverage_2X)
            myWriterFile.write("Fraction of official target region covered with at least 10X\t%.3f%%\n" % picard_pct_coverage_10X)
            myWriterFile.write("Fraction of official target region covered with at least 20X\t%.3f%%\n" % picard_pct_coverage_20X)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog="create_QCtable.py", formatter_class=argparse.RawDescriptionHelpFormatter,
    description=textwrap.dedent("""\
    Create QCtable for WGS or WES standard analysis from multi quality metrics 
    --------------------------------------------------------------------------
    include :
      - raw data quality metrics (fqstat - ${prefix}_dataqc.xls - both) 
      - alignment summary metrics (picard - ${prefix}.alignment_summary_metrics - both)
      - captrue metrics (picard - ${prefix}.capture_metrics - WES only)
      - raw wgs metrics (picard - ${prefix}.raw_wgs_metrics - WGS only)
      - flagstat metrics (sambamba - ${prefix}_flagstat.txt - both)
      - insert size metrics (picard - ${prefix}.insert_size_metrics - both)
      - capture efficiency on target region (sambamba - ${prefix}_TargetReadRegion.txt - WES only)
      - capture efficiency on target region and flank 150bp region (sambamba - ${prefix}_TargetReadRegion150.txt - WES only)
      - capture efficiency on target region and flank 500bp region (sambamba - ${prefix}_TargetReadRegion500.txt - WES only)
    --------------------------------------------------------------------------
    Examples :
	for WGS :
        python creat_QCtable.py \\
                -r ${prefix}_dataqc.xls \\
                -a ${prefix}.alignment_summary_metrics \\
                -c ${prefix}.raw_wgs_metrics \\
                -f ${prefix}_flagstat.txt \\
                -i ${prefix}.insert_size_metrics \\
                -o ${prefix} \\
                -w WGS

	for WES :
	python creat_QCtable.py \\
		-r ${prefix}_dataqc.xls \\
		-a ${prefix}.alignment_summary_metrics \\
		-c ${prefix}.capture_metrics \\
		-f ${prefix}_flagstat.txt \\
		-i ${prefix}.insert_size_metrics \\
		-o ${prefix} \\
		-p ${prefix}_TargetReadRegion.txt \\
		-p ${prefix}_TargetReadRegion150.txt \\
		-p ${prefix}_TargetReadRegion500.txt \\
		-w WES
    """))
    parser.add_argument("-r", "--rawdataMetrcis", dest="rawdataMetrcis", help="raw data Metrcis", required=True)
    parser.add_argument("-a", "--alignmentMetrics", dest="alignmentMetrics", help="alignment metrics", required=True)
    parser.add_argument("-c", "--coverageMetrics", dest="coverageMetrics", help="coverage metrics", required=True)
    parser.add_argument("-f", "--flagstat", dest="flagstat", help="flagstat metrics", required=True)
    parser.add_argument("-i", "--insertsize", dest="insertsize", help="insert size metrics", required=True)
    parser.add_argument("-o", "--ouput", dest="ouput", help="prefix of output file", required=True)
    parser.add_argument("-p", "--capture", dest="capture", help="capture efficience at target and flank region", action='append')
    parser.add_argument("-w", "--type", dest="type", help="analysis type", required=True)

    ### parse arguments
    args = parser.parse_args()
    raw_data_metrics = args.rawdataMetrcis
    alignment_summary_metrics = args.alignmentMetrics
    
    flagstat = args.flagstat
    insert_size_metrics = args.insertsize
    prefix = args.ouput
    analysis_type = args.type

    ## check if the parameters are valid and start the process
    if analysis_type == "WGS":
        raw_wgs_metrics = args.coverageMetrics
	main(raw_data_metrics, alignment_summary_metrics, raw_wgs_metrics, flagstat, insert_size_metrics, prefix, analysis_type)
    elif analysis_type == "WES":
	capture_metrics = args.coverageMetrics
	capture_efficiency = args.capture
	capture_efficiency.sort()
	if len(capture_efficiency) == 3:
	    main(raw_data_metrics, alignment_summary_metrics, capture_metrics, flagstat, insert_size_metrics, prefix, analysis_type, capture_efficiency[0], capture_efficiency[1], capture_efficiency[2])
	else:
	    print "Three TargetReadRegion metrics files should be provide!"
	    sys.exit(1)
    else:
	print "Invalid analysis type! Only WGS or WES is allowed!"
	sys.exit(1)
