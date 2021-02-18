#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on 30/06/2018

@author: lierhan
"""

import codecs
import os
import sys
import re
import argparse
import time
import logging
import subprocess
import textwrap
from yaml import load

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------
def initLogging(logFilename):
    """
    Init logging
    """
    logging.basicConfig(level=logging.DEBUG,format='%(asctime)s %(name)-6s %(levelname)-6s %(message)s',datefmt='%m/%d/%Y %H:%M:%S',filename=logFilename,filemode='w')
    # define a Handler which writes WARNING messages or higher to the sys.stderr
    console = logging.StreamHandler()
    console.setLevel(logging.WARNING)
    # set a format which is simpler for console use
    formatter = logging.Formatter('%(asctime)s %(name)-6s %(levelname)-6s %(message)s')
    # tell the handler to use this format
    console.setFormatter(formatter)
    logging.getLogger('').addHandler(console)

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------
def execute_command(cmd):
    """
    execute bash command
    """
    logging.info("CommandLine: %s" % cmd)
    s = subprocess.Popen(str(cmd), stderr = subprocess.PIPE, stdout = subprocess.PIPE, shell=True)
    stderrinfo, stdoutinfo = s.communicate()
    logging.info(stderrinfo)
    logging.info(stdoutinfo)

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------
def initAnnoStat(rtg, inputVcf, outputDir, prefix, vcfType):
    """
    Init anno stat
    variant locus : variant site in input vcf
    variant allele : alt alleles in one variant site
    Note: The 'rtg' tool only count the "PASS" variants, "clustered_events" variants need to be converted to "PASS" before run 'rtg'.
    """
    inputVcfPrefix = os.path.splitext(os.path.split(inputVcf)[-1])[0]
    rtg_stats = os.path.join(outputDir, inputVcfPrefix + ".rtg.stats")
    anno_stat_dict = {'prefix':prefix, 'Total variants locus':0, 'Total variants allele':0, 'Total SNVs':'', 'Total InDels':'', 'Ti/Tv':'', 'Het/Hom':'', 'dbSNP':0, '1000G':0, 'exonic':0, 'intronic':0, 'splicing':0, 'exonic;splicing':0, 'UTR5':0, 'UTR3':0, 'UTR5;UTR3':0, 'ncRNA':0, 'ncRNA_exonic':0, 'ncRNA_intronic':0, 'ncRNA_splicing':0, 'ncRNA_exonic;splicing':0, 'ncRNA_UTR5':0, 'ncRNA_UTR3':0, 'ncRNA_UTR5;ncRNA_UTR3':0, 'upstream':0, 'downstream':0, 'upstream;downstream':0, 'intergenic':0, 'synonymous_SNV':0, 'nonsynonymous_SNV':0, 'nonframeshift_substitution':0, 'nonframeshift_insertion':0, 'nonframeshift_deletion':0, 'frameshift_insertion':0, 'frameshift_deletion':0, 'stopgain':0, 'stoploss':0}
    
    if vcfType == 'joint':
	rtgCMD = '{rtg} vcfstats {inputVcf} > {rtg_stats}'.format(rtg=rtg, inputVcf=inputVcf, rtg_stats=rtg_stats)
	execute_command(rtgCMD)
	anno_stat_dict['Ti/Tv'] = 'NA'
	anno_stat_dict['Het/Hom'] = 'NA'
	with open(rtg_stats, 'r') as readfile:
	    lines = readfile.readlines()
	    for line in lines:
		line = line.strip()
		if line.startswith('Passed Filters'):
		    anno_stat_dict['Total variants locus'] = line.split(':')[-1].replace(' ', '')
		elif line.startswith('SNPs'):
		    anno_stat_dict['Total SNVs'] = line.split(':')[-1].replace(' ', '')
	anno_stat_dict['Total InDels'] = int(anno_stat_dict['Total variants locus']) - int(anno_stat_dict['Total SNVs'])	
    elif vcfType == 'germline':
	rtgCMD = '{rtg} vcfstats {inputVcf} > {rtg_stats}'.format(rtg=rtg, inputVcf=inputVcf, rtg_stats=rtg_stats)
	execute_command(rtgCMD)
	with open(rtg_stats, 'r') as readfile:
	    lines = readfile.readlines()
	    for line in lines:
		line = line.strip()
		if line.startswith('Passed Filters'):
		    anno_stat_dict['Total variants locus'] = line.split(':')[-1].replace(' ', '')
		elif line.startswith('SNPs'):
		    anno_stat_dict['Total SNVs'] = line.split(':')[-1].replace(' ', '')
		elif line.startswith('SNP Transitions/Transversions:'):
		    anno_stat_dict['Ti/Tv'] = line.split(':')[-1].replace(' ', '')
                elif line.startswith('Total Het/Hom ratio'):
                    anno_stat_dict['Het/Hom'] = line.split(':')[-1].replace(' ', '')
	anno_stat_dict['Total InDels'] = int(anno_stat_dict['Total variants locus']) - int(anno_stat_dict['Total SNVs'])
    elif vcfType == 'somatic':
	tumor_prefix = prefix.split('_')[0]
	rtgCMD = 'awk -F \'\t\' \'BEGIN{FS=OFS="\t"}{if ($1!~/^#/) $7="PASS"}1\' %s | %s vcfstats --sample=%s - > %s' % (inputVcf, rtg, tumor_prefix, rtg_stats)
	execute_command(rtgCMD)
        with open(rtg_stats, 'r') as readfile:
            lines = readfile.readlines()
            for line in lines:
                line = line.strip()
                if line.startswith('Passed Filters'):
                    anno_stat_dict['Total variants locus'] = line.split(':')[-1].replace(' ', '')
                elif line.startswith('SNPs'):
                    anno_stat_dict['Total SNVs'] = line.split(':')[-1].replace(' ', '')
                elif line.startswith('SNP Transitions/Transversions:'):
                    anno_stat_dict['Ti/Tv'] = line.split(':')[-1].replace(' ', '')
                elif line.startswith('Total Het/Hom ratio'):
                    anno_stat_dict['Het/Hom'] = line.split(':')[-1].replace(' ', '')
        anno_stat_dict['Total InDels'] = int(anno_stat_dict['Total variants locus']) - int(anno_stat_dict['Total SNVs'])

    return anno_stat_dict

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------
def load_config(configFile):
    """
    load and parse annotation database config file
    """
    with open(configFile,'r') as readFile:
        config = readFile.read()
        configParse = load(config)
        db_version = configParse['version']
        db_path = configParse['path']
        db_buildver = configParse['buildver']
        db_gene_items = configParse['gene'].keys()
        db_region_items = configParse['region'].keys()
        db_filter_items = configParse['filter'].keys()

    return db_version,db_path,db_buildver,db_gene_items,db_region_items,db_filter_items

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------
def check_annotaion_items(anno_gene_items,anno_filter_items,anno_region_items,db_gene_items,db_region_items,db_filter_items):
    """
    check if the query annotation items is valid
    """
    logging.info("Start checking annotation items......")
    if len(anno_gene_items) + len(anno_region_items) + len(anno_filter_items) == 0:
        logging.error("At least one annotation item need to be specified!!!")
        sys.exit()
    else:
        if len(anno_gene_items) != 0:
            for item in anno_gene_items:
                if item not in db_gene_items:
                    logging.error("Invalid gene annotation item %s!!!" % item)
                    sys.exit()
            logging.info("gene annotation items : OK")
        else:
            logging.info("gene annotation items : None")
        if len(anno_region_items) != 0:
            for item in anno_region_items:
                if item not in db_region_items:
                    logging.error("Invalid region annotation item %s!!!" % item)
                    sys.exit()
            logging.info("region annotation items : OK")
        else:
            logging.info("region annotation items : None")
        if len(anno_filter_items) != 0:
            for item in anno_filter_items:
                if item not in db_filter_items:
                    logging.error("Invalid filter annotaion item %s!!!" % item)
                    sys.exit()
            logging.info("filter annotation items : OK")
        else:
            logging.info("filter annotation items : None")

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------
def break_multi_alt_allele(inputVcf, inputVcf_dup, inputVcf_annovar, vcfType):
    """
    break multiple alt allele
    """
    logging.debug("break multiple alt allele")
    finputVcf_dup = open(inputVcf_dup, 'w')
    finputVcf_annovar = open(inputVcf_annovar, 'w')
    with open(inputVcf, 'r') as myreadfile:
        lines = myreadfile.readlines()
        for line in lines:
            if line.startswith("#"):
                print >> finputVcf_dup, line.strip()
                print >> finputVcf_annovar, line.strip()
            else:
                alt_list = []
                line_list = line.strip().split('\t')
                raw_alt_list = line_list[4].split(',')
                if vcfType == 'joint':
                    alt_list = raw_alt_list
                else:
                    gt = line_list[9].split(':')[0]
                    gt_list = list(set(gt.split('/')))
                    gt_list.sort()
                    ## remove '0' genotype
                    condition = lambda t:t != "0"
                    final_gt_list = list(filter(condition, gt_list))
                    for i in final_gt_list:
                        alt_list.append(raw_alt_list[int(i)-1])
                
                for alt in alt_list:
                    print >> finputVcf_dup, line.strip()
                    print >> finputVcf_annovar, '\t'.join(line_list[:4]) + "\t" + alt + "\t" + '\t'.join(line_list[5:])

    finputVcf_dup.close()
    finputVcf_annovar.close()

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------
def combine_germline_anno_results(multianno_vcf, anno_stat_dict):
    """
    combine anno results, and summarize anno results
    """
    anno_column = []
    final_results = []
    sample_column = []
    with open(multianno_vcf, 'r') as readfile:
        lines = readfile.readlines()
        for line in lines:
        	line = line.strip()
        	if line.startswith("##"):
        	    if "ANNOVAR" in line:
        		pattern = re.compile(r'ID=[^,]+')
        		ID = re.findall(pattern, line)[0].split('=')[-1]
        		if ID != "ANNOVAR_DATE" and ID != "ALLELE_END" and ID not in anno_column:
        		    anno_column.append(ID)
			if ID == "ANNOVAR_DATE":
			    final_results.append(line)
		    else:
			final_results.append(line)
        	elif line.startswith("#CHROM"):
        	    line_list = line.split('\t')
        	    for sample in line_list[9:]:
                        sample_column += [sample+"_Geno", sample+"_ALfreq",sample+"_DP"]
        	    header_line = line + '\t' + '\t'.join(anno_column) + '\t' + '\t'.join(sample_column)
        	    final_results.append(header_line)
        	else:
        	    line_list = line.split('\t')
        	    ## alt column
        	    alt = line_list[4]
        	    ## alt alleles
        	    alt_list = alt.split(',')
        	    ## info column
        	    info = line_list[7]
        	    ## extract simple info
        	    simple_info_pattern = re.compile(r'.+?ANNOVAR_DATE=.+?;')
        	    simple_info = re.findall(simple_info_pattern, info)[0][:-1]
        	    ## extract annotation info
        	    anno_info_pattern = re.compile(r'ANNOVAR_DATE.+?ALLELE_END')
        	    anno_info_list = re.findall(anno_info_pattern, info)
        	    ## format column 
        	    FORMAT = line_list[8]
        	    FORMAT_list = FORMAT.split(':')
        	    ## extract gt
        	    gt = line_list[9].split(':')[FORMAT_list.index('GT')]
        	    gt_list = list(set(gt.split('/')))
        	    gt_list.sort()
                    condition = lambda t:t != "0"
                    final_gt_list = list(filter(condition, gt_list))
        	    ## extract ad
        	    ad = line_list[9].split(':')[FORMAT_list.index('AD')]
        	    ad_list = ad.split(',')
        	    ## extract dp
        	    dp = line_list[9].split(':')[FORMAT_list.index('DP')]
        	    for gt_value in final_gt_list:
        		anno_stat_dict['Total variants allele'] += 1
        		## extract anno info
        		anno_info = []
        		anno_info_list_list = anno_info_list[int(gt_value) - 1].split(';')
        	        for anno_info_item in anno_info_list_list:
        	            if anno_info_item.startswith('ANNOVAR_DATE') or anno_info_item.startswith('ALLELE_END'):
        	                continue
        	            else:
        	                key = anno_info_item.split('=')[0]
        	                value = anno_info_item.split('=')[-1].decode('unicode_escape')
        	                anno_info.append(value)
        	                if key == 'Func_refGene' and value in anno_stat_dict.keys():
        	                    anno_stat_dict[value] += 1
        	                elif key == 'ExonicFunc_refGene' and value in anno_stat_dict.keys():
        	                    anno_stat_dict[value] += 1
        	                elif re.search('^(av){0,1}snp\d{3}$', key) and value != '.':
				    anno_stat_dict['dbSNP'] += 1
				## 1000g2015aug_all 1000G_ALL
        	                elif re.search('^1000(g|G).*_(ALL|all)$',key) and value != '.':   
        	                    anno_stat_dict['1000G'] += 1
        
        		## extract gt, af, dp info
            		if dp == '0':
            		    af = 0
            		    gt_af_dp = [gt, str(af), str(dp)]
            		else:
            		    af = format(float(ad_list[int(gt_value)])/float(dp), '.3f')
            		    gt_af_dp = [gt, str(af), str(dp)]
            		final_results.append('\t'.join(line_list[:7])+'\t'+simple_info+'\t'+'\t'.join(line_list[8:])+'\t'+'\t'.join(anno_info)+'\t'+'\t'.join(gt_af_dp))
    	
    return final_results, anno_stat_dict

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------
def combine_somatic_anno_results(multianno_vcf, anno_stat_dict):
    """
    combine anno results, and summarize anno results
    """
    anno_column = []
    final_results = []
    sample_column = []
    with open(multianno_vcf, 'r') as readfile:
        lines = readfile.readlines()
        for line in lines:
        	line = line.strip()
        	if line.startswith("##"):
        	    if "ANNOVAR" in line:
        		pattern = re.compile(r'ID=[^,]+')
        		ID = re.findall(pattern, line)[0].split('=')[-1]
        		if ID != "ANNOVAR_DATE" and ID != "ALLELE_END" and ID not in anno_column:
        		    anno_column.append(ID)
			if ID == "ANNOVAR_DATE":
			    final_results.append(line)
		    else:
			final_results.append(line)
        	elif line.startswith("#CHROM"):
        	    line_list = line.split('\t')
        	    for sample in line_list[9:]:
                        sample_column += [sample+"_Geno", sample+"_ALfreq",sample+"_DP"]
        	    header_line = line + '\t' + '\t'.join(anno_column) + '\t' + '\t'.join(sample_column)
        	    final_results.append(header_line)
        	else:
        	    line_list = line.split('\t')
        	    ## alt column
        	    alt = line_list[4]
        	    ## alt alleles
        	    alt_list = alt.split(',')
        	    ## info column
        	    info = line_list[7]
        	    ## extract simple info
        	    simple_info_pattern = re.compile(r'.+?ANNOVAR_DATE=.+?;')
        	    simple_info = re.findall(simple_info_pattern, info)[0][:-1]
        	    ## extract annotation info
        	    anno_info_pattern = re.compile(r'ANNOVAR_DATE.+?ALLELE_END')
        	    anno_info_list = re.findall(anno_info_pattern, info)
        	    ## format column 
        	    FORMAT = line_list[8]
        	    FORMAT_list = FORMAT.split(':')
                    ## extract tumor sample gt value
		    tumor_gt = line_list[9].split(':')[FORMAT_list.index('GT')]
                    if tumor_gt == "0/0":
                        tumor_gt = line_list[10].split(':')[FORMAT_list.index('GT')] 
		    tumor_gt_list = list(set(tumor_gt.split('/')))
		    tumor_gt_list.sort()
		    condition = lambda t:t != "0"
		    final_tumor_gt_list = list(filter(condition, tumor_gt_list))
		    ## split multi alt allele line into multi single lines
		    for tumor_gt_value in final_tumor_gt_list:
			anno_stat_dict['Total variants allele'] += 1
                        ## extract sample info
                        sample_info = []
                        for sample in line_list[9:]:
                            gt = sample.split(':')[FORMAT_list.index('GT')]
                            ad = sample.split(':')[FORMAT_list.index('AD')]
                            ad_list = ad.split(',')
                            af = sample.split(':')[FORMAT_list.index('AF')]
			    af_list = af.split(',')
			    dp = sum(map(int, ad_list))
                            sample_info += [gt, str(af_list[0]), str(dp)]
                        ## extract anno info
                        anno_info = []
                        anno_info_list_list = anno_info_list[int(tumor_gt_value)-1].split(';')
                        for anno_info_item in anno_info_list_list:
                            if anno_info_item.startswith('ANNOVAR_DATE') or anno_info_item.startswith('ALLELE_END'):
                                continue
                            else:
                                key = anno_info_item.split('=')[0]
                                value = anno_info_item.split('=')[-1].decode('unicode_escape')
                                anno_info.append(value)
                                if key == 'Func_refGene' and value in anno_stat_dict.keys():
                                    anno_stat_dict[value] += 1
                                elif key == 'ExonicFunc_refGene' and value in anno_stat_dict.keys():
                                    anno_stat_dict[value] += 1
                                elif re.search('^(av){0,1}snp\d{3}$', key) and value != '.':
                                    anno_stat_dict['dbSNP'] += 1
                                ## 1000g2015aug_all 1000G_ALL
                                elif re.search('^1000(g|G).*_(ALL|all)$',key) and value != '.':
                                    anno_stat_dict['1000G'] += 1

                        final_results.append('\t'.join(line_list[:7])+'\t'+simple_info+'\t'+'\t'.join(line_list[8:])+'\t'+'\t'.join(anno_info)+'\t'+'\t'.join(sample_info))

    return final_results, anno_stat_dict

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------
def combine_joint_anno_results(multianno_vcf, anno_stat_dict):
    """
    combine anno results, and summarize anno results
    """
    anno_column = []
    final_results = []
    sample_column = []
    with open(multianno_vcf, 'r') as readfile:
        lines = readfile.readlines()
        for line in lines:
        	line = line.strip()
        	if line.startswith("##"):
        	    if "ANNOVAR" in line:
        		pattern = re.compile(r'ID=[^,]+')
        		ID = re.findall(pattern, line)[0].split('=')[-1]
        		if ID != "ANNOVAR_DATE" and ID != "ALLELE_END" and ID not in anno_column:
        		    anno_column.append(ID)
			if ID == "ANNOVAR_DATE":
			    final_results.append(line)
		    else:
			final_results.append(line)
        	elif line.startswith("#CHROM"):
        	    line_list = line.split('\t')
        	    for sample in line_list[9:]:
                        sample_column += [sample+"_Geno", sample+"_ALfreq",sample+"_DP"]
        	    header_line = line + '\t' + '\t'.join(anno_column) + '\t' + '\t'.join(sample_column)
        	    final_results.append(header_line)
        	else:
        	    line_list = line.split('\t')
        	    ## alt column
        	    alt = line_list[4]
        	    ## alt alleles
        	    alt_list = alt.split(',')
        	    ## info column
        	    info = line_list[7]
        	    ## extract simple info
        	    simple_info_pattern = re.compile(r'.+?ANNOVAR_DATE=.+?;')
        	    simple_info = re.findall(simple_info_pattern, info)[0][:-1]
        	    ## extract annotation info
        	    anno_info_pattern = re.compile(r'ANNOVAR_DATE.+?ALLELE_END')
        	    anno_info_list = re.findall(anno_info_pattern, info)
        	    ## format column 
        	    FORMAT = line_list[8]
        	    FORMAT_list = FORMAT.split(':')
		    ## split multi alt allele line into multi single lines
		    for allele_index in range(len(alt_list)):
			anno_stat_dict['Total variants allele'] += 1
			## extract sample info
			sample_info = []
			for sample in line_list[9:]:
			    gt = sample.split(':')[FORMAT_list.index('GT')]
			    ad = sample.split(':')[FORMAT_list.index('AD')]
			    ad_list = ad.split(',')
			    dp = sample.split(':')[FORMAT_list.index('DP')]
			    if dp == '0':
				sample_info += [gt, '0', str(dp)]
			    else:
				af = format(float(ad_list[int(allele_index)+1])/float(dp), '.3f')
				sample_info += [gt, str(af), str(dp)]
        		## extract anno info
        		anno_info = []
        		anno_info_list_list = anno_info_list[int(allele_index)].split(';')
        	        for anno_info_item in anno_info_list_list:
        	            if anno_info_item.startswith('ANNOVAR_DATE') or anno_info_item.startswith('ALLELE_END'):
        	                continue
        	            else:
        	                key = anno_info_item.split('=')[0]
        	                value = anno_info_item.split('=')[-1].decode('unicode_escape')
        	                anno_info.append(value)
        	                if key == 'Func_refGene' and value in anno_stat_dict.keys():
        	                    anno_stat_dict[value] += 1
        	                elif key == 'ExonicFunc_refGene' and value in anno_stat_dict.keys():
        	                    anno_stat_dict[value] += 1
        	                elif re.search('^(av){0,1}snp\d{3}$', key) and value != '.':
				    anno_stat_dict['dbSNP'] += 1
				## 1000g2015aug_all 1000G_ALL
        	                elif re.search('^1000(g|G).*_(ALL|all)$',key) and value != '.':   
        	                    anno_stat_dict['1000G'] += 1

			final_results.append('\t'.join(line_list[:7])+'\t'+simple_info+'\t'+'\t'.join(line_list[8:])+'\t'+'\t'.join(anno_info)+'\t'+'\t'.join(sample_info))
    	
    return final_results, anno_stat_dict

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------
def write_anno_stat(anno_stat_dict, anno_stat_file, vcfType):
    dbSNP_percentage = anno_stat_dict['dbSNP']*100.0/anno_stat_dict['Total variants allele']
    oneKG_percentage = anno_stat_dict['1000G']*100.0/anno_stat_dict['Total variants allele']
    with open(anno_stat_file, 'w') as writeFile:
        writeFile.write("Variant Statistics Items\t%s\n" % anno_stat_dict['prefix'])
        writeFile.write("Total SNVs/InDels\t%s\n" % anno_stat_dict['Total variants locus'])
        writeFile.write("Total SNVs\t%s\n" % anno_stat_dict['Total SNVs'])
        writeFile.write("Total InDels\t%d\n" % anno_stat_dict['Total InDels'])
	if vcfType == "germline":
            writeFile.write("Ti/Tv\t%s\n" % anno_stat_dict['Ti/Tv'])
            writeFile.write("Het/Hom\t%s\n" % anno_stat_dict['Het/Hom'])
	writeFile.write("Total SNVs/InDels allele\t%d\n" % anno_stat_dict['Total variants allele'])
        writeFile.write("dbSNP\t%.2f%%(%d/%d)\n" % (dbSNP_percentage, anno_stat_dict['dbSNP'], anno_stat_dict['Total variants allele']))
        writeFile.write("1000G\t%.2f%%(%d/%d)\n" % (oneKG_percentage, anno_stat_dict['1000G'], anno_stat_dict['Total variants allele']))
        writeFile.write("exonic\t%d\n" % anno_stat_dict['exonic'])
        writeFile.write("intronic\t%d\n" % anno_stat_dict['intronic'])
        writeFile.write("splicing\t%d\n" % anno_stat_dict['splicing'])
        writeFile.write("exonic;splicing\t%d\n" % anno_stat_dict['exonic;splicing'])
        writeFile.write("UTR5\t%d\n" % anno_stat_dict['UTR5'])
        writeFile.write("UTR3\t%d\n" % anno_stat_dict['UTR3'])
        writeFile.write("UTR5;UTR3\t%d\n" % anno_stat_dict['UTR5;UTR3'])
        writeFile.write("ncRNA\t%d\n" % anno_stat_dict['ncRNA'])
        writeFile.write("ncRNA_exonic\t%d\n" % anno_stat_dict['ncRNA_exonic'])
        writeFile.write("ncRNA_intronic\t%d\n" % anno_stat_dict['ncRNA_intronic'])
        writeFile.write("ncRNA_splicing\t%d\n" % anno_stat_dict['ncRNA_splicing'])
        writeFile.write("ncRNA_exonic;splicing\t%d\n" % anno_stat_dict['ncRNA_exonic;splicing'])
        writeFile.write("ncRNA_UTR5\t%d\n" % anno_stat_dict['ncRNA_UTR5'])
        writeFile.write("ncRNA_UTR3\t%d\n" % anno_stat_dict['ncRNA_UTR3'])
        writeFile.write("ncRNA_UTR5;ncRNA_UTR3\t%d\n" % anno_stat_dict['ncRNA_UTR5;ncRNA_UTR3'])
        writeFile.write("upstream\t%d\n" % anno_stat_dict['upstream'])
        writeFile.write("downstream\t%d\n" % anno_stat_dict['downstream'])
        writeFile.write("upstream;downstream\t%d\n" % anno_stat_dict['upstream;downstream'])
        writeFile.write("intergenic\t%d\n" % anno_stat_dict['intergenic'])
        writeFile.write("synonymous SNV\t%d\n" % anno_stat_dict['synonymous_SNV'])
        writeFile.write("nonsynonymous SNV\t%d\n" % anno_stat_dict['nonsynonymous_SNV'])
        writeFile.write("nonframeshift substitution\t%d\n" % anno_stat_dict['nonframeshift_substitution'])
        writeFile.write("nonframeshift insertion\t%d\n" % anno_stat_dict['nonframeshift_insertion'])
        writeFile.write("nonframeshift deletion\t%d\n" % anno_stat_dict['nonframeshift_deletion'])
        writeFile.write("frameshift insertion\t%d\n" % anno_stat_dict['frameshift_insertion'])
        writeFile.write("frameshift deletion\t%d\n" % anno_stat_dict['frameshift_deletion'])
        writeFile.write("stopgain\t%d\n" % anno_stat_dict['stopgain'])
        writeFile.write("stoploss\t%d\n" % anno_stat_dict['stoploss'])

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(prog='vcf2anno.py',formatter_class=argparse.RawDescriptionHelpFormatter,
    description=textwrap.dedent("""\
Examples:
    For germline vcf:
	python vcf2anno.py \\
		-a /software_databases/software/annovar/20180416/ \\
		-s /software_databases/software/rtg-tools/rtg-tools-3.9.1/rtg \\
		-c /software_databases/database_NGS/hg38/humandb/hg38_avdblist.yaml \\
		-i R18044580LD01-16F0526_SNP_INDEL_recal.PASS.vcf \\
		-o . \\
		-p R18044580LD01-16F0526 \\
		-g refGene \\
		-f 1000g2015aug_all,avsnp150,clinvar_20180603,cosmic70,dbnsfp33a \\
		-r cytoband,rmsk,simpleRepeat \\
		-t 20 \\
		-v germline

    For somatic vcf:
        python vcf2anno.py \\
                -a /software_databases/software/annovar/20180416/ \\
  		-s /software_databases/software/rtg-tools/rtg-tools-3.9.1/rtg \\
                -c /software_databases/database_NGS/hg38/humandb/hg38_avdblist.yaml \\
                -i HCC1187C_vs_HCC1187BL_TN.PASS.vcf \\
                -o . \\
                -p HCC1187C \\
                -g refGene \\
                -f 1000g2015aug_all,avsnp150,clinvar_20180603,cosmic70,dbnsfp33a \\
                -r cytoband,rmsk,simpleRepeat \\
                -t 20 \\
                -v somatic

    For joint calling vcf:
        python vcf2anno.py \\
                -a /software_databases/software/annovar/20180416/ \\
		-s /software_databases/software/rtg-tools/rtg-tools-3.9.1/rtg \\
                -c /software_databases/database_NGS/hg38/humandb/hg38_avdblist.yaml \\
                -i 16F0526_SNP_INDEL_recal.PASS.vcf \\
                -o . \\
                -p 16F0526 \\
                -g refGene \\
                -f 1000g2015aug_all,avsnp150,clinvar_20180603,cosmic70,dbnsfp33a \\
                -r cytoband,rmsk,simpleRepeat \\
                -t 20 \\
                -v joint
    """))
    parser.add_argument("-a","--annovarDir",dest="annovarDir",required=True,help="Annovar tool")
    parser.add_argument("-s","--rtg",dest="rtg",required=True,help="rtg tool")
    parser.add_argument("-c","--config",dest="config",required=True,help="Annotation database config file")
    parser.add_argument("-i","--inputVcf",dest="inputVcf",required=True,help="vcf file to annotate")
    parser.add_argument("-o","--outputDir",dest="outputDir",required=True,help="output dir")
    parser.add_argument("-p","--prefix",dest="prefix",required=True,help="prefix")
    parser.add_argument("-g","--geneAnno",dest="geneAnno",help="gene annotation items")
    parser.add_argument("-f","--filterAnno",dest="filterAnno",help="filter annotation items")
    parser.add_argument("-r","--regionAnno",dest="regionAnno",help="region annotation items")
    parser.add_argument("-t","--threads",dest="threads",help="threads to use",default=6)
    parser.add_argument("-v","--vcfType",dest="vcfType",required=True,help="vcf type",choices=['germline','somatic','joint'])

    # parse arguments
    args = parser.parse_args()
    annovarDir = args.annovarDir
    rtg = args.rtg
    db_config = args.config
    inputVcf = args.inputVcf
    outputDir = args.outputDir
    prefix = args.prefix
    anno_gene_items = args.geneAnno.split(',')
    anno_region_items = args.regionAnno.split(',')
    anno_filter_items = args.filterAnno.split(',')
    threads = args.threads
    vcfType = args.vcfType

    # init logging
    inputVcfPrefix = os.path.splitext(os.path.split(inputVcf)[-1])[0]
    initLogging(os.path.join(outputDir, inputVcfPrefix + ".vcf2anno.log"))
    logging.debug("CommandLine: vcf2anno.py -a {annovarDir} -c {db_config} -i {inputVcf} -o {outputDir} -p {prefix} -g {anno_gene_items} -f {anno_filter_items} -r {anno_region_items} -t {threads} -v {vcfType}".format(annovarDir=annovarDir, db_config=db_config, inputVcf=inputVcf, outputDir=outputDir, prefix=prefix, anno_gene_items=','.join(anno_gene_items), anno_filter_items=','.join(anno_filter_items), anno_region_items=','.join(anno_region_items), threads=threads, vcfType=vcfType))

    # load annotation database config
    db_version,db_path,db_buildver,db_gene_items,db_region_items,db_filter_items = load_config(db_config)

    # check if arguments are valid
    check_annotaion_items(anno_gene_items,anno_filter_items,anno_region_items,db_gene_items,db_region_items,db_filter_items)

    # step0 : break multiple alt allele
#    inputVcf_dup = os.path.join(outputDir, inputVcfPrefix + "-dup.vcf")
#    inputVcf_annovar = os.path.join(outputDir, inputVcfPrefix + "-annovar.vcf")
#    break_multi_alt_allele(inputVcf, inputVcf_dup, inputVcf_annovar, vcfType)

    # step1 : perform annotation process
    protocol = ",".join(anno_gene_items) + "," + ",".join(anno_region_items) + "," + ",".join(anno_filter_items)
    operation = ",".join(["g"]*len(anno_gene_items)) + "," + ",".join(["r"]*len(anno_region_items)) +  "," + ",".join(["f"]*len(anno_filter_items))
    arguments = ",".join(["--hgvs"]*len(anno_gene_items)) + ","*(len(anno_region_items) + len(anno_filter_items))
    table_annovarCMD = "perl {table_annovar} {inputVcf} {db_path} -buildver {db_buildver} -out {out_prefix} -remove -protocol {protocol} -operation {operation} --argument {arguments} -nastring . -vcfinput --thread {threads} --dot2underline".format(table_annovar=os.path.join(annovarDir, "table_annovar.pl"), inputVcf=inputVcf, db_path=db_path, db_buildver=db_buildver, out_prefix=os.path.join(outputDir, inputVcfPrefix), protocol=protocol, operation=operation, arguments=arguments, threads=threads)
    execute_command(table_annovarCMD)

    # step2 : combine raw vcf and annotation results
    ## step2-1 : init anno stat dict
    logging.debug("Init annotation stat dict......")
    anno_stat_dict = initAnnoStat(rtg, inputVcf, outputDir, prefix, vcfType)
    ## step2-2 :
    logging.debug("Combine and summarize annotation results......")
    if vcfType == 'germline':
        final_results, anno_stat_dict_final = combine_germline_anno_results(os.path.join(outputDir, inputVcfPrefix+'.'+db_buildver+'_multianno.vcf'), anno_stat_dict)
    elif vcfType == 'somatic':
	final_results, anno_stat_dict_final = combine_somatic_anno_results(os.path.join(outputDir, inputVcfPrefix+'.'+db_buildver+'_multianno.vcf'), anno_stat_dict)
    elif vcfType == 'joint':
	final_results, anno_stat_dict_final = combine_joint_anno_results(os.path.join(outputDir, inputVcfPrefix+'.'+db_buildver+'_multianno.vcf'), anno_stat_dict)
    ## step2-3 : write final anno results
    anno_file = os.path.join(outputDir, inputVcfPrefix+'.anno.vcf')
    with codecs.open(anno_file, 'w', encoding='utf-8') as mywritefile:
	for line in final_results:
	    mywritefile.write(line+'\n')
    ## step2-4 : write final anno stat results
    anno_stat_file = os.path.join(outputDir, inputVcfPrefix+'.anno.vcf.stat')
    write_anno_stat(anno_stat_dict_final, anno_stat_file, vcfType)
    logging.debug("Annotation Processes Successfully Finished.")

if __name__ == "__main__":
    main()
