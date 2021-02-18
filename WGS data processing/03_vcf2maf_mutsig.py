#! /usr/bin/env python
import sys
import os
from glob import glob

'''
Func.refGene	exonic                  maf
##mutect
downstream      NA                      downstream
exonic          nonsynonymous SNV       Missense_Mutation
exonic;splicing nonsynonymous SNV       Missense_Mutation
exonic;splicing synonymous SNV          Synonymous
exonic;splicing unknown
exonic          stopgain                Nonsense_Mutation
exonic          stoploss                Nonstop_Mutation
exonic          synonymous SNV          Silent
exonic          unknown
intergenic      NA                      IGR
intronic        NA                      Intron
ncRNA_exonic    NA
ncRNA_intronic  NA
ncRNA_splicing  NA
splicing        NA                      Splice_Site
upstream;downstream    NA               5'Flank
upstream        NA                      5'Flank
UTR3            NA                      3'UTR
UTR5            NA                      5'UTR
UTR5;UTR3       NA
##varscan
downstream      NA                      3'Flank
exonic          frameshift deletion     Frame_Shift_Del
exonic          frameshift insertion    Frame_Shift_Ins
exonic          nonframeshift deletion  In_Frame_Del
exonic          nonframeshift insertion In_Frame_Ins
exonic;splicing nonframeshift insertion In_Frame_Ins
exonic          stopgain                Nonsense_Mutation
exonic          stoploss                Nonstop_Mutation
exonic          unknown
intergenic      NA                      IGR
intronic        NA                      Intron
NA              NA
ncRNA_exonic    NA
ncRNA_intronic  NA
ncRNA_splicing  NA
ncRNA_UTR5      NA
splicing        NA                      Splice_Site
upstream;downstream     NA              upstream;downstream
upstream        NA                      upstream
UTR3            NA                      3'UTR
UTR5            NA                      5'UTR
UTR5;UTR3       NA                      
'''

'''
Variant_Classification    effect        ANNOVAR output
Silent                    silent        
Synonymous                silent	synonymous SNV
Missense                  nonsilent
Missense_Mutation         nonsilent     nonsynonymous SNV
Nonsense                  null
Nonsense_Mutation         null          stopgain
Nonstop_Mutation          null          stoploss
Read-through              null
Frame_Shift_Del           null          frameshift deletion
Frame_Shift_Ins           null          frameshift insertion
In_Frame_Del              null          nonframeshift deletion
In_Frame_Ins              null          nonframeshift insertion
Splice                    null
Splice_Region             null
Splice_Site               null          splicing
Splice_Site_Del           null
Splice_Site_DNP           null
Splice_Site_Ins           null
Splice_Site_ONP           null
Splice_Site_SNP           null
Start_Codon_Del           null
Start_Codon_DNP           null
Start_Codon_Ins           null
Start_Codon_ONP           null
Stop_Codon_Del            null
Stop_Codon_DNP            null
Stop_Codon_Ins            null
Translation_Start_Site    null          initiator_codon_variant|start_lost
De_novo_Start             null
De_novo_Start_InFrame     null
De_novo_Start_OutOfFrame  null
IGR                       noncoding     intergenic
Intron                    noncoding     intronic
3'Flank                   noncoding
3'Promoter                noncoding
3'UTR                     noncoding     UTR3
3'-UTR                    noncoding
5'Flank                   noncoding     upstream_gene_variant
5'-Flank                  noncoding
5'Promoter                noncoding
5'UTR                     noncoding     UTR5
5'-UTR                    noncoding
downstream                noncoding     downstream
miRNA                     noncoding
NCSD                      noncoding
Non-coding_Transcript     noncoding
Promoter                  noncoding
RNA                       noncoding      
upstream                  noncoding     upstream
upstream;downstream       noncoding     upstream;downstream
'''

annovar2maf = {
'nonsynonymous_SNV':'Missense_Mutation',
'synonymous_SNV':'Silent',
'stopgain':'Nonsense_Mutation',
'stoploss':'Nonstop_Mutation',
'intergenic':'IGR',
'intronic':'Intron',
'splicing':'Splice_Site',
'upstream;downstream':"5'Flan",
'upstream':"5'Flank",
'downstream':"3'Flank",
'UTR3':"3'UTR",
'UTR5':"5'UTR",
#'UTR5;UTR3':'UTR5;UTR3',
'frameshift_deletion':'Frame_Shift_Del',
'frameshift_insertion':'Frame_Shift_Ins',
'nonframeshift_deletion':'In_Frame_Del',
'nonframeshift_insertion':'In_Frame_Ins'
#'ncRNA_exonic':'ncRNA_exonic',
#'ncRNA_intronic':'ncRNA_intronic',
#'ncRNA_splicing':'ncRNA_splicing',
#'ncRNA_UTR5':'ncRNA_UTR5'
}

def check_variant_class(func,exonic,variant_type):
	key = (exonic!="." and exonic) or func
	try:
		vc = annovar2maf[key]
	except KeyError:
		vc = "NA/unknown"
	return vc


def process_each(infile, tumor, normal):
	basename = os.path.basename(infile)
	outfile = basename.replace(".vcf",".maf")
	outf = open(outfile,'w')
	print >>outf, "Hugo_Symbol\tChromosome\tStart_Position\tEnd_Position\tVariant_Classification\tVariant_Type\tReference_Allele\tTumor_Seq_Allele1\tTumor_Seq_Allele2\tTumor_Sample_Barcode\tMatched_Norm_Sample_Barcode"
	head_names = ['Gene','chrom','pos','end',]
	indices = []
	normal_name = basename.split("_")[2]
	tumor_name = basename.split("_")[0]
	for line in open(infile,'r'):
		if line.startswith("#CHROM"):
			cols = line[:-1].split("\t")
			i_gene = cols.index("Gene_refGene")
			i_func = cols.index("Func_refGene")
			i_exonic = cols.index("ExonicFunc_refGene")
			i_tumor = cols.index(tumor)
			i_normal = cols.index(normal)
			i_format = cols.index("FORMAT")
		elif line.startswith("#"):
			continue
		else:
			cols = line[:-1].split("\t")
			chrom = cols[0]; pos = cols[1]; ref = cols[3]; alt = cols[4]
			gene = cols[i_gene]; func = cols[i_func]; exonic = cols[i_exonic]
			format = cols[i_format]
			i = format.split(":").index("GT")
			gt_tumor = cols[i_tumor].split(":")[i]
			gt_normal = cols[i_normal].split(":")[i]
			delta = len(ref) - len(alt)
			variant_type = (delta==0 and "SNP") or ((delta<0 and "INS") or "DEL")
			variant_class = check_variant_class(func,exonic,variant_type)
			#SNP
			if delta == 0:
				start = end = pos
				ref_allele = ref
				tumor_allele2 = alt
			#INS
			elif delta < 0:
				start = pos
				end = str(int(start) + len(alt)-len(ref))
				ref_allele = "-"
				tumor_allele2 = alt[len(ref):]
			#DEL
			else:
				start = str(int(pos)+1)
				end = str(int(start) + (len(ref)-len(alt)-1))
				ref_allele = ref[len(alt):]
				tumor_allele2 = "-"
			if gt_tumor == '1/1':
				tumor_allele1 = tumor_allele2
			else:
				tumor_allele1 = ref_allele
                        if variant_class != "NA/unknown":
    		        	output = [gene,chrom,start,end,variant_class,variant_type,ref_allele,tumor_allele1,tumor_allele2,tumor_name,normal_name]
				print >>outf, "\t".join(output)
	#END



def main(argv=sys.argv):
	if len(argv) not in (2,4):
		print """
OBJ
    to convert annotated vcf file to maf format for MutSig
    only the following columns will be saved in maf file
    'Hugo_Symbol Chromosome Start_Position End_Position Variant_Classification 
     Variant_Type Reference_Allele Tumor_Seq_Allele1 Tumor_Seq_Allele2 
     Tumor_Sample_Barcode Matched_Norm_Sample_Barcode'

Usage
    %s [input.list | in.vcf tumor normal]

Attention
1. each line of input.list should be 'sample_name file_path tumor normal'
   where, tumor and normal are indicators in vcf file
2. output will be named as ./{`basename in.vcf`%%vcf}maf
"""%argv[0]
		sys.exit(1)

	
	if len(argv) == 2:
		for line in open(argv[1],'r'):
			cols = line.split()
			process_each(cols[0],cols[1],cols[2])
	elif len(argv) == 4:
		process_each(argv[1],argv[2],argv[3])
	else:
		print "Error: invalid number of arguments"
		sys.exit(1)


if __name__ == '__main__':
	main()

