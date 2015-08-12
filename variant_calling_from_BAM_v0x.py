#!/usr/bin/env python

##print __name__

import os
import optparse
import re

usage_line = """
Variant_calling_from_BAM.py

Version 1.0 (20 May, 2015)
License: GNU GPLv2
To report bugs or errors, please contact Daren Card (dcard@uta.edu).
This script is provided as-is, with no support and no guarantee of proper or desirable functioning.

This script produces both an mpileup BCF file and a variants VCF file using a set of BAM mapping files \
as input. The user simple specifies the directory containing the BAM mapping files and a tab-delimited \
sample sheet specifying which samples are to be included in the mpileup and variants files. This sample \
sheet enables the user to subset data, if desired, and orders the samples in the mpileup and variants files \
so that downstream scripts will work properly. The user must also specify the faidx-indexed reference \
(FASTA) and can specify a prefix for the output files. The user also must specify the path or the executable \
name for SAMtools and BCFtools v.0.X.XX (tested with v.0.1.19). Some further capability is also provided:
	1. The option to exclude INDELS from the variants file
	2. The option to set the amount of missing data (i.e., samples without data at a locus)
	3. The option to set the p-value for the variant calling model (see SAMtools documentation)
	4. The option to run either the mpileup or the variant calling steps individually
	5. If running the mpileup and variant calling separately, a mpileup file can be specified from a \
previous run
	
This script produces two output files:
	1. Output mpileup of all loci in BCF format: <prefix>.mpileup.bcf
	2. Output variants in VCF format: <prefix>.variants.vcf
	
Dependencies include the version 0 flavors of both SAMtools and BCFtools. The user may specify the \
executable name (if the program is in the path) or may feed the entire path. This is necessary due \
to the recent release of SAMtools/BCFtools version 1, which is substantially different, and allows \
the user to specify an executable with a name other than 'samtools' or 'bcftools', in case he or \
she has both version 0 and 1 on the computer.

python Variant_calling_from_BAM.py --samplsheet <samplesheet.txt> --dir <dir_with_BAMs> --prefix <out_prefix> \
--samtools <path_to_samtools> --bcftools <path_to_bcftools> --ref <path_to_reference> [--indels --miss <0.XX> \
--pval <0.XX> --mpileup <mpileup.bcf> --exe <1,2>]
"""


#################################################
###           Parse command options           ###
#################################################

usage = usage_line
                        
parser = optparse.OptionParser(usage = usage)
parser.add_option("--samplesheet", action= "store", type = "string", dest = "sheet", help = "sample sheet containing samples being processed")
parser.add_option("--dir", action = "store", type = "string", dest = "dir", help = "directory containing sorted, indexed BAM mapping files")
parser.add_option("--prefix", action = "store", dest = "prefix", help = "prefix for output files [out]", default = "out")
parser.add_option("--samtools", action = "store", dest = "samtools", help = "path to SAMtools v0.X.XX, or just provide the program if it is in your path (this option is needed because of the multiple versions of SAMtools that exist)")
parser.add_option("--bcftools", action = "store", dest = "bcftools", help = "path to BCFtools v0.X.XX, or just provide the program if it is in your path (this option is needed because of the multiple versions of BCFtools that exist)")
parser.add_option("--ref", action = "store", dest = "ref", help = "path to the faidx-indexed reference file in FASTA format")
parser.add_option("--indels", action = "store_true", dest = "indels", help = "do not perform INDEL calling [TRUE]", default = True)
parser.add_option("--miss", action = "store", dest = "miss", help = "only include variants where fraction of samples covered by reads is above given FLOAT threshold (0-1) [0.50]", default = "0.50")
parser.add_option("--pval", action = "store", dest = "pval", help = "p-value threshold for variant calling model (i.e., if P(ref|data)<FLOAT) [0.05]", default = "0.05")
parser.add_option("--mpileup", action = "store", dest = "mpileup", help = "a mpileup file to use for variant calling (i.e., if it was already generated previously) [NA]")
parser.add_option("--exe", action = "store", dest = "exe", help = "processes to run, separated by comma (1 or 2 or 1,2): 1 = generate mpileup; 2 = call variants", default = "1,2")

options, args = parser.parse_args()


#################################################
### Create list of sample BAM files for input ###
#################################################

def make_sample_list():
	## Initialize empty sample list
	sample_list = ""
	
	## For non-commented line in sample sheet, split and focus on first column (BAM files)
	for line in open(options.sheet, "r"):
		if not line.strip().startswith("#"):
			foo = line.rstrip().split("\t")
			
	## Paste together input directory, "/", and BAM file for each line and then	append this to already existing sample list to iteratively build sample list	
			bar = options.dir+"/"+foo[0]+" "
			sample_list = sample_list + bar
	return sample_list


#################################################
###        		   Main Program               ###
#################################################

def main():
	## Make directory for BCF/VCF output
	os.system("mkdir vcf")
	
	## Create and gather sample list for command
	sample_list = make_sample_list()
	
	## Set whether to ignore INDELs
	if options.indels is True:
		indels = "-I "
	else:
		indels = ""
	
	## If user wanted to create mpileup, create command and then run it
	if "1" in options.exe:
		mpileup = options.samtools+" mpileup -P ILLUMINA -u -g -f "+options.ref+" "+sample_list+" > ./vcf/"+options.prefix+".mpileup.bcf"
		print mpileup
		os.system(mpileup)
	
	## If user wanted to create variants VCF, create command and then run it
	if "2" in options.exe:
                if options.mpileup is not None:
                        mpilein = options.mpileup
                else:
                        mpilein = "./vcf/"+options.prefix+".mpileup.bcf"
		variants = options.bcftools+" view -N -c -e -g -v -P full -t 0.001 "+indels+"-d "+options.miss+" -p "+options.pval+" "+mpilein+".mpileup.bcf > ./vcf/"+options.prefix+".variants.d"+options.miss+".p"+options.pval+".vcf"
		print variants
		os.system(variants)
	

#################################################
###        	Call Main Program             ###
#################################################

main()
        	
