#!/usr/bin/env python

##print __name__

import os
import optparse
import re

usage_line = """
genotypes_from_VCF.py

Version 1.1 (16 August, 2015)
License: GNU GPLv2
To report bugs or errors, please contact Daren Card (dcard@uta.edu).
This script is provided as-is, with no support and no guarantee of proper or desirable functioning.

Script that takes VCF file and parses it to produce input files for downstream analysis. Will produce:
	1. An output VCF based on the user defined filtering values
	2. A genotype matrix that is customizable for various downstream programs
	3. A FASTA nucleotide alignment (with IUPAC ambiguities) for phylogenetic analysis (i.e., RAxML)
	4. A FASTA trinary genotype alignment for phylogenetic analysis (i.e., SNAPP)

The script uses a sample sheet to correctly parse the desired samples, which is a tab-delimited text file \
with four columns: (1) BAM input file name, (2) Sample name, (3) Population ID, and (4) Location. \
User can designate various filtering criteria to eliminate unconfident or inappropriate variants. \
User can specify any combination of the three outputs by using the appropriate flag. For the nucleotide \
and trinary alignments, a genotype quality threshold is needed so that unreliable sites can be coded \
as missing data (?). There is also the option to thin the number of SNPs by only taking 1 SNP per 10 kb, \
so as not to violate the assumptions of many models that dictate SNPs should be independent (i.e., not \
linked). The user specifies a naming prefix that will be used for naming the output files created. \
The suffixes for the different file types are as follows:
	1. Output VCF filtered by MAF, missing data, and other options and possibly thinned.
	2. Genotype matrix output customizable for various downstream programs: .genotype
	3. Nucleotide FASTA: .nucl.fasta
	4. Trinary FASTA: .tri.fasta
	5. Log files: .maf<#>.log (needs to be saved if specifying filtered VCF)
	
Dependencies include the latest versions of R, with the package MASS installed, and VCFtools, all \
included in the user's $PATH. This VCF must include the GP, GL, and GQ format/genotype flags. \

python genotypes_from_VCF.py --samplsheet <samplesheet.txt> --vcf <in.vcf> --prefix <out_prefix> \
[--maf <0-3> --miss <0-1> --gq <PHRED_genotype_quality> --qual <PHRED_variant_quality> --thin <#> \
--biallelic --nucl --trinary --genotype --locinfo <T/F> --refalt <T/F> --headers <0-4> --delimit <1/2> \
--filvcf <file.vcf>]
"""


#################################################
###           Parse command options           ###
#################################################

usage = usage_line
                        
parser = optparse.OptionParser(usage = usage)
parser.add_option("--samplesheet", action= "store", type = "string", dest = "sheet", help = "sample sheet containing samples being processed")
parser.add_option("--vcf", action = "store", type = "string", dest = "vcf", help = "VCF input file")
parser.add_option("--prefix", action = "store", dest = "prefix", help = "prefix for output files [out]", default = "out")
parser.add_option("--nucl", action = "store_true", dest = "nucl", help = "create nucleotide FASTA alignment with IUPAC ambiguities for heterozygous sites [FALSE]", default = False)
parser.add_option("--trinary", action = "store_true", dest = "tri", help = "create trinary FASTA alignment with 0, 1, 2 genotype codes [FALSE]", default = False)
parser.add_option("--structure", action = "store_true", dest = "structure", help = "create genotype matrix input for Structure [FALSE]", default = False)
parser.add_option("--genotype", action = "store", dest = "genotype", help = "type of genotype likelihood output: 0 = Option is off (no matrix output), 1 = PHRED, 2 = -Log10, 3 = Standardized, 4 = Genotype Uncertainty [0]", default = "0")
parser.add_option("--gq", action = "store", dest = "gq", help = "threshold genotype PHRED quality score for reporting individual genotype [20]", default = "20")
parser.add_option("--thin", action = "store", dest = "thin", help = "window size to use for thinning in bp (keeps first SNP it finds and ignores others) [10000]", default = "10000")
parser.add_option("--maf", action = "store", dest = "maf", help = "the minor allele frequency range desired: 0 (all MAF), 1 (MAF >= 0.050), 2 (0.010 <= MAF < 0.050), 3 (MAF < 0.050) [1]", default = "1")
parser.add_option("--miss", action = "store", dest = "miss", help = "the proportion of missing data allowed (0 = all missing data, 1 = no missing data) [0.5]", default = "0.5") 
parser.add_option("--qual", action = "store", dest = "qual", help = "threshold variant PHRED quality score for preserving a variant (otherwise filtered out) [20]", default = "20")
parser.add_option("--biallelic", action = "store_true", dest = "biallelic", help = "only keep biallelic variants [TRUE]", default = True)
parser.add_option("--locinfo", action = "store_true", dest = "locinfo", help = "include locus positional information (format = chromosome_position) in genotype matrix [FALSE]", default = False)
parser.add_option("--refalt", action = "store_true", dest = "refalt", help = "include reference and alternative alleles in genotype matrix [FALSE]", default = False)
parser.add_option("--headers", action = "store", dest = "headers", help = "specify which type of header to include in the genotype matrix (comma separated): 0 = none, 1 = matrix dimensions, 2 = sample IDs, 3 = population IDs, 4 = position and reference/alternative headers [1,2,3,4]", default = "1,2,3,4")
parser.add_option("--delimit", action = "store", dest = "delimit", help = "specify which delimiter to use for the genotype matrix: 1 = space, 2 = tab [1]", default = "1")
parser.add_option("--filvcf", action = "store", type = "string", dest = "filvcf", help ="specify a filtered VCF for genotyping (e.g., re-running a script) - bipasses creating new VCF [N/A]", default = "")

options, args = parser.parse_args()


#################################################
###       Filter VCF using user input         ###
#################################################

## Determine user input for MAF and thinning and use it to construct VCFtools, then run command
def vcf_filter():
	## MAF routine
	if options.maf == "0":
		vcf_maf = ""
		print "\n\n***VCF will not be filtered by MAF***\n\n"
	elif options.maf == "1":
		vcf_maf = "--maf 0.0500"
		print "\n\n***Filtering VCF to MAF >= 0.05***\n\n"
	elif options.maf == "2":
		vcf_maf = "--maf 0.0100000 --max-maf 0.0499999"
		print "\n\n***Filtering VCF to 0.01 <= MAF < 0.05***\n\n"
	elif options.maf == "3":
		vcf_maf = "--maf 0 --max-maf 0.0499999"
		print "\n\n***Filtering VCF to MAF < 0.05***\n\n"
	else:
		print "\n\n***Error: a minor allele range needs to be specified!***\n\n"
	
	## Biallelic routine
	if options.biallelic is True:
		biallelic = "--min-alleles 2 --max-alleles 2"
	else:
		biallelic = ""

	## construct genotype quality, MAF, missing data, and biallelic filtering command and run it
	command = "vcftools --vcf "+str(options.vcf)+" --max-non-ref-af 0.99 "+str(vcf_maf)+" --minQ "+options.qual+" --minGQ "+options.gq+" --max-missing "+options.miss+" "+biallelic+" --recode --recode-INFO-all --out "+str(options.prefix)+".maf"+str(options.maf)+".miss"+str(options.miss)
	print "\n\n###Using the following command with VCFtools to produce MAF filtered VCF###\n\n"
	print command
	os.system(command)

	## Thinning routine (if applicable)
	if options.thin is not None:
		vcf_thin = options.thin
		print "\n\n***Thinning to one SNP per "+options.thin+" bp using the following command***\n\n"
		command = "vcftools --vcf "+str(options.prefix)+".maf"+str(options.maf)+".miss"+options.miss+".recode.vcf --thin "+str(vcf_thin)+" --recode --recode-INFO-all --out "+str(options.prefix)+".maf"+str(options.maf)+".miss"+options.miss+".thin"+options.thin
		print command
		os.system(command)
#		os.system("mv "+options.prefix+".thin.recode.vcf "+options.prefix+".maf"+options.maf+".recode.vcf")
	else:
		vcf_thin = ""
		print "\n\n***No thinning will be performed***\n\n"
	print "\n\n###The filtered VCF is named "+options.prefix+".maf"+options.maf+".recode.vcf###\n\n"


#################################################
###      Creating Genotype matrix output      ###
#################################################

## Create input file for Entropy program using sample sheet and VCF
def geno_matrix(PL, GT, filtered_vcf, delimiter):
	## Initialize output file
	genomatrix_out = open(options.prefix+".genomatrix", "w")
	
	## Get matrix dimensions (samples x loci) from VCFtools log
	if "1" in options.headers:
		[samples, loci] = get_vcf_dims()
		genomatrix_out.write(str(samples)+str(delimiter)+str(loci)+str(delimiter)+"1\n")
	
	sample_total = file_len(options.sheet)
	
	## Output line of sample names from second column of sample sheet
	if "2" in options.headers:
		if options.locinfo is True:
			if "4" in options.headers:
				genomatrix_out.write("Marker"+delimiter)
		if options.refalt is True:
			if "4" in options.headers:
				genomatrix_out.write("Ref."+delimiter+"Alt."+delimiter)
		for sline in open(options.sheet, "r"):
			if not sline.strip().startswith("#"):
				bar = sline.rstrip().split("\t")
				if "4" in options.genotype:
					l1out = bar[1]+str(delimiter)
				elif "1" or "2" or "3" in options.genotype:
					l1out = bar[1]+str(delimiter)+bar[1]+str(delimiter)+bar[1]+str(delimiter)
				genomatrix_out.write(l1out)
		genomatrix_out.write("\n")
	
	## Output line of sample populations from third column of sample sheet
	if "3" in options.headers:
		if options.locinfo is True:
			if "4" in options.headers:
				genomatrix_out.write("Marker"+delimiter)
		if options.refalt is True:
			if "4" in options.headers:
				genomatrix_out.write("Ref."+delimiter+"Alt."+delimiter)
		for sline in open(options.sheet, "r"):
			if not sline.strip().startswith("#"):
				bar = sline.rstrip().split("\t")
				if "4" in options.genotype:
					l2out = bar[2]+str(delimiter)
				elif "1" or "2" or "3" in options.genotype:
					l2out = bar[2]+str(delimiter)+bar[2]+str(delimiter)+bar[2]+str(delimiter)
				genomatrix_out.write(l2out)
		genomatrix_out.write("\n")
	
	## Output genotypes for each sample from VCF (begin at column 10)
	for vline in open(filtered_vcf, "r"):
		if not vline.strip().startswith("#"):
			bar = vline.rstrip().split("\t")
			if options.locinfo is True:
				genomatrix_out.write(bar[0]+"_"+bar[1]+str(delimiter))
			if options.refalt is True:
				genomatrix_out.write(bar[3]+str(delimiter)+bar[4]+str(delimiter))
			for sample in range(9, sample_total+9):
				vcfchunks = bar[sample].split(":")
				if vcfchunks[GT] == "./.":
					geno_out = recode_gl(genomatrix_out, "0,0,0", delimiter)
				else:
					geno_out = recode_gl(genomatrix_out, vcfchunks[PL], delimiter)		# recode genotype likelihoods user choice
				genomatrix_out.write(geno_out+str(delimiter))
			genomatrix_out.write("\n")
	
	genomatrix_out.close()
	print "\n\n###The genotype likelihood matrix can be found in "+options.prefix+".genomatrix###\n\n"


#################################################
###    Creating nucleotide alignment fasta    ###
#################################################

## Create fasta alignment of nucleotides based on genotypes from each sample in VCF
def nucl_fasta(GT, GQ, filtered_vcf):
	counter = 0
	## Initalize output file
	nucl_out = open(options.prefix+".nucl.fasta", "w")
	
	## For each individual in sample sheet
	for line in open(options.sheet, "r"):
	    if not line.strip().startswith("#"):

			## Write out fasta header with sample ID and population ID
			nucl_out.write(">"+line.split("\t")[1]+"_"+line.split("\t")[2]+"_"+line.split("\t")[3])
			
			## For each line (locus) in VCF
			for vline in open(filtered_vcf, "r"):
				if not vline.strip().startswith("#"):
					bar = vline.rstrip().split("\t")
					
					## For each individual in VCF
					target = bar[counter + 9]
					vcfchunks = target.split(":")
					
					## Only write genotypes for loci with genotype quality greater than threshold
#					if int(vcfchunks[GQ]) >= int(options.gq):

#					print vcfchunks[GT]						
					if vcfchunks[GT] == "0/0":				# homozygous reference (4th column)
#						print "homozygous ref"
						nucl_out.write(str(bar[3]))
					elif vcfchunks[GT] == "1/1":
#						print "homozygous alt"
						nucl_out.write(str(bar[4]))			# homozygous alternative (5th column)
					elif vcfchunks[GT] == "1/0":
#                                               print "heterozygous"
                                                nucl_out.write(str(get_amb(bar[3], bar[4])))    # heterozygous (use column 4/5 and subroutine to get ambiguity) 
					elif vcfchunks[GT] == "0/1":
#						print "heterozygous"
						nucl_out.write(str(get_amb(bar[3], bar[4])))	# heterozygous (use column 4/5 and subroutine to get ambiguity)
					elif vcfchunks[GT] == """./.""":
#						print "missing"
						nucl_out.write("?")
					else:
						print "Error! Unknown genotype!"

					## Else write missing data (?)
#					else:
#						nucl_out.write("?")
			nucl_out.write("\n")
			counter += 1
	
	nucl_out.close()
	print "\n\n###Nucleotide genotype alignment can be found in "+options.prefix+".nucl.fasta###\n\n"
	

#################################################
###      Creating trinary alignment fasta     ###
#################################################

## Create trinary fasta alignment based on genotypes from VCF... much the same as nucleotide function
def tri_fasta(GT, GQ, filtered_vcf):
	counter = 0
	## Initialize output file
	tri_out = open(options.prefix+".tri.fasta", "w")
	
	## For each individual in sample sheet
	for line in open(options.sheet, "r"):
	    if not line.strip().startswith("#"):
			
			## Write out fasta header with sample ID and population ID
			tri_out.write(">"+line.split()[1]+"_"+line.split()[2]+"_"+line.split()[3]+"\n")
			
			## for each line (locus) in VCF
			for vline in open(filtered_vcf, "r"):
				if not vline.strip().startswith("#"):
					bar = vline.rstrip().split("\t")
					
					## For each individual in VCF
					target = bar[counter + 9]
					vcfchunks = target.split(":")

					## Only write genotypes for loci with genotype quality greater than threshold
#					if int(vcfchunks[GQ]) >= int(options.gq):
						
					if vcfchunks[GT] == "0/0":			# homozygous reference = 0
						tri_out.write("0")
					elif vcfchunks[GT] == "1/1":			# homozygous alternative = 2
						tri_out.write("2")
					elif vcfchunks [GT] == "0/1":			# heterozygous = 1
                                                tri_out.write("1")
					elif vcfchunks [GT] == "1/0":			# heterozygous = 1
						tri_out.write("1")
					elif vcfchunks [GT] == "./.":
						tri_out.write("?")
					else:
						print "Error! Unknown genotype!"

					## Else write missing data (?)
#					else:
#						tri_out.write("?")
			tri_out.write("\n")
			counter += 1
	
	tri_out.close()
	print "\n\n###Trinary genotype alignment can be found in "+options.prefix+".tri.fasta###\n\n"
	

#################################################
###      Creating Structure input matrix      ###
#################################################

## Create genotype matrix suitable for Structure based on genotypes from each sample in VCF
def structure(GT, GQ, filtered_vcf):
	counter = 0
	## Initalize output file
	struct_out = open(options.prefix+".structure", "w")
	
	## For each individual in sample sheet
	for line in open(options.sheet, "r"):
	    if not line.strip().startswith("#"):

			## Write out first column with sample ID
			struct_out.write(line.split("\t")[1]+"\t")
			
			## For each line (locus) in VCF
			for vline in open(filtered_vcf, "r"):
				if not vline.strip().startswith("#"):
					bar = vline.rstrip().split("\t")
					
					## For each individual in VCF
					target = bar[counter + 9]
					vcfchunks = target.split(":")
					
					## Only write genotypes for loci with genotype quality greater than threshold
#					if int(vcfchunks[GQ]) >= int(options.gq):

#					print vcfchunks[GT]						
					if vcfchunks[GT] == "0/0":				# homozygous reference (use column 4 and subroutine to get code)
#						print "homozygous ref"
						struct_out.write(str(str_amb(bar[3]))+"\t")
					elif vcfchunks[GT] == "1/1":
#						print "homozygous alt"
						struct_out.write(str(str_amb(bar[4]))+"\t")			# homozygous alternative (use column 5 and subroutine to get code)
					elif vcfchunks[GT] == "1/0":
#						print "heterozygous"
						struct_out.write(str(str_amb(bar[4]))+"\t")    # heterozygous (use column 4 and subroutine to get code)
					elif vcfchunks[GT] == "0/1":
#						print "heterozygous"
						struct_out.write(str(str_amb(bar[3]))+"\t")	# heterozygous (use column 5 and subroutine to get code)
					elif vcfchunks[GT] == """./.""":
#						print "missing"
						struct_out.write("-9\t")				# Else write missing data
					else:
						print "Error! Unknown genotype!"

			struct_out.write("\n")

			## Write out first column with sample ID
			struct_out.write(line.split("\t")[1]+"\t")
			
			for vline in open(filtered_vcf, "r"):
				if not vline.strip().startswith("#"):
					bar = vline.rstrip().split("\t")
					
					## For each individual in VCF
					target = bar[counter + 9]
					vcfchunks = target.split(":")
					
					## Only write genotypes for loci with genotype quality greater than threshold
#					if int(vcfchunks[GQ]) >= int(options.gq):

#					print vcfchunks[GT]						
					if vcfchunks[GT] == "0/0":				# homozygous reference (use column 4 and subroutine to get code)
#						print "homozygous ref"
						struct_out.write(str(str_amb(bar[3]))+"\t")
					elif vcfchunks[GT] == "1/1":
#						print "homozygous alt"
						struct_out.write(str(str_amb(bar[4]))+"\t")			# homozygous alternative (use column 5 and subroutine to get code)
					elif vcfchunks[GT] == "1/0":
#						print "heterozygous"
						struct_out.write(str(str_amb(bar[3]))+"\t")    # heterozygous (use column 4 and subroutine to get code) 
					elif vcfchunks[GT] == "0/1":
#						print "heterozygous"
						struct_out.write(str(str_amb(bar[4]))+"\t")	# heterozygous (use column 5 and subroutine to get code)
					elif vcfchunks[GT] == """./.""":
#						print "missing"
						struct_out.write("-9\t")				# Else write missing data
					else:
						print "Error! Unknown genotype!"
			
			struct_out.write("\n")
						
			counter += 1
	
	struct_out.close()
	print "\n\n###Nucleotide genotype alignment can be found in "+options.prefix+".structure###\n\n"
	

#################################################
###      Subroutines for above functions      ###
#################################################	
					
## Convert PHRED genotype likelihoods to absolute genotypes (which account for uncertainty)
## Convert to likelihood (for each number alternative alleles): = 10 ^ PHRED/-10
## Standardize the likelihood (for each number alternative alleles): = likelihood/sum(all likelihoods)
## Multiple standardized likelihoods by number of alternative alles: = standardized likelihoods * # alternative alleles
## Sum to produce absolute genotype on 0 (homozygous reference) to 2 (homozygous alternative) scale
def recode_gl(outfile, genochunk, delimiter):
	bar = genochunk.split(",")
	p0 = float(10 ** (int(bar[0])/-10))
	p1 = float(10 ** (int(bar[1])/-10))
	p2 = float(10 ** (int(bar[2])/-10))
	psum = float(p0 + p1 + p2)
	g0 = float(p0/psum)
	g1 = float(p1/psum)
	g2 = float(p2/psum)
	gsum = float(float(g0*0) + float(g1*1) + float(g2*2))
	if options.genotype == "1":
		return bar[0]+delimiter+bar[1]+delimiter+bar[2]
	elif options.genotype == "2":
		return '{:.3f}'.format(p0)+delimiter+'{:.3f}'.format(p1)+delimiter+'{:.3f}'.format(p2)
	elif options.genotype == "3":
		return '{:.3f}'.format(g0)+delimiter+'{:.3f}'.format(g1)+delimiter+'{:.3f}'.format(g2)
	elif options.genotype == "4":
		return '{:.5f}'.format(gsum)
	else:
		print "\n\n***Specify the output genotype format for genotype matrix***\n\n"
		
## Determines which portion of the FORMAT column contains the genotype, PHRED genotype likelihood, and genotype quality
## Uses sample line from end of output VCF
## Splits by tabs to isolate FORMAT column (#9) and then splits by : to separate tags
## Returns the tag number of three values
def get_stat(filtered_VCF):
	os.system("tail -1 "+filtered_VCF+" > sample_VCF_line.txt")
	sample_line = open("sample_VCF_line.txt", "r").readline()
	chunks = sample_line.rstrip().split("\t")
	target = chunks[8]
	foo = target.split(":")
	GT = foo.index("GT")
	PL = foo.index("PL")
	GQ = foo.index("GQ")
	return GT, PL, GQ
	
## Determines the dimensions of the genotype matrix (loci X individuals)
## Calculates the number of uncommented (#) rows (=# loci) and the number of columns - 9 (=# individuals)
def get_vcf_dims():
	out = []
	if options.filvcf is not "":
		file = options.filvcf
	elif options.vcf is not "":
		file = options.vcf
	else:
		print "Error! Either a raw VCF or a filtered VCF must be specified!"
	snps = 0
	for line in open(file, "r"):
		if not line.strip().startswith("#"):
			bar = line.rstrip().split("\t")
			samples = len(bar) - 9
			snps += 1
	out.append(samples)
	out.append(snps)

	return out

## Determine proper ambiguity code at a locus for heterozygous individuals
def get_amb(major, minor):
	if major == "A" and minor == "G":
		return "R"
	elif major == "G" and minor == "A":
		return "R"
	elif major == "C" and minor == "T":
		return "Y"
	elif major == "T" and minor == "C":
		return "Y"
	elif major == "G" and minor == "C":
		return "S"
	elif major == "C" and minor == "G":
		return "S"
	elif major == "A" and minor == "T":
		return "W"
	elif major == "T" and minor == "A":
		return "W"
	elif major == "G" and minor == "T":
		return "K"
	elif major == "T" and minor == "G":
		return "K"
	elif major == "A" and minor == "C":
		return "M"
	elif major == "C" and minor == "A":
		return "M"
	else:
		return "?"

def file_len(fname):
    count = 0
    with open(fname) as f:
        for line in f:
        	if not line.startswith("#"):
				count += 1
    return count
    
## Determine proper a locus for Structure
def str_amb(nucleotide):
	if nucleotide == "A":
		return "1"
	elif nucleotide == "C":
		return "2"
	elif nucleotide == "G":
		return "3"
	elif nucleotide == "T":
		return "4"
	else:
#		print "Error! Unknown genotype!"
		return "-9"


#################################################
###   		   Main Program               ###
#################################################

def main():
	## If previously filtered VCF is specified, use that, otherwise filter based on user input
	if options.filvcf == "":
		print "\n\n***Producing new VCF based on MAF and SNP independence settings***\n\n"
		vcf_filter()
		if options.thin is None:
			filtered_vcf = options.prefix+".maf"+options.maf+".miss"+options.miss+".recode.vcf"
		else:
			filtered_vcf = options.prefix+".maf"+options.maf+".miss"+options.miss+".thin"+options.thin+".recode.vcf"
	else:
		print "\n\n***Working from previously filtered VCF***\n\n"
		filtered_vcf = options.filvcf
	
	## Retrieve location of FORMAT flags so we can extract the values we want
	(GT, PL, GQ) = get_stat(filtered_vcf)
	
	## If user specified genotype likelihood output, give it to them
	if options.genotype is not "0":
		print "\n\n***Creating a genotype likelihood matrix***\n\n"
		if options.delimit == "1":
			geno_matrix(PL, GT, filtered_vcf, " ")
		elif options.delimit == "2":
			geno_matrix(PL, GT, filtered_vcf, "\t")
		else:
			print "\n\n***Specify a delimiter for the genotype matrix!***\n\n"
	else:
		print "\n\n***Not creating a genotype likelihood matrix***\n\n"
	
	## If user specified nucleotide fasta output, give it to them
	if options.nucl is True:
		print "\n\n***Creating nucleotide SNP genotype alignment***\n\n"
		nucl_fasta(GT, GQ, filtered_vcf)
	else:
		print "\n\n***Not creating a nucleotide SNP genotype alignment***\n\n"
	
	## If user specified trinary fasta output, give it to them
	if options.tri is True:
		print "\n\n***Creating trinary SNP genotype alignment***\n\n"
		tri_fasta(GT, GQ, filtered_vcf)
	else:
		print "\n\n***Not creating a trinary SNP genotype alignment***\n\n"

	## If user specified trinary fasta output, give it to them
	if options.structure is True:
		print "\n\n***Creating genotype matrix input for Structure***\n\n"
		structure(GT, GQ, filtered_vcf)
	else:
		print "\n\n***Not creating genotype matrix input for Structure***\n\n"
	
	os.system("rm -f sample_VCF_line.txt")  # clean up
	
	print "\n\n###Command has finished###\n\n"


#################################################
###        	 Call Main Program            ###
#################################################

main()
        	
