#!/usr/local/env python

##print __name__

import os
import optparse
import re

usage_line = """
genotypes_from_VCF.py

Version 1.0 (12 May, 2015)
License: GNU GPLv2
To report bugs or errors, please contact Daren Card (dcard@uta.edu).
This script is provided as-is, with no support and no guarantee of proper or desirable functioning.

Script that takes VCF file and parses it to produce input files for downstream analysis. Will produce:
	1. An output VCF based on the user defined MAF ranges
	2. An input file for Entropy (Gompert & Buerkle), as well as appropriate starting points for the MCMC
	3. A FASTA nucleotide alignment (with IUPAC ambiguities) for phylogenetic analysis (i.e., RAxML)
	4. A FASTA trinary genotype alignment for phylogenetic analysis (i.e., SNAPP)

The script uses a sample sheet to correctly parse the desired samples, which is a tab-delimited text file \
with four columns: (1) BAM input file name, (2) Sample name, (3) Population ID, and (4) Location. \
User must also designate the MAF range desired for the run (this represents minimal output). \
User can specify any combination of the three outputs by using the appropriate flag. In the case of \
Entropy, starting and ended K values should be provided for creating the MCMC chain starting points. \
For the nucleotide and trinary alignments, a genotype quality threshold is needed so that unreliable \
sites can be coded as missing data (?). There is also the option to thin the number of SNPs by \
only taking 1 SNP per 10 kb, so as not to violate the assumptions of many models that dictate SNPs \
should be independent (i.e., not linked). The user specifies a naming prefix that will be used for \
naming the output files created. The suffixes for the different file types are as follows:
	1. Output VCF filtered by MAF and possibly thinned: .maf<#>.recode.vcf
	2. Entropy input: .entropy
	3. Entropy MCMC initialization (for each K): .entropy.startK<#>
	4. IN DEVELOPMENT: Various R plots from K-means clustering: .kmean.K<#>.plot
	5. IN DEVELOPMENT: Various R plots from discriminant analysis: .dapc.plot
	6. Nucleotide FASTA: .nucl.fasta
	7. Trinary FASTA: .tri.fasta
	8. Log file: .maf<#>.log (needs to be saved if specifying filtered VCF)
	
Dependencies include the latest versions of R, with the package MASS installed, and VCFtools, all \
included in the user's $PATH. The user should provide an input VCF that has already been filtered \
based upon factors like base quality, mapping quality, and missing data. This VCF must include the \
GP, GL, and GQ format/genotype flags.

python genotypes_from_VCF.py --samplsheet <samplesheet.txt> --vcf <in.vcf> --prefix <out_prefix> \
--maf <1-3> [--entropy --startK <#> --endK <#> --nucl --trinary --gq <PHRED_genotype_quality> --ind]
"""


#################################################
###           Parse command options           ###
#################################################

usage = usage_line
                        
parser = optparse.OptionParser(usage = usage)
parser.add_option("--samplesheet", action= "store", type = "string", dest = "sheet", help = "sample sheet containing samples being processed")
parser.add_option("--vcf", action = "store", type = "string", dest = "vcf", help = "VCF input file")
parser.add_option("--prefix", action = "store", dest = "prefix", help = "prefix for output files [out]", default = "out")
parser.add_option("--entropy", action = "store_true", dest = "entropy", help = "create entropy input file and proper starting points [TRUE]", default = "True")
parser.add_option("--startK", action = "store", dest = "sK", help = "starting (lower) K [1]", default = "1")
parser.add_option("--endK", action = "store", dest = "eK", help = "ending (higher) K [5]", default = "2")
parser.add_option("--nucl", action = "store_true", dest = "nucl", help = "create nucleotide FASTA alignment with IUPAC ambiguities for heterozygous sites [TRUE]", default = "True")
parser.add_option("--trinary", action = "store_true", dest = "tri", help = "create trinary FASTA alignment with 0, 1, 2 genotype codes [TRUE]", default = "True")
parser.add_option("--gq", action = "store", dest = "gq", help = "threshold genotype quality for reporting individual genotype (otherwise coded as missing - ?) [20]", default = "20")
parser.add_option("--ind", action = "store_true", dest = "ind", help = "thin dataset to only include 1 SNP per 10 kb, so as to reduce chance of linked SNPs [TRUE]", default = "True")
parser.add_option("--thin", action = "store", dest = "thin", help = "window size to use for thinning in bp (keeps first SNP it finds and ignores others) [10000]", default = "10000")
parser.add_option("--maf", action = "store", dest = "maf", help = "the minor allele frequency range desired: 0 (all MAF), 1 (MAF >= 0.050), 2 (0.010 <= MAF < 0.050), 3 (MAF < 0.050) [1]", default = "1") 
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
		vcf_maf = "--maf 0.0000001 --max-maf 0.0499999"
		print "\n\n***Filtering VCF to MAF < 0.05***\n\n"
	else:
		print "\n\n***Error: a minor allele range needs to be specified!***\n\n"
	
	## construct MAF filtering command and run it
	command = "vcftools --vcf "+options.vcf+" "+vcf_maf+" --recode --recode-INFO-all --out "+options.prefix+".maf"+options.maf
	print "\n\n###Using the following command with VCFtools to produce filtered VCF###\n\n"
	print command
	os.system(command)
	print "\n\n###The filtered VCF is named "+options.prefix+".maf"+options.maf+".recode.vcf###\n\n"

	## Thinning routine (if applicable)
	if options.ind is True:
		vcf_thin = options.thin
		print "\n\n***Thinning to one SNP per 10 kb using the following command***\n\n"
		command = "vcftools --vcf "+options.prefix+".maf"+options.maf+".recode.vcf"+" --thin "+vcf_thin+" --recode --recode-INFO-all --out "+options.prefix+".thin"
		print command
		os.system(command)
		os.system("mv "+options.prefix+".thin "+options.prefix+".maf"+options.maf+".recode.vcf")
	else:
		vcf_thin = ""
		print "\n\n***No thinning will be performed***\n\n"

#################################################
###          Creating Entropy output          ###
#################################################

## Create input file for Entropy program using sample sheet and VCF
def entropy(PL, filtered_vcf):
	## Initialize output file
	entropy_out = open(options.prefix+".entropy", "w")
	
	## Get matrix dimensions (samples x loci) from VCFtools log
	[samples, loci] = get_vcf_dims()
	entropy_out.write(samples+" "+loci+" "+"1\n")
	
	sample_total = 0
	
	## Output line of sample names from second column of sample sheet
	for sline in open(options.sheet, "r"):
		if not sline.strip().startswith("#"):
			bar = sline.rstrip().split("\t")
			l1out = bar[1]+" "
			entropy_out.write(l1out)
			sample_total += 1
	
	## Output line of sample populations from third column of sample sheet
	entropy_out.write("\n")
	for sline in open(options.sheet, "r"):
		if not sline.strip().startswith("#"):
			bar = sline.rstrip().split("\t")
			l2out = bar[2]+" "
			entropy_out.write(l2out)
	
	## Output genotypes for each sample from VCF (begin at column 10)
	entropy_out.write("\n")
	for vline in open(filtered_vcf, "r"):
		if not vline.strip().startswith("#"):
			bar = vline.rstrip().split("\t")
			for sample in range(9, sample_total+9):
				vcfchunks = bar[sample].split(":")
				geno_out = recode_gl(entropy_out, vcfchunks[PL])		# recode genotype likelihoods to genotypes
				entropy_out.write('{:.5f}'.format(geno_out)+" ")		# to go 5 decimal points
			entropy_out.write("\n")
	
	entropy_out.close()
	print "\n\n###Entropy results can be found in "+options.prefix+".entropy###\n\n"


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
					if int(vcfchunks[GQ]) >= int(options.gq):
						
						if vcfchunks[GT] == "0/0":					# homozygous reference (4th column)
							nucl_out.write(str(bar[3]))
						elif vcfchunks[GT] == "1/1":
							nucl_out.write(str(bar[4]))				# homozygous alternative (5th column)
						else:
							nucl_out.write(str(get_amb(bar[3], bar[4])))	# heterozygous (use column 4/5 and subroutine to get ambiguity)
					
					## Else write missing data (?)
					else:
						nucl_out.write("?")
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

					## Knly write genotypes for loci with genotype quality greater than threshold
					if int(vcfchunks[GQ]) >= int(options.gq):
						
						if vcfchunks[GT] == "0/0":				# homozygous reference = 0
							tri_out.write("0")
						elif vcfchunks[GT] == "1/1":			# homozygous alternative = 2
							tri_out.write("2")
						else:									# heterozygous = 1
							tri_out.write("1")
					
					## Else write missing data (?)
					else:
						tri_out.write("?")
			tri_out.write("\n")
			counter += 1
	
	tri_out.close()
	print "\n\n###Trinary genotype alignment can be found in "+options.prefix+".tri.fasta###\n\n"
	

#################################################
###      Subroutines for above functions      ###
#################################################	
					
## Convert PHRED genotype likelihoods to absolute genotypes (which account for uncertainty)
## Convert to likelihood (for each number alternative alleles): = 10 ^ PHRED/-10
## Standardize the likelihood (for each number alternative alleles): = likelihood/sum(all likelihoods)
## Multiple standardized likelihoods by number of alternative alles: = standardized likelihoods * # alternative alleles
## Sum to produce absolute genotype on 0 (homozygous reference) to 2 (homozygous alternative) scale
def recode_gl(outfile, genochunk):
	out = outfile
	foo = genochunk
	bar = foo.split(",")
	p0 = float(10 ** (int(bar[0])/-10))
	p1 = float(10 ** (int(bar[1])/-10))
	p2 = float(10 ** (int(bar[2])/-10))
	psum = float(p0 + p1 + p2)
	g0 = float((p0/psum)*0)
	g1 = float((p1/psum)*1)
	g2 = float((p2/psum)*2)
	gsum = float(round((g0 + g1 + g2),4))
	return gsum
		
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
	
## Determines the dimensions of the genotype matrix (individuals x loci)
## Uses a regular expression search of the VCFtools log
def get_vcf_dims():
	out = []
	for line in open(options.prefix+".maf"+options.maf+".log", "r"):
		if line.strip().startswith("After"):
			query = re.compile('kept (.*?) out')
			match = query.search(line)
			if match:
				element = match.group(1)
				out.append(element)
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


#################################################
###        		   Main Program               ###
#################################################

def main():
	## If previously filtered VCF is specified, use that, otherwise filter based on user input
	if options.filvcf == "":
		vcf_filter()
		filtered_vcf = options.prefix+".maf"+options.maf+".recode.vcf"
		print "\n\n***Producing new VCF based on MAF and SNP independence settings***\n\n"
	else:
		filtered_vcf = options.filvcf
		print "\n\n***Working from previously filtered VCF***\n\n"
	
	## Retrieve location of FORMAT flags so we can extract the values we want
	(GT, PL, GQ) = get_stat(filtered_vcf)
	
	## If user specified Entropy output, give it to them
	if options.entropy == True:
		entropy(PL, filtered_vcf)
		print "\n\n***Creating genotype input for Entropy***\n\n"
	else:
		print "\n\n***Not creating any input for Entropy***\n\n"
	
	## If user specified nucleotide fasta output, give it to them
	if options.nucl == True:
		nucl_fasta(GT, GQ, filtered_vcf)
		print "\n\n***Creating nucleotide SNP genotype alignment***\n\n"
	else:
		print "\n\n***Not creating a nucleotide SNP genotype alignment***\n\n"
	
	## If user specified trinary fasta output, give it to them
	if options.tri == True:
		tri_fasta(GT, GQ, filtered_vcf)
		print "\n\n***Creating trinary SNP genotype alignment***\n\n"
	else:
		print "\n\n***Not creating a trinary SNP genotype alignment***\n\n"
	
	os.system("rm -f sample_VCF_line.txt")  # clean up
	
	print "\n\n###Command has finished###\n\n"


#################################################
###        	  Call Main Program               ###
#################################################

main()
        	
