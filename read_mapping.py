#!/usr/local/env python

#print __name__

import optparse
import os

usage_line = """
read_mapping.py

Version 1.0 (5 September, 2014)
License: GNU GPLv2
To report bugs or errors, please contact Daren Card (dcard@uta.edu).
This script is provided as-is, with no support and no guarantee of proper or desirable functioning.
Note: I acknowledge that there is redundancy in my code such that commands are being run multiple times. However,
given that it works as desired and that I am a Python novice, I decided not to spend more time fixing the code.
If anyone does use this script and decide to improve the code to eliminate this redundancy, I would appreciate it
if an updated copy could be sent back to me so that I can learn how to improve.

Script that takes read files (paired or unpaired) and maps them to a reference sequence (genome) using bwa.
Note that this script only uses the 'mem' bwa mapping algorithm and can't be used with other algorithms without
modification. Input is a reference sequence in fasta format and the directory containing the read files to be
mapped. Naming conventions for read files are as follows:
1. The file extension (normally .fastq) must match that passed to the command with the '--ext' flag.
2. Prior to the extension, there must be an indication of the read type (single or paired) as follows:
	a. P1 and P2 for paired reads, with P1 designated for the single-end reads and P2 designated for the paired-end reads
	b. S1 and S2 for paired reads in which the pairs were broken during quality trimming, with S1 designated
	for the single-end broken reads and S2 designated for the paired-end broken reads.
	c. SE for non-paired, single-end reads
3. A file root with the name of the sample
Example: ID1234_Loc1_ACTTAG-GTACAG.P1.fastq = single-end reads of paired reads of sample ID1234_Loc1_ACTTAG-GTACAG
Note: Reads that do not have this file name formatting will probably not be run correctly.

One can also pass the number of threads that can be used for mapping and any of the bwa mapping flags (as one text string).
The script will create an indexed reference, map all reads from each sample to the reference to create a SAM mapping file,
convert the SAM mapping file to a BAM mapping file, merge BAM mapping files where necessary, and sort and index each sample
BAM file. It will also remove the larger SAM files to save on space unless the '--keep_sams' flag is passed. All mapping
files are outputted to the 'mapping' directory that is created in the working directory by the script.

python read_mapping.py --reference <reference.fasta> --read_dir <directory_of_reads> --ext <file_ext> [--threads <#threads>
--bwa_opts "<options_string>" --keep_sams]							
"""

#################################################
###           Parse command options           ###
#################################################

usage = usage_line
                        
parser = optparse.OptionParser(usage=usage)
parser.add_option("--reference", action = "store", type = "string", dest = "reference", help = "the reference genome you will be mapping to")
parser.add_option("--read_dir", action = "store", type = "string", dest = "directory", help = "the directory containing your reads")
parser.add_option("--threads", action = "store", type = "string", dest = "threads", help = "number of threads/cores to use [1]", default = "1")
parser.add_option("--ext", action = "store", dest = "ext", help = "the file extension of the read files")
parser.add_option("--bwa_opts", action = "store", dest = "bwa", help = "all additional bwa mapping options as a text string, in quotes")
parser.add_option("--keep_sams", action = "store_true", dest = "sams", help = "keep the SAM mapping files", default = "False")

options, args = parser.parse_args()

PE_dict = {}
SE_dict = {}

def make_PE_dict(PE_dict):
	for root,dirs,files in os.walk(options.directory):
		print "making PE_dict"
#	for file in files:
		if files.endswith("."+options.ext):
			name = file.split(os.extsep)
			if name[1] == "P1":
				root = name[0]
				print root
				read = name[1]
				print name
				ext = name[2]
				print ext
				key = str(root)+"."+str(read)+"."+str(ext)
				value = str(root)+".P2."+str(ext)
				if key not in PE_dict.keys():
					PE_dict[key] = value
		print PE_dict
		

def make_SE_dict(SE_dict):
	for root,dirs,files in os.walk(options.directory):
		print "Making SE_dict"
#	for file in files:
		if files.endswith("."+options.ext):
			name = file.split(os.extsep)
			if name[1] == "S1":
				print "Making SE_dict"
				root = name[0]
				read = name[1]
				ext = name[2]
				key = str(root)+"."+str(read)+"."+str(ext)
				value = str(root)+".S2."+str(ext)
				if key not in SE_dict.keys():
					SE_dict[key] = value
				cat_SE(SE_dict)
			if name[1] == "SE":
				print "Making SE_dict"
				root = name[0]
				read = name[1]
				ext = name[2]
				key = str(root)+"."+str(read)+"."+str(ext)
				value = "SE"
				if key not in SE_dict.keys():
					SE_dict[key] = value
	print SE_dict
			
			
def cat_SE(SE_dict):
	print "Concatenating single-end reads"
	for key in SE_dict.keys():
		foo = key.split(".")
		file = foo[0]+".SE."+options.ext
		value = SE_dict[key]
		print "cat ./"+options.directory+"/"+key+" ./"+options.directory+"/"+value+" > ./"+options.directory+"/"+file
		os.system("cat ./"+options.directory+"/"+key+" ./"+options.directory+"/"+value+" > ./"+options.directory+"/"+file)	
		
				
def index_ref():
	os.system("bwa index "+options.reference)	
	
	
def PE_map(PE_dict):
	for key in PE_dict.keys():
		print "Mapping PE reads"
		foo = key.split(".")
		print foo[0]
		file = foo[0]+".PE.sam"
		value = PE_dict[key]
		if options.bwa == None:
			params = ""
		else:
			params = options.bwa
		print "bwa mem -t "+str(options.threads)+" "+str(options.bwa)+" ./"+options.reference+" ./"+options.directory+"/"+key+" ./"+options.directory+"/"+value+" > ./mapping/"+file
		os.system("bwa mem -t "+str(options.threads)+" "+str(options.bwa)+" ./"+options.reference+" ./"+options.directory+"/"+key+" ./"+options.directory+"/"+value+" > ./mapping/"+file)


	
def SE_map(SE_dict):
	for root,dirs,files in os.walk(options.directory):
		print "Mapping SE reads"
#	for file in files:
		if files.endswith(".SE."+options.ext):
			name = file.split(os.extsep)
			if name[1] == "SE":
				input = str(name[0])+"."+str(name[1])+"."+str(name[2])
				output = name[0]+".SE.sam"
				if options.bwa == None:
					params = ""
				else:
					params = options.bwa
				print "bwa mem -t "+str(options.threads)+" "+str(params)+" ./"+options.reference+" ./"+options.directory+"/"+input+" > ./mapping/"+output
				os.system("bwa mem -t "+str(options.threads)+" "+str(params)+" ./"+options.reference+" ./"+options.directory+"/"+input+" > ./mapping/"+output)

def sam2bam():
	for root,dirs,files in os.walk("mapping"):
		print "Converting to BAM"
#	for file in files:
		if files.endswith(".sam"):
			name = file.split(os.extsep)
			input = name[0]+"."+name[1]+".sam"
			output = name[0]+"."+name[1]+".bam"
			print "samtools view -bS ./mapping/"+input+" > ./mapping/"+output
			os.system("samtools view -bS ./mapping/"+input+" > ./mapping/"+output)
				

def PE_bam_process():
	for root,dirs,files in os.walk("mapping"):
		print "Processing PE BAMs"
#	for file in files:
		if files.endswith(".bam"):
			name = file.split(os.extsep)
			sample = name[0]
			PEin = name[0]+".PE.bam"
			SEin = name[0]+".SE.bam"
			Merge_out = name[0]+".merge.bam"
			Sort_out = name[0]+".merge.sort"
			print "samtools merge -f ./mapping/"+Merge_out+" ./mapping/"+PEin+" ./mapping/"+SEin
			os.system("samtools merge -f ./mapping/"+Merge_out+" ./mapping/"+PEin+" ./mapping/"+SEin)
			print "samtools sort ./mapping/"+Merge_out+" ./mapping/"+Sort_out
			os.system("samtools sort ./mapping/"+Merge_out+" ./mapping/"+Sort_out)
			print "samtools index ./mapping/"+Sort_out+".bam"
			os.system("samtools index ./mapping/"+Sort_out+".bam")


def SE_bam_process():
	for root,dirs,files in os.walk("mapping"):
		print "Processing SE BAMs"
#	for file in files:
		if files.endswith(".bam"):
			name = file.split(os.extsep)
			sample = name[0]
			SEin = name[0]+".SE.bam"
			Sort_out = name[0]+".sort"
			print "samtools sort ./mapping/"+SEin+" ./mapping/"+Sort_out
			os.system("samtools sort ./mapping/"+SEin+" ./mapping/"+Sort_out)
			print "samtools index ./mapping/"+Sort_out+".bam"
			os.system("samtools index ./mapping/"+Sort_out+".bam")

		
def main():
#	os.system("mkdir mapping")
#	index_ref()
	for root,dirs,files in os.walk(options.directory):
		print "Running pipeline"
#	for file in files:
		if files.endswith("."+options.ext):
			name = file.split(os.extsep)
			if name[1] == "SE":
				make_SE_dict(SE_dict)
				SE_map(SE_dict)
				sam2bam()
				SE_bam_process()
			else:
				make_PE_dict(PE_dict)
				make_SE_dict(SE_dict)
				PE_map(PE_dict)
				SE_map(SE_dict)
				sam2bam()
				PE_bam_process()
	if options.sams == True:
		print "SAM output will be saved"
	else:
		os.system("rm -f ./mapping/*.sam")
		
main()