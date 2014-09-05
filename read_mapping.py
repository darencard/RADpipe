#!/usr/local/env python

#print __name__

import optparse
import os

usage_line = """
annotate_fasta.py

Version 1.0 (28 August, 2014)
License: GNU GPLv2
To report bugs or errors, please contact Daren Card (dcard@uta.edu).
This script is provided as-is, with no support and no guarantee of proper or desirable functioning.

Script that annotates assembly contigs (fasta format) using a pre-made dictionary that is based \
upon reciprocal best blast and one-way blast results using an well-annotated genome (e.g., Ensembl), which \
indicates homology. The 'make_annotation_dictionary.py' script must be run first to build the annotation \
dictionary. Input is an assembly in fasta format and the output dictionary from the 'make_annotation_dictionary.py' \
script (in json format). Output is an annotated fasta assembly with transcript annotation IDs \
(e.g., Ensembl IDs) indicated in the contig headers. Also outputs the percentage of contigs that were annotated to \
STOUT.

python annotate_fasta.py -d <dictionary> -i <input_fasta> -o <output_fasta>"""

#################################################
###           Parse command options           ###
#################################################

usage = usage_line
                        
parser = optparse.OptionParser(usage=usage)
parser.add_option("--reference", action = "store", type = "string", dest = "reference", help = "the reference genome you will be mapping to")
parser.add_option("--read_dir", action = "store", type = "string", dest = "directory", help = "the directory containing your reads")
parser.add_option("--threads", action = "store", type = "string", dest = "threads", help = "number of threads/cores to use")
parser.add_option("--single-end", action = "store_true", dest = "single", help = "pass this flag if you are working with only single-end data", default = "False")

options, args = parser.parse_args()

def make_PE_dict():

	PE_dict = {}
	
	for root,dirs,files in os.walk(options.directory):
		for file in files:
			if file.endswith(".qtrim"):
				print file
				name = file.split(os.extseq)
				if name[1] == "P1":
				print name
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
		

def make_SE_dict():

	SE_dict = {}
	
	for root,dirs,files in os.walk(options.directory):
		for file in files:
			if file.endswith(".qtrim"):
				print file
				name = file.split(os.extseq)
				if name[1] == "S1":
				print name
					root = name[0]
					print root
					read = name[1]
					print name
					ext = name[2]
					print ext
					key = str(root)+"."+str(read)+"."+str(ext)
					value = str(root)+".S2."+str(ext)
					if key not in PE_dict.keys():
						PE_dict[key] = value
		print PE_dict
		
		
def single():

	single = {}
	
	for root,dirs,files in os.walk(options.directory):
		for file in files:
			if file.endswith(".qtrim"):
				print file
				name = file.split(os.extseq)
				if name[1] == "SE":
				print name
					root = name[0]
					print root
					read = name[1]
					print name
					ext = name[2]
					print ext
					key = str(root)+"."+str(read)+"."+str(ext)
					value = 1
					if key not in PE_dict.keys():
						single[key] = value
		print single


		
def cat_SE():
	for key in SE_dict.keys():
		foo = key.split(".")
		print foo[0]
		file = foo[0]+".SE.qtrim"
		value = SE_dict[key]
		os.system("cat "+key+" "+value+" > "+file)
		
		
				
def index_ref():
	os.system("bwa index "+options.reference)
	
	
	
	
def PE_map():
	for key in PE_dict.keys():
		foo = key.split(".")
		print foo[0]
		file = foo[0]+".PE.sam"
		value = PE_dict[key]
		os.system("bwa mem "+options.reference+" "+key+" "+value+" > "+file)


	
def SE_map():
	for root,dirs,files in os.walk(options.directory):
	for file in files:
		if file.endswith(".SE.qtrim"):
				name = file.split(os.extseq)
				if name[1] == "SE":
					in = str(name[0])+"."+str(name[1])+"."+str(name[2])
					out = name[0]+".SE.sam"
					os.system("bwa mem "+options.reference+" "+in+" > "+out)

def sam2bam():
	for root,dirs,files in os.walk(options.directory):
	for file in files:
		if file.endswith(".sam"):
				name = file.split(os.extseq)
				in = name[0]+"."+name[1]+".sam"
				out = name[0]+"."+name[1]+".bam"
				os.system("samtools view -bS "+in+" > "+out)
				

def bam_process():
	for root,dirs,files in os.walk(options.directory):
	for file in files:
		if file.endswith(".bam"):
				name = file.split(os.extseq)
				sample = name[0]
				PEin = name[0]+".PE.bam"
				SEin = name[0]+".SE.bam"
				Merge_out = name[0]+".merge.bam"
				Sort_out = name[0]+".merge.sort.bam"
				os.system("samtools merge "+Merge_out+" "+PEin+" "+SEin)
				os.system("samtools sort "+Merge_out+" "+Sort_out)
				os.system("samtools index "+Sort_out)
				os.system("rm -f *.sam")

		
def main():
	os.system("mkdir mapping")
	if options.single == True:
		single
		index_ref
		SE_map
		sam2bam
		bam_process
	else:		
		PE_dict
		SE_dict
		cat_SE
		index_ref
		PE_map
		SE_map
		sam2bam
#		bam_process