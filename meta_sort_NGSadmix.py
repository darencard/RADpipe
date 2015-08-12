#!/usr/bin/env python

##print __name__

import optparse
import os
import subprocess

usage_line = """
read_mapping.py

Version 1.0 (11 June, 2015)
License: GNU GPLv2
To report bugs or errors, please contact Daren Card (dcard@uta.edu).
This script is provided as-is, with no support and no guarantee of proper or desirable functioning.

This script takes raw output from NGSadmix, adds metadata from the user sample sheet, and orders the \
rows based on either a designated column or a user-inputted list. User must specify the metadata \
(project samplesheet) and the admixture proportions from NGSadmix (.qopt). The user can sort by a column \
(1-4) alphanumerically (which can be reversed using the --rev flag). The user can also input a list of \
sample IDs (corresponding to column 2 in the sample sheet) and the data will be sorted based on the order \
in this list. The user then specifies a prefix for the two output files, which are a full file with \
metadata and admixture data put together and a similar file with these values sorted appropriately. \

python meta_sort_NGSadmix.py --metadata <samplesheet.txt> --admixture <project.qopt> --prefix <output_prefix> \
[--col <#> --rev --user <order_list.txt>]
"""

#################################################
###           Parse command options           ###
#################################################

usage = usage_line
                        
parser = optparse.OptionParser(usage=usage)
parser.add_option("--metadata", action = "store", type = "string", dest = "metadata", help = "file containing sample metadata (e.g., project samplesheet")
parser.add_option("--admixture", action = "store", type = "string", dest = "admixture", help = "file containing admixture proportions")
parser.add_option("--col", action = "store", type = "string", dest = "col", help = "column of metadata to use for sorting (1, 2, ..., N")
parser.add_option("--user", action = "store", type = "string", dest = "user", help = "file containing sample order to sort by based on sample names (column 2 of metadata), one sample per line")
parser.add_option("--prefix", action = "store", type = "string", dest = "prefix", help = "prefix for all output files produced by script [meta_sort.out]", default = "meta_sort.out")
parser.add_option("--rev", action = "store_true", dest = "rev", help = "reverse the column sorting (reverse alphanumeric order) [FALSE]", default = False)

options, args = parser.parse_args()


#################################################
###   Concatenate Metadata & Admixture Dat    ###
#################################################

def cat_meta_admix():
	# get number of columns in admixture data, so we can make proper headings
	with open(options.admixture) as indata:
		cols = indata.readline().split()
		num_cols = len(cols)
	
	# create proper headings for population numbers (1...N) and stitch together header line
	head_line = ""
	for i in range(0, num_cols):
		head_line += "Pop"+str(i+1)+"\t"
	
	# replace the spaces from NGS admix with tabs, so that all spacing is consistent
	with open(options.admixture, "r") as admix_data:
		plaintext = admix_data.read()
		plaintext = plaintext.replace(' ', '\t')
	with open(options.admixture+".tsv", "w") as out_admix_data:
		out_admix_data.write(plaintext)
	
	# call unix command that concatenates created header with tab-delimited admixture proportions
	command1 = "echo \""+head_line+"\" | cat - "+options.admixture+".tsv > temp && mv temp "+options.admixture+".reformat.tsv"
	print command1
	os.system(command1)
	
	# call unix command to paste together the metadata columns and the admixture data columns
	command2 = "paste "+options.metadata+" "+options.admixture+".reformat.tsv > "+options.prefix+".meta.tsv"
	print command2
	os.system(command2)


#################################################
###        	   Sort data by column            ###
#################################################

def sort_col():
	# store column from user
	col = options.col
	
	# subroutine to determine if data needs to be sorted in reverse
	if options.rev is True:
		rev = "-r"
	else:
		rev = ""
		
	# unix command to read header into new file	
	command1 = "head -1 "+options.prefix+".meta.tsv > "+options.prefix+".meta.sort.tsv"
	os.system(command1)
	
	# unix command to sort everything below header line by designated column and then append the sorted data to the file with waiting header
	command2 = "tail -n +2 "+options.prefix+".meta.tsv | sort -f "+rev+" -k "+col+" - >> "+options.prefix+".meta.sort.tsv"
	print command2
	os.system(command2)


#################################################
###    	       Sort data by user list         ###
#################################################

def sort_user():
	# create ordered list that reflects the order from the user's input file
	sample_list = []
	for line in open(options.user, "r"):
		if not line.strip().startswith("#"):
			foo = line.rstrip().split("\t")
			sample_list.append(foo[0])
			
	# initialize output file
	outfile = open(options.prefix+".meta.sort.tsv", "w")
	
	# read full input created under cat_meta_admix()
	full_admix = open(options.prefix+".meta.tsv", "r")
	
	# read first line from full data (see above) and output as header for new file
	head = full_admix.readline()
	outfile.write(head)
	
	# iterate through ordered list of samples and when we get a match to column 2 of full data, write that entire line to output file
	for sample in sample_list:
		for line in open(options.prefix+".meta.tsv", "r"):
			if not line.strip().startswith("#"):
				bar = line.rstrip().split("\t")
				if bar[1] == sample:
					outfile.write(line)
	outfile.close()



#################################################
###            	   Full Program               ###
#################################################
		
def main():
	if options.metadata is None:					# User didn't specify metadata (samplesheet)
		print "\n***Error: specify file containing metadata (i.e., working sample sheet)!***\n"
	if options.admixture is None:					# User didn't specify admixture data
		print "\n***Error: specify file containing admixture proportions!***\n"
	cat_meta_admix()								# concatenate metadata and admixture data
	if options.user is None:						# if user-specified order isn't present
		if options.col is None:						# check to make sure there is a column to sort by
			print "\n***Error: must specify either sorting by a column or by an inputted list!***\n"
		else:										
			sort_col()								# call column sorting function
	if options.col is None:							# if column for sorting isn't specified
		if options.user is None:					# make sure there is a user list to sort with
			print "\n***Error: must specify either sorting by a column or by an inputted list!***\n"
		else:
			sort_user()								# call user sorting function


#################################################
###              Run Full Program             ###
#################################################
		
main()