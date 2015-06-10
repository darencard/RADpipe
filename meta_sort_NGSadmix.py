#!/usr/bin/env python

##print __name__

import optparse
import os
import subprocess

usage_line = """
read_mapping.py

Version 1.0 (8 September, 2014)
License: GNU GPLv2
To report bugs or errors, please contact Daren Card (dcard@uta.edu).
This script is provided as-is, with no support and no guarantee of proper or desirable functioning.


"""

#################################################
###           Parse command options           ###
#################################################

usage = usage_line
                        
parser = optparse.OptionParser(usage=usage)
parser.add_option("--metadata", action = "store", type = "string", dest = "metadata", help = "the reference genome you will be mapping to")
parser.add_option("--admixture", action = "store", type = "string", dest = "admixture", help = "the directory containing your reads")
parser.add_option("--col", action = "store", type = "string", dest = "col", help = "number of threads/cores to use [1]", default = "1")
parser.add_option("--user", action = "store", type = "string", dest = "user", help = "the file extension of the read files")
parser.add_option("--prefix", action = "store", type = "string", dest = "prefix", help = "number of threads/cores to use [1]", default = "1")
parser.add_option("--rev", action = "store", type = "string", dest = "rev", help = "number of threads/cores to use [1]", default = "1")

options, args = parser.parse_args()


#################################################
###   Concatenate Metadata & Admixture Dat    ###
#################################################

def cat_meta_admix():
	with open(options.admixture) as indata:
		cols = indata.readline().split()
		num_cols = len(cols)
	head_line = ""
	for i in range(0, num_cols):
		head_line += "Pop"+str(i+1)+"\t"
	with open(options.admixture, "r") as admix_data:
		plaintext = admix_data.read()
		plaintext = plaintext.replace(' ', '\t')
	with open(options.admixture+".tsv", "w") as out_admix_data:
		out_admix_data.write(plaintext)
	command1 = "echo \""+head_line+"\" | cat - "+options.admixture+".tsv > temp && mv temp "+options.admixture+".tsv"
	print command1
	os.system(command1)
	command2 = "paste "+options.metadata+" "+options.admixture+".tsv > "+options.prefix+".meta.tsv"
	print command2
	os.system(command2)


#################################################
###        	   Sort data by column            ###
#################################################

def sort_col():
	col = options.col
	if options.rev is True:
		rev = "-r"
	else:
		rev = ""
	command1 = "head -1 "+options.prefix+".meta.tsv > "+options.prefix+".meta.sort.tsv"
	command2 = "tail -n +2 "+options.prefix+".meta.tsv | sort -f "+rev+" -k "+col+" - >> "+options.prefix+".meta.sort.tsv"
	print command2
	os.system(command1)
	os.system(command2)


#################################################
###    	       Sort data by user list         ###
#################################################

def sort_user():
	sample_list = []
	for line in open(options.user, "r"):
		if not line.strip().startswith("#"):
			foo = line.rstrip().split("\t")
			sample_list.append(foo[0])
	outfile = open(options.prefix+".meta.sort.tsv", "w")
	for sample in sample_list:
		for line in open(options.prefix+".meta.tsv", "r"):
			if line.strip().startswith("#"):
				outfile.write(line)
			elif not line.strip().startswith("#"):
				bar = line.rstrip().split("\t")
				if bar[1] == sample:
					outfile.write(line)
	outfile.close()


#################################################
###            	   Full Program               ###
#################################################
		
def main():
	if options.metadata is None:
		print "\n***Error: specify file containing metadata (i.e., working sample sheet)!***\n"
	if options.admixture is None:
		print "\n***Error: specify file containing admixture proportions!***\n"
	cat_meta_admix()
	if options.user is None:
		if options.col is None:
			print "\n***Error: must specify either sorting by a column or by an inputted list!***\n"
		else:
			sort_col()
	if options.col is None:
		if options.user is None:
			print "\n***Error: must specify either sorting by a column or by an inputted list!***\n"
		else:
			sort_user()


#################################################
###              Run Full Program             ###
#################################################
		
main()