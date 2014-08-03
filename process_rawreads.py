#!/usr/local/env python

#print __name__

import re
import sys
import os
import optparse
import subprocess

#main()

#################################################
###           Parse command options           ###
#################################################

usage = """Usage: filter_parse.py [options]"""
                        
parser = optparse.OptionParser(usage=usage)
parser.add_option("-t", action="store", type = "string", dest = "threads", help = "threads")
parser.add_option("-s", action="store", type = "string", dest = "sheet", help = "Sample sheet file")
parser.add_option("-p", action="store_true", dest = "paired", help = "paired reads")
parser.add_option("-c", action="store_true", dest = "clean", help = "quality trim reads using Stacks")
parser.add_option("-q", action="store_true", dest = "quality", help = "quality trim reads using Trimmomatic")
parser.add_option("-r", action="store_true", dest = "rescue", help = "rescue barcodes/restriction sites in Stacks")
parser.add_option("-1", action="store", type = "string", dest = "read1", help = "read 1 (single end)")
parser.add_option("-2", action="store", type = "string", dest = "read2", help = "read 2 (paired end)")
parser.add_option("--renz1", action="store", type = "string", dest = "renz1", help = "restriction enzyme 1 (common cutter)")
parser.add_option("--renz2", action="store", type = "string", dest = "renz2", help = "restriction enzyme 2 (rare cutter)")
parser.add_option("--bar_loc", action="store", type = "string", dest = "bar_loc", help = "location of barcode & index")
parser.add_option("-x", action="store", type = "string", dest = "run", help = "processes to run (1-X)")

options, args = parser.parse_args()

#################################################
###           Setup the environment           ###
#################################################

print "Setting up the command environment"
def setup():
### Create output directories ###
	os.system("mkdir clone_filtered")
	os.system("mkdir lead_trimmed")
	os.system("mkdir parsed")
	os.system("mkdir cleaned")

### Split file name and extension ###	
#	if options.paired == True:
#		r1nm, r1ext = os.path.splitext(options.read1)
#		r2nm, r2ext = os.path.splitext(options.read2)
#	else:
#		r1nm, r1ext = os.path.splitext(options.read1)	

#################################################
###             Clone filter reads            ###
#################################################

print "Filtering PCR duplicates"
def clone_filter():
	if options.paired == True:
	#	print "True"
		os.system("clone_filter -1 "+options.read1+" -2 "+options.read2+" -o ./clone_filtered/")
	else:
		os.system("clone_filter -1 "+options.read1+" -2 "+options.read1+" -o ./clone_filtered/")


### Trim first 8bp of each read ###
print "Trimming away leading 8bp unique molecular identifiers"
def lead_trim():
	if options.paired == True:
		r1nm, r1ext = os.path.splitext(options.read1)
		r2nm, r2ext = os.path.splitext(options.read2)
		os.system("fastx_trimmer -Q 33 -f 9 -i ./clone_filtered/"+r1nm+".fil.fq_1 -o ./lead_trimmed/"+r1nm+".1.clone.trim.fastq")
		os.system("fastx_trimmer -Q 33 -f 9 -i ./clone_filtered/"+r2nm+".fil.fq_2 -o ./lead_trimmed/"+r2nm+".2.clone.trim.fastq")
	else:
		r1nm, r1ext = os.path.splitext(options.read1)
		os.system("fastx_trimmer -Q 33 -f 9 -i ./clone_filtered/"+r1nm+".fil.fq_1 -o ./lead_trimmed/"+r1nm+".1.clone.trim.fastq")
			
#################################################
###               Parse samples               ###
#################################################
print "Parsing reads by sample"			
def parse_sample_sheet():
### Parse the sample sheet to create barcodes file ###
    barcodes = open("barcodes.txt", "w")
    for foo in open(options.sheet).read().splitlines():
        bar = foo.split()
#      	print bar[0], bar[1], bar[2], bar[3], bar[4]
        out = bar[3] + "\t" + bar[4] + "\n"
#	print out
        barcodes.write(out)
    barcodes.close()

### process_radtags subroutine ###
def sample_parser():
	if options.rescue == True:
		if options.clean == True:
			print "Quality-trimming reads using Stacks"
			if options.paired == True:
				r1nm, r1ext = os.path.splitext(options.read1)
				r2nm, r2ext = os.path.splitext(options.read2)
				os.system("process_radtags -r -c -q -b barcodes.txt -o ./parsed --inline_index --renz_1 "+options.renz1+" --renz_2 "+options.renz2+" -1 ./lead_trimmed/"+r1nm+".1.clone.trim.fastq -2 ./lead_trimmed/"+r2nm+".2.clone.trim.fastq")
			else:
				r1nm, r1ext = os.path.splitext(options.read1)
				os.system("process_radtags -r -c -q -b barcodes.txt -o ./parsed --inline_null --renz_1 "+options.renz1+" --renz_2 "+options.renz2+" -f ./lead_trimmed/"+r1nm+".1.clone.trim.fastq")
		else:
			print "Quality-trimming reads using Trimmomatic"
			if options.paired == True:
                                r1nm, r1ext = os.path.splitext(options.read1)
                                r2nm, r2ext = os.path.splitext(options.read2)
                                os.system("process_radtags -r -b barcodes.txt -o ./parsed --inline_index --renz_1 "+options.renz1+" --renz_2 "+options.renz2+" -1 ./lead_trimmed/"+r1nm+".1.clone.trim.fastq -2 ./lead_trimmed/"+r2nm+".2.clone.trim.fastq")
                        else:
                                r1nm, r1ext = os.path.splitext(options.read1)
                                os.system("process_radtags -r -b barcodes.txt -o ./parsed --inline_null --renz_1 "+options.renz1+" --renz_2 "+options.renz2+" -f ./lead_trimmed/"+r1nm+".1.clone.trim.fastq")

	else:
                if options.clean == True:
                        print "Quality-trimming reads using Stacks"
                        if options.paired == True:
                                r1nm, r1ext = os.path.splitext(options.read1)
                                r2nm, r2ext = os.path.splitext(options.read2)
                                os.system("process_radtags -c -q -b barcodes.txt -o ./parsed --inline_index --renz_1 "+options.renz1+" --renz_2 "+options.renz2+" -1 ./lead_trimmed/"+r1nm+".1.clone.trim.fastq -2 ./lead_trimmed/"+r2nm+".2.clone.trim.fastq")
                        else:
                                r1nm, r1ext = os.path.splitext(options.read1)
                                os.system("process_radtags -c -q -b barcodes.txt -o ./parsed --inline_null --renz_1 "+options.renz1+" --renz_2 "+options.renz2+" -f ./lead_trimmed/"+r1nm+".1.clone.trim.fastq")
                else:
                        print "Quality-trimming reads using Trimmomatic"
			if options.paired == True:
                                r1nm, r1ext = os.path.splitext(options.read1)
                                r2nm, r2ext = os.path.splitext(options.read2)
                                os.system("process_radtags -b barcodes.txt -o ./parsed --inline_index --renz_1 "+options.renz1+" --renz_2 "+options.renz2+" -1 ./lead_trimmed/"+r1nm+".1.clone.trim.fastq -2 ./lead_trimmed/"+r2nm+".2.clone.trim.fastq")
                        else:
                                r1nm, r1ext = os.path.splitext(options.read1)
                                os.system("process_radtags -b barcodes.txt -o ./parsed --inline_null --renz_1 "+options.renz1+" --renz_2 "+options.renz2+" -f ./lead_trimmed/"+r1nm+".1.clone.trim.fastq")


### file renaming subroutine ###
def sample_rename():
	for foo in open(options.sheet).read().splitlines():
        	bar = foo.split()
		if options.paired == True:
			parsep1_rename = "mv ./parsed/sample_"+bar[3]+"-"+bar[4]+".1.fq ./parsed/"+bar[0]+"_"+bar[3]+"-"+bar[4]+".1.fq"
			parsep2_rename = "mv ./parsed/sample_"+bar[3]+"-"+bar[4]+".2.fq ./parsed/"+bar[0]+"_"+bar[3]+"-"+bar[4]+".2.fq"
			remp1_rename = "mv ./parsed/sample_"+bar[3]+"-"+bar[4]+".rem.1.fq ./parsed/"+bar[0]+"_"+bar[3]+"-"+bar[4]+".rem.1.fq"
			remp2_rename = "mv ./parsed/sample_"+bar[3]+"-"+bar[4]+".rem.2.fq ./parsed/"+bar[0]+"_"+bar[3]+"-"+bar[4]+".rem.2.fq"
			os.system(parsep1_rename)
			os.system(parsep2_rename)
			os.system(remp1_rename)
			os.system(remp2_rename)
#			os.system("mv ./parsed/sample_"+bar[3]+"-"+bar[4]+".1.fq ./parsed/"+bar[0]+"_"+bar[3]+"-"+bar[4]+".1.fq")
#			os.system("mv ./parsed/sample_"+bar[3]+"-"+bar[4]+".2.fq ./parsed/"+bar[0]+"_"+bar[3]+"-"+bar[4]+".2.fq")
#			os.system("mv ./parsed/sample_"+bar[3]+"-"+bar[4]+".rem.1.fq ./parsed/"+bar[0]+"_"+bar[3]+"-"+bar[4]+".rem.1.fq")
#			os.system("mv ./parsed/sample_"+bar[3]+"-"+bar[4]+".rem.2.fq ./parsed/"+bar[0]+"_"+bar[3]+"-"+bar[4]+".rem.2.fq")	
		else:
			parse_single = "mv ./parsed/sample_"+bar[3]+"-"+bar[4]+".fq ./parsed/"+bar[0]+"_"+bar[3]+"-"+bar[4]+".1.fq"
			os.system(parse_single)
#			os.system("mv ./parsed/sample_"+bar[3]+"-"+bar[4]+".fq ./parsed/"+bar[0]+"_"+bar[3]+"-"+bar[4]+".1.fq")


#################################################
###     	Quality-trim samples	      ###
#################################################

def quality_trim():
	if options.quality == True:
		for foo in open(options.sheet).read().splitlines():
        		bar = foo.split()
#			print bar
#			print bar[0], bar[1], bar[2], bar[3], bar[4]
			handle = bar[0]+"_"+bar[3]+"-"+bar[4]
#			print handle
			threads = options.threads
			if options.paired == True:
				PEclean = "$jar/trimmomatic-0.32.jar PE -threads "+threads+" -trimlog ./cleaned/"+handle+".qtrim.log ./parsed/"+handle+".1.fq ./parsed/"+handle+".2.fq ./cleaned/"+handle+".P1.qtrim ./cleaned/"+handle+".S1.qtrim ./cleaned/"+handle+".P2.qtrim ./cleaned/"+handle+".S2.qtrim LEADING:10 TRAILING:10 SLIDINGWINDOW:4:15 MINLEN:36 TOPHRED33"
				os.system(str(PEclean))
			else:
				SEclean = "$jar/trimmomatic-0.32.jar SE -threads "+threads+" -trimlog ./cleaned/"+handle+".qtrim.log ./parsed/"+handle+".1.fq ./cleaned"+handle+".S1.qtrim LEADING:10 TRAILING:10 SLIDINGWINDOW:4:15 MINLEN:36 TOPHRED33"
				os.system(str(SEclean))


#setup()
#clone_filter()
#lead_trim()
#sample_sheet_parse()
#parse_samples()
#sample_rename()
#quality_trim()

#################################################
###		Specify processes	      ###
#################################################	

def main():
	if "1" in options.run:
		setup()
	if "2" in options.run:
		clone_filter()
	if "3" in options.run:
		lead_trim()
	if "4" in options.run:
		parse_sample_sheet()
		sample_parser()
		sample_rename()	
	if "5" in options.run:
		quality_trim()
#	if "6" in option.run
#		finalize()

main()
