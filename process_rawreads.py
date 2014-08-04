#!/usr/local/env python

#print __name__

import re
import sys
import os
import optparse
import subprocess


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

def setup():
	print "\n***Setting up the command environment***\n"
### Create output directories ###
	os.system("mkdir clone_filtered")
	os.system("mkdir lead_trimmed")
	os.system("mkdir parsed")
	os.system("mkdir cleaned")


#################################################
###             Clone filter reads            ###
#################################################

def clone_filter():
	print "\n***Filtering PCR duplicates***\n"
	if options.paired == True:
	#	print "True"
		os.system("clone_filter -1 "+options.read1+" -2 "+options.read2+" -o ./clone_filtered/")
	else:
		os.system("clone_filter -1 "+options.read1+" -2 "+options.read1+" -o ./clone_filtered/")


#################################################
###             Trim leading 8bp UMI          ###
#################################################

def lead_trim():
	print "\n***Trimming away leading 8bp unique molecular identifiers***\n"
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

def parse_sample_sheet():
	print "\n***Parsing reads by sample***\n"

### Parse the sample sheet to create barcodes file ###
	barcodes = open("barcodes.txt", "w")
	for foo in open(options.sheet).read().splitlines():
        	bar = foo.split()
		if options.paired == True:
#			print bar[0], bar[1], bar[2], bar[3], bar[4]
        		out = bar[3] + "\t" + bar[4] + "\n"
#			print out
        		barcodes.write(out)
		else:
#               	print bar[0], bar[1], bar[2], bar[3]
                	out = bar[3] + "\n"
#               	print out
                	barcodes.write(out)
    	barcodes.close()

### process_radtags subroutine ###
def sample_parser():
	if options.rescue == True:
		if options.clean == True:
			alert = open("./cleaned/ATTENTION", "w")
			line = "You elected to quality-trim your reads using Stacks. This trimming was done simultaneously with parsing. See the 'parsed' folder for your trimmed reads."
			alert.write(line)
			alert.close()
			if options.paired == True:
				r1nm, r1ext = os.path.splitext(options.read1)
				r2nm, r2ext = os.path.splitext(options.read2)
				os.system("process_radtags -r -c -q -b barcodes.txt -o ./parsed --inline_index --renz_1 "+options.renz1+" --renz_2 "+options.renz2+" -1 ./lead_trimmed/"+r1nm+".1.clone.trim.fastq -2 ./lead_trimmed/"+r2nm+".2.clone.trim.fastq")
			else:
				r1nm, r1ext = os.path.splitext(options.read1)
				os.system("process_radtags -r -c -q -b barcodes.txt -o ./parsed --inline_null --renz_1 "+options.renz1+" --renz_2 "+options.renz2+" -f ./lead_trimmed/"+r1nm+".1.clone.trim.fastq")
                        print "\n***Quality-trimming reads using Stacks***\n"
		else:
			if options.paired == True:
                                r1nm, r1ext = os.path.splitext(options.read1)
                                r2nm, r2ext = os.path.splitext(options.read2)
                                os.system("process_radtags -r -b barcodes.txt -o ./parsed --inline_index --renz_1 "+options.renz1+" --renz_2 "+options.renz2+" -1 ./lead_trimmed/"+r1nm+".1.clone.trim.fastq -2 ./lead_trimmed/"+r2nm+".2.clone.trim.fastq")
                        else:
                                r1nm, r1ext = os.path.splitext(options.read1)
                                os.system("process_radtags -r -b barcodes.txt -o ./parsed --inline_null --renz_1 "+options.renz1+" --renz_2 "+options.renz2+" -f ./lead_trimmed/"+r1nm+".1.clone.trim.fastq")
                        print "\n***Quality-trimming reads using Trimmomatic***\n"

	else:
                if options.clean == True:
                        alert = open("./cleaned/ATTENTION", "w")
                        line = "You elected to quality-trim your reads using Stacks. This trimming was done simultaneously with parsing. See the 'parsed' folder for your trimmed reads."
                        alert.write(line)
                        alert.close()
                        if options.paired == True:
                                r1nm, r1ext = os.path.splitext(options.read1)
                                r2nm, r2ext = os.path.splitext(options.read2)
                                os.system("process_radtags -c -q -b barcodes.txt -o ./parsed --inline_index --renz_1 "+options.renz1+" --renz_2 "+options.renz2+" -1 ./lead_trimmed/"+r1nm+".1.clone.trim.fastq -2 ./lead_trimmed/"+r2nm+".2.clone.trim.fastq")
                        else:
                                r1nm, r1ext = os.path.splitext(options.read1)
                                os.system("process_radtags -c -q -b barcodes.txt -o ./parsed --inline_null --renz_1 "+options.renz1+" --renz_2 "+options.renz2+" -f ./lead_trimmed/"+r1nm+".1.clone.trim.fastq")
                        print "\n***Quality-trimming reads using Stacks***\n"
                else:
			if options.paired == True:
                                r1nm, r1ext = os.path.splitext(options.read1)
                                r2nm, r2ext = os.path.splitext(options.read2)
                                os.system("process_radtags -b barcodes.txt -o ./parsed --inline_index --renz_1 "+options.renz1+" --renz_2 "+options.renz2+" -1 ./lead_trimmed/"+r1nm+".1.clone.trim.fastq -2 ./lead_trimmed/"+r2nm+".2.clone.trim.fastq")
                        else:
                                r1nm, r1ext = os.path.splitext(options.read1)
                                os.system("process_radtags -b barcodes.txt -o ./parsed --inline_null --renz_1 "+options.renz1+" --renz_2 "+options.renz2+" -f ./lead_trimmed/"+r1nm+".1.clone.trim.fastq")
                        print "\n***Quality-trimming reads using Trimmomatic***\n"

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
### Place restriction site trimming routine here ###
		else:
			parse_single = "mv ./parsed/sample_"+bar[3]+".fq ./parsed/"+bar[0]+"_"+bar[3]+".1.fq"
			os.system(parse_single)
### Place restriction site trimming routine here ###


#################################################
###     	Quality-trim samples	      ###
#################################################

def quality_trim():
	if options.quality == True:
		for foo in open(options.sheet).read().splitlines():
        		bar = foo.split()
#			print bar
#			print bar[0], bar[1], bar[2], bar[3], bar[4]
			threads = options.threads
			if options.paired == True:
### Put command to trim away restriction site here and below else for Trimmomatic option ###
				handle = bar[0]+"_"+bar[3]+"-"+bar[4]
#                       	print handle
				PEclean = "$jar/trimmomatic-0.32.jar PE -threads "+threads+" -trimlog ./cleaned/"+handle+".qtrim.log ./parsed/"+handle+".1.fq ./parsed/"+handle+".2.fq ./cleaned/"+handle+".P1.qtrim ./cleaned/"+handle+".S1.qtrim ./cleaned/"+handle+".P2.qtrim ./cleaned/"+handle+".S2.qtrim LEADING:10 TRAILING:10 SLIDINGWINDOW:4:15 MINLEN:36 TOPHRED33"
				os.system(str(PEclean))
			else:
                        	handle = bar[0]+"_"+bar[3]
#                       	print handle
				SEclean = "$jar/trimmomatic-0.32.jar SE -threads "+threads+" -trimlog ./cleaned/"+handle+".qtrim.log ./parsed/"+handle+".1.fq ./cleaned"+handle+".S1.qtrim LEADING:10 TRAILING:10 SLIDINGWINDOW:4:15 MINLEN:36 TOPHRED33"
				os.system(str(SEclean))


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
