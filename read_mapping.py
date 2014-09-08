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
parser.add_option("-p", action = "store_true", dest = "paired", help = "paired-end reads")
parser.add_option("-s", action = "store_true", dest = "single", help = "single-end reads only")

options, args = parser.parse_args()


		
def main():
#	os.system("mkdir mapping")
#	index_ref()
	for root,dirs,files in os.walk(options.directory):
		print "Running pipeline"
		print files
	for file in files:
		names = {}
		if file.endswith("."+options.ext):
			foo = file.split(os.extsep)
			name = foo[0]
			if name not in names.keys():
				names[name] = 1
		print name
#		if options.single == True:
#			SE_map(name)
#			sam2bam(files)
#			SE_bam_process(files)
#		elif options.paired == True:
#			PE_map(name)
#			SE_map(name)
#	if options.sams == True:
#		print "SAM output will be saved"
#	else:
#		os.system("rm -f ./mapping/*.sam")
		
main()