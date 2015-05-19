#!/usr/bin/env python

import numpy as np
import numpy.random as npr
import pylab
import optparse

usage_line = """
A script to determine a significance threshold by bootstrapping over a specified statistic. \
User must specify an input file, the column of data to use with bootstrapping, the number of \
bootstrap reps, the value of alpha (0.05 by default), and whether to report the 1-tailed (e.g., \
useful for absolute values of statistics) or 2-tailed threshold. Note: header line must begin with \
a '#'.
"""

usage = usage_line

parser = optparse.OptionParser(usage=usage)
parser.add_option("--input", action= "store", type= "string", dest="input", help="""The input file""")
parser.add_option("--permutations", action="store", type= "string", dest="perms", help="""The number of bootstrap reps/permutations""")
parser.add_option("--alpha", action="store", type= "string", dest= "alpha", help="""The value of alpha for the threshold [0.05]""", default = "0.05")
parser.add_option("--tails", action="store", type="string", dest="tails", help="""Specify whether this is one-tailed (1) or two-tailed (2)""")
parser.add_option("--column", action="store", type="string", dest="column", help="The column of data to be bootstrapped""")
options, args = parser.parse_args()

def bootstrap(data, num_samples, statistic, alpha_raw):
    n = len(data)
    idx = npr.randint(0, n, (num_samples,n))
    samples = data[idx]
    if options.tails == "1":
    	alpha = float(alpha_raw)*100
    	beta = 100-alpha
    	low = np.amin(samples)
    	high = np.sort(statistic(samples, beta))
    elif options.tails == "2":
        alpha = float((alpha_raw*100)/2)
    	beta = 100-alpha
    	low = np.sort(statistic(samples, alpha))
    	high = np.sort(statistic(samples, beta))
    else:
    	print "\n\n**Error: Specify whether you want one-tailed or two-tailed thresholds!**\n\n"
    return (np.median(low),
    		np.median(high))

if __name__ == '__main__':
	pop = []
	for line in open(options.input, "r"):
	    if not line.strip().startswith("#"):
	    	record = line.rstrip().split("\t")
	    	num = record[int(options.column)-1]
	    	pop.append(float(num))
	x = np.array(pop)

	low, high = bootstrap(x, int(options.perms), np.percentile, float(options.alpha))
   
	ciline = "The "+options.tails+"-tailed bootstrapping significance threshold after "+options.perms+" bootstraps with alpha = "+options.alpha+" is "+str(low)+" - "+str(high)+"."
	print "\n\n"+ciline+"\n\n"