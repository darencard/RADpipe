#!/usr/bin/env python

import numpy as np
import numpy.random as npr
import pylab
import optparse

usage_line = """
A script to bootstrap over Fst data from a pairwise comparison from Stacks and calculate \
a confidence interval. User must specify the <batch>.fst_<pop1>-<pop2>.tsv file, \
the number of bootstrap reps, the value of alpha (0.05 by default). \
Output to terminal is the observed average Fst in the pairwise comparison and \
a confidence interval based on the boostrapping.
"""

usage = usage_line

parser = optparse.OptionParser(usage=usage)
parser.add_option("--input", action= "store", type= "string", dest="input", help="""The input <batch>.fst_<pop1>-<pop2>.tsv file""")
parser.add_option("--permutations", action="store", type= "string", dest="perms", help="""The number of bootstrap reps/permutations""")
parser.add_option("--alpha", action="store", type= "string", dest= "alpha", help="""The value of alpha for the confidence interval width [0.05]""", default = "0.05")
options, args = parser.parse_args()

def bootstrap(data, num_samples, statistic, alpha):
    ## Returns bootstrap estimate of 100.0*(1-alpha) CI for statistic.
    n = len(data)
    idx = npr.randint(0, n, (num_samples,n))
    samples = data[idx]
    for fst in np.arange(0, 1, 0.01):
    	for row in samples:
    		#print fst
    		#print row.shape[0]
    		cond = (row >= fst).sum()
    		total = row.shape[0]
    		if float(cond)/total <= 0.05:
    			print fst
    stat = np.sort(statistic(samples, 1))
    return (stat[int((alpha/2.0)*num_samples)],
            stat[int((1-alpha/2.0)*num_samples)])

if __name__ == '__main__':
    pop = []
    for line in open(options.input, "r"):
	    if not line.strip().startswith("#"):
        	record = line.rstrip().split("\t")
        	pop.append(float(record[18]))
    observed = float(sum(pop))/len(pop)
    x = np.array(pop)

    low, high = bootstrap(x, int(options.perms), np.std, float(options.alpha))
    
    obsline = "The observation from the input dataset is "+str(observed)+"."
    ciline = "The CI after "+options.perms+" bootstraps with alpha = "+options.alpha+" is "+str(low)+" - "+str(high)+"."
    print obsline
    print ciline+"\n\n"