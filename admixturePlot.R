#!/usr/bin/env Rscript --vanilla --default-packages=base,utils,stats,RColorBrewer

## admixturePlot.R
## 
## Version 1.0 (11 June, 2015)
## License: GNU GPLv2
## To report bugs or errors, please contact Daren Card (dcard@uta.edu).
## This script is provided as-is, with no support and no guarantee of proper or desirable functioning.
## 
## This script takes formated admixture proportions with metadata (from meta_sort_NGSadmix.py) and
## creates a publication-quality bar plot with names taken from columns 2 (sample) and 4 (location)
## of the metadata (from the samplesheet). User must specify the input file, an output file name (PDF),
## and the number of populations (K). The output is a 8.5" x 11" PDF.
## 
## Rscript admixturePlot.R <input_file> <output_file> <#K>

args<-commandArgs(TRUE)

require("RColorBrewer")

admixturePlot <- function(infile, outfile, k, ...){
  admixIn <- read.table(infile, header=FALSE, sep="\t", fill=TRUE)
  print(admixIn)
  tAdmixIn <- t(as.matrix(admixIn[1:nrow(admixIn), 6:ncol(admixIn)-1]))
  print(tAdmixIn)
  mycol <- brewer.pal(as.numeric(k), "Paired")
  
  pdf(outfile, 11, 7)
  par(mar=c(12,5,2,2))
  barplot(tAdmixIn, col=mycol, space=0.2, names.arg=paste(admixIn[1:nrow(admixIn), 2], admixIn[1:nrow(admixIn), 4], sep="_"), las=2, ylab="Admixture")
  dev.off()
}

admixturePlot(args[1], args[2], args[3])
