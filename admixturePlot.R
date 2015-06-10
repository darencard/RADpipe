#!/usr/bin/env Rscript --vanilla --default-packages=base,utils,stats,RColorBrewer

args<-commandArgs(TRUE)

require("RColorBrewer")

admixturePlot <- function(infile, outfile, k, ...){
  admixIn <- read.table(infile, header=FALSE, sep="\t", fill=TRUE)
  tAdmixIn <- t(as.matrix(admixIn[,5:ncol(admixIn)]))
  mycol <- brewer.pal(k, "Paired")
  
  pdf(outfile, 11, 8.5)
  par(mar=c(12,5,2,2))
  barplot(tAdmixIn, col=mycol space=0.2, names.arg=paste(admixIn[,2], input[,4], sep="_", las=2, ylab="Admixture"))
  dev.off()
}

admixturePlot(args[1], args[2], args[3])

# 
# #input <- read.table("Boco.NGSadmix.k10.metadata.txt")
# input <- read.table("lda.k4.metadata.txt")
# admix <- t(as.matrix(input[,4:ncol(input)]))
# 
# require("RColorBrewer")
# mycol <- brewer.pal(4, "Paired")
# #pdf("Boco.NGSadmix.k10.plot.pdf", 11, 7)
# pdf("lda.k4.plot.pdf", 11, 8.5)
# par(mar=c(12,5,2,2))
# barplot(admix, col=mycol, space=0.2, names.arg=paste(input[,1], input[,3], sep="_"), las=2, ylab="Admixture")
# dev.off()
