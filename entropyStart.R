#!/usr/bin/env Rscript --vanilla --default-packages=base,utils,stats,MASS

library(MASS)

args<-commandArgs(TRUE)

entropyStart <- function(infile, startk, endk, ...){
  genoIn <- read.table(infile, header = FALSE, sep = " ", fill = TRUE)
  rmCol <- ncol(genoIn) - 1
  genoInTrim <- as.matrix(genoIn[-c(1:3), c(2:rmCol)])
  tGenoInTrim <- t(genoInTrim)
  class(tGenoInTrim) <- "numeric"

  for (k in startk:endk){
    pcaOut<-prcomp(tGenoInTrim,center=TRUE,scale=FALSE)
    kMeanOut<-kmeans(pcaOut$x[,1:5],k,iter.max=10,nstart=10,algorithm="Hartigan-Wong")
    ldaOut<-lda(x=pcaOut$x[,1:5],grouping=kMeanOut$cluster,CV=TRUE)
    write.table(round(ldaOut$posterior,7),file=paste0("lda.k",k,".out"),quote=F,row.names=F,col.names=F)
  }
}

entropyStart(args[1], args[2], args[3])
