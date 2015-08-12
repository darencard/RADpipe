#!/usr/bin/env Rscript --vanilla --default-packages=base,utils,stats,MASS

## entropyStart.R
## 
## Version 1.0 (11 June, 2015)
## License: GNU GPLv2
## To report bugs or errors, please contact Daren Card (dcard@uta.edu).
## This script is provided as-is, with no support and no guarantee of proper or desirable functioning.
## 
## This script takes a genotype matrix with genotype uncertainty values (0-2), produced by
## genotype_from_VCF.py (format 4) and produces input files for MCMC chain initialization in Entropy.
## These Entropy files are produced by performing a Discriminant Analysis of Principle Components
## (DAPC, Jombart et al. 2010) in the genotype values, which is the approach adopted in 
## Gompert et al. (2014) to initialize Entropy MCMC chains for more rapid convergence. User must specify an  
## input genotype matrix file, a starting K value, and an ending K value (script will produce Entropy input 
## for all values of K designated).
## 
## Citations:
## Gompert, Z., L.K. Lucas, C.A. Buerkle, M.L. Forister, J.A. Fordyce, & C.C. Nice. 2014.
## Admixture and the organization of gneetic diversity in a butterfly species complex revealed 
## through common and rare genetic variants. Molecular Ecology 23 (18): 4555-4573.
## doi:10.1111/mec.12811.
## 
## Jombart, T., S. Devellard, & F. Balloux. 2010 Discriminant analysis of principle components: 
## a new method for the analysis of genetically structured populations. BMC Genetics 11 (94). 
## doi:10.1186/1471-2156-11-94.
## 
## Rscript entropyStart.R <input_genotypes> <startK> <endK>

require(MASS)

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
