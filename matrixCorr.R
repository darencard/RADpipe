#!/usr/bin/env Rscript --vanilla --default-packages=base,utils,stats

## matrixCorr.R
## 
## Version 1.0 (18 November, 2015)
## License: GNU GPLv2
## To report bugs or errors, please contact Daren Card (dcard@uta.edu).
## This script is provided as-is, with no support and no guarantee of proper or desirable functioning.
## 
## This script calculates the correlation between all values two data matrices
## and report the coefficient and p-value to the user. The user simple specifies the 
## two input data matrices as two arguments to the script. It is assumed that each 
## represents comma-separated values (CSV) and that there is a single header row and 
## column describing the data. This is useful for comparing genotype matrices from 
## separate analyses and was essentially written for this purpose.
## 
## Rscript matrixCorrt.R <input_file1> <input_file2>

args <- commandArgs(TRUE)
cat("Comparison between: ")
cat(args, sep=" & ")

genoCorr <- function(infile1, infile2, ...){
    genotypes1 <- read.csv(as.character(infile1))
    genotypes2 <- read.csv(as.character(infile2))

    genoCorrelation <- cor.test(c(as.matrix(genotypes1[-1,-1])), 
                                c(as.matrix(genotypes2[-1,-1])))
    cat("\n", "Correlation coefficient = ", genoCorrelation$estimate, sep="")
    cat("\n", "P-value = ", genoCorrelation$p.value, "\n\n\n", sep="")
}

genoCorr(args[1], args[2])
