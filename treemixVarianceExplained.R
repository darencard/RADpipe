#!/usr/bin/env Rscript --vanilla --default-packages=base,utils,stats

## treemixVarianceExplained.R
## 
## Version 1.0 (24 October, 2015)
## License: GNU GPLv2
## To report bugs or errors, please contact Daren Card (dcard@uta.edu).
## This script is provided as-is, with no support and no guarantee of proper 
## or desirable functioning.
## 
## This script calculated the total standard error in the covariance matrix 
## estimated from the data in Treemix and also estimates the amount of variance 
## in population relatedness explained by the Treemix model being summarized.
## The script assumes the user is in the directory with the Treemix output files
## and simply requires the stem/prefix name for the output files. Therefore,
## the .cov.gz, .modelcov.gz, and .covse.gz must be present in the working
## directory. Optional: one can specify a second 'v' arguement and then use 
## ' > out.txt' to pipe the stdout to an output file for downstream parsing.
## 
## Rscript treemixVarianceExplained.R <stem> [v]

args <- commandArgs(TRUE)

calcVarExplain <- function(stem){
  # read in two covariance matrices
  WHat = read.table(gzfile(paste(stem, ".cov.gz", sep = "")), as.is = T, 
                    head = T, quote = "", comment.char = "")
  W = read.table(gzfile(paste(stem, ".modelcov.gz", sep = "")), as.is = T, 
                 head = T, quote = "", comment.char = "")

  # set row and column names and order properly
  names(WHat) = rownames(WHat)
  names(W) = rownames(W)
  WHat = WHat[order(names(WHat)), order(names(WHat))]
  W = W[order(names(W)), order(names(W))]

  # calculate standard error of covariance matrix
  se = read.table(gzfile(paste(stem, ".covse.gz", sep = "")), as.is = T, 
                  head = T, quote = "", comment.char = "")
  seBar = apply(se, 1, mean)
  seBar = mean(seBar)

  # calculate residuals
  R = WHat - W

  # set residual lower matrix triangle to NA
  R[lower.tri(R, diag = TRUE)] <- NA

  # calculate mean residual across upper matrix triangle
  RBar = mean(as.matrix(R), na.rm = TRUE)

  # set data-estimated W(hat) matrix lower matrix triangle to NA
  WHat[lower.tri(WHat, diag = TRUE)] <- NA

  # calculate mean W(hat)
  WHatBar = mean(as.matrix(WHat), na.rm = TRUE)

  # calculate numerator of equation 30
  fNum = sum((R - RBar)^2, na.rm = TRUE)

  # calculate denominator of equation 30
  fDen = sum((WHat - WHatBar)^2, na.rm = TRUE)

  # calculate f in equation 30
  f = 1 - (fNum/fDen)

  outList = list("StdErr" = seBar, "VarExplain" = f)
  return(outList)
}

output <- calcVarExplain(args[1])

if (length(args) == 2){
  write(paste(paste0("Standard error for all entries in the covariance",
                     "matrix estimated from the data"), 
              output$StdErr, sep = "\t"), file = stderr())
  
  write(paste(paste0("Variance of relatedness between populations",
                     "explained by the model"),
              output$VarExplain, sep = "\t"), file = stderr())
}

write(paste(paste0("Standard error for all entries in the covariance ",
                   "matrix estimated from the data"), 
                    output$StdErr, sep = "\t"), file = stdout())

write(paste(paste0("Variance of relatedness between populations ",
                   "explained by the model"),
                   output$VarExplain, sep = "\t"), file = stdout())