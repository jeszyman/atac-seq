
## - Rscript

#!/usr/bin/env Rscript

#######################################
###   Make A Txdb From A Gtf File   ###
#######################################

# Command line arguements
args = commandArgs(trailingOnly = TRUE)
ensembl_gtf = args[1]

library(GenomicFeatures)

txdb = makeTxDbFromGFF(ensembl_gtf)

saveDb(txdb, file = "/tmp/db")
