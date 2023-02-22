
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)
gtf_file = args[1]
txdb_file = args[2]

library(GenomicFeatures)

txdb = makeTxDbFromGFF(gtf_file,
                       format = "gtf")

saveDb(txdb, file = txdb_file)
