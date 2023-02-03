#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)
gtf_file = args[1]
txdb_file = args[2]

library(GenomicFeatures)

txdb = makeTxDbFromGFF(gtf_file,
                       format = "gtf")

saveDb(txdb, file = txdb_file)

#########1#########2#########3#########4#########5#########6#########7#########8

library(TxDb.Mmusculus.UCSC.mm10.ensGene)

# Create gene-transcript index as tx2gene object from a txdb
make_index = function(txdb){
  k <- keys(txdb, keytype = "TXNAME")
  tx2gene <- AnnotationDbi::select(txdb, k, "GENEID", "TXNAME")
}

tx2gene = make_index(TxDb.Mmusculus.UCSC.mm10.ensGene)

saveRDS(tx2gene, file = snakemake@output[["txdb"]])
