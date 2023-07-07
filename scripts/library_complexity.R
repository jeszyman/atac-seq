#!/usr/bin/env Rscript
#########1#########2#########3#########4#########5#########6#########7#########8
###
###   R Script to assess ATAC-seq library complexity by fragment length   ###
###

args = commandArgs(trailingOnly = TRUE)
bam = args[1]
rds = args[2]

library(preseqR)
library(ATACseqQC)
library(Rsamtools)
library(TxDb.Mmusculus.UCSC.mm10.ensGene)

libCompWrap = function(dup_bam){
  estimateLibComplexity(readsDupFreq(dup_bam))
}

complex = libCompWrap("/mnt/ris/jschwarz/cardiac-radiobiology/inputs/Rentschler_s5469_MGI2048/lc-08.TCGTGATCAG-ACACTACGTA/lc-08.TCGTGATCAG-ACACTACGTA.genome_accepted_hits.bam")

saveRDS(object = complex,
        file = rds)
