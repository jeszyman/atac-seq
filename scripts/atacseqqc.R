
## - [[file:workflow/scripts/atac-seq_qc.R][Rscript]]

# RUNS BUT DOES NOT SAVE frag.len <- fragSizeDist(proc_bam, "label")
#!/usr/bin/env Rscript

#############################
###   Atacseqqc Wrapper   ###
#############################


# Load required packages
args = commandArgs(trailingOnly = TRUE)
in_bam_dup = args[1]
in_bam_dedup = args[2]
out_rda = args[3]

# Load necessary packages
library(ATACseqQC)
library(tidyverse)
library(AnnotationDbi)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

# Load data
txdb = TxDb.Hsapiens.UCSC.hg38.knownGene

# estimateLibComplexity uses preseqR to generate "a data frame of
# 3 columns: relative sequence depth, number of distinct
# fragments, number of putative sequenced reads"
freq = readsDupFreq(in_bam_dup)
libcomp = estimateLibComplexity(freq)

# tsse_df returns the plot values of TSSEscore with the score itself being
# https://www.encodeproject.org/data-standards/terms/#enrichment
txs = transcripts(txdb)
gal = readBamFile(proc_bam, bigFile=TRUE)
tsse_list = TSSEscore(gal, txs)
tsse_df = data.frame(
  tsse = tsse_list[1],
  distance = 100*(-9:10-.5)
)
tsse = tsse_list[2]

atacqc = function(dup_bam, proc_bam, txdb){
  freq = readsDupFreq(dup_bam)
  libcomp = estimateLibComplexity(freq)
  txs = transcripts(txdb)
  gal = readBamFile(proc_bam)
  tsse_list = TSSEscore(gal, txs)
  tsse_df = data.frame(
    tsse = tsse_list[1],
    distance = 100*(-9:10-.5)
  )
  tsse = tsse_list[2]
  atac = list(libcomp, tsse, tsse_df)
  names(atac) = c("libcomp_df", "tsse", "tsse_df")
  return(atac)
}

atac_qc_out = mapply(atacqc, dup_bam_vect, proc_bam_vect, MoreArgs = list(txdb = txdb))

save(atac_qc_out, file = atac_qc_file)
