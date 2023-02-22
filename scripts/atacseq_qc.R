
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)
dup_bam_str = args[1]
proc_bam_str = args[2]
txdb_file = args[3]
atac_qc_file = args[4]

library(ATACseqQC)
library(tidyverse)
library(AnnotationDbi)

txdb = loadDb(txdb_file)

split_filename_str = function(filename_str){
  vect = strsplit(filename_str, " ")[[1]]
  return(vect)
}

dup_bam_vect = split_filename_str(dup_bam_str)
proc_bam_vect = split_filename_str(proc_bam_str)
bam_vect = data.frame(
  dup = dup_bam_vect,
  proc = proc_bam_vect
)

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
