args = commandArgs(trailingOnly = TRUE)
in_peaks_bed = args[1]
sample_sheet_tsv = args[2]
atac_macs2_dir = args[3]
out_peaks_bed = args[4]

library(tidyverse)

make_naive = function(sample_sheet,in_peaks_bed,atac_macs2_dir,out_peaks_bed){
  libraries = read_tsv(sample_sheet)
  library_filt = substr(gsub("^.*lib","lib",in_peaks_bed), 1, 6)
  group = libraries %>% filter(library == library_filt) %>% slice_head(n = 1L) %>% pull(group)
  consensus = paste0(atac_macs2_dir,"/",group,"_consensus.bed")
  system(paste("bedops --element-of 50%", in_peaks_bed, consensus, ">", out_peaks_bed))
}

make_naive(sample_sheet_tsv, in_peaks_bed, atac_macs2_dir, out_peaks_bed)
