#############################################################################
###            Script for csaw ATAC-seq local peak calling                ###
#############################################################################

# Arguements for integration testing
## library_ids_str = "/home/jeszyman/repos/atac-seq/test/bam/atac1_open_tn5.bam /home/jeszyman/repos/atac-seq/test/bam/atac2_open_tn5.bam /home/jeszyman/repos/atac-seq/test/bam/atac3_open_tn5.bam /home/jeszyman/repos/atac-seq/test/bam/atac4_open_tn5.bam"
## threads = 4
## out_background_rds = "/tmp/background.rds"
## out_rse_rds = "/tmp/rse.rds"

args = commandArgs(trailingOnly = TRUE)
library_ids_str = args[1]
threads = args[2]
out_background_rds = args[3]
out_rse_rds = args[4]

# Load packages
library(BiocParallel)
library(csaw)
library(edgeR)
library(tidyverse)

## Script-local functions
csaw_choose_window = function(bam_list){
  # Choose window width by fragment size distribution
  frag_size = lapply(bam_list, getPESizes)
  all_bam_frag_vect = frag_size %>% map(1) %>% as_vector()
  frag_vect_summary = summary(all_bam_frag_vect)
  thirdq = frag_vect_summary[[5]]
  return(thirdq)
}

# Specify params
surrounds = 2000
standard_chr <- paste0("chr", c(1:19)) # only use standard chromosomes
param = readParam(max.frag=1000, pe="both", restrict=standard_chr)

# Make bam file list
bam_list = unlist(strsplit(library_ids_str, " "))
names(bam_list) = gsub("^.*/","",bam_list)

# Filter per Reske JJ, et al. 2021. https://doi.org/10.1186/s13072-020-00342-y
(window = csaw_choose_window(bam_list))
counts = windowCounts(bam_list,
                      width = window,
                      param = param,
                      BPPARAM = MulticoreParam(workers=threads))
neighbor = suppressWarnings(resize(rowRanges(counts),
                                   surrounds, fix = "center"))
wider = regionCounts(bam_list,
                     regions = neighbor,
                     param = param,
                     BPPARAM = MulticoreParam(workers=threads))
# Remove dimnames to avoid SummarizedExperiment error in window filtering
dimnames(wider) = NULL
dimnames(counts) = NULL

filter_stat = filterWindowsLocal(counts, wider, assay.data = "counts")
filtered_counts = counts[filter_stat$filter > log2(3),]
background = windowCounts(bam_list,
                          bin = TRUE,
                          width = 10000,
                          param = param,
                          BPPARAM = MulticoreParam(workers=threads))

# Save outputs
saveRDS(object = filtered_counts, 
file = out_rse_rds)
saveRDS(object = background, 
file = out_background_rds)
