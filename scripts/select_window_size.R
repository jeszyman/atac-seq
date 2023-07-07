
#############################################################################
###            Script for csaw ATAC-seq local peak calling                ###
#############################################################################

# Setup

## Test arguements
## groups_str = "ir48h ir48h sham sham"
## library_ids_str = "/home/jeszyman/repos/atac-seq/test/bam/atac1_open_tn5.bam /home/jeszyman/repos/atac-seq/test/bam/atac2_open_tn5.bam /home/jeszyman/repos/atac-seq/test/bam/atac3_open_tn5.bam /home/jeszyman/repos/atac-seq/test/bam/atac4_open_tn5.bam"
## out_rse_rds = "/home/jeszyman/repos/atac-seq/test/csaw/norm_counts_rse.rds"
## out_dge_rds = "/home/jeszyman/repos/atac-seq/test/csaw/dge.rds"
## threads = 4

## Command line arguements
args = commandArgs(trailingOnly = TRUE)
library_ids_str = args[1]
threads = args[2]
out_rse_rds = args[3]
out_dge_rds = args[4]
groups_str = args[5]

## Load packages
library(BiocParallel)
library(csaw)
library(edgeR)
library(tidyverse)

# Specify csaw window parameters
surrounds = 2000
autosomes <- paste0("chr", c(1:19)) # only use autosomes
param = readParam(max.frag=1000, pe="both", restrict=autosomes)

# Make bam file list
bam_list = unlist(strsplit(library_ids_str, " "))
names(bam_list) = gsub("^.*/","",bam_list)

# Filter per Reske JJ, et al. 2021. https://doi.org/10.1186/s13072-020-00342-y

## Choose window width by fragment size distribution
csaw_choose_window = function(bam_list){
  frag_size = lapply(bam_list, getPESizes)
  all_bam_frag_vect = frag_size %>% map(1) %>% as_vector()
  frag_vect_summary = summary(all_bam_frag_vect)
  thirdq = frag_vect_summary[[5]]
  return(thirdq)
}

window = csaw_choose_window(bam_list)

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

# Return library names
colnames(filtered_counts) = names(bam_list)
colnames(background) = names(bam_list)

filtered_counts = normFactors(background, se.out = filtered_counts)
