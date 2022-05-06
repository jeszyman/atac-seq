#############################################################################
###              Script for csaw ATAC-seq local peak calling 
#############################################################################

args = commandArgs(trailingOnly = TRUE)
bam_dir = args[1]
bam_pattern = args[2]
filt_libs_str = args[3]
threads = args[4]
background_rds = args[5]
counts_rds = args[6]

filt_libs = unlist(strsplit(filt_libs_str, " "))

library(BiocParallel)
library(csaw)
library(edgeR)
library(tidyverse)

surrounds = 2000
standard_chr <- paste0("chr", c(1:19)) # only use standard chromosomes
param = readParam(max.frag=1000, pe="both", restrict=standard_chr)

bam_list = list.files(path = bam_dir,
                      pattern = bam_pattern,
                      full.names = TRUE)
names(bam_list) = gsub(bam_pattern, "", list.files(path = bam_dir,
                                                   pattern = bam_pattern,
                                                   full.names = FALSE))
bam_list = bam_list[names(bam_list) %in% filt_libs]

## Script-local functions
csaw_choose_window = function(bam_list){
  # Choose window width by fragment size distribution
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

dimnames(wider) = c()
dimnames(counts) = c()

filter_stat = filterWindowsLocal(counts, wider)

filtered_counts = counts[filter_stat$filter > log2(3),]

background = windowCounts(bam_list,
                          bin = TRUE,
                          width = 10000,
                          param = param,
                          BPPARAM = MulticoreParam(workers=threads))

colnames(filtered_counts) = names(bam_list)

saveRDS(object = filtered_counts, 
file = counts_rds)

colnames(background) = names(bam_list)

saveRDS(object = background, 
file = background_rds)
