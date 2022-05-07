#############################################################################
###            Script for csaw ATAC-seq local peak calling                ###
#############################################################################

# Arguements for testing
bam_dir = "~/repos/atac-seq/test/bam"
bam_pattern = "regfilt"
filt_libs_str = "atac1 atac2 atac4 atac4"
threads = 4
out_background_rds = "/tmp/background.rds"
out_rse_rds = "/tmp/rse.rds"

## args = commandArgs(trailingOnly = TRUE)
## bam_dir = args[1]

#REMOVE THIS STEP
# Split the filtered libraries string
(filt_libs = unlist(strsplit(filt_libs_str, " ")))

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

#EXPORT TO SMK
# Specify params
surrounds = 2000
standard_chr <- paste0("chr", c(1:19)) # only use standard chromosomes
param = readParam(max.frag=1000, pe="both", restrict=standard_chr)

paste0(bam_pattern,"$")

bam_list = list.files(path = bam_dir,
                       pattern = paste0(bam_pattern, ".bam$"),
                       full.names = TRUE)
names(bam_list) = gsub(bam_pattern, "", list.files(path = bam_dir,
                                                   pattern = paste0(bam_pattern, ".bam$"),                                                   
                                                   full.names = FALSE))
(bam_list = bam_list[names(bam_list) %in% filt_libs])

filt_libs
names(bam_list)

#########1#########2#########3#########4#########5#########6#########7#########8

(bam_list = list.files(path = bam_dir,
                       pattern = "_regfilt_tn5.bam$",
                       full.names = TRUE))

(names(bam_list) = gsub("_regfilt_tn5.bam$","",
                       list.files(path = bam_dir,
                                  pattern = "_regfilt_tn5.bam$",
                                  full.names = FALSE)))

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

filter_stat = filterWindowsLocal(counts, wider)

filtered_counts = counts[filter_stat$filter > log2(3),]

background = windowCounts(bam_list,
                          bin = TRUE,
                          width = 10000,
                          param = param,
                          BPPARAM = MulticoreParam(workers=threads))

saveRDS(object = filtered_counts, 
file = counts_rds)

saveRDS(object = background, 
file = background_rds)

#!/usr/bin/env Rscript

#############################################################################
###              Script for csaw ATAC-seq local peak calling 
#############################################################################

args = commandArgs(trailingOnly = TRUE)
bam_dir = args[1]
bam_pattern = args[2]
filt_libs_str = args[3]
threads = args[4]
background_rds = args[5]
rse_rds = args[6]

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

filter_stat = filterWindowsLocal(counts, wider)

filtered_counts = counts[filter_stat$filter > log2(3),]

background = windowCounts(bam_list,
                          bin = TRUE,
                          width = 10000,
                          param = param,
                          BPPARAM = MulticoreParam(workers=threads))

saveRDS(object = filtered_counts, 
file = counts_rds)

saveRDS(object = background, 
file = background_rds)

#!/usr/bin/env Rscript

#############################################################################
###              Script for csaw ATAC-seq local peak calling 
#############################################################################

args = commandArgs(trailingOnly = TRUE)
bam_dir = args[1]
bam_pattern = args[2]
filt_libs_str = args[3]
threads = args[4]
background_rds = args[5]
rse_rds = args[6]

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

filter_stat = filterWindowsLocal(counts, wider)

filtered_counts = counts[filter_stat$filter > log2(3),]

background = windowCounts(bam_list,
                          bin = TRUE,
                          width = 10000,
                          param = param,
                          BPPARAM = MulticoreParam(workers=threads))

saveRDS(object = filtered_counts, 
file = counts_rds)

saveRDS(object = background, 
file = background_rds)
