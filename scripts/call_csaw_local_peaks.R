
#########1#########2#########3#########4#########5#########6#########7#########8

#######################################################################
###    Script to call ATAC-seq peaks using local windows in csaw    ###
#######################################################################

# Setup
##
## Snakemake
args = commandArgs(trailingOnly = TRUE)
bam_dir = args[1]
bam_pattern = args[2]
rse = args[3]
dge = args[4]

## Libraries
library(csaw)
library(edgeR)
library(tidyverse)

## Script-local variables
surrounds = 2000
standard_chr <- paste0("chr", c(1:19)) # only use standard chromosomes
param = readParam(max.frag=1000, pe="both", restrict=standard_chr)

bam_list = list.files(path = bam_dir,
                  pattern = bam_pattern,
                  full.names = TRUE)

names(bam_list) = gsub(bam_pattern, "", list.files(path = bam_dir,
                                                   pattern = bam_pattern,
                                                   full.names = FALSE))

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

counts = windowCounts(bam_list, width = window, param = param)

neighbor = suppressWarnings(resize(rowRanges(counts),
                                   surrounds, fix = "center"))

wider = regionCounts(bam_list, regions = neighbor, param = param)

filter_stat = filterWindowsLocal(counts, wider)

filtered_counts = counts[filter_stat$filter > log2(3),]

background = windowCounts(bam_list, bin=TRUE, width=10000, param = param)

filtered_counts = normFactors(background, se.out = filtered_counts)

edger_input <- asDGEList(filtered_counts)

colnames(edger_input$counts) = colnames(filtered_counts)
rownames(edger_input$samples) = colnames(filtered_counts)

saveRDS(object = filtered_counts,
        file = rse)
saveRDS(object = edger_input,
        file = dge)
