#!/usr/bin/env Rscript

################################
###   Make Homer Bed Files   ###
################################

# Command line arguements
args = commandArgs(trailingOnly = TRUE)
dca_tsv = args[1]
anno_tsv = args[2]
qval_cut = args[3]
upbed_out = args[4]
downbed_out = args[5]
upbed_p_out = args[6]
downbed_p_out = args[7]
upbed_e_out = args[8]
downbed_e_out = args[9]

# Load required packages, data, and functions
library(tidyverse)
dca = read_tsv(dca_tsv)
anno = read_tsv(anno_tsv)

anno = anno %>% select(!c("start", "end","strand"))

# Main
up_tib =
  dca %>% filter(qval < qval_cut &
               logfc >0)

down_tib =
  dca %>% filter(qval < qval_cut &
               logfc < 0)

make_homer_bed = function(tibble){
  bed =
    tibble %>% mutate(chr = gsub(":.*$","",coordinate)) %>%
    mutate(start = gsub("^chr\\d+:(\\d+)-\\d+$", "\\1", coordinate)) %>%
    mutate(end = gsub("^chr\\d+:\\d+-(\\d+)$", "\\1", coordinate)) %>%
    mutate(col5 = "") %>%
    mutate(strand = "+") %>%
    select(chr, start, end, coordinate, col5, strand)
}

up_bed = make_homer_bed(up_tib)
up_bed %>% write_tsv(upbed_out, col_names = FALSE)
down_bed = make_homer_bed(down_tib)
down_bed %>% write_tsv(downbed_out, col_names = FALSE)


subset <- function(peaks_tibble, anno_tibble, annotation_string){
  result <- peaks_tibble %>%
    left_join(anno_tibble, "coordinate") %>%
    filter(grepl(annotation_string, annotation)) %>%
    select(chr, start, end, coordinate, col5, strand)

  return(result)
}

subset(up_bed, anno, "Promoter") %>% write_tsv(upbed_p_out, col_names = FALSE)
subset(down_bed, anno, "Promoter") %>% write_tsv(downbed_p_out, col_names = FALSE)
subset(up_bed, anno, "Distal") %>% write_tsv(upbed_e_out, col_names = FALSE)
subset(down_bed, anno, "Distal") %>% write_tsv(downbed_e_out, col_names = FALSE)
