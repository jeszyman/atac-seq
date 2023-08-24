#!/usr/bin/env Rscript

################################
###   Macs2 Peak Filtering   ###
################################

args = commandArgs(trailingOnly = TRUE)
chrs_tsv = args[1]
libraries_full_rds = args[2]
peak_file_str = args[3]
corces_min = args[4]
all_peaks_bed = args[5]
cluster_bed = args[6]
lib_peaks_min = args[7]
keep_bed = args[8]

# Load required packages and data
library(tidyverse)

chrs = read_tsv(chrs_tsv, col_names = c("chr")) %>% pull(chr)
libraries_full = readRDS(libraries_full_rds)
libraries_full = libraries_full %>% dplyr::select(!end)
peak_files = unlist(strsplit(peak_file_str, " "))
corces_min = as.numeric(corces_min)

# Create single peak file by library

ingest_macs2 <- function(peak) {
  col_names <- c("chr",
                 "start",
                 "end",
                 "width",
                 "strand",
                 "name",
                 "score",
                 "signalValue",
                 "pValue",
                 "qValue",
                 "peak",
                 "annotation",
                 "geneChr",
                 "geneStart",
                 "geneEnd",
                 "geneLength",
                 "geneStrand",
                 "geneId",
                 "transcriptId",
                 "distanceToTSS",
                 "simple")
  macs2peak <- read_tsv(peak, col_names = col_names)
  return(macs2peak)
}

peak_list = lapply(peak_files, ingest_macs2)
names(peak_list) = substr(gsub("^.*lib", "lib", peak_files), 1, 6)
peaks = bind_rows(peak_list, .id = "library")
peaks = peaks %>% left_join(libraries_full, by = "library")


corces_peaks =
  peaks %>%
  mutate(summit = start + peak) %>%
  mutate(start = summit - 250) %>%
  mutate(end = summit + 250) %>%
  # Remove sex chromosome and mitochondrial peaks here
  filter(chr %in% chrs) %>%
  group_by(library) %>%
  mutate(corces = pValue/sum(pValue/1000000)) %>% ungroup() %>%
  filter(corces > corces_min) %>% select(chr, start, end, library, corces, name, simple, annotation, geneId)

write_tsv(corces_peaks, file = all_peaks_bed, col_names = F)

system(paste0("bedtools sort -i ", all_peaks_bed, " | bedtools cluster -i - > ", cluster_bed))

clust = read_tsv(cluster_bed, col_names = c("chr", "start", "end", "library", "corces", "name", "simple", "annotation", "geneId", "clust"))

max =
  clust %>%
  group_by(clust) %>%
  filter(n() > 2) %>%
  slice_max(corces)

keep_libs =
  max %>%
  group_by(library) %>%
  summarize(sum = n()) %>%
  filter(sum > lib_peaks_min) %>%
  pull(library)

keep = max %>% filter(library %in% keep_libs) %>% write_tsv(keep_bed, col_names = F)

write_tsv(keep, file = keep_bed, col_names = FALSE)
