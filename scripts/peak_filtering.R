## peak_file_str = "~/cards/analysis/atac/mouse/peaks/lib051_mm10_ds9_multi_peaks.narrowPeak ~/cards/analysis/atac/mouse/peaks/lib122_mm10_ds9_multi_peaks.narrowPeak"
## libraries_full_rds = "~/cards/data-model/libraries_full.rds"
## corces_min = 5
## all_peaks_bed = "/tmp/test.bed"
## cluster_bed = "/tmp/clust.bed"
## max_bed = "/tmp/max.bed"


#!/usr/bin/env Rscript

################################
###   Macs2 Peak Filtering   ###
################################

args = commandArgs(trailingOnly = TRUE)
chrs_tsv = args[1]
libraries_full_rds = args[2]
peak_file_str = args[3]
all_peaks_bed = args[4]
corces_min = args[5]
lib_peaks_min = args[6]


# Load required packages and data
library(tidyverse)
chrs = read_tsv(chrs_tsv, col_names = c("chr")) %>% pull(chr)
libraries_full = readRDS(libraries_full_rds)
peak_files = unlist(strsplit(peak_file_str, " "))

# Create single peak file by library

ingest_macs2 = function(peak){
  macs2peak = read_tsv(peak,
                        col_names = c("chr","start","end","peak","score","strand","signal","neg_l10_pval","neg_l10_qval", "dsummit"))
}

peak_list = lapply(peak_files, ingest_macs2)
names(peak_list) = substr(gsub("^.*lib", "lib", peak_files), 1, 6)
peaks = bind_rows(peak_list, .id = "library")
peaks = peaks %>% left_join(libraries_full, by = "library")

print(peaks, n = 10)


#
corces_peaks =
  peaks %>%
  mutate(summit = start + dsummit) %>%
  mutate(start = summit - 250) %>%
  mutate(end = summit + 250) %>%
  # Remove sex chromosome and mitochondrial peaks here
  filter(chr %in% chrs) %>%
  group_by(library) %>%
  mutate(corces = neg_l10_pval/sum(neg_l10_pval/1000000)) %>% ungroup() %>%
  filter(corces > corces_min) %>% select(chr, start, end, library, corces, peak)

print(corces_peaks, n = 10)

write_tsv(corces_peaks, file = all_peaks_bed, col_names = F)

#system(paste0("bedtools sort -i ", all_peaks_bed, " | bedtools cluster -i - > ", #cluster_bed))

## clust = read_tsv(cluster_bed, col_names = c("chr","start","end","library","corces","peak","clust"))

## max =
##   clust %>%
##   group_by(clust) %>%
##   filter(n() > 2) %>%
##   slice_max(corces) %>%
##   write_tsv(max_bed, col_names = F)


## keep_libs =
##   max %>%
##   group_by(library) %>%
##   summarize(sum = n()) %>%
##   filter(sum > lib_peaks_min) %>%
##   pull(library)

## keep = max %>% filter(library %in% keep_libs) %>% write_tsv(keep_bed, col_names = F)
