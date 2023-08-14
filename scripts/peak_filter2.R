#!/usr/bin/env Rscript

# Command line arguements
args = commandArgs(trailingOnly = TRUE)
cluster_bed = args[1]
lib_peaks_min = args[2]
keep_bed = args[3]


library(tidyverse)

clust = read_tsv(cluster_bed, col_names = c("chr","start","end","library","corces","peak","clust"))

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
