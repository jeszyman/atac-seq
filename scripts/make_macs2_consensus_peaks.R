#!/usr/bin/env Rscript
#########1#########2#########3#########4#########5#########6#########7#########8
###
###   Makes consensus peak sets from GRanges MACS2 peaks   ###
###

peaks = list.files(path = "/mnt/ris/jschwarz/cardiac-radiobiology/atac/macs2",
                   pattern = "regfilt_narrow.*rds$",
                   full.names = TRUE)
names(peaks) = gsub(".rds","",list.files(path = "/mnt/ris/jschwarz/cardiac-radiobiology/atac/macs2",
                   pattern = "regfilt_narrow.*rds$",
                   full.names = FALSE))

load("/mnt/ris/jschwarz/cardiac-radiobiology/data_model/data_model.RData")

library(tidyverse)

sham_libs =
 libraries_full %>%
 filter(lib_typ == "atac") %>%
 filter(cohort_id == "sham") %>%
 pull(library_id)
sham_peaks_list = grep(paste(sham_libs, collapse = "|"), peaks, value = TRUE)

sham_peaks_list = lapply(sham_peaks_list, readRDS)

library(GenomicRanges)

library(GenomicAlignments)

gr1 = sham_peaks_list[[1]]
gr2 = sham_peaks_list[[2]]

test = subsetByOverlaps(gr1, gr2, ignore.strand = TRUE)

test


test = GRangesList(unlist(sham_peaks_list))

test2 = summarizeOverlaps(gr1, test, mode = )

head(assays(test2)$counts, n = 100)

sham_peaks = GenomicRanges::union(unlist(sham_peaks_list))

test = c(sham_peaks_list[[1]], sham_peaks_list[[2]])

test = sham_peaks_list[1][sham_peaks_list[1] %over% sham_peaks_list[2]]

sham_peaks_list[[1]] %over% sham_peaks_list[[2]] %over% sham_peaks_list[[3]]

union(sham_peaks_list[[1]]

sham_peaks_list
sham_peaks = GenomicRanges::union(unlist(sham_peaks_list))

test = unlist(sham_peaks_list)

class(test)


head(test
     )
head(unlist(sham_peaks_list))

test = union(sham_peaks_list[[1]],sham_peaks_list[[2]])

test = for (i in sham_peaks_list) {union(i)}



test = for (i in 1:length(sham_peaks_list)) {union (sham_peaks_list[[i]])}

union(i )
length(sham_peaks_list)

ir6w_libs =
  libraries_full %>%
  filter(lib_typ == "atac") %>%
  filter(cohort_id == "ir6w") %>%
  pull(library_id)
ir6w_peaks = grep(paste(ir6w_libs, collapse = "|"), peaks, value = TRUE)

ir48h_libs =
  libraries_full %>%
  filter(lib_typ == "atac") %>%
  filter(cohort_id == "ir48h") %>%
  pull(library_id)
ir48h_peaks = grep(paste(ir48h_libs, collapse = "|"), peaks, value = TRUE)
