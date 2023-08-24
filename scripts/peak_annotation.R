#!/usr/bin/env Rscript

#################################
###   Macs2 Peak Annotation   ###
#################################

# Command line arguements
args = commandArgs(trailingOnly = TRUE)
in_peak_bed = args[1]
txdb = args[2]
out_peak_bed = args[3]

# Load required packages, data, and functions
library(tidyverse)
library(ChIPpeakAnno)
library(rtracklayer)
library(ChIPseeker)
library(txdb, character.only = T)

peaks = rtracklayer::import(in_peak_bed)
anno = annotatePeak(peaks, TxDb = get(txdb))

anno = as_tibble(anno)

if (!"peak" %in% names(anno)) {
  anno <- anno %>%
    mutate(peak = end - start)
}

anno =
  anno %>% mutate(simple = case_when(
                    grepl("Promoter", annotation) ~ "promoter",
                    grepl("Exon", annotation) ~ "exon",
                    grepl("Intron", annotation) ~ "intron",
                    grepl("3' UTR", annotation) ~ "utr3",
                    grepl("5' UTR", annotation) ~ "utr5",
                    grepl("Distal Intergenic", annotation) ~ "intergenic",
                    grepl("Downstream", annotation) ~ "downstream",
                    TRUE ~ annotation
                  ))

anno =
  anno %>% dplyr::select(seqnames, start, end, width, strand, name, score, signalValue, pValue, qValue, peak, everything())

write_tsv(anno, out_peak_bed, col_names = F)
