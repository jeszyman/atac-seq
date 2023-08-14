#!/usr/bin/env Rscript

####################################
###   Bamscale Peak Annotation   ###
####################################

# Command line arguements
args = commandArgs(trailingOnly = TRUE)
cov_bed = args[1]
bmart_dataset = args[2]
txdb = args[3]
out_tsv = args[4]

# Load required packages, data, and functions

library(ChIPseeker)
library(txdb,character.only = TRUE)
library(biomaRt)
library(GenomicRanges)
library(tidyverse)

peaks = read_tsv(cov_bed, col_names = c("chr",
                                        "start",
                                        "end",
                                        "library")) %>%
  dplyr::select(chr, start, end, library)

# Create a GRanges object from the data frame
gr <- makeGRangesFromDataFrame(peaks, start.field = "start", end.field = "end",
                               seqnames.field = "chr")

anno = annotatePeak(gr, TxDb = eval(parse(text=txdb)))
anno = as_tibble(anno)

anno = anno %>% mutate(entrezgene_id = as.numeric(geneId))

geneIds = anno %>% pull(geneId) %>% unique()

mart <- useMart("ensembl")
mart <- useDataset(bmart_dataset, mart)

names = getBM(
  filters = "entrezgene_id",
  attributes=c("ensembl_gene_id",
               "entrezgene_id",
               "description",
               "external_gene_name",
               "gene_biotype"),
  values = geneIds,
  mart = mart)
names = as_tibble(names)

final = anno %>% left_join(names, by = c("entrezgene_id")) %>%
  mutate(coordinate = paste0(seqnames, ":",start,"-",end))

final <- distinct(final, coordinate, .keep_all = TRUE)
write_tsv(final, file = out_tsv)
