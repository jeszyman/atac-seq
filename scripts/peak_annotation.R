#!/usr/bin/env Rscript

####################################
###   Bamscale Peak Annotation   ###
####################################

# Command line arguements
args = commandArgs(trailingOnly = TRUE)
cov_tsv = args[1]
bmart_dataset = args[2]
txdb = args[3]
out_tsv = args[4]

# Load required packages, data, and functions

#cov_tsv = "~/cards/analysis/atac/mouse/dca/mouse_mm10.raw_coverages.tsv"
#txdb = "TxDb.Mmusculus.UCSC.mm10.knownGene"
#bmart_dataset = "mmusculus_gene_ensembl"
#out_tsv = "/tmp/test.tsv"

library(ChIPseeker)
library(txdb,character.only = TRUE)
library(biomaRt)
library(GenomicRanges)
library(tidyverse)

peaks = read_tsv(cov_tsv) %>%
  mutate(chr = gsub(":.*$","",coordinate)) %>%
  mutate(start = gsub("-.*$","", gsub("^.*:","",coordinate))) %>%
  mutate(end = gsub("^.*-","",coordinate)) %>%
  dplyr::select(coordinate, chr, start, end)

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

#########1#########2#########3#########4#########5#########6#########7#########8
###
###   Script to annotate csaw peaks   ###
###

args = commandArgs(trailingOnly = TRUE)
peaks_rds = args[1]
annotation_file = args[2]

peaks = readRDS(peaks_rds)

library(ChIPseeker)
library(csaw)
library(TxDb.Mmusculus.UCSC.mm10.ensGene)
library(tidyverse)

txdb = TxDb.Mmusculus.UCSC.mm10.ensGene

peak_loc = peaks

chipseek = annotatePeak(peak_loc, TxDb = txdb, annoDb = "org.Mm.eg.db")

annotation = as_tibble(as.data.frame(chipseek))

write.csv(annotation, row.names = F, file = annotation_file)

#########1#########2#########3#########4#########5#########6#########7#########8
###
###   Script to annotate csaw peaks   ###
###

args = commandArgs(trailingOnly = TRUE)
peaks_rds = args[1]
annotation_file = args[2]

peaks = readRDS(peaks_rds)

library(ChIPseeker)
library(csaw)
library(TxDb.Mmusculus.UCSC.mm10.ensGene)
library(tidyverse)

txdb = TxDb.Mmusculus.UCSC.mm10.ensGene

peak_loc = peaks

chipseek = annotatePeak(peak_loc, TxDb = txdb, annoDb = "org.Mm.eg.db")

annotation = as_tibble(as.data.frame(chipseek))

write.csv(annotation, row.names = F, file = annotation_file)

#########1#########2#########3#########4#########5#########6#########7#########8
###
###   Script to annotate csaw peaks   ###
###

args = commandArgs(trailingOnly = TRUE)
peaks_rds = args[1]
annotation_file = args[2]

peaks = readRDS(peaks_rds)

library(ChIPseeker)
library(csaw)
library(TxDb.Mmusculus.UCSC.mm10.ensGene)
library(tidyverse)

txdb = TxDb.Mmusculus.UCSC.mm10.ensGene

peak_loc = peaks

chipseek = annotatePeak(peak_loc, TxDb = txdb, annoDb = "org.Mm.eg.db")

annotation = as_tibble(as.data.frame(chipseek))

write.csv(annotation, row.names = F, file = annotation_file)
