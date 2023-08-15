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
write_tsv(anno, out_peak_bed)

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
