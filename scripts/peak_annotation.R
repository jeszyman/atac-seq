
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
