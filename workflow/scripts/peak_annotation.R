# Arguements for testing
## granges_rds = "/home/jeszyman/repos/atac-seq/test/dca/dca_granges_regfilt.rds"
## annotation_csv = "/home/jeszyman/repos/atac-seq/test/dca/regfilt_dca.csv"
## chipseek_file = "/home/jeszyman/repos/atac-seq/test/dca/regfilt_chipseek.rds"

# Arguements for command line input
args = commandArgs(trailingOnly = TRUE)
granges_rds = args[1]
annotation_csv = args[2]
chipseek_file = args[3]

peaks = readRDS(granges_rds)

library(ChIPseeker)
library(csaw)
library(TxDb.Mmusculus.UCSC.mm10.ensGene)
library(tidyverse)

txdb = TxDb.Mmusculus.UCSC.mm10.ensGene

chipseek = annotatePeak(peaks, TxDb = txdb, annoDb = "org.Mm.eg.db")

annotation = as_tibble(as.data.frame(chipseek)) 

write.csv(annotation, row.names = F, file = annotation_csv)

saveRDS(object = chipseek,
        file = chipseek_file)

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
