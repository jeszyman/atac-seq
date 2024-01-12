
## - [[file:./scripts/peak_union.R][Rscript]]

#!/usr/bin/env Rscript

########################################
###   Make Atac Peak Union Bedfile   ###
########################################

# Command line arguements
args = commandArgs(trailingOnly = TRUE)
macs2_str = args[1]
union_bed = args[2]

# Load required packages
library(BiocGenerics)
library(ChIPpeakAnno)
library(rtracklayer)

macs2 = unlist(strsplit(macs2_str, " "))
names(macs2) = substr(gsub("^.*lib","lib",macs2),1,6)

granges = lapply(macs2, toGRanges, format = "MACS2.broad")

all.peaks = Reduce(union, granges)

# export as a BED file
rtracklayer::export.bed(all.peaks, con = union_bed)
