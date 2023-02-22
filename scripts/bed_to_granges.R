
#!/usr/bin/env Rscript
#########1#########2#########3#########4#########5#########6#########7#########8

##################################################
###    Converts BED files to GRanges objects   ###
##################################################

args = commandArgs(trailingOnly = TRUE)
bed = args[1]
granges_file = args[2]

library(GenomicRanges)

peaks = read.table(bed, sep = "\t")[,1:3]

colnames(peaks) = c("chrom", "start", "end")

granges = GRanges(peaks)

saveRDS(object = granges, file = granges_file)
