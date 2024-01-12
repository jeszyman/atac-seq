
## - [[file:./scripts/make_background_bins.R][Base script]]

#!/usr/bin/env Rscript
#########1#########2#########3#########4#########5#########6#########7#########8
###
###   R Script to make background bins for TMM normalization   ###
###

args = commandArgs(trailingOnly = TRUE)
bam_dir = args[1]
bam_pattern = args[2]
rse_file = args[3]

library(csaw)

standard_chr <- paste0("chr", c(1:19)) # only use standard chromosomes
param = readParam(max.frag=1000, pe="both", restrict=standard_chr)


bam_list = list.files(path = bam_dir,
                  pattern = bam_pattern,
                  full.names = TRUE)

names(bam_list) = gsub(bam_pattern, "", list.files(path = bam_dir,
                                                   pattern = bam_pattern,
                                                   full.names = FALSE))


binned = windowCounts(bam_list, bin=TRUE, width=10000, param=param)

saveRDS(object = binned,
        file = rse_file)
