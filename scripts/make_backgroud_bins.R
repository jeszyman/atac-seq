
#!/usr/bin/env Rscript
#########1#########2#########3#########4#########5#########6#########7#########8
###
###   Script to make background bins for csaw TMM normalization   ###
###

# Setup
##
## Snakemake
args = commandArgs(trailingOnly = TRUE)
window_file =  args[1]
bam_dir = args[2]
bam_pattern = args[3]
filt_libs_str = args[4]
rse = args[4]
bk = args[5]

filt_libs = unlist(strsplit(filt_libs_str, " "))

## Libraries
library(csaw)
library(edgeR)
library(tidyverse)

## Script-local variables
surrounds = 2000
standard_chr <- paste0("chr", c(1:19)) # only use standard chromosomes
param = readParam(max.frag=1000, pe="both", restrict=standard_chr)

bam_list = list.files(path = bam_dir,
                  pattern = bam_pattern,
                  full.names = TRUE)

names(bam_list) = gsub(bam_pattern, "", list.files(path = bam_dir,
                                                   pattern = bam_pattern,
                                                   full.names = FALSE))

bam_list = bam_list[names(bam_list) %in% filt_libs]

binned = windowCounts(bam_list, bin=TRUE, width=10000, param=param)

#!/usr/bin/env Rscript
#########1#########2#########3#########4#########5#########6#########7#########8
###
###   Script to make background bins for csaw TMM normalization   ###
###

# Setup
##
## Snakemake
args = commandArgs(trailingOnly = TRUE)
window_file =  args[1]
bam_dir = args[2]
bam_pattern = args[3]
filt_libs_str = args[4]
rse = args[4]
bk = args[5]

filt_libs = unlist(strsplit(filt_libs_str, " "))

## Libraries
library(csaw)
library(edgeR)
library(tidyverse)

## Script-local variables
surrounds = 2000
standard_chr <- paste0("chr", c(1:19)) # only use standard chromosomes
param = readParam(max.frag=1000, pe="both", restrict=standard_chr)

bam_list = list.files(path = bam_dir,
                  pattern = bam_pattern,
                  full.names = TRUE)

names(bam_list) = gsub(bam_pattern, "", list.files(path = bam_dir,
                                                   pattern = bam_pattern,
                                                   full.names = FALSE))

bam_list = bam_list[names(bam_list) %in% filt_libs]

binned = windowCounts(bam_list, bin=TRUE, width=10000, param=param)
