#!/usr/bin/env Rscript
#########1#########2#########3#########4#########5#########6#########7#########8

#############################################################################
###   Counts reads overlapping MACS2 consenses peaks as GRanges objects   ###
#############################################################################

args = commandArgs(trailingOnly = TRUE)
bam_dir = args[1]
bam_pattern = args[2]
peaks = args[3]
rse = args[4]
dge = args[5]

## Libraries
library(csaw)
library(edgeR)
library(tidyverse)

## Script-local variables
standard_chr <- paste0("chr", c(1:19)) # only use standard chromosomes
param = readParam(max.frag=1000, pe="both", restrict=standard_chr)

peaks = readRDS(peaks)

bam_list = list.files(path = bam_dir,
                      pattern = bam_pattern,
                      full.names = TRUE)

names(bam_list) = gsub(bam_pattern, "", list.files(path = bam_dir,
                                                   pattern = bam_pattern,
                                                   full.names = FALSE))

unfilt = regionCounts(bam_list, peaks, param = param)

abundance = aveLogCPM(asDGEList(unfilt))

counts = unfilt[abundance > -3, ]

edger_input = asDGEList(counts)

saveRDS(object = counts,
        file = rse)
saveRDS(object = edger_input,
        file = dge)

#!/usr/bin/env Rscript
#########1#########2#########3#########4#########5#########6#########7#########8

#############################################################################
###   Counts reads overlapping MACS2 consenses peaks as GRanges objects   ###
#############################################################################

args = commandArgs(trailingOnly = TRUE)
bam_dir = args[1]
bam_pattern = args[2]
peaks = args[3]
rse = args[4]
dge = args[5]

## Libraries
library(csaw)
library(edgeR)
library(tidyverse)

## Script-local variables
standard_chr <- paste0("chr", c(1:19)) # only use standard chromosomes
param = readParam(max.frag=1000, pe="both", restrict=standard_chr)

peaks = readRDS(peaks)

bam_list = list.files(path = bam_dir,
                      pattern = bam_pattern,
                      full.names = TRUE)

names(bam_list) = gsub(bam_pattern, "", list.files(path = bam_dir,
                                                   pattern = bam_pattern,
                                                   full.names = FALSE))

unfilt = regionCounts(bam_list, peaks, param = param)

abundance = aveLogCPM(asDGEList(unfilt))

counts = unfilt[abundance > -3, ]

edger_input = asDGEList(counts)

saveRDS(object = counts,
        file = rse)
saveRDS(object = edger_input,
        file = dge)

#!/usr/bin/env Rscript
#########1#########2#########3#########4#########5#########6#########7#########8

#############################################################################
###   Counts reads overlapping MACS2 consenses peaks as GRanges objects   ###
#############################################################################

args = commandArgs(trailingOnly = TRUE)
bam_dir = args[1]
bam_pattern = args[2]
peaks = args[3]
rse = args[4]
dge = args[5]

## Libraries
library(csaw)
library(edgeR)
library(tidyverse)

## Script-local variables
standard_chr <- paste0("chr", c(1:19)) # only use standard chromosomes
param = readParam(max.frag=1000, pe="both", restrict=standard_chr)

peaks = readRDS(peaks)

bam_list = list.files(path = bam_dir,
                      pattern = bam_pattern,
                      full.names = TRUE)

names(bam_list) = gsub(bam_pattern, "", list.files(path = bam_dir,
                                                   pattern = bam_pattern,
                                                   full.names = FALSE))

unfilt = regionCounts(bam_list, peaks, param = param)

abundance = aveLogCPM(asDGEList(unfilt))

counts = unfilt[abundance > -3, ]

edger_input = asDGEList(counts)

saveRDS(object = counts,
        file = rse)
saveRDS(object = edger_input,
        file = dge)

#!/usr/bin/env Rscript
#########1#########2#########3#########4#########5#########6#########7#########8

#############################################################################
###   Counts reads overlapping MACS2 consenses peaks as GRanges objects   ###
#############################################################################

args = commandArgs(trailingOnly = TRUE)
bam_dir = args[1]
bam_pattern = args[2]
peaks = args[3]
rse = args[4]
dge = args[5]

## Libraries
library(csaw)
library(edgeR)
library(tidyverse)

## Script-local variables
standard_chr <- paste0("chr", c(1:19)) # only use standard chromosomes
param = readParam(max.frag=1000, pe="both", restrict=standard_chr)

peaks = readRDS(peaks)

bam_list = list.files(path = bam_dir,
                      pattern = bam_pattern,
                      full.names = TRUE)

names(bam_list) = gsub(bam_pattern, "", list.files(path = bam_dir,
                                                   pattern = bam_pattern,
                                                   full.names = FALSE))

unfilt = regionCounts(bam_list, peaks, param = param)

abundance = aveLogCPM(asDGEList(unfilt))

counts = unfilt[abundance > -3, ]

edger_input = asDGEList(counts)

saveRDS(object = counts,
        file = rse)
saveRDS(object = edger_input,
        file = dge)

#!/usr/bin/env Rscript
#########1#########2#########3#########4#########5#########6#########7#########8

#############################################################################
###   Counts reads overlapping MACS2 consenses peaks as GRanges objects   ###
#############################################################################

args = commandArgs(trailingOnly = TRUE)
bam_dir = args[1]
bam_pattern = args[2]
peaks = args[3]
rse = args[4]
dge = args[5]

## Libraries
library(csaw)
library(edgeR)
library(tidyverse)

## Script-local variables
standard_chr <- paste0("chr", c(1:19)) # only use standard chromosomes
param = readParam(max.frag=1000, pe="both", restrict=standard_chr)

peaks = readRDS(peaks)

bam_list = list.files(path = bam_dir,
                      pattern = bam_pattern,
                      full.names = TRUE)

names(bam_list) = gsub(bam_pattern, "", list.files(path = bam_dir,
                                                   pattern = bam_pattern,
                                                   full.names = FALSE))

unfilt = regionCounts(bam_list, peaks, param = param)

abundance = aveLogCPM(asDGEList(unfilt))

counts = unfilt[abundance > -3, ]

edger_input = asDGEList(counts)

saveRDS(object = counts,
        file = rse)
saveRDS(object = edger_input,
        file = dge)

#!/usr/bin/env Rscript
#########1#########2#########3#########4#########5#########6#########7#########8

#############################################################################
###   Counts reads overlapping MACS2 consenses peaks as GRanges objects   ###
#############################################################################

args = commandArgs(trailingOnly = TRUE)
bam_dir = args[1]
bam_pattern = args[2]
peaks = args[3]
rse = args[4]
dge = args[5]

## Libraries
library(csaw)
library(edgeR)
library(tidyverse)

## Script-local variables
standard_chr <- paste0("chr", c(1:19)) # only use standard chromosomes
param = readParam(max.frag=1000, pe="both", restrict=standard_chr)

peaks = readRDS(peaks)

bam_list = list.files(path = bam_dir,
                      pattern = bam_pattern,
                      full.names = TRUE)

names(bam_list) = gsub(bam_pattern, "", list.files(path = bam_dir,
                                                   pattern = bam_pattern,
                                                   full.names = FALSE))

unfilt = regionCounts(bam_list, peaks, param = param)

abundance = aveLogCPM(asDGEList(unfilt))

counts = unfilt[abundance > -3, ]

edger_input = asDGEList(counts)

saveRDS(object = counts,
        file = rse)
saveRDS(object = edger_input,
        file = dge)
