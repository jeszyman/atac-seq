#!/usr/bin/env Rscript

# Command line arguements
args = commandArgs(trailingOnly = TRUE)
design_rds = args[1]
fit_rds = args[2]
contrast_string = args[3]
res_tsv = args[4]

# Load required packages, data, and functions
library(edgeR)
library(tidyverse)

design = readRDS(design_rds)
fit = readRDS(fit_rds)

contrast <- makeContrasts(eval(parse(text = contrast_string)), levels=design)
lrt = glmLRT(fit, contrast = contrast)

res =
  as.data.frame(topTags(lrt, n = Inf)) %>%
  rownames_to_column(var = "coordinate") %>%
  as_tibble()

write_tsv(res, file = res_tsv)
