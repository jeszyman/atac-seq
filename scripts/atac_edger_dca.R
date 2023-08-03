#!/usr/bin/env Rscript

# Command line arguements
args = commandArgs(trailingOnly = TRUE)
design_rds = args[1]
fit_rds = args[2]
annotation_tsv = args[3]
cohorts_str = args[4]
res_tsv = args[5]

# Load required packages, data, and functions
library(edgeR)
library(tidyverse)

design = readRDS(design_rds)
fit = readRDS(fit_rds)
annotation = read_tsv(annotation_tsv)

cohorts_vec = strsplit(cohorts_str, " ")[[1]]
contrast_string <- paste(cohorts_vec[[1]], "-", cohorts_vec[[2]])

contrast <- makeContrasts(eval(parse(text = contrast_string)), levels=design)

qlf = glmQLFTest(fit, contrast = contrast)

res =
  as.data.frame(topTags(qlf, n = Inf)) %>%
  rownames_to_column(var = "ensembl_gene_id") %>%
  as_tibble() %>%
  left_join(annotation, by = "ensembl_gene_id") %>%
  mutate(sign = sign(logFC)) %>%
  mutate(score = sign * -log10(PValue)) %>%
  mutate(rank = rank(-score, ties.method = "random"))

write_tsv(res, file = res_tsv)
