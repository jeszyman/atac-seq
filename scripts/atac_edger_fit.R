#!/usr/bin/env Rscript

# Command line arguements
args = commandArgs(trailingOnly = TRUE)
counts_tsv = args[1]
design_rds = args[2]
libs_rds = args[3]
dge_rds = args[4]
fit_rds = args[5]

library(tidyverse)
library(edgeR)

counts = read_tsv(counts_tsv)
design = readRDS(design_rds)

mat = as.matrix(counts[,-1])
row.names(mat) = counts$coordinate
colnames(mat) = substr(colnames(mat),1,6)

dge = DGEList(counts = mat)
keep = filterByExpr(dge, design)
dge = dge[keep,]
dge = calcNormFactors(dge)
dge = estimateGLMCommonDisp(dge, design)
dge = estimateGLMTrendedDisp(dge, design)
dge = estimateGLMTagwiseDisp(dge, design)
fit = glmFit(dge, design)

saveRDS(dge, dge_rds)
saveRDS(fit, fit_rds)
