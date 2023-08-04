#!/usr/bin/env Rscript

##############################
###   Human Dca With Rvu   ###
##############################

# Command line arguements
args = commandArgs(trailingOnly = TRUE)
counts_tsv = args[1]
libraries_full_rds = args[2]
design_rds = args[3]
ruv_k = args[4]
ruv_counts_rds = args[5]
fit_rds = args[6]

# Load required packages, data, and functions
library(RUVSeq)
library(tidyverse)
libraries_full = readRDS(libraries_full_rds)
counts = read_tsv(counts_tsv)
design = readRDS(design_rds)

# Setup data objects
mat = as.matrix(counts[,-1])
row.names(mat) = counts$coordinate
colnames(mat) = substr(colnames(mat),1,6)

model_df = as.data.frame(design)
model_df$'(Intercept)' <- NULL

set = newSeqExpressionSet(mat,
                          phenoData = AnnotatedDataFrame(model_df))

dge <- DGEList(counts = counts(set))
y <- DGEList(counts = counts(set))
y <- calcNormFactors(y, method = "upperquartile")
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
fit <- glmFit(y, design)
dev <- residuals(fit, type="deviance")
peaks = row.names(dev)
adjust <- RUVr(set, peaks, k = as.numeric(ruv_k), dev)
adjust_counts = normCounts(adjust)

saveRDS(adjust_counts, ruv_counts_rds)

adjust_counts = adjust_counts + 1
dge = DGEList(counts = adjust_counts)
dge = calcNormFactors(dge)
dge = estimateGLMCommonDisp(dge, design)
dge = estimateGLMTrendedDisp(dge, design)
dge = estimateGLMTagwiseDisp(dge, design)
fit = glmFit(dge, design)

saveRDS(fit, fit_rds)
