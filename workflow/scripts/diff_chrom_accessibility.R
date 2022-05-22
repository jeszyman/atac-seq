#########1#########2#########3#########4#########5#########6#########7#########8
###                                                                          ###
###   Script to generate differential accessibility model with EdgeR         ###                
###                                                                          ###
#########1#########2#########3#########4#########5#########6#########7#########8

# Setup

## Arguements for testing
dge_rds = "~/repos/atac-seq/test/csaw/dge_regfilt.rds"
groups_str = "ir48h ir48h sham sham"
contrast = "ir48h-sham"
norm_counts_rds = "~/repos/atac-seq/test/csaw/norm_counts_rse_regfilt.rds"
dca_granges_rds = "/tmp/test.rds"

args = commandArgs(trailingOnly = TRUE)
dge_rds = args[1]
groups_str =args [2]
contrast = args[3]
norm_counts_rds = args[4]
dca_granges_rds = args[5]

library(csaw)
library(edgeR)
library(tidyverse)

y = readRDS(dge_rds)
groups = as.factor(unlist(strsplit(groups_str, " ")))

design = model.matrix(~0 + groups, data=y$samples)
colnames(design) = levels(groups)
y = estimateDisp(y, design)
fit = glmQLFit(y, design, robust=TRUE)
results = glmQLFTest(fit, contrast=makeContrasts(contrast, levels=design))

# combine GRanges rowdata with DA statistics
counts = readRDS(norm_counts_rds)
rowData(counts) = cbind(rowData(counts), results$table)

# merge nearby windows
# up to "tol" distance apart: 500 bp in this case
# max merged window width: 5000 bp
merged.peaks <- mergeWindows(rowRanges(counts), tol=500L, max.width=5000L)

# use most significant window as statistical representation for p-value and FDR for merged windows
tab.best = getBestTest(merged.peaks$id, results$table)

# combine merged peaks window range with statistics
final.merged.peaks = merged.peaks$region
final.merged.peaks@elementMetadata = cbind(final.merged.peaks@elementMetadata, tab.best[,-1])
final.merged.peaks = final.merged.peaks[order(final.merged.peaks@elementMetadata$FDR), ] # sort by FDR

saveRDS(object = final.merged.peaks, 
        file = dca_granges_rds)
