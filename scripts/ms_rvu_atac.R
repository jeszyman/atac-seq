#!/usr/bin/env Rscript

##############################
###   Mouse Dca With Rvu   ###
##############################

# Command line arguements
args = commandArgs(trailingOnly = TRUE)
libraries_full_rds = args[1]
counts_tsv = args[2]
keep_tsv = args[3]
ir2d_sham_tsv = args[4]
ir6w_sham_tsv = args[5]
ir6w_ir2d_tsv = args[6]
rvu_k = args[7]

# Load required packages, data, and functions
library(RUVSeq)
library(tidyverse)
libraries_full = readRDS(libraries_full_rds)
counts = read_tsv(counts_tsv)
keep = read_tsv(keep_tsv)

# Setup data objects
mat = as.matrix(counts[,-1])
row.names(mat) = counts$coordinate
colnames(mat) = substr(colnames(mat),1,6)
mat = mat[, colnames(mat) %in% keep$library]
(group = data.frame(library = colnames(mat)) %>% left_join(libraries_full) %>% pull(cohort) %>% droplevels())
set = newSeqExpressionSet(mat,
                          phenoData = data.frame(group, row.names = colnames(mat)))

print("Setup complete")

# Perform RUV adjustment by k unwanted factors
design <- model.matrix(~0+group, data = pData(set))
y <- DGEList(counts = counts(set), group = group)
y <- calcNormFactors(y, method = "upperquartile")
#y <- calcNormFactors(y)
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
fit <- glmFit(y, design)
dev <- residuals(fit, type="deviance")
peaks = row.names(dev)
adjust <- RUVr(set, peaks, k = as.numeric(rvu_k), dev)

print("\nModified counts created\n")

# Perform DCA by EdgeR on adjusted gene expression
if (rvu_k != 0) {
  modcounts = normCounts(adjust)
} else {
  modcounts = mat
}

design = model.matrix(~ 0 + group)
dge = DGEList(counts = modcounts)
dge = calcNormFactors(dge)
dge = estimateGLMCommonDisp(dge, design)
dge = estimateGLMTrendedDisp(dge, design)
dge = estimateGLMTagwiseDisp(dge, design)
fit = glmFit(dge, design)

name_mapping= data.frame(newname = c("logfc", "logcpm", "lr", "pval", "qval"),
                 oldname = c("logFC", "logCPM", "LR", "PValue", "FDR"))

print("\nDGE performed with modified counts\n")

# Generate specific comparisons
ir2d_sham_cont = makeContrasts(groupir2d - groupsham, levels = colnames(design))
ir2d_sham_lrt = glmLRT(fit, contrast = ir2d_sham_cont)
ir2d_sham = data.frame(topTags(ir2d_sham_lrt, n = nrow(dge))$table) %>% rownames_to_column(var = "coordinate") %>% as_tibble() %>%
  rename(!!!setNames(as.list(name_mapping$oldname), name_mapping$newname)) %>%
  write_tsv(ir2d_sham_tsv)

ir6w_sham_cont = makeContrasts(groupir6w - groupsham, levels = colnames(design))
ir6w_sham_lrt = glmLRT(fit, contrast = ir6w_sham_cont)
ir6w_sham = data.frame(topTags(ir6w_sham_lrt, n = nrow(dge))$table) %>% rownames_to_column(var = "coordinate") %>% as_tibble() %>%
  rename(!!!setNames(as.list(name_mapping$oldname), name_mapping$newname)) %>%
  write_tsv(ir6w_sham_tsv)

ir6w_ir2d_cont = makeContrasts(groupir6w - groupir2d, levels = colnames(design))
ir6w_ir2d_lrt = glmLRT(fit, contrast = ir6w_ir2d_cont)
ir6w_ir2d = data.frame(topTags(ir6w_ir2d_lrt, n = nrow(dge))$table) %>% rownames_to_column(var = "coordinate") %>% as_tibble() %>%
  rename(!!!setNames(as.list(name_mapping$oldname), name_mapping$newname)) %>%
  write_tsv(ir6w_ir2d_tsv)
