#!/usr/bin/env Rscript
#########1#########2#########3#########4#########5#########6#########7#########8
###
###   Script to generate differential accessibility model with EdgeR   ###                
###

args = commandArgs(trailingOnly = TRUE)
counts_rds = args[1]
background_rds = args[2]
data_model = args[3]
dca_granges_file = args[4]

library(csaw)
library(edgeR)
library(tidyverse)

counts = readRDS(counts_rds)
load(data_model)
background = readRDS(background_rds)

counts = normFactors(background, se.out = counts)

y <- asDGEList(counts)
colnames(y$counts) <- colnames(counts)
rownames(y$samples) <- colnames(counts)

groups =
  data.frame(library_id = rownames(y$samples)) %>%
  left_join(libraries_full, by = "library_id") %>%  
  droplevels() %>%
  pull(cohort_id)
groups = fct_relevel(groups, "sham", "ir48h") 

y$samples$group = groups

design <- model.matrix(~0 + groups, data=y$samples)

colnames(design) = levels(groups)

# stabilize dispersion estimates with empirical bayes
y <- estimateDisp(y, design)

fit <- glmQLFit(y, design, robust=TRUE)

# testing for differentially-accessible windows
results <- glmQLFTest(fit, contrast=makeContrasts(ir48h-sham, levels=design))


# combine GRanges rowdata with DA statistics
rowData(working.windows) <- cbind(rowData(working.windows), results$table)
#working.windows@rowRanges

# merge nearby windows
# up to "tol" distance apart: 500 bp in this case
# max merged window width: 5000 bp
merged.peaks <- mergeWindows(rowRanges(working.windows), tol=500L, max.width=5000L)
merged.peaks <- mergeWindows(rowRanges(filtered_counts), tol=500L, max.width=5000L)

# summary(width(merged.peaks$region))
# should merge some peaks; change as desired

# use most significant window as statistical representation for p-value and FDR for merged windows
tab.best <- getBestTest(merged.peaks$id, results$table)


# combine merged peaks window range with statistics
final.merged.peaks <- merged.peaks$region
final.merged.peaks@elementMetadata <- cbind(final.merged.peaks@elementMetadata, tab.best[,-1])
final.merged.peaks <- final.merged.peaks[order(final.merged.peaks@elementMetadata$FDR), ] # sort by FDR

saveRDS(object = final.merged.peaks, 
        file = dca_grange_file)

#!/usr/bin/env Rscript
#########1#########2#########3#########4#########5#########6#########7#########8
###
###   Script to generate differential accessibility model with EdgeR   ###                
###

args = commandArgs(trailingOnly = TRUE)
counts_rds = args[1]
background_rds = args[2]
data_model = args[3]
dca_granges_file = args[4]

library(csaw)
library(edgeR)
library(tidyverse)

counts = readRDS(counts_rds)
load(data_model)
background = readRDS(background_rds)

counts = normFactors(background, se.out = counts)

y <- asDGEList(counts)
colnames(y$counts) <- colnames(counts)
rownames(y$samples) <- colnames(counts)

groups =
  data.frame(library_id = rownames(y$samples)) %>%
  left_join(libraries_full, by = "library_id") %>%  
  droplevels() %>%
  pull(cohort_id)
groups = fct_relevel(groups, "sham", "ir48h") 

y$samples$group = groups

design <- model.matrix(~0 + groups, data=y$samples)

colnames(design) = levels(groups)

# stabilize dispersion estimates with empirical bayes
y <- estimateDisp(y, design)

fit <- glmQLFit(y, design, robust=TRUE)

# testing for differentially-accessible windows
results <- glmQLFTest(fit, contrast=makeContrasts(ir48h-sham, levels=design))


# combine GRanges rowdata with DA statistics
rowData(working.windows) <- cbind(rowData(working.windows), results$table)
#working.windows@rowRanges

# merge nearby windows
# up to "tol" distance apart: 500 bp in this case
# max merged window width: 5000 bp
merged.peaks <- mergeWindows(rowRanges(working.windows), tol=500L, max.width=5000L)
merged.peaks <- mergeWindows(rowRanges(filtered_counts), tol=500L, max.width=5000L)

# summary(width(merged.peaks$region))
# should merge some peaks; change as desired

# use most significant window as statistical representation for p-value and FDR for merged windows
tab.best <- getBestTest(merged.peaks$id, results$table)


# combine merged peaks window range with statistics
final.merged.peaks <- merged.peaks$region
final.merged.peaks@elementMetadata <- cbind(final.merged.peaks@elementMetadata, tab.best[,-1])
final.merged.peaks <- final.merged.peaks[order(final.merged.peaks@elementMetadata$FDR), ] # sort by FDR

saveRDS(object = final.merged.peaks, 
        file = dca_grange_file)

#!/usr/bin/env Rscript
#########1#########2#########3#########4#########5#########6#########7#########8
###
###   Script to generate differential accessibility model with EdgeR   ###                
###

args = commandArgs(trailingOnly = TRUE)
counts_rds = args[1]
background_rds = args[2]
data_model = args[3]
dca_granges_file = args[4]

library(csaw)
library(edgeR)
library(tidyverse)

counts = readRDS(counts_rds)
load(data_model)
background = readRDS(background_rds)

counts = normFactors(background, se.out = counts)

y <- asDGEList(counts)
colnames(y$counts) <- colnames(counts)
rownames(y$samples) <- colnames(counts)

groups =
  data.frame(library_id = rownames(y$samples)) %>%
  left_join(libraries_full, by = "library_id") %>%  
  droplevels() %>%
  pull(cohort_id)
groups = fct_relevel(groups, "sham", "ir48h") 

y$samples$group = groups

design <- model.matrix(~0 + groups, data=y$samples)

colnames(design) = levels(groups)

# stabilize dispersion estimates with empirical bayes
y <- estimateDisp(y, design)

fit <- glmQLFit(y, design, robust=TRUE)

# testing for differentially-accessible windows
results <- glmQLFTest(fit, contrast=makeContrasts(ir48h-sham, levels=design))


# combine GRanges rowdata with DA statistics
rowData(working.windows) <- cbind(rowData(working.windows), results$table)
#working.windows@rowRanges

# merge nearby windows
# up to "tol" distance apart: 500 bp in this case
# max merged window width: 5000 bp
merged.peaks <- mergeWindows(rowRanges(working.windows), tol=500L, max.width=5000L)
merged.peaks <- mergeWindows(rowRanges(filtered_counts), tol=500L, max.width=5000L)

# summary(width(merged.peaks$region))
# should merge some peaks; change as desired

# use most significant window as statistical representation for p-value and FDR for merged windows
tab.best <- getBestTest(merged.peaks$id, results$table)


# combine merged peaks window range with statistics
final.merged.peaks <- merged.peaks$region
final.merged.peaks@elementMetadata <- cbind(final.merged.peaks@elementMetadata, tab.best[,-1])
final.merged.peaks <- final.merged.peaks[order(final.merged.peaks@elementMetadata$FDR), ] # sort by FDR

saveRDS(object = final.merged.peaks, 
        file = dca_grange_file)

#!/usr/bin/env Rscript
#########1#########2#########3#########4#########5#########6#########7#########8
###
###   Script to generate differential accessibility model with EdgeR   ###                
###

args = commandArgs(trailingOnly = TRUE)
counts_rds = args[1]
background_rds = args[2]
data_model = args[3]
dca_granges_file = args[4]

library(csaw)
library(edgeR)
library(tidyverse)

counts = readRDS(counts_rds)
load(data_model)
background = readRDS(background_rds)

counts = normFactors(background, se.out = counts)

y <- asDGEList(counts)
colnames(y$counts) <- colnames(counts)
rownames(y$samples) <- colnames(counts)

groups =
  data.frame(library_id = rownames(y$samples)) %>%
  left_join(libraries_full, by = "library_id") %>%  
  droplevels() %>%
  pull(cohort_id)
groups = fct_relevel(groups, "sham", "ir48h") 

y$samples$group = groups

design <- model.matrix(~0 + groups, data=y$samples)

colnames(design) = levels(groups)

# stabilize dispersion estimates with empirical bayes
y <- estimateDisp(y, design)

fit <- glmQLFit(y, design, robust=TRUE)

# testing for differentially-accessible windows
results <- glmQLFTest(fit, contrast=makeContrasts(ir48h-sham, levels=design))


# combine GRanges rowdata with DA statistics
rowData(working.windows) <- cbind(rowData(working.windows), results$table)
#working.windows@rowRanges

# merge nearby windows
# up to "tol" distance apart: 500 bp in this case
# max merged window width: 5000 bp
merged.peaks <- mergeWindows(rowRanges(working.windows), tol=500L, max.width=5000L)
merged.peaks <- mergeWindows(rowRanges(filtered_counts), tol=500L, max.width=5000L)

# summary(width(merged.peaks$region))
# should merge some peaks; change as desired

# use most significant window as statistical representation for p-value and FDR for merged windows
tab.best <- getBestTest(merged.peaks$id, results$table)


# combine merged peaks window range with statistics
final.merged.peaks <- merged.peaks$region
final.merged.peaks@elementMetadata <- cbind(final.merged.peaks@elementMetadata, tab.best[,-1])
final.merged.peaks <- final.merged.peaks[order(final.merged.peaks@elementMetadata$FDR), ] # sort by FDR

saveRDS(object = final.merged.peaks, 
        file = dca_grange_file)
