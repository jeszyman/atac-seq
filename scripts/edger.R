#!/usr/bin/env Rscript
#########1#########2#########3#########4#########5#########6#########7#########8
###
###   Do differential expression of ATAC-seq peaks through edgeR   ###
###

args = commandArgs(trailingOnly = TRUE)
= args[1]

library(csaw)
library(DESeq2)
library(edgeR)
library(tidyverse)

# Load counts as DGE list
counts = readRDS(input)

counts = readRDS("/mnt/ris/jschwarz/cardiac-radiobiology/atac/norm/macs2_all_union_open_narrow_tmm_rse.rds")
load("/mnt/ris/jschwarz/cardiac-radiobiology/data_model/data_model.RData")

# setup design matrix
# see edgeR manual for more information
y <- asDGEList(counts)
colnames(y$counts) =
rownames(y$samples) = names(counts$bam.files)
groups =
  data.frame(library_id = names(counts$bam.files)) %>%
  left_join(libraries_full, by = "library_id") %>%
  droplevels() %>%
  pull(cohort_id)
groups = fct_relevel(groups, "sham", "ir48h")
y$samples$group = groups
colors = as.character(factor(y$samples$group, levels = c("sham", "ir48h", "ir6w"), labels = c("darkgreen", "red", "blue")))

plotMDS(y, col = colors, gene.selection = "common", top = 1000000)

test

test = subset(counts, select = !(colnames(counts) %in% c("lib001", "lib002","lib007","lib005", "lib006", "lib003", "lib004")))
test = subset(counts, select = !(colnames(counts) %in% c("lib001", "lib002","lib007","lib005", "lib006", "lib003", "lib004","lib013","lib018")))
counts = test

# setup design matrix
# see edgeR manual for more information
y <- asDGEList(counts)
colnames(y$counts) = rownames(y$samples) = names(counts$bam.files)
groups =
  data.frame(library_id = names(counts$bam.files)) %>%
  left_join(libraries_full, by = "library_id") %>%
  droplevels() %>%
  pull(cohort_id)
groups = fct_relevel(groups, "sham", "ir48h")
groups
y$samples$group = groups
colors = as.character(factor(y$samples$group, levels = c("sham", "ir48h", "ir6w"), labels = c("darkgreen", "red", "blue")))

pdf("/tmp/pca.pdf")
plotMDS(y, col = colors, gene.selection = "common", top = 80)
dev.off()

plotMDS(y, col = colors, top = 100)

design <- model.matrix(~group, data=y$samples)
colnames(design) = levels(groups)


# stabilize dispersion estimates with empirical bayes
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design, robust=TRUE)

summary(fit$df.prior)

fit <- glmQLFit(y, design)

class(design)
# testing for differentially-accessible windows
results <- glmQLFTest(fit, contrast=makeContrasts(sham-ir6w, levels=design))
# head(results$table)

topTags(results)

# combine GRanges rowdata with DA statistics
rowData(counts) <- cbind(rowData(counts), results$table)

res = as.data.frame(topTags(results, n = Inf))

ggplot(res, aes(x = logFC)) + geom_density()
summary(as.data.frame(topTags(results, n = Inf))$FDR)

test = as_tibble(as.data.frame(topTags(results, n = Inf)))

max(test$FDR)

summary(results$table$PValue)

head(results$table$PValue)

fit = glmFit(y, design, contrast = makeContrasts(ir48h-sham, levels = design))

fit
lrt = glmLRT(fit, contrast = makeContrasts(ir48h-sham, levels = design))
test=as.data.frame(topTags(lrt, n = 10000))
class(test)
summary(test$FDR)
lrt
head(lrt$table)
et = exactTest(y)
topTags(et)

# merge nearby windows
# up to "tol" distance apart: 500 bp in this case
# max merged window width: 5000 bp
merged.peaks <- mergeWindows(rowRanges(counts), tol=500L, max.width=5000L)
# summary(width(merged.peaks$region))
# should merge some peaks; change as desired

# use most significant window as statistical representation for p-value and FDR for merged windows
tab.best <- getBestTest(merged.peaks$id, results$table)
head(tab.best)
min(tab.best$PValue)
min(tab.best$FDR)

# combine merged peaks window range with statistics
final.merged.peaks <- merged.peaks$region
final.merged.peaks@elementMetadata <- cbind(final.merged.peaks@elementMetadata, tab.best[,-1])
final.merged.peaks <- final.merged.peaks[order(final.merged.peaks@elementMetadata$FDR), ] # sort by FDR
final.merged.peaks # all windows

# filter by FDR threshold
FDR.thresh <- 0.05 # set as desired
final.merged.peaks.sig <- final.merged.peaks[final.merged.peaks@elementMetadata$FDR < FDR.thresh, ]
final.merged.peaks.sig # significant differentially-accessible windows




colnames(design) = levels(counts$samples$group)

test = rlog(assays(counts)$counts)
rld = test

class(rld)
mat = t(rld)
pca = prcomp(mat)

summary(pca)

head(counts$counts)
rownames(counts$counts)

class(working.windows)

working.windows

# stabilize dispersion estimates with empirical bayes
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design, robust=TRUE)

# testing for differentially-accessible windows
#results <- glmQLFTest(fit, contrast=makeContrasts(treat-control, levels=design))

results <- glmQLFTest(fit, contrast=makeContrasts(ir48h-sham, levels=design))
# head(results$table)

test = results$table
min(test$PValue)

class(working.windows)

test = working.windows[,8:15]


# combine GRanges rowdata with DA statistics
#rowData(working.windows) <- cbind(rowData(working.windows), results$table)
rowData(test) = cbind(rowData(test), results$table)

test@rowRanges
working.windows = test

# merge nearby windows
# up to "tol" distance apart: 500 bp in this case
# max merged window width: 5000 bp
merged.peaks <- mergeWindows(rowRanges(working.windows), tol=500L, max.width=5000L)
summary(width(merged.peaks$region))
# should merge some peaks; change as desired

# use most significant window as statistical representation for p-value and FDR for merged windows
tab.best <- getBestTest(merged.peaks$id, results$table)
head(tab.best)
# combine merged peaks window range with statistics
final.merged.peaks <- merged.peaks$region
final.merged.peaks@elementMetadata <- cbind(final.merged.peaks@elementMetadata, tab.best[,-1])
final.merged.peaks <- final.merged.peaks[order(final.merged.peaks@elementMetadata$FDR), ] # sort by FDR
final.merged.peaks # all windows

# filter by FDR threshold
#FDR.thresh <- 0.05 # set as desired
#final.merged.peaks.sig <- final.merged.peaks[final.merged.peaks@elementMetadata$FDR < FDR.thresh, ]
#final.merged.peaks.sig # significant differentially-accessible windows



#########1#########2#########3#########4#########5#########6#########7#########8

library(DESeq2)

test = subset(counts, select = !(colnames(counts) %in% c("lib001", "lib002","lib007","lib005", "lib006", "lib003", "lib004","lib013","lib018", "lib023", "lib014")))
counts = test


test = rlog(assays(counts)$counts)
rld = test

rld = vst(assays(counts)$counts)
mat = t(rld)
pca = prcomp(mat)

summary(pca)

pca_plot = as.data.frame(pca$x) %>%
  rownames_to_column(var = "library_id") %>%
  left_join(libraries_full, by = "library_id") %>%
  ggplot(., aes(x = PC1, y = PC2, color = cohort_id)) +
  geom_point(size = 4)
pca_plot



#lowdate = as.character(data.frame(library_id = colnames(y)) %>% left_join(libraries_full, by = "library_id") %>% pull(flow_date))

#########1#########2#########3#########4#########5#########6#########7#########8

#!/usr/bin/env Rscript
#########1#########2#########3#########4#########5#########6#########7#########8
###
###   Do differential expression of ATAC-seq peaks through edgeR   ###
###

args = commandArgs(trailingOnly = TRUE)
= args[1]

library(csaw)
library(DESeq2)
library(edgeR)
library(tidyverse)

# Load counts as DGE list
counts = readRDS(input)

counts = readRDS("/mnt/ris/jschwarz/cardiac-radiobiology/atac/norm/macs2_all_union_open_narrow_tmm_rse.rds")
load("/mnt/ris/jschwarz/cardiac-radiobiology/data_model/data_model.RData")

# setup design matrix
# see edgeR manual for more information
y <- asDGEList(counts)
colnames(y$counts) =
rownames(y$samples) = names(counts$bam.files)
groups =
  data.frame(library_id = names(counts$bam.files)) %>%
  left_join(libraries_full, by = "library_id") %>%
  droplevels() %>%
  pull(cohort_id)
groups = fct_relevel(groups, "sham", "ir48h")
y$samples$group = groups
colors = as.character(factor(y$samples$group, levels = c("sham", "ir48h", "ir6w"), labels = c("darkgreen", "red", "blue")))

plotMDS(y, col = colors, gene.selection = "common", top = 1000000)

test

test = subset(counts, select = !(colnames(counts) %in% c("lib001", "lib002","lib007","lib005", "lib006", "lib003", "lib004")))
test = subset(counts, select = !(colnames(counts) %in% c("lib001", "lib002","lib007","lib005", "lib006", "lib003", "lib004","lib013","lib018")))
counts = test

# setup design matrix
# see edgeR manual for more information
y <- asDGEList(counts)
colnames(y$counts) = rownames(y$samples) = names(counts$bam.files)
groups =
  data.frame(library_id = names(counts$bam.files)) %>%
  left_join(libraries_full, by = "library_id") %>%
  droplevels() %>%
  pull(cohort_id)
groups = fct_relevel(groups, "sham", "ir48h")
groups
y$samples$group = groups
colors = as.character(factor(y$samples$group, levels = c("sham", "ir48h", "ir6w"), labels = c("darkgreen", "red", "blue")))

pdf("/tmp/pca.pdf")
plotMDS(y, col = colors, gene.selection = "common", top = 80)
dev.off()

plotMDS(y, col = colors, top = 100)

design <- model.matrix(~group, data=y$samples)
colnames(design) = levels(groups)


# stabilize dispersion estimates with empirical bayes
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design, robust=TRUE)

summary(fit$df.prior)

fit <- glmQLFit(y, design)

class(design)
# testing for differentially-accessible windows
results <- glmQLFTest(fit, contrast=makeContrasts(sham-ir6w, levels=design))
# head(results$table)

topTags(results)

# combine GRanges rowdata with DA statistics
rowData(counts) <- cbind(rowData(counts), results$table)

res = as.data.frame(topTags(results, n = Inf))

ggplot(res, aes(x = logFC)) + geom_density()
summary(as.data.frame(topTags(results, n = Inf))$FDR)

test = as_tibble(as.data.frame(topTags(results, n = Inf)))

max(test$FDR)

summary(results$table$PValue)

head(results$table$PValue)

fit = glmFit(y, design, contrast = makeContrasts(ir48h-sham, levels = design))

fit
lrt = glmLRT(fit, contrast = makeContrasts(ir48h-sham, levels = design))
test=as.data.frame(topTags(lrt, n = 10000))
class(test)
summary(test$FDR)
lrt
head(lrt$table)
et = exactTest(y)
topTags(et)

# merge nearby windows
# up to "tol" distance apart: 500 bp in this case
# max merged window width: 5000 bp
merged.peaks <- mergeWindows(rowRanges(counts), tol=500L, max.width=5000L)
# summary(width(merged.peaks$region))
# should merge some peaks; change as desired

# use most significant window as statistical representation for p-value and FDR for merged windows
tab.best <- getBestTest(merged.peaks$id, results$table)
head(tab.best)
min(tab.best$PValue)
min(tab.best$FDR)

# combine merged peaks window range with statistics
final.merged.peaks <- merged.peaks$region
final.merged.peaks@elementMetadata <- cbind(final.merged.peaks@elementMetadata, tab.best[,-1])
final.merged.peaks <- final.merged.peaks[order(final.merged.peaks@elementMetadata$FDR), ] # sort by FDR
final.merged.peaks # all windows

# filter by FDR threshold
FDR.thresh <- 0.05 # set as desired
final.merged.peaks.sig <- final.merged.peaks[final.merged.peaks@elementMetadata$FDR < FDR.thresh, ]
final.merged.peaks.sig # significant differentially-accessible windows




colnames(design) = levels(counts$samples$group)

test = rlog(assays(counts)$counts)
rld = test

class(rld)
mat = t(rld)
pca = prcomp(mat)

summary(pca)

head(counts$counts)
rownames(counts$counts)

class(working.windows)

working.windows

# stabilize dispersion estimates with empirical bayes
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design, robust=TRUE)

# testing for differentially-accessible windows
#results <- glmQLFTest(fit, contrast=makeContrasts(treat-control, levels=design))

results <- glmQLFTest(fit, contrast=makeContrasts(ir48h-sham, levels=design))
# head(results$table)

test = results$table
min(test$PValue)

class(working.windows)

test = working.windows[,8:15]


# combine GRanges rowdata with DA statistics
#rowData(working.windows) <- cbind(rowData(working.windows), results$table)
rowData(test) = cbind(rowData(test), results$table)

test@rowRanges
working.windows = test

# merge nearby windows
# up to "tol" distance apart: 500 bp in this case
# max merged window width: 5000 bp
merged.peaks <- mergeWindows(rowRanges(working.windows), tol=500L, max.width=5000L)
summary(width(merged.peaks$region))
# should merge some peaks; change as desired

# use most significant window as statistical representation for p-value and FDR for merged windows
tab.best <- getBestTest(merged.peaks$id, results$table)
head(tab.best)
# combine merged peaks window range with statistics
final.merged.peaks <- merged.peaks$region
final.merged.peaks@elementMetadata <- cbind(final.merged.peaks@elementMetadata, tab.best[,-1])
final.merged.peaks <- final.merged.peaks[order(final.merged.peaks@elementMetadata$FDR), ] # sort by FDR
final.merged.peaks # all windows

# filter by FDR threshold
#FDR.thresh <- 0.05 # set as desired
#final.merged.peaks.sig <- final.merged.peaks[final.merged.peaks@elementMetadata$FDR < FDR.thresh, ]
#final.merged.peaks.sig # significant differentially-accessible windows



#########1#########2#########3#########4#########5#########6#########7#########8

library(DESeq2)

test = subset(counts, select = !(colnames(counts) %in% c("lib001", "lib002","lib007","lib005", "lib006", "lib003", "lib004","lib013","lib018", "lib023", "lib014")))
counts = test


test = rlog(assays(counts)$counts)
rld = test

rld = vst(assays(counts)$counts)
mat = t(rld)
pca = prcomp(mat)

summary(pca)

pca_plot = as.data.frame(pca$x) %>%
  rownames_to_column(var = "library_id") %>%
  left_join(libraries_full, by = "library_id") %>%
  ggplot(., aes(x = PC1, y = PC2, color = cohort_id)) +
  geom_point(size = 4)
pca_plot



#lowdate = as.character(data.frame(library_id = colnames(y)) %>% left_join(libraries_full, by = "library_id") %>% pull(flow_date))

#########1#########2#########3#########4#########5#########6#########7#########8

#!/usr/bin/env Rscript
#########1#########2#########3#########4#########5#########6#########7#########8
###
###   Do differential expression of ATAC-seq peaks through edgeR   ###
###

args = commandArgs(trailingOnly = TRUE)
= args[1]

library(csaw)
library(DESeq2)
library(edgeR)
library(tidyverse)

# Load counts as DGE list
counts = readRDS(input)

counts = readRDS("/mnt/ris/jschwarz/cardiac-radiobiology/atac/norm/macs2_all_union_open_narrow_tmm_rse.rds")
load("/mnt/ris/jschwarz/cardiac-radiobiology/data_model/data_model.RData")

# setup design matrix
# see edgeR manual for more information
y <- asDGEList(counts)
colnames(y$counts) =
rownames(y$samples) = names(counts$bam.files)
groups =
  data.frame(library_id = names(counts$bam.files)) %>%
  left_join(libraries_full, by = "library_id") %>%
  droplevels() %>%
  pull(cohort_id)
groups = fct_relevel(groups, "sham", "ir48h")
y$samples$group = groups
colors = as.character(factor(y$samples$group, levels = c("sham", "ir48h", "ir6w"), labels = c("darkgreen", "red", "blue")))

plotMDS(y, col = colors, gene.selection = "common", top = 1000000)

test

test = subset(counts, select = !(colnames(counts) %in% c("lib001", "lib002","lib007","lib005", "lib006", "lib003", "lib004")))
test = subset(counts, select = !(colnames(counts) %in% c("lib001", "lib002","lib007","lib005", "lib006", "lib003", "lib004","lib013","lib018")))
counts = test

# setup design matrix
# see edgeR manual for more information
y <- asDGEList(counts)
colnames(y$counts) = rownames(y$samples) = names(counts$bam.files)
groups =
  data.frame(library_id = names(counts$bam.files)) %>%
  left_join(libraries_full, by = "library_id") %>%
  droplevels() %>%
  pull(cohort_id)
groups = fct_relevel(groups, "sham", "ir48h")
groups
y$samples$group = groups
colors = as.character(factor(y$samples$group, levels = c("sham", "ir48h", "ir6w"), labels = c("darkgreen", "red", "blue")))

pdf("/tmp/pca.pdf")
plotMDS(y, col = colors, gene.selection = "common", top = 80)
dev.off()

plotMDS(y, col = colors, top = 100)

design <- model.matrix(~group, data=y$samples)
colnames(design) = levels(groups)


# stabilize dispersion estimates with empirical bayes
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design, robust=TRUE)

summary(fit$df.prior)

fit <- glmQLFit(y, design)

class(design)
# testing for differentially-accessible windows
results <- glmQLFTest(fit, contrast=makeContrasts(sham-ir6w, levels=design))
# head(results$table)

topTags(results)

# combine GRanges rowdata with DA statistics
rowData(counts) <- cbind(rowData(counts), results$table)

res = as.data.frame(topTags(results, n = Inf))

ggplot(res, aes(x = logFC)) + geom_density()
summary(as.data.frame(topTags(results, n = Inf))$FDR)

test = as_tibble(as.data.frame(topTags(results, n = Inf)))

max(test$FDR)

summary(results$table$PValue)

head(results$table$PValue)

fit = glmFit(y, design, contrast = makeContrasts(ir48h-sham, levels = design))

fit
lrt = glmLRT(fit, contrast = makeContrasts(ir48h-sham, levels = design))
test=as.data.frame(topTags(lrt, n = 10000))
class(test)
summary(test$FDR)
lrt
head(lrt$table)
et = exactTest(y)
topTags(et)

# merge nearby windows
# up to "tol" distance apart: 500 bp in this case
# max merged window width: 5000 bp
merged.peaks <- mergeWindows(rowRanges(counts), tol=500L, max.width=5000L)
# summary(width(merged.peaks$region))
# should merge some peaks; change as desired

# use most significant window as statistical representation for p-value and FDR for merged windows
tab.best <- getBestTest(merged.peaks$id, results$table)
head(tab.best)
min(tab.best$PValue)
min(tab.best$FDR)

# combine merged peaks window range with statistics
final.merged.peaks <- merged.peaks$region
final.merged.peaks@elementMetadata <- cbind(final.merged.peaks@elementMetadata, tab.best[,-1])
final.merged.peaks <- final.merged.peaks[order(final.merged.peaks@elementMetadata$FDR), ] # sort by FDR
final.merged.peaks # all windows

# filter by FDR threshold
FDR.thresh <- 0.05 # set as desired
final.merged.peaks.sig <- final.merged.peaks[final.merged.peaks@elementMetadata$FDR < FDR.thresh, ]
final.merged.peaks.sig # significant differentially-accessible windows




colnames(design) = levels(counts$samples$group)

test = rlog(assays(counts)$counts)
rld = test

class(rld)
mat = t(rld)
pca = prcomp(mat)

summary(pca)

head(counts$counts)
rownames(counts$counts)

class(working.windows)

working.windows

# stabilize dispersion estimates with empirical bayes
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design, robust=TRUE)

# testing for differentially-accessible windows
#results <- glmQLFTest(fit, contrast=makeContrasts(treat-control, levels=design))

results <- glmQLFTest(fit, contrast=makeContrasts(ir48h-sham, levels=design))
# head(results$table)

test = results$table
min(test$PValue)

class(working.windows)

test = working.windows[,8:15]


# combine GRanges rowdata with DA statistics
#rowData(working.windows) <- cbind(rowData(working.windows), results$table)
rowData(test) = cbind(rowData(test), results$table)

test@rowRanges
working.windows = test

# merge nearby windows
# up to "tol" distance apart: 500 bp in this case
# max merged window width: 5000 bp
merged.peaks <- mergeWindows(rowRanges(working.windows), tol=500L, max.width=5000L)
summary(width(merged.peaks$region))
# should merge some peaks; change as desired

# use most significant window as statistical representation for p-value and FDR for merged windows
tab.best <- getBestTest(merged.peaks$id, results$table)
head(tab.best)
# combine merged peaks window range with statistics
final.merged.peaks <- merged.peaks$region
final.merged.peaks@elementMetadata <- cbind(final.merged.peaks@elementMetadata, tab.best[,-1])
final.merged.peaks <- final.merged.peaks[order(final.merged.peaks@elementMetadata$FDR), ] # sort by FDR
final.merged.peaks # all windows

# filter by FDR threshold
#FDR.thresh <- 0.05 # set as desired
#final.merged.peaks.sig <- final.merged.peaks[final.merged.peaks@elementMetadata$FDR < FDR.thresh, ]
#final.merged.peaks.sig # significant differentially-accessible windows



#########1#########2#########3#########4#########5#########6#########7#########8

library(DESeq2)

test = subset(counts, select = !(colnames(counts) %in% c("lib001", "lib002","lib007","lib005", "lib006", "lib003", "lib004","lib013","lib018", "lib023", "lib014")))
counts = test


test = rlog(assays(counts)$counts)
rld = test

rld = vst(assays(counts)$counts)
mat = t(rld)
pca = prcomp(mat)

summary(pca)

pca_plot = as.data.frame(pca$x) %>%
  rownames_to_column(var = "library_id") %>%
  left_join(libraries_full, by = "library_id") %>%
  ggplot(., aes(x = PC1, y = PC2, color = cohort_id)) +
  geom_point(size = 4)
pca_plot



#lowdate = as.character(data.frame(library_id = colnames(y)) %>% left_join(libraries_full, by = "library_id") %>% pull(flow_date))

#########1#########2#########3#########4#########5#########6#########7#########8

#!/usr/bin/env Rscript
#########1#########2#########3#########4#########5#########6#########7#########8
###
###   Do differential expression of ATAC-seq peaks through edgeR   ###
###

args = commandArgs(trailingOnly = TRUE)
= args[1]

library(csaw)
library(DESeq2)
library(edgeR)
library(tidyverse)

# Load counts as DGE list
counts = readRDS(input)

counts = readRDS("/mnt/ris/jschwarz/cardiac-radiobiology/atac/norm/macs2_all_union_open_narrow_tmm_rse.rds")
load("/mnt/ris/jschwarz/cardiac-radiobiology/data_model/data_model.RData")

# setup design matrix
# see edgeR manual for more information
y <- asDGEList(counts)
colnames(y$counts) =
rownames(y$samples) = names(counts$bam.files)
groups =
  data.frame(library_id = names(counts$bam.files)) %>%
  left_join(libraries_full, by = "library_id") %>%
  droplevels() %>%
  pull(cohort_id)
groups = fct_relevel(groups, "sham", "ir48h")
y$samples$group = groups
colors = as.character(factor(y$samples$group, levels = c("sham", "ir48h", "ir6w"), labels = c("darkgreen", "red", "blue")))

plotMDS(y, col = colors, gene.selection = "common", top = 1000000)

test

test = subset(counts, select = !(colnames(counts) %in% c("lib001", "lib002","lib007","lib005", "lib006", "lib003", "lib004")))
test = subset(counts, select = !(colnames(counts) %in% c("lib001", "lib002","lib007","lib005", "lib006", "lib003", "lib004","lib013","lib018")))
counts = test

# setup design matrix
# see edgeR manual for more information
y <- asDGEList(counts)
colnames(y$counts) = rownames(y$samples) = names(counts$bam.files)
groups =
  data.frame(library_id = names(counts$bam.files)) %>%
  left_join(libraries_full, by = "library_id") %>%
  droplevels() %>%
  pull(cohort_id)
groups = fct_relevel(groups, "sham", "ir48h")
groups
y$samples$group = groups
colors = as.character(factor(y$samples$group, levels = c("sham", "ir48h", "ir6w"), labels = c("darkgreen", "red", "blue")))

pdf("/tmp/pca.pdf")
plotMDS(y, col = colors, gene.selection = "common", top = 80)
dev.off()

plotMDS(y, col = colors, top = 100)

design <- model.matrix(~group, data=y$samples)
colnames(design) = levels(groups)


# stabilize dispersion estimates with empirical bayes
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design, robust=TRUE)

summary(fit$df.prior)

fit <- glmQLFit(y, design)

class(design)
# testing for differentially-accessible windows
results <- glmQLFTest(fit, contrast=makeContrasts(sham-ir6w, levels=design))
# head(results$table)

topTags(results)

# combine GRanges rowdata with DA statistics
rowData(counts) <- cbind(rowData(counts), results$table)

res = as.data.frame(topTags(results, n = Inf))

ggplot(res, aes(x = logFC)) + geom_density()
summary(as.data.frame(topTags(results, n = Inf))$FDR)

test = as_tibble(as.data.frame(topTags(results, n = Inf)))

max(test$FDR)

summary(results$table$PValue)

head(results$table$PValue)

fit = glmFit(y, design, contrast = makeContrasts(ir48h-sham, levels = design))

fit
lrt = glmLRT(fit, contrast = makeContrasts(ir48h-sham, levels = design))
test=as.data.frame(topTags(lrt, n = 10000))
class(test)
summary(test$FDR)
lrt
head(lrt$table)
et = exactTest(y)
topTags(et)

# merge nearby windows
# up to "tol" distance apart: 500 bp in this case
# max merged window width: 5000 bp
merged.peaks <- mergeWindows(rowRanges(counts), tol=500L, max.width=5000L)
# summary(width(merged.peaks$region))
# should merge some peaks; change as desired

# use most significant window as statistical representation for p-value and FDR for merged windows
tab.best <- getBestTest(merged.peaks$id, results$table)
head(tab.best)
min(tab.best$PValue)
min(tab.best$FDR)

# combine merged peaks window range with statistics
final.merged.peaks <- merged.peaks$region
final.merged.peaks@elementMetadata <- cbind(final.merged.peaks@elementMetadata, tab.best[,-1])
final.merged.peaks <- final.merged.peaks[order(final.merged.peaks@elementMetadata$FDR), ] # sort by FDR
final.merged.peaks # all windows

# filter by FDR threshold
FDR.thresh <- 0.05 # set as desired
final.merged.peaks.sig <- final.merged.peaks[final.merged.peaks@elementMetadata$FDR < FDR.thresh, ]
final.merged.peaks.sig # significant differentially-accessible windows




colnames(design) = levels(counts$samples$group)

test = rlog(assays(counts)$counts)
rld = test

class(rld)
mat = t(rld)
pca = prcomp(mat)

summary(pca)

head(counts$counts)
rownames(counts$counts)

class(working.windows)

working.windows

# stabilize dispersion estimates with empirical bayes
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design, robust=TRUE)

# testing for differentially-accessible windows
#results <- glmQLFTest(fit, contrast=makeContrasts(treat-control, levels=design))

results <- glmQLFTest(fit, contrast=makeContrasts(ir48h-sham, levels=design))
# head(results$table)

test = results$table
min(test$PValue)

class(working.windows)

test = working.windows[,8:15]


# combine GRanges rowdata with DA statistics
#rowData(working.windows) <- cbind(rowData(working.windows), results$table)
rowData(test) = cbind(rowData(test), results$table)

test@rowRanges
working.windows = test

# merge nearby windows
# up to "tol" distance apart: 500 bp in this case
# max merged window width: 5000 bp
merged.peaks <- mergeWindows(rowRanges(working.windows), tol=500L, max.width=5000L)
summary(width(merged.peaks$region))
# should merge some peaks; change as desired

# use most significant window as statistical representation for p-value and FDR for merged windows
tab.best <- getBestTest(merged.peaks$id, results$table)
head(tab.best)
# combine merged peaks window range with statistics
final.merged.peaks <- merged.peaks$region
final.merged.peaks@elementMetadata <- cbind(final.merged.peaks@elementMetadata, tab.best[,-1])
final.merged.peaks <- final.merged.peaks[order(final.merged.peaks@elementMetadata$FDR), ] # sort by FDR
final.merged.peaks # all windows

# filter by FDR threshold
#FDR.thresh <- 0.05 # set as desired
#final.merged.peaks.sig <- final.merged.peaks[final.merged.peaks@elementMetadata$FDR < FDR.thresh, ]
#final.merged.peaks.sig # significant differentially-accessible windows



#########1#########2#########3#########4#########5#########6#########7#########8

library(DESeq2)

test = subset(counts, select = !(colnames(counts) %in% c("lib001", "lib002","lib007","lib005", "lib006", "lib003", "lib004","lib013","lib018", "lib023", "lib014")))
counts = test


test = rlog(assays(counts)$counts)
rld = test

rld = vst(assays(counts)$counts)
mat = t(rld)
pca = prcomp(mat)

summary(pca)

pca_plot = as.data.frame(pca$x) %>%
  rownames_to_column(var = "library_id") %>%
  left_join(libraries_full, by = "library_id") %>%
  ggplot(., aes(x = PC1, y = PC2, color = cohort_id)) +
  geom_point(size = 4)
pca_plot



#lowdate = as.character(data.frame(library_id = colnames(y)) %>% left_join(libraries_full, by = "library_id") %>% pull(flow_date))

#########1#########2#########3#########4#########5#########6#########7#########8

#!/usr/bin/env Rscript
#########1#########2#########3#########4#########5#########6#########7#########8
###
###   Do differential expression of ATAC-seq peaks through edgeR   ###
###

args = commandArgs(trailingOnly = TRUE)
= args[1]

library(csaw)
library(DESeq2)
library(edgeR)
library(tidyverse)

# Load counts as DGE list
counts = readRDS(input)

counts = readRDS("/mnt/ris/jschwarz/cardiac-radiobiology/atac/norm/macs2_all_union_open_narrow_tmm_rse.rds")
load("/mnt/ris/jschwarz/cardiac-radiobiology/data_model/data_model.RData")

# setup design matrix
# see edgeR manual for more information
y <- asDGEList(counts)
colnames(y$counts) =
rownames(y$samples) = names(counts$bam.files)
groups =
  data.frame(library_id = names(counts$bam.files)) %>%
  left_join(libraries_full, by = "library_id") %>%
  droplevels() %>%
  pull(cohort_id)
groups = fct_relevel(groups, "sham", "ir48h")
y$samples$group = groups
colors = as.character(factor(y$samples$group, levels = c("sham", "ir48h", "ir6w"), labels = c("darkgreen", "red", "blue")))

plotMDS(y, col = colors, gene.selection = "common", top = 1000000)

test

test = subset(counts, select = !(colnames(counts) %in% c("lib001", "lib002","lib007","lib005", "lib006", "lib003", "lib004")))
test = subset(counts, select = !(colnames(counts) %in% c("lib001", "lib002","lib007","lib005", "lib006", "lib003", "lib004","lib013","lib018")))
counts = test

# setup design matrix
# see edgeR manual for more information
y <- asDGEList(counts)
colnames(y$counts) = rownames(y$samples) = names(counts$bam.files)
groups =
  data.frame(library_id = names(counts$bam.files)) %>%
  left_join(libraries_full, by = "library_id") %>%
  droplevels() %>%
  pull(cohort_id)
groups = fct_relevel(groups, "sham", "ir48h")
groups
y$samples$group = groups
colors = as.character(factor(y$samples$group, levels = c("sham", "ir48h", "ir6w"), labels = c("darkgreen", "red", "blue")))

pdf("/tmp/pca.pdf")
plotMDS(y, col = colors, gene.selection = "common", top = 80)
dev.off()

plotMDS(y, col = colors, top = 100)

design <- model.matrix(~group, data=y$samples)
colnames(design) = levels(groups)


# stabilize dispersion estimates with empirical bayes
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design, robust=TRUE)

summary(fit$df.prior)

fit <- glmQLFit(y, design)

class(design)
# testing for differentially-accessible windows
results <- glmQLFTest(fit, contrast=makeContrasts(sham-ir6w, levels=design))
# head(results$table)

topTags(results)

# combine GRanges rowdata with DA statistics
rowData(counts) <- cbind(rowData(counts), results$table)

res = as.data.frame(topTags(results, n = Inf))

ggplot(res, aes(x = logFC)) + geom_density()
summary(as.data.frame(topTags(results, n = Inf))$FDR)

test = as_tibble(as.data.frame(topTags(results, n = Inf)))

max(test$FDR)

summary(results$table$PValue)

head(results$table$PValue)

fit = glmFit(y, design, contrast = makeContrasts(ir48h-sham, levels = design))

fit
lrt = glmLRT(fit, contrast = makeContrasts(ir48h-sham, levels = design))
test=as.data.frame(topTags(lrt, n = 10000))
class(test)
summary(test$FDR)
lrt
head(lrt$table)
et = exactTest(y)
topTags(et)

# merge nearby windows
# up to "tol" distance apart: 500 bp in this case
# max merged window width: 5000 bp
merged.peaks <- mergeWindows(rowRanges(counts), tol=500L, max.width=5000L)
# summary(width(merged.peaks$region))
# should merge some peaks; change as desired

# use most significant window as statistical representation for p-value and FDR for merged windows
tab.best <- getBestTest(merged.peaks$id, results$table)
head(tab.best)
min(tab.best$PValue)
min(tab.best$FDR)

# combine merged peaks window range with statistics
final.merged.peaks <- merged.peaks$region
final.merged.peaks@elementMetadata <- cbind(final.merged.peaks@elementMetadata, tab.best[,-1])
final.merged.peaks <- final.merged.peaks[order(final.merged.peaks@elementMetadata$FDR), ] # sort by FDR
final.merged.peaks # all windows

# filter by FDR threshold
FDR.thresh <- 0.05 # set as desired
final.merged.peaks.sig <- final.merged.peaks[final.merged.peaks@elementMetadata$FDR < FDR.thresh, ]
final.merged.peaks.sig # significant differentially-accessible windows




colnames(design) = levels(counts$samples$group)

test = rlog(assays(counts)$counts)
rld = test

class(rld)
mat = t(rld)
pca = prcomp(mat)

summary(pca)

head(counts$counts)
rownames(counts$counts)

class(working.windows)

working.windows

# stabilize dispersion estimates with empirical bayes
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design, robust=TRUE)

# testing for differentially-accessible windows
#results <- glmQLFTest(fit, contrast=makeContrasts(treat-control, levels=design))

results <- glmQLFTest(fit, contrast=makeContrasts(ir48h-sham, levels=design))
# head(results$table)

test = results$table
min(test$PValue)

class(working.windows)

test = working.windows[,8:15]


# combine GRanges rowdata with DA statistics
#rowData(working.windows) <- cbind(rowData(working.windows), results$table)
rowData(test) = cbind(rowData(test), results$table)

test@rowRanges
working.windows = test

# merge nearby windows
# up to "tol" distance apart: 500 bp in this case
# max merged window width: 5000 bp
merged.peaks <- mergeWindows(rowRanges(working.windows), tol=500L, max.width=5000L)
summary(width(merged.peaks$region))
# should merge some peaks; change as desired

# use most significant window as statistical representation for p-value and FDR for merged windows
tab.best <- getBestTest(merged.peaks$id, results$table)
head(tab.best)
# combine merged peaks window range with statistics
final.merged.peaks <- merged.peaks$region
final.merged.peaks@elementMetadata <- cbind(final.merged.peaks@elementMetadata, tab.best[,-1])
final.merged.peaks <- final.merged.peaks[order(final.merged.peaks@elementMetadata$FDR), ] # sort by FDR
final.merged.peaks # all windows

# filter by FDR threshold
#FDR.thresh <- 0.05 # set as desired
#final.merged.peaks.sig <- final.merged.peaks[final.merged.peaks@elementMetadata$FDR < FDR.thresh, ]
#final.merged.peaks.sig # significant differentially-accessible windows



#########1#########2#########3#########4#########5#########6#########7#########8

library(DESeq2)

test = subset(counts, select = !(colnames(counts) %in% c("lib001", "lib002","lib007","lib005", "lib006", "lib003", "lib004","lib013","lib018", "lib023", "lib014")))
counts = test


test = rlog(assays(counts)$counts)
rld = test

rld = vst(assays(counts)$counts)
mat = t(rld)
pca = prcomp(mat)

summary(pca)

pca_plot = as.data.frame(pca$x) %>%
  rownames_to_column(var = "library_id") %>%
  left_join(libraries_full, by = "library_id") %>%
  ggplot(., aes(x = PC1, y = PC2, color = cohort_id)) +
  geom_point(size = 4)
pca_plot



#lowdate = as.character(data.frame(library_id = colnames(y)) %>% left_join(libraries_full, by = "library_id") %>% pull(flow_date))

#########1#########2#########3#########4#########5#########6#########7#########8
