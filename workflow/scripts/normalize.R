#!/usr/bin/env Rscript
#########1#########2#########3#########4#########5#########6#########7#########8

######################################
###   Normalize csaw peak counts   ###
######################################

args = commandArgs(trailingOnly = TRUE)
rse_file = args[1]
bk_filt = args[2]
tmm_file = args[3]
loess_file = args[4]

rse = readRDS(rse_file)
bk = readRDS(bk_filt)

library(csaw)
library(edgeR)

tmm = normFactors(bk, se.out = rse)
loess = normOffsets(rse, se.out = TRUE)

saveRDS(object = tmm,
        file = tmm_file)
saveRDS(object = loess,
        file = loess_file)

#!/usr/bin/env Rscript
#########1#########2#########3#########4#########5#########6#########7#########8

######################################
###   Normalize csaw peak counts   ###
######################################

args = commandArgs(trailingOnly = TRUE)
rse_file = args[1]
bk_filt = args[2]
tmm_file = args[3]
loess_file = args[4]

rse = readRDS(rse_file)
bk = readRDS(bk_filt)

library(csaw)
library(edgeR)

tmm = normFactors(bk, se.out = rse)
loess = normOffsets(rse, se.out = TRUE)

saveRDS(object = tmm,
        file = tmm_file)
saveRDS(object = loess,
        file = loess_file)

#!/usr/bin/env Rscript
#########1#########2#########3#########4#########5#########6#########7#########8

######################################
###   Normalize csaw peak counts   ###
######################################

args = commandArgs(trailingOnly = TRUE)
rse_file = args[1]
bk_filt = args[2]
tmm_file = args[3]
loess_file = args[4]

rse = readRDS(rse_file)
bk = readRDS(bk_filt)

library(csaw)
library(edgeR)

tmm = normFactors(bk, se.out = rse)
loess = normOffsets(rse, se.out = TRUE)

saveRDS(object = tmm,
        file = tmm_file)
saveRDS(object = loess,
        file = loess_file)

#!/usr/bin/env Rscript
#########1#########2#########3#########4#########5#########6#########7#########8
###
###   Normalize peak counts   ###
###

args = commandArgs(trailingOnly = TRUE)
rse_file = args[1]
bk_filt = args[2]
tmm_file = args[3]
loess_file = args[4]

rse = readRDS(rse_file)
bk = readRDS(bk_filt)

library(csaw)
library(edgeR)

tmm = normFactors(bk, se.out = rse)
loess = normOffsets(rse, se.out = TRUE)

# Make logCPM counts of normalized data
make_logcpm = function(in_norm){
  dge = asDGEList(in_norm)
  colnames(dge) = colnames(in_norm)
  log_cpm = cpm(dge, normalized.lib.sizes = TRUE, log = TRUE, prior.count = 2)
  return(log_cpm)
}

in_counts = readRDS("/mnt/ris/jschwarz/cardiac-radiobiology/atac/counts/macs2_all_intersect_open_broad_peaks_rse.rds")

test =make_logcpm(in_counts)

in_norm = readRDS("/mnt/ris/jschwarz/cardiac-radiobiology/atac/norm/macs2_all_intersect_open_broad_tmm_rse.rds")

dge = asDGEList(in_norm)

in_norm

head(assays(in_norm)$counts)

colnames(dge) = colnames(in_norm)

# offsets for a log-link GLM.
normMat <- sweep(normMat, 2, eff.lib, "*")
normMat <- log(normMat)

# Creating a DGEList object for use in edgeR.
y <- DGEList(cts, group = factor(c(rep("ms_6wk",6), rep("ms_sham", 5))))

y <- scaleOffset(y, normMat)

# filtering
keep <- filterByExpr(y)
## Warning in filterByExpr.DGEList(y): All samples appear to belong to the same
## group.
y <- y[keep, ]
# y is now ready for estimate dispersion functions see edgeR User's Guide

y = calcNormFactors(y)

design <- model.matrix(~0+group, data=y$samples)

y = estimateDisp(y, design, robust = T)



tmm_logcpm = make_logcpm(tmm)



loess_logcpm = make_logcpm(loess)

head(tmm_logcpm)

colnames(tmm)
colnames(test2) = colnames(norm)

pca = prcomp(t(test2))

summary(pca)


saveRDS(object = tmm,
        file = tmm_file)
saveRDS(object = loess,
        file = loess_file)

#!/usr/bin/env Rscript
#########1#########2#########3#########4#########5#########6#########7#########8

######################################
###   Normalize csaw peak counts   ###
######################################

args = commandArgs(trailingOnly = TRUE)
rse_file = args[1]
bk_filt = args[2]
tmm_file = args[3]
loess_file = args[4]

rse = readRDS(rse_file)
bk = readRDS(bk_filt)

library(csaw)
library(edgeR)

tmm = normFactors(bk, se.out = rse)
loess = normOffsets(rse, se.out = TRUE)

saveRDS(object = tmm,
        file = tmm_file)
saveRDS(object = loess,
        file = loess_file)

#!/usr/bin/env Rscript
#########1#########2#########3#########4#########5#########6#########7#########8
###
###   Normalize peak counts   ###
###

args = commandArgs(trailingOnly = TRUE)
rse_file = args[1]
bk_filt = args[2]
tmm_file = args[3]
loess_file = args[4]

rse = readRDS(rse_file)
bk = readRDS(bk_filt)

library(csaw)
library(edgeR)

tmm = normFactors(bk, se.out = rse)
loess = normOffsets(rse, se.out = TRUE)

# Make logCPM counts of normalized data
make_logcpm = function(in_norm){
  dge = asDGEList(in_norm)
  colnames(dge) = colnames(in_norm)
  log_cpm = cpm(dge, normalized.lib.sizes = TRUE, log = TRUE, prior.count = 2)
  return(log_cpm)
}

in_counts = readRDS("/mnt/ris/jschwarz/cardiac-radiobiology/atac/counts/macs2_all_intersect_open_broad_peaks_rse.rds")

test =make_logcpm(in_counts)

in_norm = readRDS("/mnt/ris/jschwarz/cardiac-radiobiology/atac/norm/macs2_all_intersect_open_broad_tmm_rse.rds")

dge = asDGEList(in_norm)

in_norm

head(assays(in_norm)$counts)

colnames(dge) = colnames(in_norm)

# offsets for a log-link GLM.
normMat <- sweep(normMat, 2, eff.lib, "*")
normMat <- log(normMat)

# Creating a DGEList object for use in edgeR.
y <- DGEList(cts, group = factor(c(rep("ms_6wk",6), rep("ms_sham", 5))))

y <- scaleOffset(y, normMat)

# filtering
keep <- filterByExpr(y)
## Warning in filterByExpr.DGEList(y): All samples appear to belong to the same
## group.
y <- y[keep, ]
# y is now ready for estimate dispersion functions see edgeR User's Guide

y = calcNormFactors(y)

design <- model.matrix(~0+group, data=y$samples)

y = estimateDisp(y, design, robust = T)



tmm_logcpm = make_logcpm(tmm)



loess_logcpm = make_logcpm(loess)

head(tmm_logcpm)

colnames(tmm)
colnames(test2) = colnames(norm)

pca = prcomp(t(test2))

summary(pca)


saveRDS(object = tmm,
        file = tmm_file)
saveRDS(object = loess,
        file = loess_file)

#!/usr/bin/env Rscript
#########1#########2#########3#########4#########5#########6#########7#########8
###
###   Normalize peak counts   ###
###

args = commandArgs(trailingOnly = TRUE)
rse_file = args[1]
bk_filt = args[2]
tmm_file = args[3]
loess_file = args[4]

rse = readRDS(rse_file)
bk = readRDS(bk_filt)

library(csaw)
library(edgeR)

tmm = normFactors(bk, se.out = rse)
loess = normOffsets(rse, se.out = TRUE)

# Make logCPM counts of normalized data
make_logcpm = function(in_norm){
  dge = asDGEList(in_norm)
  colnames(dge) = colnames(in_norm)
  log_cpm = cpm(dge, normalized.lib.sizes = TRUE, log = TRUE, prior.count = 2)
  return(log_cpm)
}

in_counts = readRDS("/mnt/ris/jschwarz/cardiac-radiobiology/atac/counts/macs2_all_intersect_open_broad_peaks_rse.rds")

test =make_logcpm(in_counts)

in_norm = readRDS("/mnt/ris/jschwarz/cardiac-radiobiology/atac/norm/macs2_all_intersect_open_broad_tmm_rse.rds")

dge = asDGEList(in_norm)

in_norm

head(assays(in_norm)$counts)

colnames(dge) = colnames(in_norm)

# offsets for a log-link GLM.
normMat <- sweep(normMat, 2, eff.lib, "*")
normMat <- log(normMat)

# Creating a DGEList object for use in edgeR.
y <- DGEList(cts, group = factor(c(rep("ms_6wk",6), rep("ms_sham", 5))))

y <- scaleOffset(y, normMat)

# filtering
keep <- filterByExpr(y)
## Warning in filterByExpr.DGEList(y): All samples appear to belong to the same
## group.
y <- y[keep, ]
# y is now ready for estimate dispersion functions see edgeR User's Guide

y = calcNormFactors(y)

design <- model.matrix(~0+group, data=y$samples)

y = estimateDisp(y, design, robust = T)



tmm_logcpm = make_logcpm(tmm)



loess_logcpm = make_logcpm(loess)

head(tmm_logcpm)

colnames(tmm)
colnames(test2) = colnames(norm)

pca = prcomp(t(test2))

summary(pca)


saveRDS(object = tmm,
        file = tmm_file)
saveRDS(object = loess,
        file = loess_file)

#!/usr/bin/env Rscript
#########1#########2#########3#########4#########5#########6#########7#########8

######################################
###   Normalize csaw peak counts   ###
######################################

args = commandArgs(trailingOnly = TRUE)
rse_file = args[1]
bk_filt = args[2]
tmm_file = args[3]
loess_file = args[4]

rse = readRDS(rse_file)
bk = readRDS(bk_filt)

library(csaw)
library(edgeR)

tmm = normFactors(bk, se.out = rse)
loess = normOffsets(rse, se.out = TRUE)

saveRDS(object = tmm,
        file = tmm_file)
saveRDS(object = loess,
        file = loess_file)

#!/usr/bin/env Rscript
#########1#########2#########3#########4#########5#########6#########7#########8

######################################
###   Normalize csaw peak counts   ###
######################################

args = commandArgs(trailingOnly = TRUE)
rse_file = args[1]
bk_filt = args[2]
tmm_file = args[3]
loess_file = args[4]

rse = readRDS(rse_file)
bk = readRDS(bk_filt)

library(csaw)
library(edgeR)

tmm = normFactors(bk, se.out = rse)
loess = normOffsets(rse, se.out = TRUE)

saveRDS(object = tmm,
        file = tmm_file)
saveRDS(object = loess,
        file = loess_file)

#!/usr/bin/env Rscript
#########1#########2#########3#########4#########5#########6#########7#########8
###
###   Normalize peak counts   ###
###

args = commandArgs(trailingOnly = TRUE)
rse_file = args[1]
bk_filt = args[2]
tmm_file = args[3]
loess_file = args[4]

rse = readRDS(rse_file)
bk = readRDS(bk_filt)

library(csaw)
library(edgeR)

tmm = normFactors(bk, se.out = rse)
loess = normOffsets(rse, se.out = TRUE)

# Make logCPM counts of normalized data
make_logcpm = function(in_norm){
  dge = asDGEList(in_norm)
  colnames(dge) = colnames(in_norm)
  log_cpm = cpm(dge, normalized.lib.sizes = TRUE, log = TRUE, prior.count = 2)
  return(log_cpm)
}

in_counts = readRDS("/mnt/ris/jschwarz/cardiac-radiobiology/atac/counts/macs2_all_intersect_open_broad_peaks_rse.rds")

test =make_logcpm(in_counts)

in_norm = readRDS("/mnt/ris/jschwarz/cardiac-radiobiology/atac/norm/macs2_all_intersect_open_broad_tmm_rse.rds")

dge = asDGEList(in_norm)

in_norm

head(assays(in_norm)$counts)

colnames(dge) = colnames(in_norm)

# offsets for a log-link GLM.
normMat <- sweep(normMat, 2, eff.lib, "*")
normMat <- log(normMat)

# Creating a DGEList object for use in edgeR.
y <- DGEList(cts, group = factor(c(rep("ms_6wk",6), rep("ms_sham", 5))))

y <- scaleOffset(y, normMat)

# filtering
keep <- filterByExpr(y)
## Warning in filterByExpr.DGEList(y): All samples appear to belong to the same
## group.
y <- y[keep, ]
# y is now ready for estimate dispersion functions see edgeR User's Guide

y = calcNormFactors(y)

design <- model.matrix(~0+group, data=y$samples)

y = estimateDisp(y, design, robust = T)



tmm_logcpm = make_logcpm(tmm)



loess_logcpm = make_logcpm(loess)

head(tmm_logcpm)

colnames(tmm)
colnames(test2) = colnames(norm)

pca = prcomp(t(test2))

summary(pca)


saveRDS(object = tmm,
        file = tmm_file)
saveRDS(object = loess,
        file = loess_file)
