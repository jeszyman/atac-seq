
## :LOGBOOK:
## - State "WAITING"    from "TODO"       [2022-03-31 Thu 14:24]
## :END:

## - [X] Functions and test

args = commandArgs(trailingOnly = TRUE)
counts_rds = args[1]
background_rds = args[2]
logcpm_file = args[3]

background = readRDS(background_rds)
counts = readRDS(counts_rds)

library(csaw)
library(edgeR)
library(tidyverse)

counts = normFactors(background, se.out = counts)

make_logcpm = function(in_norm){
  dge = asDGEList(in_norm)
  colnames(dge) = colnames(in_norm)
  log_cpm = cpm(dge, normalized.lib.sizes = TRUE, log = TRUE, prior.count = 2)
  return(log_cpm)
}

logcpm = make_logcpm(counts)

saveRDS(object = logcpm,
        file = logcpm_file)
