#!/usr/bin/env Rscript

# ---   Setup   --- #
# ----------------- #

## ---   Load Packages   --- ##
## ------------------------- ##

library(optparse)
# Load last
library(tidyverse)

## ---   Load Inputs   --- ##
## ----------------------- ##

option_list <- list(
  make_option(c("--agg_samstats_tsv"),
              type = "character",
              default = "~/cards/analysis/atac/qc/human_filt/agg_samstats.tsv"),
  make_option(c("--libraries_full_rds"),
              type = "character",
              default = "~/cards/data-model/lists/libraries_full.rds"),
  make_option(c("--mil_reads"),
              type = "numeric",
              default = "9"),
  make_option(c("--out_tsv"),
              type = "character",
              default = "~/cards/analysis/atac/qc/human_filt/keep.tsv")

)

opts <- parse_args(OptionParser(option_list = option_list))

list_of_options <- names(opts)
for (opt_name in list_of_options) {
  assign(opt_name, opts[[opt_name]], envir = .GlobalEnv)
}

libraries_full <- readRDS(libraries_full_rds)
samstats = read_tsv(agg_samstats_tsv)

keep = samstats %>% left_join(libraries_full) %>%
  mutate(adequate_reads = total_reads > (mil_reads * 1000000)) %>%
  mutate(atac_qc = ifelse(adequate_reads == TRUE,
                          "PASS", "FAIL")) %>%
  select(library, atac_qc, adequate_reads)

keep

write_tsv(keep, file = out_tsv)
