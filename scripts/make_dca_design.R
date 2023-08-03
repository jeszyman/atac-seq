#!/usr/bin/env Rscript

###############################
###   Make Rna-Seq Design   ###
###############################

# ---   Command Line Arguements   --- #
# ----------------------------------- #

args = commandArgs(trailingOnly = TRUE)
libraries_full_rds = args[1]
formula = args[2]
libs_str = args[3]
design_rds = args[4]

# ---   Load   --- #
# ---------------- #

library(tidyverse)
libraries_full = readRDS(libraries_full_rds)
libs_vect = strsplit(libs_str, " ")[[1]]

# ---   Run   --- #
# --------------- #

libs =
  data.frame(library = libs_vect) %>%
  left_join(libraries_full) %>%
  mutate(across(where(is.factor), droplevels))

design = model.matrix(as.formula(formula), data = libs)


# Splitting the formula to isolate the relevant part ("cohort")
split_formula <- strsplit(formula, " ")[[1]]
for (part in split_formula) {
  # Skipping over "+" as it is an operator
  if (part != "+") {
    colnames(design) <- gsub(part, "", colnames(design))
  }
}

rownames(design) = libs$library


saveRDS(object = design,
        file = design_rds)
