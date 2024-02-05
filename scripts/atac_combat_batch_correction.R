#!/usr/bin/env Rscript

# ---   Load Packages   --- #
# ------------------------- #

# List of packages to load
packages <- c(
  "cowplot",
  "optparse",
  "edgeR",
  "ggrepel",
  "gridExtra",
  "patchwork",
  "sva",
  "tidyverse")

sapply(packages, require, character.only = TRUE, quietly = TRUE)

# ---   Load Inputs   --- #
# ----------------------- #

option_list <- list(
  make_option(c("--batch_var"), type = "character", default = "run"),
  make_option(c("--counts_tsv"), type = "character", default = "~/cards/analysis/atac/models/unadjusted/filt_mouse/bamscale/raw_coverages.tsv"),
  make_option(c("--covars"), type = "character", default = "cohort"),
  make_option(c("--design_rds"), type = "character", default = "~/cards/analysis/atac/models/unadjusted/filt_mouse/design.rds"),
  make_option(c("--libraries_full_rds"), type = "character", default = "~/cards/data-model/lists/libraries_full.rds"),
  make_option(c("--out_dir"), type = "character", default = "/tmp/")
)

opts <- parse_args(OptionParser(option_list = option_list))

batch_var = opts$batch_var
design <- readRDS(opts$design_rds)
libraries_full <- readRDS(opts$libraries_full_rds)
counts <- read_tsv(opts$counts_tsv)
covars_int <- unlist(strsplit(opts$covars, " "))

# ---   Main   --- #
# ---------------- #

main = function(counts, batch_var, libraries_full, covars_int, design){
  counts_mat = format_counts_matrix(counts)
  batch_int = make_batch_int(counts_mat, libraries_full, batch_var)
  adjusted = do_combat(counts_mat, batch_int, covars_int, libraries_full)

  edger_unadjusted = run_edger(counts_mat, design)
  edger_adjusted = run_edger(adjusted, design)

  pca_unadjusted = make_pca(edger_unadjusted$logcpm, libraries_full, covars_int)
  pca_adjusted = make_pca(edger_adjusted$logcpm, libraries_full, covars_int)

  return(list(dge_adjusted = edger_adjusted$dge,
              fit_adjusted = edger_adjusted$fit,
              logcpm_adjusted = edger_adjusted$logcpm,
              pca_adjusted = pca_adjusted,
              dge_unadjusted = edger_unadjusted$dge,
              fit_unadjusted = edger_unadjusted$fit,
              logcpm_unadjusted = edger_unadjusted$logcpm,
              pca_unadjusted = pca_unadjusted))
}

# ---   Functions   --- #
# --------------------- #

format_counts_matrix = function(counts){
  counts2 <- counts %>%
    rename_with(~ifelse(is.na(str_extract(.x, "lib\\d+")), .x, str_extract(.x, "lib\\d+")))

  counts2 = counts2 %>%
    column_to_rownames(var = names(counts2)[1])

  counts2 = counts2[-1,]

  counts_mat = as.matrix(counts2)

  return(counts_mat)
}

make_batch_int = function(counts_mat, libraries_full, batch_var){
  batch <- data.frame(library = colnames(counts_mat)) %>%
    left_join(libraries_full, by = "library") %>%
    pull(!!sym(batch_var)) %>%
    as.factor() %>%
    as.integer()

  return(batch)
}

do_combat = function(counts_mat, batch_int, covars_int, libraries_full){

  # Check if there is only one covariate
  if (length(covars_int) == 1) {
    covar <- data.frame(library = colnames(counts_mat)) %>%
      left_join(libraries_full, by = "library") %>%
      pull(!!sym(covars_int[1])) %>%
      as.factor() %>%
      as.integer()
    adjusted <- ComBat_seq(counts_mat,
                           batch = batch_int,
                           group = covar)
  } else { # If there are multiple covariates
    covar_list <- list() # Initialize an empty list to store covariate vectors
    for (covar_name in covars_int) {
      covar <- data.frame(library = colnames(counts_mat)) %>%
        left_join(libraries_full, by = "library") %>%
        pull(!!sym(covar_name)) %>%
        as.factor() %>%
        as.integer()
      covar_list[[covar_name]] <- covar
    }
    # Combine all covariate vectors into a matrix
    covar_mat <- do.call(cbind, covar_list)
    adjusted <- ComBat_seq(counts_mat,
                           batch = batch_int,
                           group = NULL,
                           covar_mod = covar_mat)
  }

  return(adjusted)
}

run_edger = function(counts_mat, design){

  dge = DGEList(counts = counts_mat)
  keep = filterByExpr(dge, design)
  dge = dge[keep,]
  dge = calcNormFactors(dge)
  dge = estimateGLMCommonDisp(dge, design)
  dge = estimateGLMTrendedDisp(dge, design)
  dge = estimateGLMTagwiseDisp(dge, design)
  fit = glmFit(dge, design)
  logcpm <- edgeR::cpm(dge, normalized.lib.sizes = TRUE, log = TRUE, prior.count = 2)

   return(list(dge = dge, fit = fit, logcpm = logcpm))
}

make_pca <- function(logcpm, libraries_full, covars_int) {
  pca <- prcomp(t(as.matrix(logcpm[,-1])))
  pve_pc1 <- round(100 * summary(pca)$importance[2,1])
  pve_pc2 <- round(100 * summary(pca)$importance[2,2])
  pca_plot <- as.data.frame(pca$x) %>%
    rownames_to_column(var = "library") %>%
    left_join(libraries_full, by = "library") %>%
    ggplot(aes(x = PC1, y = PC2, label = library, color = !!sym(covars_int[[1]]))) +
    geom_point(size = 4) +
    geom_text_repel() +
    xlab(paste("PC1, ", pve_pc1, "% variance explained", sep ="")) +
    ylab(paste("PC2, ", pve_pc2, "% variance explained", sep ="")) +
    coord_fixed(ratio = 1)

  return(pca_plot)
}

# ---   Run   --- #
# --------------- #

out = main(counts, batch_var, libraries_full, covars_int, design)

# Correct iteration over the names of the list elements
for (name in names(out)) {
  item <- out[[name]]  # Access the list item by name
  file_path <- file.path(opts$out_dir, paste0(name, ".rds"))  # Construct the file path
  saveRDS(item, file_path)  # Save the item as RDS
}
