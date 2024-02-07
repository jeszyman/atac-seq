#!/usr/bin/env Rscript

# ---   Load Packages   --- #
# ------------------------- #

packages <- c("ggrepel", "optparse", "RUVSeq", "tidyverse")
sapply(packages, require, character.only = TRUE, quietly = TRUE)

# ---   Load Inputs   --- #
# ----------------------- #

option_list <- list(
  make_option(c("--counts_tsv"), type = "character"),
  make_option(c("--design_rds"), type = "character"),
  make_option(c("--libraries_full_rds"), type = "character"),
  make_option(c("--out_dir"), type = "character"),
  make_option(c("--formula"), type = "character")
  )

opts <- parse_args(OptionParser(option_list = option_list))

## opts$counts_tsv = "/home/jeszyman/cards/analysis/atac/models/unadjusted/human_icell/bamscale/raw_coverages.tsv"
## opts$design_rds = "/home/jeszyman/cards/analysis/atac/models/unadjusted/human_icell/design.rds"
## opts$libraries_full_rds = "~/cards/data-model/lists/libraries_full.rds"
## opts$ruv_k = 4
## opts$out_dir = "/tmp/test"
## opts$formula = "~ 0 + cohort"

inputs_list = function(opts){
  counts_raw = read_tsv(opts$counts_tsv)
  design = readRDS(opts$design_rds)
  formula = opts$formula
  libraries_full = readRDS(opts$libraries_full_rds)
  out_dir = opts$out_dir

  split_formula <- strsplit(formula, " ")[[1]]
  formula_vect <- split_formula[!split_formula %in% c("~", "0", "+")]

  return(list(counts_raw = counts_raw,
              design = design,
              libraries_full = libraries_full,
              formula_vect = formula_vect,
              out_dir = out_dir))

}

# ---   Main   --- #
# ---------------- #

main = function(opts){
  # Process inputs
  inputs = inputs_list(opts)

  # Adjust raw BAMscale counts using the provided ruv
  adjusted_ruv2 = make_ruv_adjustment(inputs$counts_raw, inputs$design, 2)
  adjusted_ruv3 = make_ruv_adjustment(inputs$counts_raw, inputs$design, 3)
  adjusted_ruv4 = make_ruv_adjustment(inputs$counts_raw, inputs$design, 4)

  # Perform edgeR normalization and DCA on the adjusted counts
  edger_ruv2 = run_edger(adjusted_ruv2, inputs$design)
  edger_ruv3 = run_edger(adjusted_ruv3, inputs$design)
  edger_ruv4 = run_edger(adjusted_ruv4, inputs$design)

  # Get PCA for comparison
  pca_ruv2 = make_pca(edger_ruv2$logcpm, inputs$libraries_full, inputs$formula_vect)
  pca_ruv3 = make_pca(edger_ruv3$logcpm, inputs$libraries_full, inputs$formula_vect)
  pca_ruv4 = make_pca(edger_ruv4$logcpm, inputs$libraries_full, inputs$formula_vect)

  return(list(dge_ruv2 = edger_ruv2$dge,
              fit_ruv2 = edger_ruv2$fit,
              logcpm_ruv2 = edger_ruv2$logcpm,
              pca_ruv2 = pca_ruv2,
              dge_ruv3 = edger_ruv3$dge,
              fit_ruv3 = edger_ruv3$fit,
              logcpm_ruv3 = edger_ruv3$logcpm,
              pca_ruv3 = pca_ruv3,
              dge_ruv4 = edger_ruv4$dge,
              fit_ruv4 = edger_ruv4$fit,
              logcpm_ruv4 = edger_ruv4$logcpm,
              pca_ruv4 = pca_ruv4))
}

# ---   Functions   --- #
# --------------------- #

make_ruv_adjustment = function(counts, design, ruv_k){
  mat = as.matrix(counts[-1,-1])
  row.names(mat) = counts$coordinate[-1]
  colnames(mat) = substr(colnames(mat),1,6)

  model_df = as.data.frame(design)
  model_df$'(Intercept)' <- NULL

  mat <- mat[, rownames(model_df)]

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

  colnames(adjust_counts) = row.names(model_df)

  return(adjust_counts)
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

out = main(opts)

for (name in names(out)) {
  item <- out[[name]]
  file_path <- file.path(opts$out_dir, paste0(name, ".rds"))
  saveRDS(item, file_path)
}

pdf(paste0(opts$out_dir,"/pca_ruv2.pdf"), title = "RUV2 Adjusted PCA")
print(out$pca_ruv2)
dev.off()

pdf(paste0(opts$out_dir,"/pca_ruv3.pdf"), title = "RUV3 Adjusted PCA")
print(out$pca_ruv3)
dev.off()

pdf(paste0(opts$out_dir,"/pca_ruv4.pdf"), title = "RUV4 Adjusted PCA")
print(out$pca_ruv4)
dev.off()
