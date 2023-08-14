dca_str="~/cards/analysis/atac/mouse/dca/ms_atac_k0_ir2d_sham.tsv ~/cards/analysis/atac/mouse/dca/ms_atac_k1_ir2d_sham.tsv ~/cards/analysis/atac/mouse/dca/ms_atac_k2_ir2d_sham.tsv ~/cards/analysis/atac/mouse/dca/ms_atac_k3_ir2d_sham.tsv ~/cards/analysis/atac/mouse/dca/ms_atac_k4_ir2d_sham.tsv"
libraries_full_rds="~/cards/data-model/lists/libraries_full.rds"

#!/usr/bin/env Rscript

#############################
###   Atac Rvu Spearman   ###
#############################

# Command line arguements
args = commandArgs(trailingOnly = TRUE)
dca_str = args[1]
libraries_full_rds = args[2]
out_rda = args[3]
out_pdf = args[4]

# Load required packages, data, and functions
library(tidyverse)

# Setup data
libraries_full = readRDS(libraries_full_rds)
(dca_char = strsplit(dca_str, " ")[[1]])
tibble_list = lapply(dca_char, read_tsv)
(names(tibble_list) = substr(gsub("^.*atac_", "", dca_char), 1, 2))


# Rename columns and merge tibbles
merged <- reduce(names(tibble_list), function(result, suffix) {
  tibble_i <- tibble_list[[suffix]]
  colnames(tibble_i) <- ifelse(colnames(tibble_i) == "coordinate", "coordinate", paste0(suffix, "_", colnames(tibble_i)))

  if (is.null(result)) {
    return(tibble_i)
  } else {
    return(left_join(result, tibble_i, by = "coordinate"))
  }
}, .init = NULL)

plot_merged <- function(merged) {
  # Create a data frame to store the R² values and p-values
  stats_values <- merged %>%
    select(coordinate, ends_with("logfc")) %>%
    pivot_longer(cols = !c("coordinate", "k0_logfc"), names_to = "rvu_k", values_to = "logfc") %>%
    mutate(dif = abs(k0_logfc - logfc)) %>% group_by(rvu_k) %>%
    summarize(r2 = cor(k0_logfc, logfc)^2, p_value_diff = wilcox.test(k0_logfc, logfc)$p.value)
  stats_values

  exceeds <- merged %>% select(coordinate, ends_with("qval")) %>%
    pivot_longer(cols = !coordinate, names_to = "rvu_k", values_to = "qval")  %>%
    mutate(qsig = ifelse(qval < 0.05, "yes", "no")) %>%
    group_by(rvu_k) %>%
    summarize(qsum = sum(qsig == "yes")) %>%
    mutate(rvu_k = paste0(substr(rvu_k, 1, 2), "_logfc")) %>%
    full_join(stats_values) %>%
    filter(rvu_k != "k0_logfc")

  plot <- merged %>%
    select(coordinate, ends_with("logfc")) %>%
    pivot_longer(cols = !c("coordinate", "k0_logfc"), names_to = "rvu_k", values_to = "logfc") %>%
    mutate(dif = abs(k0_logfc - logfc)) %>% group_by(rvu_k) %>%
    mutate(sd2 = abs(mean(dif, na.rm = T)) + 3 * abs(sd(dif, na.rm = T))) %>%
    filter(dif > sd2) %>%
    ggplot(., aes(x = k0_logfc, y = logfc)) + geom_point() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "blue") +
    geom_text(data = exceeds, aes(x = -Inf, y = Inf, label = paste0("R² = ", round(r2, 3))),
              hjust = 0, nudge_x = -1, vjust = 2, size = 4) +
    geom_text(data = exceeds, aes(x = -Inf, y = Inf, label = paste0("p(diff) = ", format.pval(p_value_diff, digits = 3))),
              hjust = 0, nudge_x = -1, vjust = 4, size = 4) +
    geom_text(data = exceeds, aes(x = -Inf, y = Inf, label = paste0(qsum, " regions with FDR < 0.05")),
              hjust = 0, nudge_x = -1, vjust = 6, size = 4) +
    facet_wrap(~rvu_k) +     theme(text = element_text(size = 14), plot.margin = margin(1, 1, 1, 1, "cm"))


  return(plot)
}

plot = plot_merged(merged)

save(merged, plot, file = out_rda)

ggsave(plot, file = out_pdf)
