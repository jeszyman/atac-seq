#########1#########2#########3#########4#########5#########6#########7#########8

###############################################################################
###    Script to make a custom DESeq2 object for QC-filtered nuclear bulk   ###
###    RNA-seq data  WITHOUT FILTERING                                      ###
###############################################################################

args = commandArgs(trailingOnly = TRUE)
data_model = args[1]
txi_counts = readRDS(args[2])
nuc_deseq = args[3]

load(data_model)
library(tidyverse)
library(DESeq2)
library(EnsDb.Mmusculus.v79)
edb = EnsDb.Mmusculus.v79

deseq_column_data =
  data.frame(library_id = colnames(txi_counts$counts)) %>%
  left_join(libraries_full, by = "library_id") %>%
  dplyr::select(run_id,flow_date,cohort_id, lib_typ) %>%
  mutate(flow_date = replace_na(as.character(flow_date), "noflow")) %>%
  mutate(sampletype = paste(run_id,flow_date,cohort_id,lib_typ, sep = "_"))

dds = DESeqDataSetFromTximport(txi = txi_counts,
                               colData = deseq_column_data,
                               design = ~ cohort_id + run_id)
dds$cohort_id = factor(dds$cohort_id, levels = c("sham", "ir24h", "ir2w", "ir6w"))

nuc = dds[ ,dds$lib_typ == "nuc_rna"]
nuc$lib_typ = droplevels(nuc$lib_typ)
nuc$cohort_id = droplevels(nuc$cohort_id)
design(nuc) = ~ cohort_id

full_nuc_vsd <- vst(nuc)
full_nuc_rld = rlog(nuc)

save(full_nuc_vsd,
     full_nuc_rld,
     file = nuc_deseq)
