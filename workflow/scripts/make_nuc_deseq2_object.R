#########1#########2#########3#########4#########5#########6#########7#########8

###############################################################################
###    Script to make a custom DESeq2 object for QC-filtered nuclear bulk   ###
###    RNA-seq data                                                         ###
###############################################################################

args = commandArgs(trailingOnly = TRUE)
data_model = args[1]
txi_counts = readRDS(args[2])
nuc_deseq = args[3]
nuc_res_csv = args[4]

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
nuc = nuc[ ,!(colnames(nuc) %in% c("lib033",
                                   "lib031",
                                   "lib030",
                                   "lib027",
                                   "lib045",
                                   "lib042"))]
nuc$lib_typ = droplevels(nuc$lib_typ)
nuc$cohort_id = droplevels(nuc$cohort_id)
design(nuc) = ~ cohort_id

nuc_vsd <- vst(nuc)
nuc_dds = DESeq(nuc)
nuc_res = results(nuc_dds, name = "cohort_id_ir6w_vs_sham")

res_genes = mapIds(edb, keys=rownames(nuc_res), column="GENENAME", keytype="GENEID", multiVals = 'first')

name_index = data.frame(name = res_genes,
                        ensembl_id = names(res_genes))

nuc_res_anno = as_tibble(as.data.frame(nuc_res)) %>%
  mutate(ensembl_id = row.names(nuc_res)) %>%
  left_join(name_index, by = "ensembl_id") %>%
  arrange(padj)

save(nuc_vsd,
     nuc_dds,
     nuc_res,
     nuc_res_anno,
     file = nuc_deseq)

write.csv(nuc_res_anno, file = nuc_res_csv, row.names = F)
