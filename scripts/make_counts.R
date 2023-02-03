#########1#########2#########3#########4#########5#########6#########7#########8

source(snakemake@config[["r_libs"]])

load("/mnt/ris/jschwarz/cardiac-radiobiology/data_model/data_model.RData")

salmon_paths_4630 = list.files("/mnt/ris/jschwarz/cardiac-radiobiology/inputs/Rentschler_s4630_MGI0042",
                                         pattern = "\\.sf$",
                                         recursive = TRUE,
                                         full.names = TRUE)

(salmon_paths_4630_alt_ids = gsub("\\..*$",
                                             "",
                                             list.files("/mnt/ris/jschwarz/cardiac-radiobiology/inputs/Rentschler_s4630_MGI0042",
                                         pattern = "\\.sf$",
                                         recursive = TRUE,
                                         full.names = FALSE)
                                         ))

salmon_paths_4730 = list.files("/mnt/ris/jschwarz/cardiac-radiobiology/inputs/Rentschler_s4730_MGI0070",
                                         pattern = "\\.sf$",
                                         recursive = TRUE,
                                         full.names = TRUE)

(salmon_paths_4730_alt_ids = gsub("\\..*$","", list.files("/mnt/ris/jschwarz/cardiac-radiobiology/inputs/Rentschler_s4730_MGI0070",
                                         pattern = "\\.sf$",
                                         recursive = TRUE,
                                         full.names = FALSE)))

salmon_paths_5469 = list.files("/mnt/ris/jschwarz/cardiac-radiobiology/inputs/Rentschler_s5469_MGI2048",
                                         pattern = "\\.sf$",
                                         recursive = TRUE,
                                         full.names = TRUE)

salmon_paths_5469_alt_ids = gsub("-","_",gsub("^","sample.",substr(list.files("/mnt/ris/jschwarz/cardiac-radiobiology/inputs/Rentschler_s5469_MGI2048",
                                         pattern = "\\.sf$",
                                         recursive = TRUE,
                                         full.names = FALSE), 1, 5)))

names(salmon_paths_5469) = as.character(data.frame(alt_lib_id = salmon_paths_5469_alt_ids) %>%
  left_join(libraries, by = "alt_lib_id") %>%
  pull(library_id))

salmon_paths_5708 = list.files("/mnt/ris/jschwarz/cardiac-radiobiology/inputs/Rentschler_s5708_MGI2548",
                                         pattern = "\\.sf$",
                                         recursive = TRUE,
                                         full.names = TRUE)

salmon_paths_5708_alt_ids = paste0("sample.",substr(list.files("/mnt/ris/jschwarz/cardiac-radiobiology/inputs/Rentschler_s5708_MGI2548",
                                         pattern = "\\.sf$",
                                         recursive = TRUE,
                                         full.names = FALSE), 1, 4))

salmon_paths = c(salmon_paths_4630,
                 salmon_paths_4730,
                 salmon_paths_5469,
                 salmon_paths_5708)

alt_ids_salmon_paths = c(salmon_paths_4630_alt_ids,
                         salmon_paths_4730_alt_ids,
                         salmon_paths_5469_alt_ids,
                         salmon_paths_5708_alt_ids)

names(salmon_paths) =
as.character(data.frame(alt_lib_id = alt_ids_salmon_paths) %>%
  left_join(libraries, by = "alt_lib_id") %>% pull(library_id))


txdb = readRDS(snakemake@input[["txdb"]])

txi_nuc = tximport(salmon_paths,
                   type = "salmon",
                   tx2gene = txdb)

saveRDS(txi_nuc, file = snakemake@output[["nuc_gene_cnts"]])
