#!/usr/bin/env Rscript

#########1#########2#########3#########4#########5#########6#########7#########8
# Script to annotate MACS2 peak files with the ChIPseeker package. Accepts
#  either broadPeak or narrowPeak files. Returns narrowPeak format (i.e. 10
#  columns including a final peak column) followed by annotation columns from
#  ChIPseeker. Annotation text is simplified for use subsequent computational
#  processing.

# ---   Load Packages   --- #
# ------------------------- #

packages <- c(
  "ChIPseeker",
  "optparse",
  "rtracklayer",
  "tidyverse")

sapply(packages, require, character.only = TRUE, quietly = TRUE)


# ---   Load Inputs   --- #
# ----------------------- #

option_list <- list(
  make_option(c("--in_peak_bed"), type = "character"),
  make_option(c("--out_anno_bed"), type = "character"),
  make_option(c("--txdb"), type = "character")
)

opts <- parse_args(OptionParser(option_list = option_list))

#opts$in_peak_bed = "~/cards/analysis/atac/peaks/lib223_hg38_ds9_peaks.broadPeak"
#opts$out_anno_bed = "/tmp/thetest.bed"
#opts$txdb = "TxDb.Hsapiens.UCSC.hg38.knownGene"

# ---   Main   --- #
# ---------------- #

main = function(opts){

  library(opts$txdb,character.only = TRUE)

  # Annotate
  annotated = annotate(opts$in_peak_bed, opts$txdb)

  # Return a list of outputs
  return(annotated)
}

# ---   Functions   --- #
# --------------------- #

annotate = function(in_peak_bed, txdb){
  peaks = rtracklayer::import(in_peak_bed)
  anno = annotatePeak(peaks, TxDb = get(txdb))
  anno = as_tibble(anno)

  if (!"peak" %in% names(anno)) {
    anno <- anno %>%
      mutate(peak = end - start)
  }

  anno =
    anno %>% mutate(simple = case_when(
                      grepl("Promoter", annotation) ~ "promoter",
                      grepl("Exon", annotation) ~ "exon",
                      grepl("Intron", annotation) ~ "intron",
                      grepl("3' UTR", annotation) ~ "utr3",
                      grepl("5' UTR", annotation) ~ "utr5",
                      grepl("Distal Intergenic", annotation) ~ "intergenic",
                      grepl("Downstream", annotation) ~ "downstream",
                      TRUE ~ annotation
                    ))

  anno =
    anno %>% dplyr::select(seqnames, start, end, width, strand, name, score, signalValue, pValue, qValue, peak, simple, geneChr, geneStart, geneEnd, geneLength, geneStrand, geneId, transcriptId, distanceToTSS, simple)
}

# ---   Run   --- #
# --------------- #

out = main(opts)

write_tsv(out,
          file = opts$out_anno_bed,
          col_names = F)
