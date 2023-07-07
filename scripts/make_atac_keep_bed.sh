#!/usr/bin/env bash

hg38_chrome_sizes="${1}"
hg38_atac_keep_bed="${2}"

cat $hg38_chrome_sizes |
    # Grep out all the non-canonical contigs and the mitochondrial reads
    grep -vE 'chrM|_|\*' |
    # Convert to bedfile format
    awk -v FS='\t' -v OFS='\t' '$2 = "1" FS $2' > $hg38_atac_keep_bed
