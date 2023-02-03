#!/usr/bin/env bash

# For unit testing
#in_bam="test/analysis/atac/bams/lib003_dedup.bam"
#out_bam="test/analysis/atac/bams/lib003_filt.bam"

inbam="${1}"
outbam="${2}"
threads="${4}"

samtools view -@ $threads -b -f 1 -h -q 20 -o $outbam $inbam
samtools index $outbam
