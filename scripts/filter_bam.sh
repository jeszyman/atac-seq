#!/usr/bin/env bash

# For unit testing
in_bam="test/analysis/atac/bams/lib003_dedup.bam"
out_bam="test/analysis/atac/bams/lib003_filt.bam"
bed="test/ref/"

inbam="${1}"
outbam="${2}"
bed="${3}"
threads="${4}"

samtools view -@ $threads -b -f 3 -h -L $bed -M -q 20 -o $outbam $inbam
samtools index $outbam
