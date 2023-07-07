#!/usr/bin/env bash

in_bam="${1}"
in_keep="${2}"
out_bam="${3}"
threads="${4}"

samtools view --bam --with-header -o $out_bam -L $in_keep --use-index --threads $threads $in_bam
samtools index -@ $threads $out_bam
