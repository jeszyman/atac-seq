#!/usr/bin/env bash

# Script variables
raw_bam="${1}"
dedup_bam="${2}"
threads="${3}"

samtools sort -@ $threads -n -o - $raw_bam |
    samtools fixmate -m - - |
    samtools sort -@ $threads -o - - |
    samtools markdup -@ $threads -r - $dedup_bam
samtools index $dedup_bam
