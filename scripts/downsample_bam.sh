#!/bin/bash

in_bam="${1}"
milreads="${2}"
threads="${3}"
out_bam="${4}"

reads=$(echo |awk -v var1=$milreads '{ print 1000000*var1 }')

## Calculate the sampling factor based on the intended number of reads:

FACTOR=$(samtools idxstats $in_bam | cut -f3 |awk -v COUNT=$reads 'BEGIN {total=0} {total += $1} END {print COUNT/total}')

samtools view -s $FACTOR -b -@ $threads $in_bam > $out_bam

samtools index -@ threads $out_bam
