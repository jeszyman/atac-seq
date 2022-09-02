#!/usr/bin/env bash

inbam="${1}"
outbam="${2}"
bed="${3}"
threads="${4}"
samtools view -@ $threads -b -f 3 -h -L $bed -M -q 20 -o $outbam $inbam
samtools index $outbam
