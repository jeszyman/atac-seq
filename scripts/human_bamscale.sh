#!/usr/bin/env bash

bam_dir="${1}"
input_bed="${2}"
out_dir="${3}"

bam_files=($(ls "$bam_dir"/*.bam))

# build the BAMscale command with the --bam flags
bams=""
for bam in "${bam_files[@]}"
do
    bams+="--bam $bam "
done

BAMscale cov --bed $input_bed \
         --prefix human_hg38 --outdir $out_dir --threads 16 \
         $bams
