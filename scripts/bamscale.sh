#!/usr/bin/env bash

bams="${1}"
bais="${2}"
bed="${3}"
tmp_dir="${4}"
atac_set="${5}"
out_dir="${6}"

rm -rf $tmp_dir
mkdir -p $tmp_dir

echo $bams $bais | tr ' ' '\n' | parallel --max-args 1 cp {} $tmp_dir

# set the directory containing the input BAM files
bam_dir=$tmp_dir
# get a list of BAM files in the directory
bam_files=($(ls "$bam_dir"/*.bam))
# build the BAMscale command with the --bam flags
bams=""
for bam in "${bam_files[@]}"
do
bams+="--bam $bam "
done

BAMscale cov --bed $bed --outdir $out_dir --threads 16 $bams
rm -rf $tmp_dir
