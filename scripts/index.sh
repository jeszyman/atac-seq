#!/usr/bin/env bash
reference=$1
bt2_index_base=$2
dir="${3}"
threads="${4}"

mkdir -p $dir

bowtie2-build \
    --threads $threads \
    $reference \
    $bt2_index_base
