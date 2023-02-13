#!/usr/bin/env bash
reference=$1
bt2_index_base=$2
output=$3

mkdir -p $output
bowtie2-build \
    $reference \
    $bt2_index_base
