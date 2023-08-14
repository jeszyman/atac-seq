#!/usr/bin/env bash

peaks_str="${1}"
ref_bed="${2}"
qval_cut="${3}"
threads="${4}"
out_tsv="${5}"


process_file() {
    file="$1"
    base=$(basename $file)
    ref_bed="${2}"
    qval_cut="${3}"
    cat "$file" |
        awk -v cut="$qval_cut" '$9 > cut' |
        sort-bed - |
        bedmap --echo --bases-uniq --delim '\t' $ref_bed - |
        awk -v base="$base" '{print $0 "\t" base}'

}

export -f process_file

parallel -j "$threads" process_file {} "$ref_bed" "$qval_cut" ::: $peaks_str > "$out_tsv"
