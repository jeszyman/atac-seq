#!/usr/bin/env bash
peaks_str="${1}"
ref_bed="${2}"
state="${3}"
qval_cut="${4}"
threads="${5}"
out_tsv="${6}"

process_file() {
    file="$1"
    base=$(basename $file)
    state="${2}"
    ref_bed="${3}"
    qval_cut="${4}"

    # Check if the state is "all"; if not, apply state filtering
    if [ "$state" = "all" ]; then
        cat "$file" |
            awk -v cut="$qval_cut" '$8 > cut' |
            sort-bed - |
            bedmap --echo --bases-uniq --delim '\t' $ref_bed - |
            awk -v base="$base" '{print $0 "\t" base}'
    else
        cat "$file" |
            awk -v state="$state" '$10 == state' |
            awk -v cut="$qval_cut" '$8 > cut' |
            sort-bed - |
            bedmap --echo --bases-uniq --delim '\t' $ref_bed - |
            awk -v base="$base" '{print $0 "\t" base}'
    fi
}

export -f process_file

parallel -j "$threads" process_file {} "$state" "$ref_bed" "$qval_cut" ::: $peaks_str > "$out_tsv"
