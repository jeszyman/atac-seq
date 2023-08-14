#!/usr/bin/env bash

input_str="${1}"
output="${2}"

#echo "library	reads	frip" > $output

IFS=' ' read -r -a array <<< "$input_str"

# Define a function to process each file
process_file() {
    file="$1"
    library=$(echo $(basename $file) | sed 's/_.*$//g')
    peaks=sed 's|/bams/lib[0-9]\{3\}_mm10_ds9.bam|/peaks/&_multi_peaks.narrowPeak|' filename
    col1=$(bedtools intersect -a $file -b /mnt/ris/jschwarz/Active/cardiac-radiobiology/analysis/atac/${species}/peaks/${library}_${build}_ds9_multi_peaks.narrowPeak -wa -f 0.2 -r | wc -l)
    col2=$(cat $file | wc -l)
    frip=$(awk 'BEGIN {print '${col1}'/'${col2}' }')
    echo "$library	$col2	$frip"
}

# Export the function to make it available to parallel
export -f process_file

# Run the loop in parallel using parallel command
parallel -j 4 process_file ::: "${array[@]}" > "$output"

sed -i '1i\library\tpeaks\tfrip' $output
