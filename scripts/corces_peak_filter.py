#!/usr/bin/env python3

"""
Script to impliment Corces et al., 2018 per-peak filtering (https://doi.org/10.1126/science.aav1898)
"""

# ---   Setup   --- #
# ----------------- #

## ---   Load Packages   --- ##
## ------------------------- ##

import argparse
import pandas as pd
import pybedtools
from intervaltree import Interval, IntervalTree

## ---   Load Inputs   --- ##
## ----------------------- ##

def load_inputs():
    # Setup argparse to handle command line arguments
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--blacklist_bed", type=str, default="/mnt/ris/jschwarz/Active/cardiac-radiobiology/ref/mm10-blacklist.v2.bed")
    parser.add_argument("--narrowPeak_file", type=str, default="/mnt/ris/jschwarz/Active/cardiac-radiobiology/analysis/atac/peaks/lib051_mm10_filt_peaks.narrowPeak")
    parser.add_argument("--chrom_size_bed", type=str, default="/mnt/ris/jschwarz/Active/cardiac-radiobiology/ref/mm10_chrom.bed")
    parser.add_argument("--out_narrowPeak", type=str, default="/tmp/test.narroPeak")
    args = parser.parse_args()

    # Load the narrowPeak file as a DataFrame
    narrowPeak_df = pd.read_csv(args.narrowPeak_file, sep='\t', header=None)
    narrowPeak_df.columns = ['chrom', 'start', 'end', 'name', 'score', 'strand', 'signalValue', 'pValue', 'qValue', 'peak']

    return args, narrowPeak_df

# ---   Main   --- #
# ---------------- #

def main():
    # Load inputs
    args, narrowPeak_df = load_inputs()

    # Standardize and filter peaks
    valid_bed = peak_std_and_filt(narrowPeak_df, args.blacklist_bed, args.chrom_size_bed)

    # Process to ensure unique peaks
    final_peaks_df = it_unique_peaks(valid_bed)

    final_peaks_df.to_csv(args.out_narrowPeak, index=False, sep='\t', header=False)


# ---   Functions   --- #
# --------------------- #

def peak_std_and_filt(narrowPeak_df, blacklist_bed, chrom_size_bed):
    """
    """
    summit_std_df = narrowPeak_df
    summit_std_df['std_start'] = narrowPeak_df['start'] + narrowPeak_df['peak'] - 250
    summit_std_df['std_end'] = narrowPeak_df['start'] + narrowPeak_df['peak'] + 250
    summit_std_str = summit_std_df[['chrom', 'std_start', 'std_end', 'name', 'score', 'strand', 'signalValue', 'pValue', 'qValue', 'peak']].to_csv(sep='\t', index=False, header=False)
    summit_bed = pybedtools.BedTool(summit_std_str, from_string=True)

    blacklist_bed = pybedtools.BedTool(blacklist_bed)
    filtered_bed = summit_bed.subtract(blacklist_bed)

    chrom_size_bed = pybedtools.BedTool(chrom_size_bed)
    valid_bed = filtered_bed.intersect(chrom_size_bed, wa=True, u=True, f=1.0)
    num_peaks_initial = len(narrowPeak_df)
    num_peaks_after_blacklist = len(filtered_bed.to_dataframe())
    num_peaks_after_valid = len(valid_bed.to_dataframe())

    print(f"Number of peaks (initial): {num_peaks_initial}")
    print(f"Number of peaks (after blacklist filtering): {num_peaks_after_blacklist}")
    print(f"Number of peaks (after validity check): {num_peaks_after_valid}")

    assert num_peaks_after_blacklist <= num_peaks_initial, "Blacklist filtering increased the number of peaks!"
    assert num_peaks_after_valid <= num_peaks_after_blacklist, "Validity check increased the number of peaks!"

    return valid_bed

def it_unique_peaks(bed):
    # Convert to DataFrame and sort
    bed_df = bed.to_dataframe(names=['chrom', 'start', 'end', 'name', 'score', 'strand', 'signalValue', 'pValue', 'qValue', 'peak'])
    sorted_df = bed_df.sort_values('pValue', ascending=False)

    # Convert to list of tuples
    peaks_list = list(sorted_df.itertuples(index=False, name=None))

    # Initialize variables for interval tree processing
    remaining_peaks = []
    interval_trees = {}

    # Process each peak
    for peak in peaks_list:
        chrom, start, end, *rest = peak
        if chrom not in interval_trees:
            interval_trees[chrom] = IntervalTree()

        overlapping = interval_trees[chrom][start:end]
        if not overlapping:
            interval_trees[chrom].addi(start, end, peak)
            remaining_peaks.append(peak)

    # Create a DataFrame from the non-overlapping peaks
    column_names = ['chrom', 'start', 'end', 'name', 'score', 'strand', 'signalValue', 'pValue', 'qValue', 'peak']
    final_peaks_df = pd.DataFrame(remaining_peaks, columns=column_names)

  # Report and check final number of peaks
    final_num_peaks = len(final_peaks_df)
    print(f"Number of unique peaks (final): {final_num_peaks}")

    # Convert the final DataFrame to a BedTool object
    final_peaks_bedtool = pybedtools.BedTool.from_dataframe(final_peaks_df)

    # Use intersect to find overlaps (with itself), with a minimum overlap of 1bp
    overlaps = final_peaks_bedtool.intersect(final_peaks_bedtool, u=True, f=1.0)

    # Assert that there are no overlaps by comparing the counts
    assert len(overlaps) == len(final_peaks_df), "There are overlapping peaks in the final set!"

    return final_peaks_df

# ---   Main Guard   --- #
# ---------------------- #

if __name__ == "__main__":
    main()
