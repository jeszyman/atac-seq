#!/usr/bin/env python3
"""
Merge and filter MACS2 peaks according to Corces et al., 2018 (https://doi.org/10.1126/science.aav1898)
"""

# ---   Setup   --- #
# ----------------- #

## ---   Load Packages   --- ##
## ------------------------- ##

import argparse
from intervaltree import IntervalTree
import pandas as pd
import pybedtools

# ---   Load Inputs   --- #
# ----------------------- #


def load_inputs():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--peak_list", nargs="+")
    parser.add_argument("--overlap", type=float)
    parser.add_argument("--count", type=int)
    parser.add_argument("--out_peakset")
    args = parser.parse_args()
    peak_df_list = [pd.read_csv(file, sep="\t", header=None) for file in args.peak_list]
    return (args, peak_df_list)


# ---   Main   --- #
# ---------------- #


def main():
    # Load inputs
    args, peak_df_list = load_inputs()

    # Merge peaks
    merged = merge_peaks(peak_df_list)

    # Iteratively filter by score per million
    score_filt = iterative_score_filter(merged)

    # Identify reproducible peaks
    reproducible = filter_by_count(score_filt, peak_df_list, args.overlap, args.count)

    # Write output
    reproducible.to_csv(args.out_peakset, index=False, sep="\t")


# ---   Functions   --- #
# --------------------- #


def merge_peaks(peak_list):
    merged_df = pd.concat(peak_list, ignore_index=True)
    column_names = [
        "chrom",
        "start",
        "end",
        "name",
        "score",
        "strand",
        "signalValue",
        "pValue",
        "qValue",
        "peak",
        "spm",
    ]
    merged_df.columns = column_names
    return merged_df


def iterative_score_filter(bed_df):
    # Convert to DataFrame and sort
    sorted_df = bed_df.sort_values("spm", ascending=False)

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
    column_names = [
        "chrom",
        "start",
        "end",
        "name",
        "score",
        "strand",
        "signalValue",
        "pValue",
        "qValue",
        "peak",
        "spm",
    ]
    filt_peaks_df = pd.DataFrame(remaining_peaks, columns=column_names)

    # Report and check final number of peaks
    final_num_peaks = len(filt_peaks_df)
    print(f"Final peak set: {final_num_peaks}")

    # Convert the final DataFrame to a BedTool object
    filt_peaks_bedtool = pybedtools.BedTool.from_dataframe(filt_peaks_df)

    # Use intersect to find overlaps (with itself), with a minimum overlap of 1bp
    overlaps = filt_peaks_bedtool.intersect(filt_peaks_bedtool, u=True)

    # Assert that there are no overlaps by comparing the counts
    assert len(overlaps) == len(
        filt_peaks_df
    ), "There are overlapping peaks in the final set!"

    return filt_peaks_df


def filter_by_count(iterative_filt_df, original_df_list, overlap, count):
    """
    For iteratively filtered peaks, counts overlaps with original peak files. Retains based on counts at a set reciprocal overlap.
    """

    # Convert to bedfiles
    iterative_filt_bed = pybedtools.BedTool.from_dataframe(iterative_filt_df)
    original_bed_list = [
        pybedtools.BedTool.from_dataframe(df) for df in original_df_list
    ]

    # Intersect, counting overlaps
    intersected = iterative_filt_bed.intersect(
        b=original_bed_list, wa=True, f=overlap, r=True, c=True
    )

    # Convert back to dataframe, label, filter by overlap count
    intersected_df = intersected.to_dataframe()
    column_names = [
        "chrom",
        "start",
        "end",
        "name",
        "score",
        "strand",
        "signalValue",
        "pValue",
        "qValue",
        "peak",
        "spm",
        "overlap_count",
    ]
    intersected_df.columns = column_names
    filtered_df = intersected_df[intersected_df["overlap_count"] >= count]

    return filtered_df


# ---   Main Guard   --- #
# ---------------------- #

if __name__ == "__main__":
    main()
