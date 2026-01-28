#!/usr/bin/env python3

"""
Summarize signal from multiple BigWig files over a set of genomic peaks.

Inputs:
    - peaks.bed: BED file with at least 3 columns: chr, start, end
    - directory of BigWig files (*.bw)

Outputs:
    - peaks_summary_all_bw.tsv: Tab-separated file with summed signal values for each peak per sample

Author: Saroja Somasundaram
Date: 2025-07
"""

import pyBigWig
import pandas as pd
import os
import numpy as np
from glob import glob

# === Configuration ===
PEAKS_FILE = "peaks_NC_mapped_filtered.bed"
BIGWIG_DIR = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/data/macaque/ATAC/Group_bigwig_TSSnorm/"
OUTPUT_FILE = "peaks_summary_all_bw.tsv"

# === Load and validate input peaks ===
print(f"Reading peaks from {PEAKS_FILE}")
peaks_df = pd.read_csv(PEAKS_FILE, sep="\t", header=None)

# Ensure it has at least 3 columns
if peaks_df.shape[1] < 3:
    raise ValueError("peaks.bed must have at least 3 columns (chr, start, end)")

# Truncate to max 6 columns
peaks_df = peaks_df.iloc[:, :6]

# Assign column names dynamically
default_cols = ["chr", "start", "end", "group", "status", "peak_id"]
assigned_cols = default_cols[:peaks_df.shape[1]]
peaks_df.columns = assigned_cols

# Ensure 'id' column exists and is complete
if "peak_id" not in peaks_df.columns:
    peaks_df["peak_id"] = [f"peak_{i}" for i in range(len(peaks_df))]
else:
    # Fill missing or empty id values
    peaks_df["peak_id"] = peaks_df["peak_id"].replace("", pd.NA)
    missing_ids = peaks_df["peak_id"].isna()
    if missing_ids.sum() > 0:
        print(f"Warning: {missing_ids.sum()} missing 'id' values found. Filling them with 'peak_<index>'.")
        peaks_df.loc[missing_ids, "peak_id"] = [f"peak_{i}" for i in peaks_df.index[missing_ids]]

# Convert start and end to integers
peaks_df["start"] = pd.to_numeric(peaks_df["start"], errors="coerce")
peaks_df["end"] = pd.to_numeric(peaks_df["end"], errors="coerce")
peaks_df = peaks_df.dropna(subset=["start", "end"])
peaks_df["start"] = peaks_df["start"].astype(int)
peaks_df["end"] = peaks_df["end"].astype(int)

# Drop duplicates
dup_mask = peaks_df.duplicated(subset=["chr", "start", "end"])
num_dups = dup_mask.sum()
if num_dups > 0:
    print(f"Warning: Found {num_dups} duplicated peak(s), removing duplicates!")
    peaks_df = peaks_df.drop_duplicates(subset=["chr", "start", "end"])

# === List BigWig files ===
bw_files = sorted(glob(os.path.join(BIGWIG_DIR, "*.bw")))
if len(bw_files) == 0:
    raise FileNotFoundError(f"No BigWig files found in directory: {BIGWIG_DIR}")

print(f"Found {len(bw_files)} BigWig file(s): {[os.path.basename(f) for f in bw_files]}")

# === Function: summarize signal over a peak ===
def summarize_peak(bw, chrom, start, end):
    try:
        values = bw.values(chrom, int(start), int(end), numpy=True)
        values = values[~np.isnan(values)]  # Remove NaNs
        return values.sum() if len(values) > 0 else 0
    except Exception as e:
        print(f"Warning: Failed to retrieve values for {chrom}:{start}-{end} — {e}")
        return 0

# === Main summarization ===
results = peaks_df[["chr", "start", "end", "group", "status", "peak_id"]].copy()

for bw_path in bw_files:
    sample_name = os.path.basename(bw_path).replace(".bw", "")
    print(f"Processing {sample_name}")
    
    with pyBigWig.open(bw_path) as bw:
        signal = peaks_df.apply(
            lambda row: summarize_peak(bw, row["chr"], row["start"], row["end"]),
            axis=1
        )
        results[sample_name] = signal

# === Save output ===
print(f"Saving summary to {OUTPUT_FILE}")
results.to_csv(OUTPUT_FILE, sep="\t", index=False)

print("Done.")

