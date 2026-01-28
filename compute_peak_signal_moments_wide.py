#!/usr/bin/env python3

"""
Compute skewness, kurtosis, and bimodality coefficient per peak across multiple BigWig files (cell types).

Outputs one row per peak, with metrics for each cell type in separate columns.

Author: Saroja Somasundaram
Modified: 2025-07 (peak_idx removed)
"""

import pyBigWig
import pandas as pd
import os
import numpy as np
from glob import glob
from scipy.stats import skew, kurtosis

PEAKS_FILE = "peaks_NC_mapped_filtered.bed"
BIGWIG_DIR = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/data/macaque/ATAC/Group_bigwig_TSSnorm/"
OUTPUT_FILE = "peaks_moments_wide.tsv"

peak_cols = ["chr", "start", "end", "group", "status", "peak_id"]

print(f"Reading peaks from {PEAKS_FILE}")
peaks_df = pd.read_csv(PEAKS_FILE, sep="\t", header=None, names=peak_cols)

# Remove duplicate peaks by genomic coordinates
dup_mask = peaks_df.duplicated(subset=["chr", "start", "end"])
if dup_mask.sum() > 0:
    print(f"Warning: Found {dup_mask.sum()} duplicate peaks, removing.")
    peaks_df = peaks_df.drop_duplicates(subset=["chr", "start", "end"])

# --- Check for peak_id completeness ---
if 'peak_id' not in peaks_df.columns:
    peaks_df["peak_id"] = [f"peak_{i}" for i in range(len(peaks_df))]
else:
    # Convert to string and fill missing or empty strings
    peaks_df["peak_id"] = peaks_df["peak_id"].fillna("").astype(str)
    missing_mask = peaks_df["peak_id"].str.strip() == ""
    peaks_df.loc[missing_mask, "peak_id"] = [f"peak_{i}" for i in peaks_df.index[missing_mask]]

bw_files = sorted(glob(os.path.join(BIGWIG_DIR, "*.bw")))
if not bw_files:
    raise FileNotFoundError(f"No BigWig files found in {BIGWIG_DIR}")

print(f"Found {len(bw_files)} BigWig files: {[os.path.basename(f) for f in bw_files]}")

def get_per_base_signal(bw, chrom, start, end):
    try:
        values = bw.values(chrom, int(start), int(end), numpy=True)
        return np.nan_to_num(values, nan=0.0)
    except Exception as e:
        print(f"Warning: Error reading {chrom}:{start}-{end}: {e}")
        return np.array([])

def bimodality_coefficient(vals):
    n = len(vals)
    if n < 4:
        return np.nan
    s = skew(vals)
    k = kurtosis(vals, fisher=False)
    denom = k + (3*(n-1)**2)/((n-2)*(n-3))
    if denom == 0:
        return np.nan
    return (s**2 + 1) / denom

all_records = []

for bw_path in bw_files:
    celltype = os.path.basename(bw_path).replace(".bw", "")
    print(f"Processing {celltype}")
    with pyBigWig.open(bw_path) as bw:
        for _, row in peaks_df.iterrows():
            chrom, start, end = row["chr"], row["start"], row["end"]
            vals = get_per_base_signal(bw, chrom, start, end)
            if vals.size == 0:
                skew_val = np.nan
                kurt_val = np.nan
                bc_val = np.nan
            elif np.all(vals == 0) or np.all(vals == vals[0]):
                skew_val = 0.0
                kurt_val = 0.0
                bc_val = np.nan
            else:
                skew_val = skew(vals)
                kurt_val = kurtosis(vals, fisher=False)
                bc_val = bimodality_coefficient(vals)
            all_records.append({
                "chr": chrom,
                "start": start,
                "end": end,
                "peak_id": row["peak_id"],
                "celltype": celltype,
                "skewness": skew_val,
                "kurtosis": kurt_val,
                "bimodality_coefficient": bc_val
            })

df_long = pd.DataFrame(all_records)

# Pivot to wide format
index_cols = ["chr", "start", "end", "peak_id"]
df_skew = df_long.pivot(index=index_cols, columns="celltype", values="skewness")
df_skew.columns = [f"{c}_skewness" for c in df_skew.columns]

df_kurt = df_long.pivot(index=index_cols, columns="celltype", values="kurtosis")
df_kurt.columns = [f"{c}_kurtosis" for c in df_kurt.columns]

df_bc = df_long.pivot(index=index_cols, columns="celltype", values="bimodality_coefficient")
df_bc.columns = [f"{c}_bimodality" for c in df_bc.columns]

# Combine metrics
df_metrics = pd.concat([df_skew, df_kurt, df_bc], axis=1).reset_index()

# Merge with original peaks (to retain group, status info)
df_final = peaks_df.merge(df_metrics, on=["chr", "start", "end", "peak_id"], how="left")

print(f"Saving wide-format output to {OUTPUT_FILE}")
df_final.to_csv(OUTPUT_FILE, sep="\t", index=False)

