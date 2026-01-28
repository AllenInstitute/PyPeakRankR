#!/usr/bin/env python3

"""
Compute GC content per peak from reference genome fasta.

Inputs:
  - peaks.bed: BED file with chr, start, end (at least 3 columns)
  - reference.fa: Reference genome fasta (indexed by pyfaidx)

Outputs:
  - peaks_gc_content.tsv: tab-separated with GC content per peak (fraction 0-1)

Author: Saroja Somasundaram
Date: 2025-07
"""

import pandas as pd
from pyfaidx import Fasta

# === Configuration ===
PEAKS_FILE = "peaks_NC_mapped_filtered.bed"
REFERENCE_FASTA = "/allen/programs/celltypes/workgroups/rnaseqanalysis/references/macaque/ncbi/mmul10/genome/fasta/genome.fa"
OUTPUT_FILE = "peaks_gc_content.tsv"

# Expected minimum columns
min_cols = ["chr", "start", "end"]

# Load peaks
peaks_df = pd.read_csv(PEAKS_FILE, sep="\t", header=None)

# Assign standard column names based on number of columns
col_count = peaks_df.shape[1]
if col_count < 3:
    raise ValueError("BED file must have at least 3 columns: chr, start, end.")

# Extend column names
all_cols = min_cols + ["group", "status", "peak_id"]
peaks_df.columns = all_cols[:col_count]

# If 'peak_id' is missing, or has missing/incomplete values, fill them
if 'peak_id' not in peaks_df.columns:
    peaks_df["peak_id"] = [f"peak_{i}" for i in range(len(peaks_df))]
else:
    peaks_df["peak_id"] = peaks_df["peak_id"].fillna("").astype(str)
    peaks_df.loc[peaks_df["peak_id"].str.strip() == "", "peak_id"] = [
        f"peak_{i}" for i in peaks_df.index[peaks_df["peak_id"].str.strip() == ""]
    ]



# Load fasta with pyfaidx
fasta = Fasta(REFERENCE_FASTA)

def calc_gc(seq):
    seq = seq.upper()
    gc_count = seq.count('G') + seq.count('C')
    return gc_count / len(seq) if len(seq) > 0 else 0

# Compute GC content
gc_list = []
for idx, row in peaks_df.iterrows():
    chrom, start, end = row["chr"], int(row["start"]), int(row["end"])
    try:
        seq = fasta[chrom][start:end].seq
        gc = calc_gc(seq)
    except KeyError:
        print(f"Warning: Chromosome {chrom} not found in fasta")
        gc = None
    except Exception as e:
        print(f"Warning: Error fetching sequence {chrom}:{start}-{end} - {e}")
        gc = None
    gc_list.append(gc)

peaks_df["GC_content"] = gc_list

# Save output
peaks_df.to_csv(OUTPUT_FILE, sep="\t", index=False)
print(f"GC content saved to {OUTPUT_FILE}")

