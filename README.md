# PyPeakRanker

**PyPeakRanker** is a Python package for extracting quantitative features from a predefined set of genomic peaks and assembling them into a reproducible, analysis-ready table.

It generates a standardized **peak × feature matrix**, enabling systematic ranking and comparison of regulatory elements across cell types, conditions, or species.

PyPeakRanker does **not perform peak calling**. Instead, it standardizes feature extraction so that peak prioritization can be performed reproducibly and transparently using downstream statistical or machine-learning approaches.

Given a fixed set of genomic peaks, PyPeakRanker:

- Extracts quantitative features per peak from signal tracks (e.g., BigWig)
- Computes sequence-based features from reference genomes
- Optionally integrates conservation tracks
- Produces a unified table where rows represent peaks and columns represent features
- Separates feature generation from downstream ranking logic

This design allows users to apply custom scoring functions, statistical tests, or predictive models after feature extraction.

---

# Statement of Need

ATAC-seq and related assays generate large sets of candidate regulatory regions. However, peak prioritization across cell types or biological conditions is often performed using ad hoc scripts with inconsistent feature definitions and aggregation strategies.

Existing tools focus primarily on:

- Peak calling  
- Differential accessibility testing  
- Genomic annotation  

They typically do not provide a standardized framework for **reproducible peak-level feature extraction**.

PyPeakRanker addresses this gap by:

- Systematically aggregating quantitative features for predefined peaks
- Producing a single, analysis-ready feature table
- Enabling transparent and reproducible peak ranking workflows

---

# Supported Feature Modules

PyPeakRanker provides modular feature extraction through CLI subcommands.

## `init`
Initialize a feature table from a peaks file (BED/TSV with chr, start, end).

Creates a clean, deduplicated peak table.

---

## `add-signal`
Summarize BigWig signal tracks over peaks.

Supported summary statistics:
- `sum`
- `mean`
- `max`

Produces one feature column per BigWig file.

---

## `add-gc`
Compute GC fraction per peak from a reference genome FASTA.

Adds:
- `GC_content`

---

## `add-phylop`
Compute mean PhyloP conservation score per peak from a PhyloP BigWig track.

Optional:
- LiftOver support via UCSC chain files for cross-assembly scoring

Adds:
- `phyloP_mean` (customizable column name)

---

## `add-moments`
Compute higher-order distribution statistics across BigWig signal per peak:

- Skewness
- Kurtosis
- Bimodality coefficient

Produces three feature columns per BigWig file.

---

# Installation

Install from source:

```bash
git clone https://github.com/AllenInstitute/PyPeakRankR
cd PyPeakRankR
pip install -e .
```

Or install directly from GitHub:

```
pip install git+https://github.com/AllenInstitute/PyPeakRankR.git
```

## Quick Example

Initialize a feature table:

```
pypeakranker init \
  --peaks peaks.bed \
  --out features.tsv

pypeakranker add-signal \
  --table features.tsv \
  --bigwig-files sample1.bw sample2.bw \
  --stat sum \
  --suffix summary \
  --out features.tsv

pypeakranker add-gc \
  --table features.tsv \
  --reference-fasta genome.fa \
  --out features.tsv

pypeakranker add-phylop \
  --table features.tsv \
  --phylop-bw phyloP.bw \
  --out features.tsv

For cross-assembly scoring (--chain), UCSC liftOver must be installed and available on PATH (or provide --liftover-exe).

pypeakranker add-moments \
  --table features.tsv \
  --bigwig-files sample1.bw sample2.bw \
  --out features.tsv

```

The resulting features.tsv will contain:

Original peak coordinates

One column per signal summary

GC_content

phyloP_mean

Skewness, kurtosis, and bimodality metrics per track

## Design Philosophy

PyPeakRanker explicitly separates:

Feature extraction → deterministic and reproducible
Peak ranking / modeling → user-defined and flexible

This ensures that ranking logic remains transparent and adaptable to specific biological questions.

## Author

Saroja Somasundaram

## License

MIT License. See LICENSE file for details.