---
title: "PyPeakRanker: Reproducible Peak-Level Feature Extraction for Regulatory Element Ranking"
tags:
  - Python
  - ATAC-seq
  - genomics
  - regulatory elements
  - bioinformatics
  - BigWig
  - conservation
authors:
  - name: Saroja Somasundaram
    affiliation: 1
affiliations:
  - name: Allen Institute for Brain Science
    index: 1
date: 2026-02-26
bibliography: references.bib
---

# Summary

High-throughput chromatin accessibility assays such as ATAC-seq generate large sets of candidate regulatory elements. Downstream analyses often require prioritizing peaks across cell types, experimental conditions, or species. However, peak ranking workflows are frequently implemented using ad hoc scripts, with inconsistent feature definitions and aggregation strategies, limiting reproducibility and cross-study comparability.

**PyPeakRanker** is a Python package that standardizes quantitative feature extraction for predefined genomic peaks. It produces a reproducible peak × feature matrix that can be used for downstream statistical analysis, ranking, or machine-learning applications. By separating feature generation from ranking logic, PyPeakRanker promotes transparency, modularity, and reproducibility in regulatory element analysis.

# Statement of Need

While many tools exist for:

- Peak calling
- Differential accessibility analysis
- Genomic annotation

there is limited support for standardized, reproducible extraction of peak-level quantitative features across heterogeneous signal tracks.

Peak prioritization typically requires combining:

- Signal intensity summaries from BigWig tracks
- Sequence-based features (e.g., GC content)
- Conservation metrics (e.g., PhyloP)
- Higher-order distribution statistics

These features are often computed using custom scripts that vary across projects and laboratories. Such variability complicates benchmarking, replication, and cross-study integration.

PyPeakRanker addresses this gap by providing:

- A modular, CLI-driven and API-accessible framework
- Deterministic aggregation of signal and sequence features
- Support for conservation scoring with optional cross-assembly liftOver
- A unified, analysis-ready feature table

This design allows researchers to focus on downstream modeling and biological interpretation while ensuring consistent and reproducible feature computation.

# Functionality

PyPeakRanker operates on a predefined set of genomic peaks and supports modular feature extraction:

## Signal Summaries

Signal tracks stored in BigWig format can be summarized over peak intervals using:

- Sum
- Mean
- Maximum

One feature column is produced per BigWig track.

## Sequence-Based Features

GC content per peak is computed directly from a reference genome FASTA file.

## Conservation Scoring

Mean PhyloP conservation scores can be computed from a PhyloP BigWig track.  
Optional support for UCSC liftOver enables scoring peaks across genome assemblies.

## Distribution Moments

Higher-order statistics describing signal distribution within peaks include:

- Skewness
- Kurtosis
- Bimodality coefficient

These metrics provide additional structural information about regulatory signal profiles.

# Design Philosophy

PyPeakRanker explicitly separates:

**Feature extraction** — deterministic, standardized, and reproducible  
**Peak ranking / modeling** — user-defined and flexible  

By decoupling these steps, PyPeakRanker enables:

- Transparent benchmarking of ranking strategies
- Integration with statistical pipelines
- Use within machine-learning workflows
- Cross-species and cross-condition comparisons

# Implementation

PyPeakRanker is implemented in Python (>=3.9) and relies on widely used scientific libraries:

- `pandas` for tabular data processing
- `numpy` for numerical computation
- `pyBigWig` for BigWig signal extraction
- `pyfaidx` for FASTA access
- `scipy` for statistical metrics

The package provides:

- A command-line interface (`pypeakranker`)
- A programmatic Python API
- Unit tests for core functionality
- Installable packaging via `pyproject.toml`

# Example Workflow

A typical workflow includes:

1. Initializing a feature table from peaks
2. Adding signal summaries from BigWig tracks
3. Computing GC content
4. Adding conservation scores
5. Computing distribution moments

The resulting feature table can then be used for:

- Ranking regulatory elements
- Cell-type specificity analysis
- Comparative genomics
- Machine-learning model input

# Availability

- Source code: https://github.com/AllenInstitute/PyPeakRankR
- License: MIT
- Python version: >=3.9

# Acknowledgements

Development was assisted by AI-based coding tools and informed by regulatory genomics workflows developed at the Allen Institute for Brain Science.

# References