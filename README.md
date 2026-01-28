# PyPeakRanker

PyPeakRankR is a Python package for collecting quantitative features for a predefined set of ATAC-seq peaks and assembling them into a reproducible, analysis-ready table.
The resulting peak × feature matrix enables systematic ranking and comparison of regulatory elements across cell types, conditions, or species.

PyPeakRankR does not perform peak calling. Instead, it standardizes feature extraction so that peak prioritization can be performed reproducibly and transparently using downstream ranking or modeling approaches.

Given a fixed set of genomic peaks, PyPeakRankR:

* Extracts multiple quantitative features per peak from ATAC-seq data
* Aggregates features consistently across cell types or groups
* Produces a single unified table where rows are peaks and columns are features
* Enables reproducible ranking of peaks across biological contexts

This design separates feature generation from ranking logic, allowing users to apply custom scoring functions, statistical tests, or machine-learning models downstream.

# Statement of Need
ATAC-seq experiments generate large sets of candidate regulatory regions, yet peak prioritization across cell types or conditions is often performed using ad hoc scripts with inconsistent feature definitions and normalization. This limits reproducibility and cross-study comparability. Existing tools focus on peak calling, differential accessibility testing, or annotation, but lack a standardized framework for reproducible peak-level feature extraction. PeakRankR addresses this gap by providing a Python package that systematically aggregates quantitative features for predefined ATAC-seq peaks into a single, analysis-ready table, enabling transparent and reproducible peak ranking and comparative analysis.

## Contents
Features collected:

- ATAC specificity
- Seq conservation - PhyloP score
- GC content
- TSS distance
- Peak skewness
- Peak kurtosis
- Peak bimodality
- Gene marker score

The framework is designed to be easily extended with additional peak-level features.
 
## Author
Saroja Somasundaram

## Acknowledgements

OpenAI was used as a development aid.


