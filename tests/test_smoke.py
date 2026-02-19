import os
import pandas as pd
from pypeakranker.summary import init_table, add_signal
from pypeakranker.gc import add_gc


def test_pipeline(tmp_path):
    # Create tiny peaks file
    peaks = tmp_path / "peaks.bed"
    with open(peaks, "w") as f:
        f.write("chr1\t0\t10\n")
        f.write("chr1\t10\t20\n")

    # Create tiny fake FASTA
    fasta = tmp_path / "genome.fa"
    with open(fasta, "w") as f:
        f.write(">chr1\n")
        f.write("ACGTACGTACGTACGTACGTACGTACGTACGT\n")

    # Initialize table
    table = tmp_path / "features.tsv"
    init_table(str(peaks), str(table))

    # Add GC
    add_gc(
        table_tsv=str(table),
        reference_fasta=str(fasta),
        out_tsv=str(table),
        allow_missing_chroms=True,
    )

    df = pd.read_csv(table, sep="\t")

    # Check expected structure
    assert "chr" in df.columns
    assert "start" in df.columns
    assert "end" in df.columns
    assert "GC_content" in df.columns
    assert len(df) == 2