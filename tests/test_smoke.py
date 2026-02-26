# tests/test_smoke.py

import pandas as pd
import pytest

from pypeakranker.summary import init_table, add_signal
from pypeakranker.gc import add_gc


def test_pipeline_init_and_gc(tmp_path):
    peaks = tmp_path / "peaks.bed"
    peaks.write_text("chr1\t0\t10\nchr1\t10\t20\n")

    fasta = tmp_path / "genome.fa"
    fasta.write_text(">chr1\nACGTACGTACGTACGTACGTACGTACGTACGT\n")

    table = tmp_path / "features.tsv"
    init_table(str(peaks), str(table))

    add_gc(
        table_tsv=str(table),
        reference_fasta=str(fasta),
        out_tsv=str(table),
        allow_missing_chroms=True,
    )

    df = pd.read_csv(table, sep="\t")

    assert {"chr", "start", "end", "GC_content"}.issubset(df.columns)
    assert len(df) == 2
    assert df["GC_content"].between(0, 1).all()
    assert df["GC_content"].notna().all()
    assert (df["GC_content"] == 0.5).all()


def test_public_api_imports():
    import pypeakranker

    assert hasattr(pypeakranker, "init_table")
    assert hasattr(pypeakranker, "add_signal")
    assert hasattr(pypeakranker, "add_gc")
    assert hasattr(pypeakranker, "add_phylop")
    assert hasattr(pypeakranker, "add_moments")


def test_add_signal_missing_bigwig_raises(tmp_path):
    table = tmp_path / "table.tsv"
    table.write_text("chr\tstart\tend\nchr1\t0\t10\n")

    # pyBigWig.open raises RuntimeError when it can't open a file
    with pytest.raises((FileNotFoundError, OSError, RuntimeError)):
        add_signal(
            table_tsv=str(table),
            bigwig_files=[str(tmp_path / "nonexistent.bw")],
            out_tsv=str(tmp_path / "out.tsv"),
        )