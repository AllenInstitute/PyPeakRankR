import argparse
from pypeakranker.summary import add_signal
from pypeakranker.gc import add_gc


def build_parser():
    p = argparse.ArgumentParser(prog="pypeakranker")
    sub = p.add_subparsers(dest="cmd", required=True)

    # init: create a feature table from a BED/TSV peaks file
    p_init = sub.add_parser("init", help="Initialize a feature table from peaks")
    p_init.add_argument("--peaks", required=True, help="Peaks BED/TSV (>=3 cols)")
    p_init.add_argument("--out", required=True, help="Output TSV (feature table)")

    # add-signal: add bigwig summaries to an existing feature table
    p_sig = sub.add_parser("add-signal", help="Add bigWig signal summaries to feature table")
    p_sig.add_argument("--table", required=True, help="Existing feature table TSV")
    p_sig.add_argument("--bigwig-files", nargs="+", required=True, help="BigWig files")
    p_sig.add_argument("--out", required=True, help="Output TSV (can overwrite same file)")
    p_sig.add_argument("--stat", default="sum", choices=["sum", "mean", "max"], help="Summary statistic")
    p_sig.add_argument("--suffix", default="summary", help="Suffix for new columns (e.g., summary)")
    p_sig.add_argument("--allow-missing-chroms", action="store_true", help="Missing chrom -> 0 instead of error")

    # add-gc: add GC_content column to an existing feature table
    p_gc = sub.add_parser("add-gc", help="Add GC content to feature table")
    p_gc.add_argument("--table", required=True, help="Existing feature table TSV")
    p_gc.add_argument("--reference-fasta", required=True, help="Reference FASTA")
    p_gc.add_argument("--out", required=True, help="Output TSV (can overwrite same file)")
    p_gc.add_argument("--allow-missing-chroms", action="store_true", help="Missing chrom -> NA instead of error")

    return p


def main():
    p = build_parser()
    args = p.parse_args()

    if args.cmd == "init":
        # just initialize by reading peaks and writing as TSV
        from pypeakranker.summary import init_table
        init_table(args.peaks, args.out)

    elif args.cmd == "add-signal":
        add_signal(
            table_tsv=args.table,
            bigwig_files=args.bigwig_files,
            out_tsv=args.out,
            stat=args.stat,
            suffix=args.suffix,
            allow_missing_chroms=args.allow_missing_chroms,
        )

    elif args.cmd == "add-gc":
        add_gc(
            table_tsv=args.table,
            reference_fasta=args.reference_fasta,
            out_tsv=args.out,
            allow_missing_chroms=args.allow_missing_chroms,
        )