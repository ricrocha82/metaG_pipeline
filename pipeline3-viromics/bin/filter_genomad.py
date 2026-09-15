#!/usr/bin/env python3
"""
Filter geNomad table by thresholds and output contig IDs.

Usage:
  filter_vs2_genomad.py --input genomad_summary.tsv --sample SAMPLE_ID [--outdir .]

Writes:
  SAMPLE_ID_contigs.txt in the output directory.
"""

import argparse
from pathlib import Path
import pandas as pd
import sys

def main():
    parser = argparse.ArgumentParser(description="Filter VirSorter2 + geNomad combined table.")
    parser.add_argument(
        "--input", "-i", required=True, type=Path,
        help="Combined TSV file from VirSorter2 + geNomad"
    )
    parser.add_argument(
        "--sample", "-s", required=True,
        help="Sample ID to use for naming the output file"
    )
    parser.add_argument(
        "--outdir", "-o", default=".",
        help="Output directory (default: current directory)"
    )
    parser.add_argument(
        "--length", type=int, default=5000,
        help="Minimum contig length threshold (default: 5000)"
    )
    parser.add_argument(
        "--virus-score", type=float, default=0.7,
        help="Minimum virus_score threshold (default: 0.7)"
    )
    args = parser.parse_args()

    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    out_file = outdir / f"{args.sample}_contigs.txt"

    # Read TSV
    try:
        df = pd.read_csv(args.input, sep="\t", dtype=str)
    except Exception as e:
        print(f"[ERROR] Could not read {args.input}: {e}", file=sys.stderr)
        sys.exit(1)

    # Normalize contig column name
    if "contig" not in df.columns:
        for cand in ("seq_name", "seqname", "sequence"):
            if cand in df.columns:
                df = df.rename(columns={cand: "contig"})
                break
    if "contig" not in df.columns:
        print(f"[ERROR] Could not find a contig column in {args.input}", file=sys.stderr)
        sys.exit(1)

    # Safely coerce to numeric
    ln = pd.to_numeric(df.get("length"), errors="coerce")
    vs = pd.to_numeric(df.get("virus_score"), errors="coerce")

    mask = (ln >= args.length) & (vs >= args.virus_score)

    contigs = df.loc[mask, "contig"].dropna()

    if contigs.empty:
        print(f"[WARN] No contigs passed the filters for {args.sample}", file=sys.stderr)

    contigs.to_csv(out_file, index=False, header=False)
    print(f"[OK] Wrote {len(contigs)} contigs to {out_file}")

if __name__ == "__main__":
    main()