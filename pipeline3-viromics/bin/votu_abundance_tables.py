#!/usr/bin/env python3
"""
Merge per-sample CoverM contig tables (<sample>.coverm.tsv, sparse format)
into vOTU x sample tables, applying the covered-fraction rule.

CoverM is run with --min-covered-fraction 0 (it refuses count/length/
reads_per_base otherwise), so the rule is applied here: when a vOTU covers
< --min-covered-fraction of its length in a sample, its depth/count values
(Mean, Trimmed Mean, Variance, Read Count, Reads per base, RPKM) are set to 0.
TPM is recomputed from the thresholded Read Count / Length so it still sums
to 1e6 per sample. Length, Covered Bases and Covered Fraction stay raw, so
you can see why a vOTU was zeroed.

Writes:
  <prefix>_<metric>.tsv   one vOTU x sample matrix per CoverM metric
  <prefix>_long.tsv       sample, vOTU, raw values of every metric, and
                          passes_covered_fraction (True/False)
"""
import argparse
import csv
from pathlib import Path

RAW_METRICS = {"Length", "Covered Bases", "Covered Fraction"}


def file_name(metric):
    return metric.lower().replace(" ", "_")


def to_float(value):
    try:
        return float(value)
    except ValueError:
        return 0.0


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--min-covered-fraction", type=float, required=True,
                        help="fraction (0.70) or percent (70)")
    parser.add_argument("--out-prefix", required=True)
    parser.add_argument("tables", nargs="+", type=Path)
    args = parser.parse_args()
    threshold = args.min_covered_fraction / 100 if args.min_covered_fraction > 1 else args.min_covered_fraction

    raw = {}        # (sample, vOTU) -> {metric: value}
    metrics = None
    samples = []
    for table in args.tables:
        # Sample name from the file (CoverM's own Sample column is <ref>/<read1>)
        sample = table.name[: -len(".coverm.tsv")]
        samples.append(sample)
        with open(table) as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            table_metrics = [c for c in reader.fieldnames if c not in ("Sample", "Contig")]
            if metrics is None:
                metrics = table_metrics
            if "Covered Fraction" not in table_metrics:
                raise SystemExit(f"{table}: needs the covered_fraction method for the covered-fraction rule")
            for row in reader:
                raw[(sample, row["Contig"])] = {m: row[m] for m in metrics}

    samples.sort()
    votus = sorted({votu for _, votu in raw})

    # Long table: raw values + whether the vOTU passes the rule
    passes = {}
    with open(f"{args.out_prefix}_long.tsv", "w") as out:
        out.write("\t".join(["sample", "vOTU"] + metrics + ["passes_covered_fraction"]) + "\n")
        for (sample, votu), values in sorted(raw.items()):
            ok = to_float(values["Covered Fraction"]) >= threshold
            passes[(sample, votu)] = ok
            out.write("\t".join([sample, votu] + [values[m] for m in metrics] + [str(ok)]) + "\n")

    # Thresholded values; vOTUs absent from a sparse table are 0 (no reads)
    def value(sample, votu, metric):
        values = raw.get((sample, votu))
        if values is None:
            return 0.0
        if metric not in RAW_METRICS and not passes[(sample, votu)]:
            return 0.0
        return to_float(values[metric])

    # TPM recomputed after thresholding: (count / length) scaled to 1e6 per sample
    tpm = {}
    if "Read Count" in metrics and "Length" in metrics:
        for sample in samples:
            rate = {}
            for votu in votus:
                length = value(sample, votu, "Length")
                rate[votu] = value(sample, votu, "Read Count") / length if length else 0.0
            total = sum(rate.values())
            for votu in votus:
                tpm[(sample, votu)] = rate[votu] / total * 1e6 if total else 0.0

    for metric in metrics:
        with open(f"{args.out_prefix}_{file_name(metric)}.tsv", "w") as out:
            out.write("vOTU\t" + "\t".join(samples) + "\n")
            for votu in votus:
                if metric == "TPM" and tpm:
                    row = [tpm[(s, votu)] for s in samples]
                else:
                    row = [value(s, votu, metric) for s in samples]
                out.write(votu + "\t" + "\t".join(f"{v:g}" for v in row) + "\n")

    n_zeroed = sum(1 for ok in passes.values() if not ok)
    print(f"{len(votus)} vOTUs x {len(samples)} samples; "
          f"{n_zeroed} sample-vOTU pairs below covered fraction {threshold} set to 0")


if __name__ == "__main__":
    main()
