#!/usr/bin/env python3

import argparse
import csv
from pathlib import Path


def read_fasta_ids(fasta):
    """Read representative gene IDs in catalog order."""
    ids = []
    seen = set()

    with open(fasta) as handle:
        for line in handle:
            if line.startswith(">"):
                gene = line[1:].split()[0]

                if gene in seen:
                    raise ValueError(f"Duplicate catalog ID: {gene}")

                seen.add(gene)
                ids.append(gene)

    if not ids:
        raise ValueError("Gene catalog is empty")

    return ids


def find_column(header, method):
    """Find the CoverM column corresponding to a method."""
    matches = [
        name for name in header[1:]
        if name.lower().split()[-1] == method.lower()
    ]

    if len(matches) != 1:
        raise ValueError(
            f"Expected one {method} column; found {matches}"
        )

    return matches[0]


def read_coverm(path):
    """Read one sample's gene counts, mean coverage and TPM."""
    result = {}

    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        if not reader.fieldnames:
            raise ValueError(f"Empty CoverM table: {path}")

        gene_col = reader.fieldnames[0]
        count_col = find_column(reader.fieldnames, "Count")
        mean_col = find_column(reader.fieldnames, "Mean")
        tpm_col = find_column(reader.fieldnames, "TPM")

        for row in reader:
            gene = row[gene_col]

            if gene in result:
                raise ValueError(
                    f"Duplicate gene {gene} in {path}"
                )

            result[gene] = (
                float(row[count_col]),
                float(row[mean_col]),
                float(row[tpm_col]),
            )

    return result


def write_matrix(path, genes, samples, values):
    """Write a gene-by-sample TSV matrix."""
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["gene_id", *samples])

        for gene in genes:
            writer.writerow([
                gene,
                *[
                    f"{values[sample][gene]:.6f}"
                    for sample in samples
                ],
            ])


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("--catalog", required=True)
    parser.add_argument("--tables", nargs="+", required=True)
    parser.add_argument("--outdir", default=".")

    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    genes = read_fasta_ids(args.catalog)
    gene_set = set(genes)

    counts = {}
    coverage = {}
    tpm = {}

    suffix = ".gene_abundance.tsv"

    for filename in args.tables:
        path = Path(filename)

        if not path.name.endswith(suffix):
            raise ValueError(f"Unexpected filename: {path.name}")

        sample = path.name[:-len(suffix)]

        if sample in counts:
            raise ValueError(f"Duplicate sample: {sample}")

        data = read_coverm(path)

        unknown = set(data) - gene_set

        if unknown:
            raise ValueError(
                f"{sample}: {len(unknown)} unknown gene IDs"
            )

        counts[sample] = {
            gene: data.get(gene, (0.0, 0.0, 0.0))[0]
            for gene in genes
        }

        coverage[sample] = {
            gene: data.get(gene, (0.0, 0.0, 0.0))[1]
            for gene in genes
        }

        tpm[sample] = {
            gene: data.get(gene, (0.0, 0.0, 0.0))[2]
            for gene in genes
        }

    if not counts:
        raise ValueError("No CoverM tables provided")

    samples = sorted(counts)

    cpm = {}

    for sample in samples:
        total = sum(counts[sample].values())

        cpm[sample] = {
            gene: (
                counts[sample][gene] * 1_000_000 / total
                if total > 0 else 0.0
            )
            for gene in genes
        }

    write_matrix(
        outdir / "gene_catalog_counts.tsv",
        genes, samples, counts
    )

    write_matrix(
        outdir / "gene_catalog_cpm.tsv",
        genes, samples, cpm
    )

    write_matrix(
        outdir / "gene_catalog_mean_coverage.tsv",
        genes, samples, coverage
    )

    write_matrix(
        outdir / "gene_catalog_tpm.tsv",
        genes, samples, tpm
    )


if __name__ == "__main__":
    main()