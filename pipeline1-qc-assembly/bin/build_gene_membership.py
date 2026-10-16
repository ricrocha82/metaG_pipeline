#!/usr/bin/env python3

import argparse
import csv
from pathlib import Path


def read_gene_coordinates(gff_files):
    """Map each Prodigal gene to its contig, sample and coordinates.

    Prodigal names genes differently in its two outputs: FASTA headers (used
    by MMseqs2) are <contig>_<n>, while the GFF ID attribute is
    <contig counter>_<n> (e.g. 1_1, repeated in every sample). Build the
    FASTA-style ID from the GFF contig column and the gene number.

    Only complete genes (partial=00) are kept, matching the genes that are
    clustered into the catalog.
    """
    genes = {}

    for filename in gff_files:
        name = Path(filename).name
        suffix = ".genes.gff"
        sample = name[:-len(suffix)] if name.endswith(suffix) else name

        with open(filename) as handle:
            for line in handle:
                if line.startswith("#"):
                    continue

                fields = line.rstrip("\n").split("\t")

                if len(fields) != 9 or fields[2] != "CDS":
                    continue

                attributes = dict(
                    field.split("=", 1)
                    for field in fields[8].split(";")
                    if "=" in field
                )

                if attributes.get("partial") != "00":
                    continue

                gff_id = attributes.get("ID")

                if not gff_id:
                    raise ValueError(
                        f"Missing gene ID in {filename}"
                    )

                gene_number = gff_id.rsplit("_", 1)[-1]
                gene_id = f"{fields[0]}_{gene_number}"

                if gene_id in genes:
                    raise ValueError(
                        f"Duplicate gene ID: {gene_id}"
                    )

                genes[gene_id] = {
                    "contig_id": fields[0],
                    "source_sample": sample,
                    "start": fields[3],
                    "end": fields[4],
                    "strand": fields[6],
                }

    return genes


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("--clusters", required=True)
    parser.add_argument("--gff", nargs="+", required=True)
    parser.add_argument("--output", required=True)

    args = parser.parse_args()

    genes = read_gene_coordinates(args.gff)

    columns = [
        "catalog_gene_id",
        "original_gene_id",
        "contig_id",
        "source_sample",
        "start",
        "end",
        "strand",
    ]

    with open(args.clusters) as infile, \
         open(args.output, "w", newline="") as outfile:

        reader = csv.reader(infile, delimiter="\t")
        writer = csv.DictWriter(
            outfile, fieldnames=columns, delimiter="\t", lineterminator="\n"
        )

        writer.writeheader()

        seen = set()

        for row in reader:
            if len(row) != 2:
                raise ValueError(
                    f"Invalid MMseqs2 cluster row: {row}"
                )

            representative, original = row

            if original in seen:
                raise ValueError(
                    f"Duplicate cluster member: {original}"
                )

            seen.add(original)

            if original not in genes:
                raise ValueError(
                    f"Gene missing from Prodigal GFF: {original}"
                )

            writer.writerow({
                "catalog_gene_id": representative,
                "original_gene_id": original,
                **genes[original],
            })

        missing = set(genes) - seen

        if missing:
            raise ValueError(
                f"{len(missing)} Prodigal genes are absent "
                "from the MMseqs2 cluster membership table"
            )


if __name__ == "__main__":
    main()