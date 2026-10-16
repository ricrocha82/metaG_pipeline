#!/usr/bin/env python3
"""
Pick the candidate viral sequences of one sample, record which tool found
each one, and carry over every score/annotation column the tools report.

  geNomad only : geNomad's viral sequences (proviruses already excised)
  VS2 only     : VS2's viral sequences (whole contigs, --keep-original-seq)
  both         : contigs found by both tools -> original contig (CheckV trims
                 host regions later); found by one tool -> that tool's sequence

Sequence IDs are made safe for BLAST and kept traceable:
  - '|' becomes '__' (geNomad 'contig|provirus_1_9000' -> 'contig__provirus_1_9000')
  - VS2 '||full' / '||lt2gene' suffixes are removed
  - IDs are prefixed with '<sample>_' unless they already are, so IDs stay
    unique once all samples are pooled

Columns are fixed (NA when a tool was not run or gave no value), so the
per-sample tables can be concatenated:
  seq_id, contig_id, sample, detected_by, length, genomad_score,
  genomad_<virus_summary columns>, genomad_agg_<whole-contig scores>,
  vs2_<final-viral-score columns>
A contig with several geNomad proviruses or VS2 rows gets their values
joined with ';' when the whole contig is used (found by both tools).
genomad_score is the score the filtering uses: the provirus score for
excised proviruses, otherwise the whole-contig score.
"""
import argparse
import csv

GENOMAD_SUMMARY = ["seq_name", "length", "topology", "coordinates", "n_genes", "genetic_code",
                   "virus_score", "fdr", "n_hallmarks", "marker_enrichment", "taxonomy"]
GENOMAD_AGG = ["chromosome_score", "plasmid_score", "virus_score"]
VS2_INFO = ["seqname", "max_score", "max_score_group", "length", "hallmark", "viral", "cellular"]


def read_fasta(path):
    name, chunks = None, []
    with open(path) as handle:
        for line in handle:
            line = line.rstrip()
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name, chunks = line[1:].split()[0], []
            elif line:
                chunks.append(line)
    if name is not None:
        yield name, "".join(chunks)


def read_rows(path, key):
    """key column -> list of row dicts (an empty file gives no rows)."""
    rows = {}
    if not path:
        return rows
    with open(path) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            rows.setdefault(row[key], []).append(row)
    return rows


def contig_of(seq_name):
    """Assembly contig a tool-reported sequence came from."""
    return seq_name.split("|")[0]


def by_contig(rows):
    grouped = {}
    for name, name_rows in rows.items():
        grouped.setdefault(contig_of(name), []).extend(name_rows)
    return grouped


def joined(rows, column):
    values = [row.get(column) or "NA" for row in rows]
    return ";".join(values) if values else "NA"


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--contigs", required=True, help="assembly contigs screened by the tools")
    parser.add_argument("--genomad-virus", help="geNomad <id>_virus.fna")
    parser.add_argument("--genomad-scores", help="geNomad <id>_aggregated_classification.tsv (every contig)")
    parser.add_argument("--genomad-summary", help="geNomad <id>_virus_summary.tsv")
    parser.add_argument("--vs2-virus", help="VirSorter2 final-viral-combined.fa")
    parser.add_argument("--vs2-scores", help="VirSorter2 final-viral-score.tsv")
    parser.add_argument("--vs2-groups", default="", help="comma-separated VS2 --include-groups (score columns)")
    parser.add_argument("--out-prefix", required=True)
    args = parser.parse_args()

    if not args.genomad_virus and not args.vs2_virus:
        parser.error("give --genomad-virus and/or --vs2-virus")

    genomad = dict(read_fasta(args.genomad_virus)) if args.genomad_virus else {}
    vs2 = {contig_of(name): seq for name, seq in read_fasta(args.vs2_virus)} if args.vs2_virus else {}

    genomad_summary = read_rows(args.genomad_summary, "seq_name")
    genomad_summary_by_contig = by_contig(genomad_summary)
    genomad_agg = read_rows(args.genomad_scores, "seq_name")
    vs2_rows = by_contig(read_rows(args.vs2_scores, "seqname"))
    vs2_groups = [g for g in args.vs2_groups.split(",") if g]

    genomad_contigs = {contig_of(name) for name in genomad}
    shared = genomad_contigs & set(vs2) if (genomad and vs2) else set()

    selected = []  # (original_name, contig_id, detected_by, sequence)
    if shared:
        for name, seq in read_fasta(args.contigs):
            if name in shared:
                selected.append((name, name, "both", seq))
    for name, seq in genomad.items():
        if contig_of(name) not in shared:
            selected.append((name, contig_of(name), "genomad", seq))
    for name, seq in vs2.items():
        if name not in shared:
            selected.append((name, name, "vs2", seq))

    vs2_columns = vs2_groups + VS2_INFO
    header = (["seq_id", "contig_id", "sample", "detected_by", "length", "genomad_score"]
              + [f"genomad_{c}" for c in GENOMAD_SUMMARY]
              + [f"genomad_agg_{c}" for c in GENOMAD_AGG]
              + [f"vs2_{c}" for c in vs2_columns])

    prefix = f"{args.sample}_"
    with open(f"{args.out_prefix}.fna", "w") as fasta, open(f"{args.out_prefix}.tsv", "w") as table:
        table.write("\t".join(header) + "\n")
        for name, contig_id, detected_by, seq in selected:
            seq_id = name.replace("|", "__")
            if not seq_id.startswith(prefix):
                seq_id = prefix + seq_id

            # geNomad summary rows describing this sequence: the exact
            # sequence for geNomad's own output, every provirus/virus row of
            # the contig when the whole contig is used.
            if detected_by == "genomad":
                summary_rows = genomad_summary.get(name, [])
            else:
                summary_rows = genomad_summary_by_contig.get(contig_id, [])
            agg_rows = genomad_agg.get(contig_id, [])

            if detected_by == "genomad" and summary_rows:
                score = summary_rows[0]["virus_score"]
            else:
                score = joined(agg_rows, "virus_score")

            row = [seq_id, contig_id, args.sample, detected_by, str(len(seq)), score]
            row += [joined(summary_rows, c) for c in GENOMAD_SUMMARY]
            row += [joined(agg_rows, c) for c in GENOMAD_AGG]
            row += [joined(vs2_rows.get(contig_id, []), c) for c in vs2_columns]
            fasta.write(f">{seq_id}\n{seq}\n")
            table.write("\t".join(row) + "\n")

    counts = {tool: sum(1 for s in selected if s[2] == tool) for tool in ("both", "genomad", "vs2")}
    print(f"{args.sample}: {len(selected)} candidates {counts}")


if __name__ == "__main__":
    main()
