#!/usr/bin/env python3
"""
Split the CheckV-processed candidates into final viral sequences
(host-trimmed proviruses + viruses) and removed ones, with the reason.

Rules:
  all modes      : final sequence length >= --min-length
  --mode both    : sequences found by only one tool must also have
                   CheckV contig_length >= --exclusive-min-length,
                   viral_genes >= --exclusive-min-viral-genes and
                   geNomad virus score >= --exclusive-min-genomad-score

Writes <prefix>.fna/.tsv and removed_<prefix>.fna/.tsv (extra 'reason' column).
CheckV names excised proviruses '<contig_id>_<n>'.
"""
import argparse
import csv
import re


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


def read_tsv(path, key):
    with open(path) as handle:
        return {row[key]: row for row in csv.DictReader(handle, delimiter="\t")}


def exclusive_rule_failures(qc, candidate, args):
    """Reasons a sequence found by only one tool fails the --mode both rule."""
    reasons = []
    if int(qc["contig_length"]) < args.exclusive_min_length:
        reasons.append(f"one_tool:contig_length<{args.exclusive_min_length}")
    if int(qc["viral_genes"]) < args.exclusive_min_viral_genes:
        reasons.append(f"one_tool:viral_genes<{args.exclusive_min_viral_genes}")
    score = candidate["genomad_score"]
    if score == "NA" or float(score) < args.exclusive_min_genomad_score:
        reasons.append(f"one_tool:genomad_score<{args.exclusive_min_genomad_score}")
    return reasons


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--quality-summary", required=True)
    parser.add_argument("--viruses", required=True)
    parser.add_argument("--proviruses", required=True)
    parser.add_argument("--candidates", required=True)
    parser.add_argument("--mode", required=True, choices=["genomad", "vs2", "both"])
    parser.add_argument("--min-length", type=int, required=True)
    parser.add_argument("--exclusive-min-length", type=int, required=True)
    parser.add_argument("--exclusive-min-viral-genes", type=int, required=True)
    parser.add_argument("--exclusive-min-genomad-score", type=float, required=True)
    parser.add_argument("--out-prefix", required=True)
    args = parser.parse_args()

    checkv = read_tsv(args.quality_summary, "contig_id")
    candidates = read_tsv(args.candidates, "seq_id")

    sequences = [(name, name, seq) for name, seq in read_fasta(args.viruses)]
    sequences += [(name, re.sub(r"_\d+$", "", name), seq) for name, seq in read_fasta(args.proviruses)]

    # Every candidate column (tool scores/annotations) + every CheckV column.
    # 'length' is the final (CheckV-trimmed) sequence; the candidate's own
    # length becomes 'candidate_length'.
    with open(args.candidates) as handle:
        candidate_columns = [c for c in next(csv.reader(handle, delimiter="\t")) if c != "seq_id"]
    with open(args.quality_summary) as handle:
        checkv_columns = [c for c in next(csv.reader(handle, delimiter="\t")) if c != "contig_id"]
    header = "\t".join(["seq_id", "candidate_id", "length"]
                       + [("candidate_length" if c == "length" else c) for c in candidate_columns]
                       + [(c if c.startswith("checkv_") else f"checkv_{c}") for c in checkv_columns])
    kept = removed = 0
    with open(f"{args.out_prefix}.fna", "w") as fasta, open(f"{args.out_prefix}.tsv", "w") as table, \
            open(f"removed_{args.out_prefix}.fna", "w") as removed_fasta, \
            open(f"removed_{args.out_prefix}.tsv", "w") as removed_table:
        table.write(header + "\n")
        removed_table.write(header + "\treason\n")
        for seq_id, candidate_id, seq in sequences:
            candidate, qc = candidates[candidate_id], checkv[candidate_id]
            reasons = []
            # Length of the final (CheckV-trimmed) sequence
            if len(seq) < args.min_length:
                reasons.append(f"length<{args.min_length}")
            if args.mode == "both" and candidate["detected_by"] != "both":
                reasons += exclusive_rule_failures(qc, candidate, args)

            row = ([seq_id, candidate_id, str(len(seq))]
                   + [candidate[c] for c in candidate_columns]
                   + [qc[c] for c in checkv_columns])
            if reasons:
                removed_fasta.write(f">{seq_id}\n{seq}\n")
                removed_table.write("\t".join(row) + "\t" + ";".join(reasons) + "\n")
                removed += 1
            else:
                fasta.write(f">{seq_id}\n{seq}\n")
                table.write("\t".join(row) + "\n")
                kept += 1

    print(f"{len(checkv)} candidates checked: {kept} viral sequences kept, {removed} removed")


if __name__ == "__main__":
    main()
