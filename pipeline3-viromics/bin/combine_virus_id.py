#!/usr/bin/env python3
"""
Merge one VirSorter2 score table and one geNomad virus summary into a single TSV.

Example:
  combine_virus_id.py \
    --vs2-score /.../2019_02_viral_fraction.nextera_mspd_asm_final-viral-score.tsv \
    --genomad-summary /.../2019_02_viral_fraction.nextera_mspd_asm_virus_summary.tsv \
    --samples /.../samples.csv \
    --outdir .

Output:
  <prefix>_VS2_GENOMAD.txt  (prefix is --prefix OR inferred from --samples OR input filenames)
"""

from pathlib import Path
import argparse, sys
import pandas as pd
from typing import List, Optional

# ---------------- utilities ----------------

def read_samples(samples_path: Path) -> List[str]:
    """Read a samples file (CSV/TSV/whitespace). Returns list of sample IDs."""
    try:
        df = pd.read_csv(samples_path, sep=None, engine="python", dtype=str, comment="#")
    except Exception:
        df = pd.read_csv(samples_path, sep=r"\s+", engine="python", dtype=str, comment="#", header=None)
    if "sample_id" in df.columns:
        s = df["sample_id"]
    else:
        s = df.iloc[:, 0]
    s = s.astype(str).str.strip()
    return [x for x in s.tolist() if x]

def infer_prefix_from_samples(vs2_path: Path, gmd_path: Path, samples: List[str]) -> Optional[str]:
    """If samples list has one ID, use it. If multiple, pick one that matches either filename stem."""
    if not samples:
        return None
    if len(samples) == 1:
        return samples[0]
    base = f"{vs2_path.stem} {gmd_path.stem}"
    matches = [sid for sid in samples if sid in base]
    if len(matches) == 1:
        return matches[0]
    return None

# ---------------- loaders ----------------

def load_vs2(tsv: Path) -> pd.DataFrame:
    df = pd.read_csv(tsv, sep="\t", dtype=str)
    # normalize to 'seqname'
    if "seqname" not in df.columns:
        for cand in ("seq_name", "sequence"):
            if cand in df.columns:
                df = df.rename(columns={cand: "seqname"})
                break
    if "seqname" not in df.columns:
        raise ValueError(f"[VS2] Could not find seqname in {tsv}")

    # contig = part before '||...' (e.g., '...||full')
    df["contig"] = df["seqname"].str.replace(r"\|\|.*$", "", regex=True)

    keep = ["contig", "max_score", "max_score_group", "hallmark", "viral", "cellular", "length"]
    keep = [c for c in keep if c in df.columns]
    return df[keep].copy()

def load_genomad(tsv: Path) -> pd.DataFrame:
    df = pd.read_csv(tsv, sep="\t", dtype=str)
    # normalize contig column name
    for cand in ("seq_name", "seqname", "sequence", "contig"):
        if cand in df.columns:
            df = df.rename(columns={cand: "contig"})
            break
    if "contig" not in df.columns:
        raise ValueError(f"[geNomad] Could not find contig column in {tsv}")

    keep = ["contig", "length", "genetic_code", "virus_score", "n_hallmarks", "marker_enrichment", "taxonomy"]
    keep = [c for c in keep if c in df.columns]
    return df[keep].copy()

# ---------------- main ----------------

def main():
    ap = argparse.ArgumentParser(description="Merge VirSorter2 +/or geNomad (file inputs).")
    ap.add_argument("--vs2-score", type=Path, default=None, help="VS2 final-viral-score.tsv")
    ap.add_argument("--genomad-summary", type=Path, default=None, help="geNomad *_virus_summary.tsv")
    ap.add_argument("--samples", type=Path, default=None, help="Optional samples file to pick output prefix")
    ap.add_argument("--prefix", default=None, help="Explicit output basename (overrides --samples inference)")
    ap.add_argument("--outdir", type=Path, default=Path("."), help="Output directory (default: .)")
    args = ap.parse_args()

    if not args.vs2_score and not args.genomad_summary:
        print("ERROR: provide at least one of --vs2-score or --genomad-summary", file=sys.stderr)
        sys.exit(2)

    vs2_path = args.vs2_score.resolve() if args.vs2_score else None
    gmd_path = args.genomad_summary.resolve() if args.genomad_summary else None

    outdir = args.outdir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    # choose prefix
    prefix = args.prefix
    sample_ids: List[str] = []
    if not prefix and args.samples:
        try:
            sample_ids = read_samples(args.samples.resolve())
        except Exception as e:
            print(f"[WARN] Could not read samples file: {e}", file=sys.stderr)
    if not prefix:
        prefix = infer_prefix_from_samples(vs2_path, gmd_path, sample_ids)
    if not prefix:
        # fall back to whichever file is present
        if vs2_path:
            prefix = vs2_path.stem.replace("_final-viral-score", "")
        if (not prefix or prefix == (vs2_path.stem if vs2_path else None)) and gmd_path:
            alt = gmd_path.stem.replace("_virus_summary", "")
            prefix = alt or prefix
    if not prefix:
        prefix = "merged"

    # load inputs
    vs2_df = None
    gmd_df = None

    if vs2_path:
        try:
            vs2_df = load_vs2(vs2_path)
        except Exception as e:
            print(f"[ERROR] VS2 parse failed: {vs2_path} :: {e}", file=sys.stderr)
            sys.exit(1)

    if gmd_path:
        try:
            gmd_df = load_genomad(gmd_path)
        except Exception as e:
            print(f"[ERROR] geNomad parse failed: {gmd_path} :: {e}", file=sys.stderr)
            sys.exit(1)

    # merge logic for 3 cases
    if vs2_df is not None and gmd_df is not None:
        merged = vs2_df.merge(gmd_df, on="contig", how="outer")
        # unify length (prefer geNomad if both present)
        for c in ("length_x", "length_y"):
            if c in merged:
                merged[c] = pd.to_numeric(merged[c], errors="coerce")
        merged["length"] = merged.get("length_y").combine_first(merged.get("length_x")).astype("Int64")
        if {"length_x","length_y"} <= set(merged.columns):
            merged["length_agree"] = (
                (merged["length_x"] == merged["length_y"]) |
                merged["length_x"].isna() | merged["length_y"].isna()
            )
        merged = merged.drop(columns=[c for c in ("length_x","length_y") if c in merged])
    elif vs2_df is not None:
        merged = vs2_df.copy()
        # ensure Int64 dtype if length exists
        if "length" in merged.columns:
            merged["length"] = pd.to_numeric(merged["length"], errors="coerce").astype("Int64")
    else:  # gmd_df only
        merged = gmd_df.copy()
        if "length" in merged.columns:
            merged["length"] = pd.to_numeric(merged["length"], errors="coerce").astype("Int64")

    # reorder a bit
    first = ["contig", "length"]
    merged = merged[[*(c for c in first if c in merged.columns), *(c for c in merged.columns if c not in first)]]

    out_path = outdir / f"{prefix}_VS2_GENOMAD.txt"
    merged.to_csv(out_path, sep="\t", index=False)
    print(f"[OK] wrote {out_path} ({len(merged)} rows)")

if __name__ == "__main__":
    main()