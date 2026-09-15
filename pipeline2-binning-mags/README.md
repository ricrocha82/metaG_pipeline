# pipeline2-binning-mags

Steps 10–14: metagenome binning (4 binners), bin refinement, quality,
taxonomy, and MAG abundance.

Deliberately decoupled from `pipeline1-qc-assembly`: this pipeline only
needs an assembly + a pair of clean reads per sample, so it can be re-run
independently (e.g. re-refining bins with different DAS_Tool parameters)
without re-running QC/assembly. Can also be run as part of the repo-root
umbrella pipeline alongside pipeline1 and pipeline3 -- see the root
`README.md`.

## Usage

```bash
nextflow run main.nf -profile singularity,slurm \
    --binning_input /path/to/pipeline1/results/pipeline2_samplesheet.csv \
    --binning_outdir results \
    -resume
```

**Note:** `--input`/`--outdir`/`--min_contig_size` were renamed to
`--binning_input`/`--binning_outdir`/`--binning_min_contig_size` so this
pipeline can run alongside pipeline1 (which has its own
`--min_contig_size`) in the root umbrella pipeline without one silently
overriding the other.

`--gtdbtk_db` already defaults to the path you gave
(`/fs/project/PAS1117/modules/GTDB-Tk/gtdbtk_data/gtdbtk_r95_data.tar.gz`).
`--checkm2_db` is left unset by default, which triggers an automatic
one-time download (see caveats below) — pass your own path if you already
have a CheckM2 database somewhere.

## Before you run this — real compatibility concerns

1. **GTDB-Tk r95 database vs. v2.1.1 binary.** This is the biggest risk in
   this pipeline. GTDB release 95 is from ~2020 and was built for GTDB-Tk
   v1.x. `GTDB-Tk-2.1.1-PAS1117.sif` is a v2.x binary, which expects
   reference package release ~207 or later — this mismatch will likely
   cause `gtdbtk classify_wf` to fail its internal compatibility check, or
   at best silently produce unreliable classifications. Before relying on
   this step:
   ```bash
   find /fs/project/PAS1117 -iname "*gtdbtk*" 2>/dev/null
   ```
   If nothing newer exists, this is worth flagging to whoever manages your
   group's shared module library — GTDB-Tk reference data is large (~100GB+)
   but a stale r95 copy paired with a v2 binary isn't really usable as-is.

2. **CheckM2 database auto-download.** `CheckM2_database` (the DIAMOND
   reference CheckM2 needs) is *not* the same as the CheckM v1
   `CheckMdata` you have — v1's marker-gene HMM database and CheckM2's
   DIAMOND-based reference are structurally different. If `--checkm2_db`
   isn't supplied, the pipeline downloads CheckM2's database on first run
   via `checkm2 database --download`, cached under
   `${binning_outdir}/checkm2_db` so it only happens once. This needs outbound
   internet access from the compute node — confirm your SLURM
   partition/queue allows it (it did for your BUSCO/Kraken2 test runs
   earlier, so likely fine).

3. **COMEBin and SemiBin2 run via conda, not containers.** No `.sif` for
   either exists in `/fs/project/PAS1117/modules/singularity/`. Both
   processes declare a `conda` directive with no `container`, and the
   config is set up so this resolves correctly even under
   `-profile singularity` (see comments in `nextflow.config`). First run
   will be slower while conda builds these environments — consider a
   dry/small-data test run first.

4. **VAMB via container, but check the CLI syntax matches your version.**
   `VAMB-3.0.2.sif` is used — VAMB's CLI has changed somewhat across major
   versions (3.x vs. newer 4.x), so if you hit an argument-parsing error,
   check `singularity exec <sif> vamb --help` against the command in
   `modules/binning.nf`.

## Pipeline flow

```
per sample:  assembly + reads
                  │
                  ▼
          CoverM mapping
                  │
      ┌───────────┼───────────┬───────────┐
  MetaBAT2     COMEBin    SemiBin2      VAMB          (4 binners, parallel)
      └───────────┼───────────┴───────────┘
                  ▼
              DAS_Tool                                 (refinement, all 4 sets)
                  │
                  ▼
              CheckM2                                  (per-sample quality scores)

  ── all 4 samples' bins + CheckM2 scores pooled together ──

                  ▼
                dRep                                   (cross-sample dereplication,
                  │                                      99% ANI, best-quality rep kept)
      ┌───────────┼───────────┐
      ▼           ▼           ▼
  CheckM2      GTDB-Tk     CoverM genome matrix
  (final)     (taxonomy)   (MAG abundance, all samples)
```

**Why dRep matters here:** without it, a strain present in multiple samples
(quite plausible across your burn timepoints) would get assembled and
binned independently in each sample and counted as multiple separate MAGs.
dRep pools everything, uses CheckM2's completeness/contamination scores to
pick the best representative among near-duplicates, and gives you one
non-redundant MAG catalog. GTDB-Tk and CoverM then run once on that catalog
instead of once per sample.

## dRep specifics

- **Container**: `dRep-3.6.2.sif` (already in your PAS1117 module library).
- **Thresholds used**: `-pa 0.9 -sa 0.99 -nc 0.30` (primary clustering ANI
  90%, secondary/dereplication ANI 99%, minimum genome coverage overlap
  30%) — dRep's commonly-used defaults for species-level dereplication.
  Adjust in `modules/drep.nf` if you want a different stringency.
- **Genome quality input**: reuses the per-sample CheckM2 scores (step 4)
  rather than re-running CheckM2 before dereplication — avoids redundant
  compute. A **second, final CheckM2 pass** (`CHECKM2_FINAL`) then runs on
  the actual dereplicated output, since those are the quality numbers that
  matter for your final MAG set.

## MAG abundance matrix specifics

`COVERM_GENOME_MATRIX` maps *all 4 samples'* reads against the dereplicated
representative genomes in one CoverM call, producing a proper MAG × sample
abundance matrix — not 4 separate per-sample tables. This relies on all
samples' clean read files matching a consistent `*_clean_1/2.fastq.gz`
naming pattern (as pipeline1 produces) and sorting into matching pairs;
the process logs the pairing it detected before running CoverM, so check
that log looks right on your first run.


