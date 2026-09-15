# pipeline1-qc-assembly

Steps 1-9 of the pig_burn metagenomics workflow: QC, trimming, host/PhiX
removal, taxonomic classification (Kraken2 + SingleM), assembly,
evaluation, mapping/abundance, and gene prediction/annotation/abundance.

No skip toggles -- every step always runs, every reference path is
required. The pipeline fails fast at startup if any is missing, rather
than silently skipping a step you didn't mean to skip.

Feeds into `pipeline2-binning-mags` (steps 10-14) -- see that pipeline's
README for how they connect. Can also be run as part of the repo-root
umbrella pipeline alongside pipeline2 and pipeline3 -- see the root
`README.md`.

## Usage

```bash
nextflow run main.nf -profile singularity,slurm \
    --qc_input assets/samplesheet.csv \
    --qc_outdir results \
    --host_fasta /path/to/pig_genome.fa \
    --phix_fasta /fs/project/PAS1117/bioinformatic_tools/bbmap_38.51/resources/phix174_ill.ref.fa.gz \
    --adapters_fasta /fs/project/PAS1117/bioinformatic_tools/bbmap_38.51/resources/adapters.fa \
    --kraken2_db /path/to/kraken2_db \
    -resume
```

**Note:** `--input`/`--outdir`/`--min_contig_size` were renamed to
`--qc_input`/`--qc_outdir`/`--qc_min_contig_size` so this pipeline can run
alongside pipeline2 (which has its own `--min_contig_size`) in the root
umbrella pipeline without one silently overriding the other.

## Required inputs -- still open

1. **`--host_fasta`**: pig reference genome fasta. Not yet provided.
2. **`--kraken2_db`**: Kraken2 database directory. Not yet provided.
3. **`--singlem_db`**: removed. SingleM now runs with no `--metapackage`
   flag at all, relying on whatever default `SingleM-0.13.2.sif` resolves
   on its own. **Unverified** whether this container actually has a
   working default configured — if `SINGLEM_PIPE` fails with something
   like "no metapackage found" on your first real run, that means it
   doesn't, and a `--metapackage` path will need adding back into
   `modules/singlem.nf`. Check before running for real:
   ```bash
   singularity exec /fs/project/PAS1117/modules/singularity/SingleM-0.13.2.sif \
       bash -c 'echo $SINGLEM_METAPACKAGE_PATH'
   singularity exec /fs/project/PAS1117/modules/singularity/SingleM-0.13.2.sif \
       find / -iname "*.smpkg*" 2>/dev/null
   ```
4. **`--phix_fasta`** and **`--adapters_fasta`**: confirmed available at
   `/fs/project/PAS1117/bioinformatic_tools/bbmap_38.51/resources/`
   (`phix174_ill.ref.fa.gz` and `adapters.fa` respectively).

## SingleM steps

Matches the official two-step workflow
(https://wwood.github.io/singlem/tools/prokaryotic_fraction) exactly:
1. `SINGLEM_PIPE` -- runs `singlem pipe` to produce a taxonomic profile.
2. `SINGLEM_PROKARYOTIC_FRACTION` -- runs `singlem prokaryotic_fraction`
   on that profile (this subcommand was previously named
   `microbial_fraction` in older SingleM docs/versions).

Both steps run with no `--metapackage` argument, relying on
`SingleM-0.13.2.sif`'s own default (unverified whether it has one --
see item 3 above).

## Output -> Pipeline 2 handoff

At the end of the run, this pipeline auto-generates
`${qc_outdir}/pipeline2_samplesheet.csv` with columns
`sample,assembly,short_reads_1,short_reads_2` -- this is the exact
samplesheet format `pipeline2-binning-mags` expects. Just point pipeline 2's
`--binning_input` at that file (or, in the umbrella pipeline, this handoff
happens in-memory automatically -- the CSV is still written for
inspection/resume purposes).
