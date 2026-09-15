# Nextflow pipeline to process metaG samples
Pipeline to process metagenomic samples (Constantly being updated)

## Structure

Three independent Nextflow pipelines, each still fully runnable on its own:

- [`pipeline1-qc-assembly`](pipeline1-qc-assembly/README.md) -- QC, trimming,
  host/PhiX removal, taxonomic classification, assembly, gene
  prediction/annotation.
- [`pipeline2-binning-mags`](pipeline2-binning-mags/README.md) -- binning (4
  binners), bin refinement, cross-sample dereplication, MAG taxonomy,
  MAG abundance. Consumes pipeline1's assemblies + clean reads.
- [`pipeline3-viromics`](pipeline3-viromics/README.md) -- virus detection,
  CheckV, clustering, abundance, microdiversity. Also consumes pipeline1's
  assemblies + clean reads, independently of pipeline2.

pipeline2 and pipeline3 both branch off pipeline1's output and don't depend
on each other.

## Running all three together

The root `main.nf` chains all three pipelines into one Nextflow run,
wiring pipeline1's assemblies directly into pipeline2 and pipeline3
in-memory (no intermediate CSV round-trip needed, though pipeline1 still
writes `pipeline2_samplesheet.csv` for inspection/resume purposes).

```bash
nextflow run main.nf -profile singularity,slurm \
    --qc_input assets/samplesheet.csv \
    --outdir results \
    --host_fasta /path/to/pig_genome.fa \
    --phix_fasta /path/to/phix174_ill.ref.fa.gz \
    --adapters_fasta /path/to/adapters.fa \
    --kraken2_db /path/to/kraken2_db \
    --bakta_db /path/to/bakta_db \
    --gtdbtk_db /path/to/gtdbtk_data.tar.gz \
    --genomad_db /path/to/genomad_db \
    -resume
```

This writes into `results/01_qc_assembly`, `results/02_binning_mags`, and
`results/03_viromics` respectively, plus one combined
`results/pipeline_info/{timeline,report}.html` and `trace.txt` for the
whole run.

Use `--skip_binning true` or `--skip_viromics true` to run only a subset of
stages (pipeline1 always runs -- it's the source of both downstream
stages' input).

Each pipeline also still runs completely standalone via its own `main.nf`
(see each subdirectory's README) -- useful for e.g. re-running pipeline2
with different DAS_Tool/dRep parameters without re-running QC/assembly.

## Why the params are pipeline-prefixed (`qc_*` / `binning_*` / `viromics_*`)

Nextflow has no per-module params namespace -- every `params.x` lives in one
flat global map. Before this umbrella existed, pipeline1 and pipeline2 both
defined `--outdir`/`--input`/`--min_contig_size` independently, which was
fine when each ran as its own isolated `nextflow run` invocation. Running
all three in one session would have made whichever config loaded last
silently win for every one of those shared names -- e.g. pipeline2's
`min_contig_size = 1500` overwriting pipeline1's `1000` and quietly
changing MEGAHIT's assembly behavior. Each pipeline's ambiguous params were
renamed to be pipeline-scoped specifically to make that impossible, while
keeping every other pipeline-specific param name (e.g. `--gtdbtk_db`,
`--checkm2_db`, `--genomad_db`, `--detect`) as-is since those were already
unique.

| Old name | New name (pipeline1) | New name (pipeline2) | New name (pipeline3) |
|---|---|---|---|
| `--input` | `--qc_input` | `--binning_input` | `--viromics_input` (was `--input_csv`) |
| `--outdir` | `--qc_outdir` | `--binning_outdir` | `--viromics_outdir` |
| `--min_contig_size` | `--qc_min_contig_size` | `--binning_min_contig_size` | n/a |

## Adding new tools/analyses later

Each pipeline is a set of modules (`modules/*.nf`) wired together in that
pipeline's `main.nf`. To add a new tool to, say, pipeline3: write a new
module file, add its container/conda definition to `conf/base.config` /
`nextflow.config`, `include` it in `pipeline3-viromics/main.nf`, and wire it
into the channel graph -- pipeline1 and pipeline2 are untouched, and the
umbrella `main.nf` doesn't need to change unless the new step's *output* is
meant to feed the other pipelines. Follow pipeline3's existing
loader/processor split (`LOAD_*` / `RUN_*` subworkflows) as the convention
for keeping new steps resumable and independently testable.
