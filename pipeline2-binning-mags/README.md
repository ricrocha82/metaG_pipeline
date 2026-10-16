# pipeline2-binning-mags

Metagenome binning, bin refinement, dereplication, MAG quality, taxonomy,
MAG abundance, and MAG functional annotation. 

Deliberately decoupled from `pipeline1-qc-assembly`: it only needs an
assembly + a pair of clean reads per sample, so binning can be re-run
without re-running QC/assembly.

## Usage

```bash
nextflow run main.nf -profile singularity,slurm \
    --input /path/to/pipeline1/results/pipeline2_samplesheet.csv \
    --outdir results \
    -ansi-log false \
    -resume
```

`--input` is the `pipeline2_samplesheet.csv` that pipeline1 writes
(`sample,assembly,short_reads_1,short_reads_2`, pointing to its published
files). Sample IDs must be unique.

## Pipeline flow

```
per sample:  assembly + clean reads
                  │
             CoverM mapping ─ samtools index
                  │
      ┌───────────┼───────────┐
  MetaBAT2     SemiBin2     COMEBin          (COMEBin skipped below comebin_min_contigs)
      └───────────┼───────────┘
               Binette                        (refinement, CheckM2 scoring)

  ── all samples' Binette bins + quality pooled ──
                  │
                dRep                          (quality filter + cross-sample dereplication)
      ┌───────────┼───────────┐
  CheckM2      GTDB-Tk     CoverM genome
  (final)      (r226)      (MAG x sample abundance)
      └─────┬─────┘
          DRAM                                (per MAG, with taxonomy + quality)
            │
     merge ─ distill
```

## Parameters

| Param | Default | Notes |
|---|---|---|
| `input` | required | pipeline1's `pipeline2_samplesheet.csv`. |
| `checkm2_db` | `CheckM2_database/uniref100.KO.1.dmnd` | Used by Binette and the final CheckM2. |
| `gtdbtk_db` | `GTDBTK_data/release226` | Extracted reference package for GTDB-Tk 2.6.1. |
| `comebin_simclr_patch` | `comebin/patched_simclr.py` | Bound over COMEBin 1.0.3's `simclr.py`; guards an empty-batch crash in the training-accuracy metric only. |
| `dram_db` | `DRAM_dbs/1.4.6` | Bound at `/databases/DRAM` in the DRAM 1.4.6 (no KEGG) image. |
| `binning_seed` | `1` | MetaBAT2 and SemiBin2. |
| `min_contig_size` | `1500` | MetaBAT2 `-m` (1500 is its minimum). |
| `semibin2_environment` | `global` | SemiBin2 pretrained model, e.g. `human_oral`, `human_gut`. |
| `comebin_min_contigs` | `1024` | Smaller assemblies skip COMEBin (its marker-gene seeding is unreliable there). |
| `drep_min_completeness` | `70` | value. MIMAG medium quality would be 50. |
| `drep_max_contamination` | `10` | Binette already caps contamination at 10. |
| `drep_secondary_ani` | `0.98` | value. Use 0.95 for species-level representatives. |
| `ignore_genome_quality` | `false` | Test runs only; disables dRep's quality filter. |
| `coverm_min_covered_fraction` | `0.70` | A MAG counts as present only if this fraction of it is covered. |

## Outputs

| Folder | Content |
|---|---|
| `00_coverm_mapping/` | Per-sample contig depth. |
| `01_binning/<sample>/` | MetaBAT2, SemiBin2 and COMEBin bins. |
| `02_binette/<sample>/` | Refined bins (`<sample>_binN.fa`), their quality, and CheckM2 quality of every input bin set. |
| `03_drep/` | `genomeInfo.csv` and dRep output; representatives in `drep_out/dereplicated_genomes`. |
| `04_checkm2_final/` | CheckM2 report for the dereplicated MAGs. |
| `05_gtdbtk/` | GTDB-Tk classification. |
| `06_mag_abundance/` | MAG x sample abundance matrix. |
| `07_dram/` | `dram_merged/` (annotations, genes, tRNAs/rRNAs for all MAGs) and `dram_distill/` (metabolism summary, product heatmap). |

## Notes

- A binner can finish with zero bins (small or low-coverage assemblies);
  Binette is only given the bin sets that contain bins. A sample with no
  bins at all yields an empty Binette result instead of stopping the run.
- COMEBin has `errorStrategy 'ignore'` and exits 0 on zero bins, so a
  failed sample is cached and refined from MetaBAT2 + SemiBin2 on `-resume`.
- DRAM runs one job per dereplicated MAG, after GTDB-Tk and the final
  CheckM2, so taxonomy and quality are included in its annotations.
  DRAM's own default drops contigs shorter than 2500 bp before calling
  genes. The image has no KEGG database, so KEGG-based distillate modules
  rely on the other databases (KOfam, Pfam, dbCAN, MEROPS, ...).
- Contig IDs are `<sample>_k141_N` (pipeline1), so contigs stay unique
  when bins are pooled across samples.
- `publishDir` paths that use `meta` are closures, as required by the
  strict parser (default from Nextflow 26.04).
