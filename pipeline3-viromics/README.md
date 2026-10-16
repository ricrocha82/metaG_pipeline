# pipeline3-viromics

Viral contig detection, CheckV QC, vOTU clustering, vOTU abundance and
host prediction from per-sample assemblies + clean paired reads.
Based on the virion3 prototype and the CoMS dsDNA virus SOP
(`5.dsDNA_virus_ID-annotate_genomad_vs2.smk` → `8.dsDNA_virus_iphop.smk`).

```
assembly -> geNomad (+/- VS2) -> candidate pool (one file)
         -> CheckV -> filtering (length, score) -> vOTU clustering (CheckV ANI)
         -> CoverM (minimap2) -> vOTU abundance -> iPHoP host prediction
```

## Usage

```bash
nextflow run main.nf -profile singularity,slurm -ansi-log false \
    --input /path/to/pipeline1/results/pipeline2_samplesheet.csv \
    --outdir results \
    -resume

# geNomad + VirSorter2, with a custom iPHoP DB built from your MAGs
nextflow run main.nf -profile singularity,slurm -ansi-log false \
    --input samplesheet.csv \
    --detect both \
    --host_mags /path/to/pipeline2/results/03b_drep/drep_out/dereplicated_genomes \
    -resume
```

Input: `sample,assembly,short_reads_1,short_reads_2` (extra columns are
ignored), so pipeline1's `pipeline2_samplesheet.csv` works as is.

All containers and databases are on the CoMS shared paths; every database
is a param with a default in `nextflow.config`.

## Steps and outputs

| Folder | Step |
|---|---|
| `02_genomad/<sample>`, `02_virsorter2/<sample>` | Virus detection on the whole assembly (no length pre-filter) |
| `03_candidates` | All samples pooled: `candidates.fna` + `candidates.tsv` (which tool found each sequence + all geNomad and VS2 scores/annotations) |
| `04_checkv` | CheckV `end_to_end` on the pool |
| `05_viral_contigs` | `viral_contigs.*`: final viral sequences (CheckV-trimmed proviruses + viruses) with source, tool and CheckV quality. `removed_viral_contigs.*`: viral sequences that failed a rule, with a `reason` column (e.g. `length<5000`) |
| `06_votus` | `vOTUs.fna` (longest member of each cluster), `vOTU_clusters.tsv` |
| `07_coverm/<sample>` | CoverM table + BAM against the vOTUs |
| `08_abundance` | `vOTU_trimmed_mean.tsv`, `vOTU_read_count.tsv`, `vOTU_covered_fraction.tsv` |
| `09_iphop` | Merged `Host_prediction_to_genus_m90.csv` / `..._genome_m90.csv`; per-chunk outputs in `chunks/` |

geNomad/VS2 screen every contig, so short viral contigs are detected and
CheckV-assessed too; they end up in `removed_viral_contigs.*` with
`length<5000` instead of being discarded. Only sequences passing all rules
(`--min_contig_length`, default 5000) become vOTUs.

## Detection and filtering (`--detect`)

| Mode | Candidates | After CheckV |
|---|---|---|
| `genomad` (default) | geNomad viruses, score >= 0.7 (`--genomad_min_score`); proviruses already excised | length >= 5 kb |
| `vs2` | VS2 viruses, score >= 0.5, whole contigs | length >= 5 kb |
| `both` | Found by both: original contig. Found by one: that tool's sequence | length >= 5 kb; one-tool sequences also need contig length >= 10 kb, CheckV viral genes >= 1 and geNomad score >= 0.1 (CoMS SOP) |

`--detect vs2` alone has no CheckV gene filter (the SOP's `vs2_sop.sh`
curation is not included), so it is more permissive than the other modes.

### Columns of `candidates.tsv`, `viral_contigs.tsv`, `removed_viral_contigs.tsv`

| Columns | Content |
|---|---|
| `seq_id`, `contig_id`, `sample`, `detected_by` | Sequence, assembly contig it came from, sample, `genomad` / `vs2` / `both` |
| `length` / `candidate_length` | Final (CheckV-trimmed) length / length before CheckV |
| `genomad_score` | Score used by the filtering (provirus score for excised proviruses, else whole-contig score) |
| `genomad_*` | Every column of geNomad's `virus_summary.tsv` (coordinates, n_genes, hallmarks, fdr, taxonomy, ...) |
| `genomad_agg_*` | geNomad chromosome / plasmid / virus scores of the whole contig |
| `vs2_*` | Every column of VS2's `final-viral-score.tsv` (per-group scores, max_score, hallmark, viral, cellular) |
| `checkv_*` | Every column of CheckV's `quality_summary.tsv` |
| `reason` | Removed table only: rules the sequence failed |

`NA` = tool not run or no value. When a contig found by both tools has
several geNomad proviruses, their values are joined with `;`.

## vOTUs and abundance

* Clustering: blastn all-vs-all → `anicalc.py` → `aniclust.py` (CheckV),
  95% ANI over 85% of the shorter sequence (MIUViG). `--votu_min_ani`,
  `--votu_min_tcov`.
* Abundance: `coverm contig --mapper minimap2-sr`, reads >= 95% identity
  over >= 75% of their length; vOTUs with < 70% of their length covered are
  set to 0 (CoMS SOP). `--coverm_min_identity`, `--coverm_min_aligned`,
  `--coverm_min_covered_fraction`.

## Host prediction (iPHoP)

vOTUs are split into chunks of `--iphop_chunk_size` (5000) and predicted
as independent jobs against `--iphop_db` (CoMS `Jun_2025_pub_rw`).

With `--host_mags <dir>` (files `*.{host_mags_ext}`, default `fa`):
GTDB-Tk `de_novo_wf` (bacteria, then archaea; the archaeal run may fail
when there are no archaeal MAGs and is then skipped) → `iphop add_to_db`.
The new DB contains the default hosts plus your MAGs, so predictions use
only this DB. It is stored in `09_iphop/custom_db/iphop_custom_db` and not
rebuilt on later runs; reuse it elsewhere with `--iphop_db`.

## Before the first real run

* CheckV 0.8.1 + `checkv-db-v1.0` are the newest on CoMS. Both work
  together but are older than the current release (CheckV 1.0.x, DB v1.5);
  completeness estimates may differ from recent papers.
* geNomad 1.7.4 is used (as in the SOP). Check the DB version with
  `cat .../databases/genomad_db/version.txt`.
* VirSorter2 is run without `--db-dir`, assuming the CoMS image includes
  its database (as in the SOP).
* Resources for iPHoP and GTDB-Tk in `conf/base.config` are starting
  points; adjust from `pipeline_info/trace.txt`.

Tested with stand-in tools (output layouts of each tool) in `genomad`,
`vs2` and `both` modes, with and without `--host_mags`; not yet on real
data.
