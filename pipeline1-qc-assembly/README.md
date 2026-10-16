# pipeline1-qc-assembly

QC, trimming, host/PhiX removal, read-level taxonomy, assembly,
evaluation, mapping/abundance, and gene prediction/annotation/abundance.

No skip toggles -- every step always runs, every reference path is
required. The pipeline fails fast at startup if any is missing, rather
than silently skipping a step you didn't mean to skip.

Feeds into `pipeline2-binning-mags` (steps 10-14) -- see that pipeline's
README for how they connect.

## Usage

```bash
nextflow run main.nf -profile singularity,slurm \
    --input assets/samplesheet.csv \
    --outdir results \
    --host_fasta /path/to/host_genome.fna \
    -resume
```

Choose read-level taxonomy tools with a comma-separated list (default
`metaphlan`); `kraken2` also runs Bracken:

```bash
    --read_taxonomy metaphlan,singlem,kraken2
```

Or (needs adjustments)

```bash
nextflow run main.nf -profile singularity,slurm \
    --input assets/samplesheet.csv \
    --outdir results \
    --host_fasta /path/to/host_genome.fa \
    --phix_fasta ../phix174_ill.ref.fa.gz \
    --adapters_fasta ../adapters.fa \
    --kraken2_db /path/to/kraken2_db \
    [...]
    -resume
```

## Parameters

| Param | Default | Notes |
|---|---|---|
| `input` | required | CSV: `sample,short_reads_1,short_reads_2`. Sample IDs must be unique. |
| `host_fasta` | required | Host genome fasta. |
| `phix_fasta` | `/fs/ess/PAS0822/references/phix174_ill.ref.fa.gz` | |
| `adapters_fasta` | `adapters` | BBDuk keyword for its bundled adapters, or a fasta path. |
| `bowtie2_index_dir` | `${outdir}/bowtie2_indexes` | Host/PhiX indexes are built once here and reused. |
| `read_taxonomy` | `metaphlan` | Any of `kraken2`, `singlem`, `metaphlan`. SingleM also reports the bacterial/archaeal read fraction. |
| `kraken2_db` |  `kraken2_pluspf` | Must contain Bracken `.kmer_distrib` files. |
| `kraken2_confidence` | `0` | Kraken2 `--confidence`. |
| `bracken_sif` | `/fs/ess/PAS0822/modules/singularity/bracken-3.1p1.sif` | Used with `kraken2`. |
| `bracken_read_length` | `150` | Must match a `database<N>mers.kmer_distrib` file. |
| `bracken_level` | `S` | |
| `singlem_metapackage` |  `S3.1.0.metapackage_20221209` | |
| `metaphlan_db` / `metaphlan_index` |  `mph4_2_db` / `mpa_vJan25_CHOCOPhlAnSGB_202503` | Index must exist locally; compute nodes have no internet. |
| `min_contig_size` | `1000` | MEGAHIT and QUAST. |


## Notes

- Contigs are renamed `<sample>_k141_N` so IDs stay unique when pooled
  across samples. Prokka rejects contig IDs longer than 37 characters, so
  keep sample IDs short.
- `00_summary/host_phix_summary.tsv` reports Bowtie2's overall alignment
  rate; the PhiX percentage is relative to host-removed reads.
- Kraken2 loads the whole database into RAM; its memory is set in
  `conf/base.config`.
- `publishDir` paths that use `meta` are wrapped in a closure
  (`publishDir { "...${meta.id}" }`), required by the strict parser that is
  the default from Nextflow 26.04; this also works on older versions.


## SingleM steps

1. `SINGLEM_PIPE` -- runs `singlem pipe` to produce an OTU table and a taxonomic profile.
2. `SINGLEM_READ_FRACTION` -- `singlem read_fraction` (the CoMS image is
   SingleM 1.0.0beta5; newer releases call this `prokaryotic_fraction`),
   giving the fraction of reads that are bacterial/archaeal.
   Works around a beta5 bug where `read_fraction` and `pipe` derive
   different sample names from the same read file.

## Gene catalog (`10_gene_catalog/`)


```
     Prodigal genes, all samples (.fna + .gff)
                      |
                      v
          MMseqs2 (95% nucleotide identity)
                      |
      +---------------+----------------------------+
      |                                            |
      v                                            v
Nonredundant gene catalog (.fna)          Cluster table (rep -> member)
      |                                            |  + Prodigal .gff
      +-------------------+                        v
      |                   |               Membership table
      v                   v               (gene -> rep, contig, sample)
   CoverM            eggNOG-mapper                 |
 (+ clean reads,     (--itype CDS)                 |
  per sample)             |                        |
      |                   v                        |
      v           Functional annotation            |
Gene x sample             |                        |
matrices                  |                        |
(counts, CPM, coverage)   |                        |
      |                   |                        |
      +-------------------+------------------------+
                          |
                          v
            Pipeline2 (+ MAGs, GTDB-Tk, CheckM2)
                          |
                          v
                MAG-linked master table
```

Follows Maranga et al. 2023 (mSystems, doi:10.1128/msystems.01178-22):

1. `MMSEQS2` -- complete Prodigal genes (`partial=00`) from all samples,
   clustered at 95% nucleotide identity (90% coverage of each member,
   greedy like CD-HIT `-c 0.95 -aS 0.9`) into a non-redundant catalog:
   `mmseqs2/gene_catalog_rep_seq.fasta` and `mmseqs2/gene_catalog_cluster.tsv`
   (representative, member). Partial genes are left out, as in the paper.
2. `EGGNOG` + `EGGNOG_MERGE` -- eggNOG-mapper 2.1.12 (DB 5.0.2) annotates the
   representatives in parallel chunks of `--eggnog_chunk_size` genes
   (default 50000), with the paper's options (`--itype CDS --translate
   --go_evidence non-electronic`, same e-value/score/coverage thresholds),
   merged into `eggnog/gene_catalog.emapper.annotations`.
3. `GENE_CATALOG_MEMBERSHIP` -- `membership/gene_catalog_membership.tsv`:
   one row per complete gene with its representative, contig, source
   sample and coordinates. Used by pipeline2 to link genes to MAGs.
4. `COVERM_GENE_CATALOG` + `GENE_CATALOG_MATRIX` -- each sample's clean
   reads mapped to the catalog (reads at >= 95% identity), merged into
   gene x sample matrices in `matrices/`: read counts, CPM (counts per
   million reads mapped to the catalog, not length-normalized), mean
   coverage, and TPM (length-weighted, comparable to the paper's CPM).

All catalog files use the same gene IDs (the representative's
`<contig>_<n>`), so annotation, membership and abundance join directly.
Scripts used by these steps are in `bin/` (must be executable).


## Output -> Pipeline 2 handoff

At the end of the run, this pipeline auto-generates
`${outdir}/pipeline2_samplesheet.csv` with columns
`sample,assembly,short_reads_1,short_reads_2` -- this is the exact
samplesheet format `pipeline2-binning-mags` expects. Just point pipeline 2's
`--input` at that file.
