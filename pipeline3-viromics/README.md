# virion3
codes to analyse data for the Virion3 project

# sequences data
We have two platforms: ONT and PacBio. Along with that we have short-reads from illumina.

# Pipeline
assembly -> virus detection (genomad +/or vs2) -> filtering -> contig pool (one file) ->  checkV -> filtering -> cluster (checkV) -> map (minimap2 + coverM) -> vOTU abundance

Each assembly type were independently searched for virus contigs using `VS2` and `geNomad`.

to run go to the directory and run
```bash
conda activate nf-core

# go to your directory
cd my/path/to/nf/results

# this will run only genomad
nextflow run main.nf --outdir results -profile singularity

# to run both genomad + vs2
nextflow run main.nf --outdir results -profile singularity -detect both
```

Resume entry workflows
If you want to stop at a specific stage:
```bash
nextflow run main.nf -entry VIRUS_DETECTION
nextflow run main.nf -entry CHECKV
nextflow run main.nf -entry CLUSTERING
nextflow run main.nf -entry MAPPING
nextflow run main.nf -entry ABUNDANCE
nextflow run main.nf -entry ENHANCED
nextflow run main.nf -entry MICRODIVERSITY   # same as FULL_PIPELINE
```

If you already have some files and want to run some steps only
```bash
# Have CheckV FASTAs, want to re-run everything after
nextflow run main.nf -entry FROM_CHECKV \
    --checkv_dir results/virus/checkv \
    --input_csv   /results/samples.csv \
    --outdir      results \
    -profile singularity

# Have cluster FASTAs, want mapping + abundance + enhanced
nextflow run main.nf -entry FROM_CLUSTERING \
    --cluster_dir results/virus/cluster \
    --input_csv   /results/samples.csv \
    --outdir      results \
    -profile singularity

# Have BAMs and cluster FASTAs, want abundance + enhanced
nextflow run main.nf -entry FROM_MAPPING \
    --bam_dir     results/virus/mapping \
    --cluster_dir results/virus/cluster \
    --input_csv   /results/samples.csv \
    --outdir      results \
    -profile singularity

# Have BAMs + clusters, want enhanced virome + microdiversity
nextflow run main.nf -entry FROM_ABUNDANCE \
    --bam_dir     results/virus/mapping \
    --cluster_dir results/virus/cluster \
    --input_csv   /results/samples.csv \
    --outdir      results \
    -profile singularity

# Have everything up to enhanced, want microdiversity only
nextflow run main.nf -entry FROM_ENHANCED \
    --bam_dir      results/virus/mapping \
    --cluster_dir  results/virus/cluster \
    --enhanced_ref results/virus/cluster/enhanced.self-blastn.clusters.fna \
    --enhanced_bam results/virus/mapping/enhanced_sorted.bam \
    --input_csv   /results/samples.csv \
    --outdir      results \
    -profile singularity
```

Filter steps:

 1. VS2: length >= 5000, maxscore >= 0.5 (if used)
 2. geNomad: length >= 5000, virus_score >=0.7
 3. checkV: contig_length >= 5000 (filter only by contig length -> more ecologically meaningful)
 4. checkV: contig_length >= 5000 and completeness >=50 (not ideal for the analysis)
 5. checkV: contig_length >= 10000, gene counts >= 13, viral_genes >=1, host_genes <2, completeness >=70 (removed too many contigs)

 After that, virus contigs were pulled and dereplicated into viral populations using `CheckV` clustering method with 70% coverage and 95% nucleotide clustering.

The next step is mapping the clustered sequences to produce an abundance table using `minimap2` and `CoverM`.

Trimmed mean coverage was calculated. Trimmed mean coverage values were then divided by total Gbp per quality filtered paired reads. 

The last process is to calculate microdiversity of each assembly.


# Samples
# WEC samples
path to WEC samples: `/fs/ess/PAS1117/virion3/03_WEC_GS`

We followed two approaches:
- full data
- subsampling (based on basedpair)

in both approaches `Mylo`, `Hfam`, `Canu`, `Flye` and `MDBG` were tested in each long-read sequences.

## 1. full assemblies

### Long-reads

#### raw reads
ONT: `/fs/project/PAS1117/virion3/03_WEC_GS/reads_ont/filt/WECR10_EXP-PBC001_barcode10_qlf.fq`

PacBio: `/fs/project/PAS1117/virion3/03_WEC_GS/reads_pcb/filt/m84039_2509_hifi_reads.bc2004_qlf.fq`

1. Full pacbio assemlbies:

    Mylo: `/fs/ess/PAS1117/virion3/03_WEC_GS/assemblies/mylo/m84039_2509_hifi_reads.bc2004/m84039_2509_hifi_reads.bc2004_mylo_asm.fa`

    Hfam: `/fs/ess/PAS1117/virion3/03_WEC_GS/assemblies/hfam/m84039_2509_hifi_reads.bc2004/m84039_2509_hifi_reads.bc2004_hfam_asm.fa`

    Hybrid: `/fs/ess/PAS1117/virion3/03_WEC_GS/assemblies/hspd/ill_WEC_2019_02_viral_fraction_pcb_m84039_2509_hifi_reads.bc2004_hspd/WEC_2019_02_viral_fraction_pcb_m84039_2509_hifi_reads.bc2004_hspd_asm.fa`

2. full nanopore:

    Mylo: `/fs/ess/PAS1117/virion3/03_WEC_GS/assemblies/mylo/WECR10_EXP-PBC001_barcode10/WECR10_EXP-PBC001_barcode10_mylo_asm.fa` (considered the best)

    canu: `assemblies/canu/WECR10_EXP-PBC001_barcode10/WECR10_EXP-PBC001_barcode10_canu_asm.fa` (it was removed)

    flye: `assemblies/flye/WECR10_EXP-PBC001_barcode10/WECR10_EXP-PBC001_barcode10_flye_asm.fa`

    mdbg: `assemblies/mdbg/WECR10_EXP-PBC001_barcode10/WECR10_EXP-PBC001_barcode10_mdbg_asm.fa`

    Hybrid: Hybrid (LR + SR): `/fs/ess/PAS1117/virion3/03_WEC_GS/assemblies/hspd/ill_WEC_2019_02_viral_fraction_ont_WECR10_EXP-PBC001_barcode10_hspd/WEC_2019_02_viral_fraction_ont_WECR10_EXP-PBC001_barcode10_hspd_asm.fa`

## #Short reads
#### Assemblies
Only Short-reads: `/fs/ess/PAS1117/virion3/03_WEC_GS/assemblies/mspd/WEC_2019_02_viral_fraction_mspd/WEC_2019_02_viral_fraction_mspd_asm.fa`

#### raw reads
R1: `/fs/project/PAS1117/virion3/03_WEC_GS/reads_ill/trim/WEC_2019_02_viral_fraction_R1_PE.fastq`

R2: `/fs/project/PAS1117/virion3/03_WEC_GS/reads_ill/trim/WEC_2019_02_viral_fraction_R2_PE.fastq`



### Subsampling 

Try to make the reads "comparable" in terms of based-pairs. The PacBio sequences is much more deep compared to ONT. After analysis, the subsampled sequences were stable in terms of quality and recovered viral contigs. So, we decided to use the sequences with the pattern `sample1Gbp-04`

#### raw reads
ONT: `/fs/ess/PAS1117/virion3/03_WEC_GS/reads_ont/filt/WECR10_EXP-PBC001_barcode10_sample1Gbp-04_qlf.fq,`

PacBio: `/fs/ess/PAS1117/virion3/03_WEC_GS/reads_pcb/filt/m84039_2509_hifi_reads.bc2004_sample1Gbp-04_qlf.fq`

#### long reads

1. subsampled pacbio assemlbies:

    Mylo: `/fs/ess/PAS1117/virion3/03_WEC_GS/assemblies/mylo/m84039_2509_hifi_reads.bc2004_sample1Gbp-04/m84039_2509_hifi_reads.bc2004_sample1Gbp-04_mylo_asm.fa`

    Hfam: `/fs/ess/PAS1117/virion3/03_WEC_GS/assemblies/hfam/m84039_2509_hifi_reads.bc2004_sample1Gbp-04/m84039_2509_hifi_reads.bc2004_sample1Gbp-04_hfam_asm.fa`

    flye: `/fs/ess/PAS1117/virion3/03_WEC_GS/assemblies/flye/m84039_2509_hifi_reads.bc2004_sample1Gbp-04/m84039_2509_hifi_reads.bc2004_sample1Gbp-04_flye_asm.fa`

    mdbg: `/fs/ess/PAS1117/virion3/03_WEC_GS/assemblies/mdbg/m84039_2509_hifi_reads.bc2004_sample1Gbp-04/m84039_2509_hifi_reads.bc2004_sample1Gbp-04_mdbg_asm.fa`

2. Subsampled nanopore:

    Mylo: `/fs/ess/PAS1117/virion3/03_WEC_GS/assemblies/mylo/WECR10_EXP-PBC001_barcode10_sample1Gbp-04/WECR10_EXP-PBC001_barcode10_sample1Gbp-04_mylo_asm.fa` (considered the best)

    flye: `/fs/ess/PAS1117/virion3/03_WEC_GS/assemblies/flye/WECR10_EXP-PBC001_barcode10_sample1Gbp-04/WECR10_EXP-PBC001_barcode10_sample1Gbp-04_flye_asm.fa`

    mdbg: `/fs/ess/PAS1117/virion3/03_WEC_GS/assemblies/mdbg/WECR10_EXP-PBC001_barcode10_sample1Gbp-04/WECR10_EXP-PBC001_barcode10_sample1Gbp-04_mdbg_asm.fa`

## Short reads
### Assemblies
Hybrid (LR + SR): `/fs/project/PAS1117/ricardo/virion3/data/fasta/hybrid/hspd/ill_WEC_2019_02_viral_fraction_sample1Gbp-04_ont_WECR10_EXP-PBC001_barcode10_sample1Gbp-04_hspd_asm.fa`

Only Short-reads: `/fs/project/PAS1117/ricardo/virion3/data/fasta/short-reads/mspd/WEC_2019_02_viral_fraction_sample1Gbp-04_mspd_asm.fa`
### raw reads
R1: `/fs/ess/PAS1117/virion3/03_WEC_GS/reads_ill/trim/WEC_2019_02_viral_fraction_sample1Gbp-04_R1_PE.fastq`

R2: `/fs/ess/PAS1117/virion3/03_WEC_GS/reads_ill/trim/WEC_2019_02_viral_fraction_sample1Gbp-04_R2_PE.fastq`


## input csv table 

read_path: `assembly file`

raw_read[1,2]: `sequence reads`

| sample_id | seq_name | read_path | raw_read1 | raw_read2 |
|---|---|---|---|---|
| `WECR10_EXP-PBC001_barcode10_canu_asm` | `WECR10_EXP-PBC001_barcode10_canu_asm.fa` | `/fs/ess/PAS1117/virion3/03_WEC_GS/assemblies/canu/WECR10_EXP-PBC001_barcode10/WECR10_EXP-PBC001_barcode10_canu_asm.fa` | `/fs/ess/PAS1117/virion3/03_WEC_GS/fastqs_filt/WECR10_EXP-PBC001_barcode10_qlf.fq` |  |
| `WECR10_EXP-PBC001_barcode10_flye_asm` | `WECR10_EXP-PBC001_barcode10_flye_asm.fa` | `/fs/ess/PAS1117/virion3/03_WEC_GS/assemblies/flye/WECR10_EXP-PBC001_barcode10/WECR10_EXP-PBC001_barcode10_flye_asm.fa` | `/fs/ess/PAS1117/virion3/03_WEC_GS/fastqs_filt/WECR10_EXP-PBC001_barcode10_qlf.fq` |  |
| `WECR10_EXP-PBC001_barcode10_mdbg_asm` | `WECR10_EXP-PBC001_barcode10_mdbg_asm.fa` | `/fs/ess/PAS1117/virion3/03_WEC_GS/assemblies/mdbg/WECR10_EXP-PBC001_barcode10/WECR10_EXP-PBC001_barcode10_mdbg_asm.fa` | `/fs/ess/PAS1117/virion3/03_WEC_GS/fastqs_filt/WECR10_EXP-PBC001_barcode10_qlf.fq` |  |
| `2019_02_viral_fraction.nextera__WECR10_EXP-PBC001_barcode10__hspd_asm` | `2019_02_viral_fraction.nextera__WECR10_EXP-PBC001_barcode10__hspd_asm.fa` | `/fs/ess/PAS1117/virion3/03_WEC_GS/assemblies/hspd/2019_02_viral_fraction.nextera__WECR10_EXP-PBC001_barcode10/2019_02_viral_fraction.nextera__WECR10_EXP-PBC001_barcode10__hspd_asm.fa` | `/fs/ess/PAS1117/virion3/00_data/short_reads/ben_WEC/trimm/2019_02_viral_fraction.nextera_pe1.fq.gz` | `/fs/ess/PAS1117/virion3/00_data/short_reads/ben_WEC/trimm/2019_02_viral_fraction.nextera_pe2.fq.gz` |
| `2019_02_viral_fraction.nextera_mspd_asm` | `2019_02_viral_fraction.nextera_mspd_asm.fa` | `/fs/ess/PAS1117/virion3/03_WEC_GS/assemblies/mspd/2019_02_viral_fraction.nextera/2019_02_viral_fraction.nextera_mspd_asm.fa` | `/fs/ess/PAS1117/virion3/00_data/short_reads/ben_WEC/trimm/2019_02_viral_fraction.nextera_pe1.fq.gz` | `/fs/ess/PAS1117/virion3/00_data/short_reads/ben_WEC/trimm/2019_02_viral_fraction.nextera_pe2.fq.gz` |



## Steps of the analysis

1. subsampling raw reads based on bp
2. assembly
3. run 2 pipelines for each platform with all the assemblers together
3. select the best assembler for each platform

    Based on the results we've chosen `Hifiasm` for PacBio and `Mylo` for ONT.
    
    inputs: `/fs/project/PAS1117/ricardo/virion3/data/samples_ont_all.csv` and `/fs/project/PAS1117/ricardo/virion3/data/samples_pcb_all.csv`

4. run the 2 pipelines with selected assemblers and get the `enhanced` virome

    inputs: `/fs/project/PAS1117/ricardo/virion3/data/samples_ont.csv` and `/fs/project/PAS1117/ricardo/virion3/data/samples_pcb.csv`



# Comparing shared vTOUs among full assemblies

### 1. ONT-Mylo vs PB-Mylo

- PB-Mylo: `m84039_2509_hifi_reads.bc2004_mylo_asm_final_viral_contigs.fasta`
- ONT-Mylo: `WECR10_EXP-PBC001_barcode10_mylo_asm_final_viral_contigs.fasta`
- output dir: `/fs/project/PAS1117/ricardo/virion3/results/analysis/cluster_mylo_full`

### 2. ONT-Mylo vs - PB-HFAM

- ONT-Mylo: `WECR10_EXP-PBC001_barcode10_mylo_asm_final_viral_contigs.fasta`
- PB-HFAM: `m84039_2509_hifi_reads.bc2004_hfam_asm_final_viral_contigs.fasta`
- output dir: `/fs/project/PAS1117/ricardo/virion3/results/analysis/cluster_ont_pb_full`

### 3. PB-HFAM vs PB-Mylo

- PB-HFAM: `m84039_2509_hifi_reads.bc2004_hfam_asm_final_viral_contigs.fasta`
- PB-Mylo: `m84039_2509_hifi_reads.bc2004_mylo_asm_final_viral_contigs.fast`
- output dir: `/fs/project/PAS1117/ricardo/virion3/results/analysis/cluster_pb_full`



# NOTES

## MetaPop SNP Linking Bug (`UnboundLocalError` / `IndexError`)

**Affected version:** MetaPop 1.0.0 (Python 3.7)  
**File:** `metapop_mine_reads.py`, function `formatted_to_valid_combos` (~line 370)

---

### The Error

When running the microdiversity module, MetaPop may crash at the SNP linking step with one of the following errors:

```
UnboundLocalError: local variable 'new' referenced before assignment
```
or
```
IndexError: list index out of range
```

---

### Root Cause

This happens in the SNP linking step (do_mine_reads → formatted_to_valid_combos). The variable new is only assigned inside a conditional block, but the code tries to append it unconditionally — a classic Python scoping bug in MetaPop's source.

This is a known MetaPop bug that's triggered when a contig has SNPs on the same codon but the co-occurrence logic hits an edge case (e.g., reads that partially span the codon, or unusual codon positions).

The function `formatted_to_valid_combos` builds a `new_keys` list by iterating over SNP codon combinations. The original code contained two bugs:

1. **`UnboundLocalError`** — Inside the `if 1 < len(key) < 3` branch, the variable `new` was only assigned when `codon_positions` matched one of three expected patterns (`[1,2]`, `[2,3]`, `[1,3]`). If none matched, `new` was never defined but was still appended, causing the crash. The three `if` statements should have been `if/elif/elif/else`.

2. **`IndexError`** — The outer loop had no `else` branch for keys where `len(key)` was not between 1 and 3. These keys were silently skipped, making `new_keys` shorter than `key_results`, causing an index mismatch at line 380 when building `matchable_dict`.

---

### Fix

The patched version of `metapop_mine_reads.py` is available at `conf/metapop/metapop_mine_reads.py`. The corrected logic is:

```python
for key in key_results:
    if 1 < len(key) < 3:
        if codon_positions == [1, 2]:
            new = key + "*"
        elif codon_positions == [2, 3]:
            new = "*" + key
        elif codon_positions == [1, 3]:
            new = key[0] + "*" + key[1]
        else:
            new = key          # fallback for unexpected codon position patterns
        new_keys.append(new)
    else:
        new_keys.append(key)   # keys with len != 2 must still be appended
```

The correct py file is located in
```bash
./bin/metapop_mine_reads.py
```
