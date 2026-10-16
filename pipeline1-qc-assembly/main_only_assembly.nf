#!/usr/bin/env nextflow
/*
 * pipeline1-qc-assembly: QC, trimming, host/PhiX removal, read-level
 * taxonomy, assembly, evaluation, mapping/abundance, gene prediction,
 * and a cross-sample gene catalog with functional annotation and abundance.
 *
 * Feeds into pipeline2-binning-mags, which is a separate pipeline by
 * design -- see README.md for how the two connect.
 *
 * Steps:
 *   1. FastQC (raw) + BBDuk              - adapter/quality trimming
 *   2. Bowtie2                           - host read removal
 *   3. Bowtie2                           - PhiX removal
 *   4. FastQC (clean) + MultiQC          - QC reporting
 *   5. Read taxonomy (--read_taxonomy)   - Kraken2+Bracken, SingleM (+ read fraction), MetaPhlAn
 *   6. MEGAHIT                           - assembly (sample-prefixed contig IDs)
 *   7. QUAST                             - assembly evaluation
 *   8. CoverM                            - mapping + contig abundance
 *   9. Prodigal                          - gene prediction (Prokka disabled)
 *  10. Gene catalog (all samples)        - complete genes, MMseqs2 95% clustering,
 *                                          eggNOG-mapper (parallel chunks), CoverM
 *                                          abundance, gene x sample matrices (counts,
 *                                          CPM, mean coverage, TPM), membership table
 *  11. featureCounts                     - per-sample gene abundance table
 */

nextflow.enable.dsl = 2

include { FASTQC_RAW; FASTQC_TRIMMED                } from './modules/fastqc.nf'
include { BBDUK_TRIM                                } from './modules/bbtools.nf'
include { BOWTIE2_BUILD as BOWTIE2_BUILD_HOST       } from './modules/bowtie2.nf'
include { BOWTIE2_BUILD as BOWTIE2_BUILD_PHIX       } from './modules/bowtie2.nf'
include { BOWTIE2_REMOVAL as BOWTIE2_HOST_REMOVAL   } from './modules/bowtie2.nf'
include { BOWTIE2_REMOVAL as BOWTIE2_PHIX_REMOVAL   } from './modules/bowtie2.nf'
include { HOST_PHIX_SUMMARY                         } from './modules/host_phix_summary.nf'
include { MULTIQC                                   } from './modules/multiqc.nf'
include { KRAKEN2; BRACKEN                          } from './modules/kraken2.nf'
include { SINGLEM_PIPE; SINGLEM_READ_FRACTION       } from './modules/singlem.nf'
include { METAPHLAN                                 } from './modules/metaphlan.nf'
include { MEGAHIT                                   } from './modules/megahit.nf'
include { QUAST                                     } from './modules/quast.nf'
include { COVERM_CONTIG; COVERM_GENE_CATALOG        } from './modules/coverm.nf'
include { SAMTOOLS_INDEX                            } from './modules/samtools_index.nf'
include { PRODIGAL                                  } from './modules/prodigal_prokka.nf'
// include { PRODIGAL; PROKKA                       } from './modules/prodigal_prokka.nf'   // Prokka disabled: annotation is done on the gene catalog (eggNOG-mapper)
include { MMSEQS2                                   } from './modules/mmseqs2.nf'
include { EGGNOG; EGGNOG_MERGE                      } from './modules/eggnog.nf'
include { GENE_CATALOG_MEMBERSHIP                   } from './modules/gene_catalog_membership.nf'
include { GENE_CATALOG_MATRIX                       } from './modules/gene_catalog_matrix.nf'
include { GENE_ABUNDANCE                            } from './modules/gene_abundance.nf'

workflow {

    // ---- 0. Validate parameters and samplesheet before anything runs ----
    if (!params.input)      error "Missing --input (samplesheet: sample,short_reads_1,short_reads_2)"
    if (!params.host_fasta) error "Missing --host_fasta (host genome fasta)"

    def valid_read_taxonomy_tools = ['kraken2', 'singlem', 'metaphlan']
    def read_taxonomy_tools = params.read_taxonomy.tokenize(',')*.trim()
    def unknown_tools = read_taxonomy_tools - valid_read_taxonomy_tools
    if (unknown_tools) {
        error "Unknown --read_taxonomy value(s): ${unknown_tools.join(', ')}. Valid: ${valid_read_taxonomy_tools.join(', ')}"
    }
    if ('kraken2' in read_taxonomy_tools && !params.bracken_sif) {
        error "--read_taxonomy includes kraken2, so --bracken_sif is required"
    }

    def sample_ids = file(params.input, checkIfExists: true).splitCsv(header: true).collect { row -> row.sample }
    def duplicated_ids = sample_ids.countBy { sample_id -> sample_id }.findAll { _sample_id, count -> count > 1 }.keySet()
    if (duplicated_ids) error "Duplicated sample IDs in ${params.input}: ${duplicated_ids.join(', ')}"

    file(params.eggnog_db, checkIfExists: true)

    ch_input = channel
        .fromPath(params.input)
        .splitCsv(header: true)
        .map { row ->
            def meta = [id: row.sample]
            [meta, file(row.short_reads_1, checkIfExists: true), file(row.short_reads_2, checkIfExists: true)]
        }

    // ---- 1. Raw QC + adapter/quality trimming ----
    FASTQC_RAW(ch_input)
    BBDUK_TRIM(ch_input, params.adapters_fasta)

    ch_kraken2_db = channel.value(file(params.kraken2_db, checkIfExists: true))
    KRAKEN2(BBDUK_TRIM.out.reads, ch_kraken2_db)

    // ---- 2-3. Host and PhiX removal (indexes built once, then reused) ----
    BOWTIE2_BUILD_HOST(channel.value(['host', file(params.host_fasta, checkIfExists: true)]))
    BOWTIE2_BUILD_PHIX(channel.value(['phix', file(params.phix_fasta, checkIfExists: true)]))

    BOWTIE2_HOST_REMOVAL(BBDUK_TRIM.out.reads, BOWTIE2_BUILD_HOST.out.index)

    BOWTIE2_PHIX_REMOVAL(BOWTIE2_HOST_REMOVAL.out.reads, BOWTIE2_BUILD_PHIX.out.index)
    ch_clean_reads = BOWTIE2_PHIX_REMOVAL.out.reads

    HOST_PHIX_SUMMARY(
        BOWTIE2_HOST_REMOVAL.out.log.collect(),
        BOWTIE2_PHIX_REMOVAL.out.log.collect()
    )

    // ---- 4. Post-cleaning QC + aggregate report ----
    FASTQC_TRIMMED(ch_clean_reads)
    ch_multiqc_input = FASTQC_RAW.out.zip
        .mix(
            FASTQC_TRIMMED.out.zip,
            BBDUK_TRIM.out.stats,
            BOWTIE2_HOST_REMOVAL.out.log,
            BOWTIE2_PHIX_REMOVAL.out.log
        )
        .collect()
    MULTIQC(ch_multiqc_input)

    // ---- 5. Read-level taxonomy (tools chosen with --read_taxonomy) ----
    if ('kraken2' in read_taxonomy_tools) {
        ch_kraken2_db = channel.value(file(params.kraken2_db, checkIfExists: true))
        KRAKEN2(ch_clean_reads, ch_kraken2_db)
        BRACKEN(KRAKEN2.out.report, ch_kraken2_db)
    }
    if ('singlem' in read_taxonomy_tools) {
        ch_singlem_metapackage = channel.value(file(params.singlem_metapackage, checkIfExists: true))
        SINGLEM_PIPE(ch_clean_reads, ch_singlem_metapackage)
        SINGLEM_READ_FRACTION(ch_clean_reads.join(SINGLEM_PIPE.out.profile), ch_singlem_metapackage)
    }
    if ('metaphlan' in read_taxonomy_tools) {
        METAPHLAN(ch_clean_reads, channel.value(file(params.metaphlan_db, checkIfExists: true)))
    }

    // ---- 6. Assembly ----
    MEGAHIT(ch_clean_reads)

    // ---- 7. Assembly evaluation ----
    QUAST(MEGAHIT.out.assembly)

    // ---- 8. Mapping + contig abundance (CoverM) ----
    ch_coverm_input = MEGAHIT.out.assembly.join(ch_clean_reads)
    COVERM_CONTIG(ch_coverm_input)
    SAMTOOLS_INDEX(COVERM_CONTIG.out.bam)

    // ---- 9. Gene prediction ----
    // PRODIGAL(MEGAHIT.out.assembly)
    // // PROKKA(MEGAHIT.out.assembly)   // disabled: see include above

    // // ---- 10. Cross-sample gene catalog ----
    // // Nucleotide genes, not proteins: the same catalog is used for read
    // // mapping (DNA) and for eggNOG-mapper (--itype CDS).
    // ch_all_genes = PRODIGAL.out.nt_fasta
    //     .map { _meta, fna -> fna }
    //     .collect()
    // MMSEQS2(ch_all_genes)

    // // Annotate the representative genes once, in parallel chunks.
    // ch_catalog_chunks = MMSEQS2.out.representatives
    //     .splitFasta(by: params.eggnog_chunk_size as Integer, file: true)
    // EGGNOG(ch_catalog_chunks)
    // EGGNOG_MERGE(
    //     EGGNOG.out.annotations.collect(),
    //     EGGNOG.out.orthologs.collect()
    // )

    // // Link every original gene to its catalog representative and contig.
    // ch_all_gff = PRODIGAL.out.gff
    //     .map { _meta, gff -> gff }
    //     .collect()
    // GENE_CATALOG_MEMBERSHIP(MMSEQS2.out.clusters, ch_all_gff)

    // // Map each sample's reads to the catalog, then merge into matrices.
    // COVERM_GENE_CATALOG(ch_clean_reads.combine(MMSEQS2.out.representatives))
    // GENE_CATALOG_MATRIX(
    //     MMSEQS2.out.representatives,
    //     COVERM_GENE_CATALOG.out.abundance.map { _meta, table -> table }.collect()
    // )

    // // ---- 11. Per-sample gene abundance table ----
    // ch_gene_input = PRODIGAL.out.gff.join(SAMTOOLS_INDEX.out.bam_bai)
    // GENE_ABUNDANCE(ch_gene_input)

    // ---- Emit a manifest for pipeline2 to consume ----
    // (sample, assembly, clean_reads_1, clean_reads_2) -- see README.
    // Points to the published copies in outdir, not to work/, so the
    // manifest stays valid after work/ is cleaned.
    def published_dir = file(params.outdir).toAbsolutePath()
    MEGAHIT.out.assembly
        .map { meta, _assembly ->
            [
                meta.id,
                "${published_dir}/06_megahit/${meta.id}/${meta.id}.contigs.fa",
                "${published_dir}/04a_clean_reads/${meta.id}_clean_1.fastq.gz",
                "${published_dir}/04a_clean_reads/${meta.id}_clean_2.fastq.gz"
            ].join(',')
        }
        .collectFile(
            name: 'pipeline2_samplesheet.csv',
            storeDir: params.outdir,
            seed: 'sample,assembly,short_reads_1,short_reads_2',
            newLine: true,
            sort: true
        )
}
