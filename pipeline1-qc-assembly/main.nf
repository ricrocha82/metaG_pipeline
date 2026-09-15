#!/usr/bin/env nextflow
/*
 * pipeline1-qc-assembly: steps 1-9 of ricardo's pig_burn metagenomics
 * workflow -- QC, trimming, host/PhiX removal, taxonomic classification,
 * assembly, evaluation, mapping/abundance, gene prediction + annotation,
 * and gene abundance.
 *
 * No skip toggles -- every step always runs. host_fasta, phix_fasta,
 * adapters_fasta, and kraken2_db are required and the pipeline fails fast
 * if any are missing. SingleM runs with no --metapackage flag, relying on
 * whatever default the SingleM-0.13.2.sif container resolves on its own
 * (unverified -- see README).
 *
 * Feeds into pipeline2-binning-mags (steps 10-14), which is a separate
 * pipeline by design -- see README.md for how the two connect.
 *
 * Steps:
 *   1. FastQC (raw) + BBDuk          - adapter/quality trimming
 *   2. BBMap                          - host read removal
 *   3. BBDuk                          - PhiX removal
 *   4. FastQC (trimmed) + MultiQC     - QC reporting
 *   5. Kraken2 + SingleM              - read taxonomic classification
 *   6. MEGAHIT                        - assembly
 *   7. QUAST                          - assembly evaluation
 *   8. CoverM                         - mapping + contig abundance
 *   9. Prodigal + Prokka              - gene prediction + annotation
 *  10. featureCounts (subread)        - gene abundance table
 */

nextflow.enable.dsl = 2

include { FASTQC_RAW; FASTQC_TRIMMED               } from './modules/fastqc.nf'
include { BBDUK_TRIM; BBMAP_HOST_REMOVAL; BBDUK_PHIX_REMOVAL } from './modules/bbtools.nf'
// include { MULTIQC                                  } from './modules/multiqc.nf'
include { KRAKEN2                                  } from './modules/kraken2.nf'
include { SINGLEM_PIPE } from './modules/singlem.nf'
// include { SINGLEM_PIPE; SINGLEM_PROKARYOTIC_FRACTION } from './modules/singlem.nf'
include { MEGAHIT                                  } from './modules/megahit.nf'
include { QUAST                                    } from './modules/quast.nf'
include { COVERM_CONTIG                            } from './modules/coverm.nf'
include { SAMTOOLS_INDEX                           } from './modules/samtools_index.nf'
include { PRODIGAL; BAKTA                         } from './modules/prodigal_prokka.nf'
include { GENE_ABUNDANCE                           } from './modules/gene_abundance.nf'
include { HOST_PHIX_SUMMARY                        } from './modules/host_phix_summary.nf'

workflow QC_ASSEMBLY {
    take:
        ch_input   // [meta, short_reads_1, short_reads_2] -- pass null to read from params.qc_input instead

    main:

    // ---- 0. Input samplesheet: sample,short_reads_1,short_reads_2 ----
    if (ch_input == null) {
        ch_input = Channel
            .fromPath(params.qc_input, checkIfExists: true)
            .splitCsv(header: true)
            .map { row ->
                def meta = [id: row.sample]
                [meta, file(row.short_reads_1), file(row.short_reads_2)]
            }
    }

    // ---- 1. Raw QC + adapter/quality trimming ----
    FASTQC_RAW(ch_input)
    ch_adapters_fasta = Channel.fromPath(params.adapters_fasta, checkIfExists: true)
    BBDUK_TRIM(ch_input, ch_adapters_fasta.first())

    // ---- 2. Host read removal ----
    ch_host_fasta = Channel.fromPath(params.host_fasta, checkIfExists: true)
    BBMAP_HOST_REMOVAL(BBDUK_TRIM.out.reads, ch_host_fasta.first())

    // ---- 3. PhiX removal ----
    ch_phix_fasta = Channel.fromPath(params.phix_fasta, checkIfExists: true)
    BBDUK_PHIX_REMOVAL(BBMAP_HOST_REMOVAL.out.reads, ch_phix_fasta.first())
    ch_clean_reads = BBDUK_PHIX_REMOVAL.out.reads

    // ---- Host/PhiX contamination % summary across all samples ----
    HOST_PHIX_SUMMARY(
        BBMAP_HOST_REMOVAL.out.stats.collect(),
        BBDUK_PHIX_REMOVAL.out.stats.collect()
    )

    // ---- 4. Post-trim QC + aggregate report ----
    FASTQC_TRIMMED(ch_clean_reads)
    ch_multiqc_input = FASTQC_RAW.out.zip
        .mix(FASTQC_TRIMMED.out.zip)
        .mix(BBDUK_TRIM.out.stats)
        .mix(BBMAP_HOST_REMOVAL.out.stats)
        .mix(BBDUK_PHIX_REMOVAL.out.stats)
        .collect()
    // MULTIQC(ch_multiqc_input)

    // ---- 5. Taxonomic classification of reads (Kraken2 + SingleM) ----
    ch_kraken2_db = Channel.fromPath(params.kraken2_db, checkIfExists: true)
    KRAKEN2(ch_clean_reads, ch_kraken2_db.first())

    SINGLEM_PIPE(ch_clean_reads)
    // SINGLEM_PROKARYOTIC_FRACTION(SINGLEM_PIPE.out.profile)

    // ---- 6. Assembly ----
    MEGAHIT(ch_clean_reads)

    // ---- 7. Assembly evaluation ----
    QUAST(MEGAHIT.out.assembly)

    // ---- 8. Mapping + contig abundance (CoverM) ----
    ch_coverm_input = MEGAHIT.out.assembly.join(ch_clean_reads)
    COVERM_CONTIG(ch_coverm_input)
    SAMTOOLS_INDEX(COVERM_CONTIG.out.bam)

    // ---- 9. Gene prediction + annotation ----
    PRODIGAL(MEGAHIT.out.assembly)

    ch_bakta_db = Channel
        .fromPath(params.bakta_db, checkIfExists: true)
    BAKTA(MEGAHIT.out.assembly, ch_bakta_db.first()
        )

    // ---- 10. Gene abundance table ----
    ch_gene_input = PRODIGAL.out.gff.join(SAMTOOLS_INDEX.out.bam_bai)
    GENE_ABUNDANCE(ch_gene_input)

    // ---- Emit a manifest for pipeline2/pipeline3 to consume ----
    // (sample, assembly, clean_reads_1, clean_reads_2) -- see README
    ch_assemblies = MEGAHIT.out.assembly.join(ch_clean_reads)   // [meta, asm, r1, r2]

    ch_manifest = ch_assemblies
        .map { meta, asm, r1, r2 ->
            "${meta.id},${asm.toAbsolutePath()},${r1.toAbsolutePath()},${r2.toAbsolutePath()}"
        }
        .collectFile(
            name: 'pipeline2_samplesheet.csv',
            storeDir: params.qc_outdir,
            seed: 'sample,assembly,short_reads_1,short_reads_2\n',
            newLine: true
        )

    emit:
        assemblies = ch_assemblies   // [meta, assembly, short_reads_1, short_reads_2] for in-memory umbrella wiring
        manifest   = ch_manifest     // path to pipeline2_samplesheet.csv, for standalone pipeline2/3 runs
}

// ---- Standalone entry point: `nextflow run pipeline1-qc-assembly/main.nf` ----
workflow {
    QC_ASSEMBLY(null)
}
