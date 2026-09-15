#!/usr/bin/env nextflow
/*
 * pipeline2-binning-mags: steps 10-14, now with cross-sample dereplication.
 *
 * Per-sample (unchanged): mapping -> binning (4 binners) -> DAS_Tool -> CheckM2
 * Cross-sample (new):     pool all samples' bins + CheckM2 scores -> dRep
 *                         -> GTDB-Tk (once) -> CoverM genome matrix (once)
 *
 * Steps:
 *   1. CoverM                 - mapping (per-sample)
 *   2. MetaBAT2/COMEBin/SemiBin2/VAMB - binning (per-sample, 4 binners) # removed for now
 *   3. DAS_Tool                - bin refinement (per-sample)
 *   4. CheckM2                  - bin quality (per-sample, feeds dRep)
 *   5. dRep                      - cross-sample dereplication
 *   6. CheckM2 (final)            - quality report on dereplicated MAG set
 *   7. GTDB-Tk                     - MAG taxonomy (once, on dereplicated set)
 *   8. CoverM genome matrix          - MAG abundance across all samples (once)
 */

nextflow.enable.dsl = 2

include { COVERM_MAPPING; SAMTOOLS_SORT_INDEX      } from './modules/coverm_mapping.nf'
include { COUNT_CONTIGS; METABAT2; COMEBIN; SEMIBIN2 } from './modules/binning.nf'
// include { METABAT2; COMEBIN; SEMIBIN2; VAMB          } from './modules/binning.nf'
include { FASTA_TO_CONTIG2BIN as FASTA_TO_CONTIG2BIN_METABAT2 } from './modules/dastool.nf'
include { FASTA_TO_CONTIG2BIN as FASTA_TO_CONTIG2BIN_COMEBIN   } from './modules/dastool.nf'
include { FASTA_TO_CONTIG2BIN as FASTA_TO_CONTIG2BIN_SEMIBIN2  } from './modules/dastool.nf'
// include { DASTOOL                                                } from './modules/dastool.nf'
include { METAWRAP_REFINE                                       } from './modules/metawrap.nf'
include { CHECKM2_DOWNLOAD_DB; CHECKM2_FINAL; CHECKM2   } from './modules/checkm2.nf'
include { POOL_BINS; POOL_CHECKM2_GENOMEINFO             } from './modules/pool_bins.nf'
include { DREP                                            } from './modules/drep.nf'
include { GTDBTK                        } from './modules/gtdbtk.nf'
// include { GTDBTK_PREPARE_DB; GTDBTK                        } from './modules/gtdbtk.nf'
include { COVERM_GENOME_MATRIX                              } from './modules/coverm_genome.nf'
include { CHECKM2_SING as CHECKM2_SING_METABAT2 } from './modules/checkm2.nf'
include { CHECKM2_SING as CHECKM2_SING_SEMIBIN2 } from './modules/checkm2.nf'
include { CHECKM2_SING as CHECKM2_SING_COMEBIN  } from './modules/checkm2.nf'

workflow BINNING_MAGS {
    take:
        ch_input   // [meta, assembly, short_reads_1, short_reads_2] -- pass null to read from params.binning_input instead

    main:

    // ---- 0. Input samplesheet: sample,assembly,short_reads_1,short_reads_2 ----
    if (ch_input == null) {
        ch_input = Channel
            .fromPath(params.binning_input, checkIfExists: true)
            .splitCsv(header: true)
            .map { row ->
                def meta = [id: row.sample]
                [meta, file(row.assembly), file(row.short_reads_1), file(row.short_reads_2)]
            }
    }

    // ---- 1. Mapping (per-sample) ----
    COVERM_MAPPING(ch_input)
    SAMTOOLS_SORT_INDEX(COVERM_MAPPING.out.bam)

    ch_assembly = ch_input.map { meta, asm, r1, r2 -> [meta, asm] }
    ch_binning_input = ch_assembly.join(SAMTOOLS_SORT_INDEX.out.bam_bai)

    // ---- 2. Binning (per-sample, 4 binners, parallel) ----
    METABAT2(ch_binning_input)
    // COMEBIN(ch_binning_input)
    SEMIBIN2(ch_binning_input)
    // VAMB(ch_binning_input)

    // COMEBin's marker-gene seeding step is unreliable below
    // params.comebin_min_contigs contigs (can yield zero detected seeds,
    // producing zero bins regardless of other tuning). Check contig count
    // up front and only submit qualifying samples to COMEBin at all --
    // avoids burning a SLURM job on a predictably-failing run.
    COUNT_CONTIGS(ch_binning_input)

    ch_comebin_branched = COUNT_CONTIGS.out.counted
        .map { meta, assembly, bam, bai, n_contigs_file ->
            def n_contigs = n_contigs_file.text.trim().toInteger()
            [meta, assembly, bam, bai, n_contigs]
        }
        .branch { meta, assembly, bam, bai, n_contigs ->
            run_comebin:  n_contigs >= params.comebin_min_contigs
            skip_comebin: true
        }

    COMEBIN(
        ch_comebin_branched.run_comebin.map { meta, assembly, bam, bai, n -> [meta, assembly, bam, bai] }
    )

    ch_comebin_skipped = ch_comebin_branched.skip_comebin.map { meta, assembly, bam, bai, n ->
        def empty_dir = file("${workDir}/comebin_skipped_${meta.id}")
        empty_dir.mkdirs()
        log.warn "Skipping COMEBin for ${meta.id}: only ${n} contigs (< ${params.comebin_min_contigs}) -- known marker-gene seeding limitation, refining from MetaBAT2 + SemiBin2 only"
        [meta, empty_dir]
    }

    ch_comebin_final = COMEBIN.out.bins.mix(ch_comebin_skipped)

    // ---- 2.1. Bin quality (per-sample) ----
    if (!params.checkm2_db) {
        CHECKM2_DOWNLOAD_DB()
        ch_checkm2_db = CHECKM2_DOWNLOAD_DB.out.db
    } else {
        ch_checkm2_db = Channel.fromPath(params.checkm2_db, checkIfExists: true)
    }
    // CHECKM2(DASTOOL.out.bins, ch_checkm2_db.first())

    ch_metabat2_raw = METABAT2.out.bins.map     { meta, dir -> [meta + [binner: 'metabat2'], dir] }
    ch_semibin2_raw = SEMIBIN2.out.bins.map     { meta, dir -> [meta + [binner: 'semibin2'], dir] }
    ch_comebin_raw  = ch_comebin_final.map      { meta, dir -> [meta + [binner: 'comebin'],  dir] }

    CHECKM2_SING_METABAT2(ch_metabat2_raw, ch_checkm2_db.first())
    CHECKM2_SING_SEMIBIN2(ch_semibin2_raw, ch_checkm2_db.first())
    CHECKM2_SING_COMEBIN(ch_comebin_raw,  ch_checkm2_db.first())

    // ---- 3. Bin refinement (per-sample, DAS_Tool consolidates all 4 sets) ----
    // FASTA_TO_CONTIG2BIN_METABAT2(METABAT2.out.bins.map { meta, dir -> [meta, dir, 'metabat2'] })
    // FASTA_TO_CONTIG2BIN_COMEBIN(COMEBIN.out.bins.map { meta, dir -> [meta, dir, 'comebin'] })
    // FASTA_TO_CONTIG2BIN_SEMIBIN2(SEMIBIN2.out.bins.map { meta, dir -> [meta, dir, 'semibin2'] })

    ch_refine_input = METABAT2.out.bins
    .join(ch_comebin_final, remainder: true)
    // .join(COMEBIN.out.bins, remainder: true)
    .join(SEMIBIN2.out.bins, remainder: true)
    .map { meta, metabat2_dir, comebin_dir, semibin2_dir ->
        def empty_dir = file("${workDir}/empty_bins_${meta.id}")
        empty_dir.mkdirs()
        [
            meta,
            metabat2_dir ?: empty_dir,
            comebin_dir  ?: empty_dir,
            semibin2_dir ?: empty_dir
        ]
    }

    // ---- 4. Bin refiment + quality (per-sample, feeds dRep's genomeInfo) ----
    METAWRAP_REFINE(ch_refine_input)
    // DASTOOL(ch_dastool_input)
    CHECKM2(METAWRAP_REFINE.out.bins, ch_checkm2_db.first())
    // CHECKM2(DASTOOL.out.bins, ch_checkm2_db.first())

    // ---- 5. Cross-sample dereplication ----
    ch_all_bin_dirs = METAWRAP_REFINE.out.bins.map { meta, dir -> dir }.collect()
    // ch_all_bin_dirs = DASTOOL.out.bins.map { meta, dir -> dir }.collect()
    ch_all_checkm2  = CHECKM2.out.results.map { meta, tsv -> tsv }.collect()

    POOL_BINS(ch_all_bin_dirs)
    POOL_CHECKM2_GENOMEINFO(ch_all_checkm2)
    DREP(POOL_BINS.out.bins, POOL_CHECKM2_GENOMEINFO.out.genome_info)

    // ---- 6. Final quality report on the dereplicated MAG set ----
    CHECKM2_FINAL(DREP.out.representatives, ch_checkm2_db.first())

    // ---- 7. MAG taxonomy (once, on dereplicated set) ----
    // ch_gtdbtk_archive = Channel.fromPath(params.gtdbtk_db, checkIfExists: true)
    // ch_gtdbtk_db_dir = Channel.fromPath(params.gtdbtk_db, checkIfExists: true, type: 'dir')
    GTDBTK(
        DREP.out.representatives.map { dir -> [[id: 'dereplicated_mags'], dir] },
        params.gtdbtk_db
    )
    // GTDBTK_PREPARE_DB(ch_gtdbtk_archive)
    // GTDBTK(
    //     DREP.out.representatives.map { dir -> [[id: 'dereplicated_mags'], dir] },
    //     GTDBTK_PREPARE_DB.out.db_dir.first()
    // )

    // ---- 8. MAG abundance matrix (once, all samples vs. dereplicated set) ----
    ch_coupled_reads = ch_input
        .map { meta, asm, r1, r2 -> "${r1} ${r2}" }
        .collect()
        .map { it.join(' ') }

    COVERM_GENOME_MATRIX(DREP.out.representatives, ch_coupled_reads)

    emit:
        representatives = DREP.out.representatives   // dereplicated MAG dir, for downstream use
        gtdbtk          = GTDBTK.out.results         // taxonomy assignments
}

// ---- Standalone entry point: `nextflow run pipeline2-binning-mags/main.nf` ----
workflow {
    BINNING_MAGS(null)
}
