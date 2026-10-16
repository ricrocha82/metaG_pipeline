#!/usr/bin/env nextflow
/*
 * pipeline2-binning-mags
 *
 * Per sample:  CoverM mapping -> MetaBAT2 / SemiBin2 / COMEBin -> Binette
 * Pooled:      all samples' Binette bins -> dRep -> CheckM2 (final)
 *              -> GTDB-Tk -> CoverM MAG x sample abundance matrix
 *              -> DRAM annotation (per MAG) -> merge -> distill
 *
 * Input is the pipeline2_samplesheet.csv written by pipeline1-qc-assembly
 * (sample,assembly,short_reads_1,short_reads_2).
 */

nextflow.enable.dsl = 2

include { COVERM_MAPPING; SAMTOOLS_INDEX             } from './modules/coverm_mapping.nf'
include { COUNT_CONTIGS; METABAT2; COMEBIN; SEMIBIN2 } from './modules/binning.nf'
include { BINETTE; MERGE_INPUT_BIN_QUALITY           } from './modules/binette.nf'
include { POOL_BINS; POOL_GENOMEINFO                 } from './modules/pool_bins.nf'
include { DREP                                       } from './modules/drep.nf'
include { CHECKM2_FINAL                              } from './modules/checkm2.nf'
include { GTDBTK                                     } from './modules/gtdbtk.nf'
include { COVERM_GENOME_MATRIX                       } from './modules/coverm_genome.nf'
include { DRAM_ANNOTATE; DRAM_MERGE; DRAM_DISTILL    } from './modules/dram.nf'

workflow {

    // ---- 0. Validate parameters and samplesheet before anything runs ----
    if (!params.input) error "Missing --input (pipeline1's pipeline2_samplesheet.csv: sample,assembly,short_reads_1,short_reads_2)"

    def sample_ids = file(params.input, checkIfExists: true).splitCsv(header: true).collect { row -> row.sample }
    def duplicated_ids = sample_ids.countBy { sample_id -> sample_id }.findAll { _sample_id, count -> count > 1 }.keySet()
    if (duplicated_ids) error "Duplicated sample IDs in ${params.input}: ${duplicated_ids.join(', ')}"

    file(params.comebin_simclr_patch, checkIfExists: true)
    file(params.dram_db, checkIfExists: true)
    ch_checkm2_db = channel.value(file(params.checkm2_db, checkIfExists: true))
    ch_gtdbtk_db  = channel.value(file(params.gtdbtk_db, checkIfExists: true))

    ch_input = channel
        .fromPath(params.input)
        .splitCsv(header: true)
        .map { row ->
            def meta = [id: row.sample]
            [
                meta,
                file(row.assembly, checkIfExists: true),
                file(row.short_reads_1, checkIfExists: true),
                file(row.short_reads_2, checkIfExists: true)
            ]
        }

    // ---- 1. Mapping (per sample) ----
    COVERM_MAPPING(ch_input)
    SAMTOOLS_INDEX(COVERM_MAPPING.out.bam)

    ch_assembly = ch_input.map { meta, assembly, _reads1, _reads2 -> [meta, assembly] }
    ch_binning_input = ch_assembly.join(SAMTOOLS_INDEX.out.bam_bai)

    // ---- 2. Binning (per sample, 3 binners in parallel) ----
    METABAT2(ch_binning_input)

    SEMIBIN2(ch_binning_input)

    // COMEBin's marker-gene seeding is unreliable below
    // params.comebin_min_contigs contigs, so small assemblies are not
    // submitted at all and are refined from MetaBAT2 + SemiBin2 only.
    COUNT_CONTIGS(ch_binning_input)

    ch_comebin_input = COUNT_CONTIGS.out.counted
        .filter { meta, _assembly, _bam, _bai, n_contigs_file ->
            def n_contigs = n_contigs_file.text.trim().toInteger()
            if (n_contigs < params.comebin_min_contigs) {
                log.warn "Skipping COMEBin for ${meta.id}: ${n_contigs} contigs (< ${params.comebin_min_contigs})"
            }
            n_contigs >= params.comebin_min_contigs
        }
        .map { meta, assembly, bam, bai, _n_contigs_file -> [meta, assembly, bam, bai] }

    COMEBIN(ch_comebin_input)

    // ---- 3. Bin refinement (per sample) ----
    // COMEBin is missing for skipped or failed samples (errorStrategy
    // 'ignore'), so it is joined with remainder and dropped when absent.
    ch_bin_sets = METABAT2.out.bins
        .join(SEMIBIN2.out.bins)
        .join(COMEBIN.out.bins, remainder: true)
        .map { meta, metabat2_bins, semibin2_bins, comebin_bins ->
            [meta, [metabat2_bins, semibin2_bins, comebin_bins].findAll { bins -> bins != null }]
        }

    BINETTE(ch_assembly.join(ch_bin_sets), ch_checkm2_db)
    MERGE_INPUT_BIN_QUALITY(BINETTE.out.input_quality)

    // ---- 4. Cross-sample dereplication ----
    POOL_BINS(BINETTE.out.bins.map { _meta, bins -> bins }.collect())
    POOL_GENOMEINFO(BINETTE.out.quality.map { _meta, report -> report }.collect())
    DREP(POOL_BINS.out.bins, POOL_GENOMEINFO.out.genome_info)

    // // ---- 5. Final quality report on the dereplicated MAG set ----
    CHECKM2_FINAL(DREP.out.representatives, ch_checkm2_db)

    // ---- 6. MAG taxonomy ----
    GTDBTK(
        DREP.out.representatives.map { mags -> [[id: 'dereplicated_mags'], mags] },
        ch_gtdbtk_db
    )

    // ---- 7. MAG abundance matrix (all samples vs. dereplicated set) ----
    ch_coupled_reads = ch_input
        .map { _meta, _assembly, reads1, reads2 -> "${reads1} ${reads2}" }
        .collect()
        .map { read_pairs -> read_pairs.join(' ') }

    COVERM_GENOME_MATRIX(DREP.out.representatives, ch_coupled_reads)

    // ---- 8. MAG functional annotation (DRAM, one job per MAG) ----
    ch_mags = DREP.out.representatives
        .flatMap { mags -> files("${mags}/*.fa") }
        .map { mag -> [[id: mag.baseName], mag] }

    DRAM_ANNOTATE(
        ch_mags,
        GTDBTK.out.summaries.collect(),
        CHECKM2_FINAL.out.results
    )
    DRAM_MERGE(DRAM_ANNOTATE.out.annotations.map { _meta, annotations -> annotations }.collect())
    DRAM_DISTILL(DRAM_MERGE.out.merged)
}
