#!/usr/bin/env nextflow
/*
 * pipeline3-viromics: viral contig detection, CheckV QC, vOTU clustering,
 * vOTU abundance and host prediction.
 *
 * Input is the same samplesheet as pipeline2 (pipeline1 writes it as
 * results/pipeline2_samplesheet.csv): sample,assembly,short_reads_1,short_reads_2
 *
 * Steps:
 *   1-2. geNomad / VS2   - virus detection on the whole assembly, --detect genomad|vs2|both (per sample)
 *   3. Candidate pool    - one FASTA + table of which tool found each sequence
 *   4. CheckV            - quality + host-region trimming (once, on the pool)
 *   5. Filtering         - length / score rules; removed sequences saved with the reason (once)
 *   6. vOTU clustering   - 95% ANI / 85% AF, CheckV method (once)
 *   7. CoverM            - minimap2 mapping + vOTU coverage (per sample)
 *   8. Abundance tables  - vOTU x sample matrices (once)
 *   9. iPHoP             - host prediction, optional custom DB from MAGs
 */

nextflow.enable.dsl = 2

include { SPLIT_VOTUS                              } from './modules/seqkit.nf'
include { GENOMAD                                  } from './modules/genomad.nf'
include { VIRSORTER2                               } from './modules/virsorter2.nf'
include { SELECT_CANDIDATES                        } from './modules/candidates.nf'
include { CHECKV; FILTER_VIRAL_CONTIGS             } from './modules/checkv.nf'
include { VOTU_CLUSTERING                          } from './modules/votu_clustering.nf'
include { COVERM_VOTU; VOTU_ABUNDANCE_TABLES       } from './modules/coverm.nf'
include { GTDBTK_DE_NOVO; IPHOP_ADD_TO_DB; IPHOP_PREDICT } from './modules/iphop.nf'

workflow {

    // ---- 0. Validate parameters and samplesheet before anything runs ----
    if (!params.input) error "Missing --input (samplesheet: sample,assembly,short_reads_1,short_reads_2)"
    if (!(params.detect in ['genomad', 'vs2', 'both'])) {
        error "Invalid --detect '${params.detect}'. Use genomad, vs2 or both"
    }

    def rows = file(params.input, checkIfExists: true).splitCsv(header: true)
    def missing_columns = ['sample', 'assembly', 'short_reads_1', 'short_reads_2'] - rows[0].keySet()
    if (missing_columns) error "Samplesheet ${params.input} is missing column(s): ${missing_columns.join(', ')}"
    def duplicated_ids = rows.countBy { row -> row.sample }.findAll { _sample_id, count -> count > 1 }.keySet()
    if (duplicated_ids) error "Duplicated sample IDs in ${params.input}: ${duplicated_ids.join(', ')}"

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

    def run_genomad = params.detect in ['genomad', 'both']
    def run_vs2     = params.detect in ['vs2', 'both']

    // ---- 1-2. Virus detection on the whole assembly (per sample) ----
    // No length pre-filter: the tools see every contig, and length/score
    // rules are applied afterwards (step 5), so short viral contigs are
    // still reported in removed_viral_contigs.*
    ch_assembly = ch_input.map { meta, assembly, _r1, _r2 -> [meta, assembly] }

    // Missing tool outputs are passed as [] so SELECT_CANDIDATES has one input shape.
    if (run_genomad) {
        GENOMAD(ch_assembly, file(params.genomad_db, checkIfExists: true))
        ch_genomad = GENOMAD.out.virus.join(GENOMAD.out.scores).join(GENOMAD.out.summary)
    } else {
        ch_genomad = ch_assembly.map { meta, _assembly -> [meta, [], [], []] }
    }

    if (run_vs2) {
        VIRSORTER2(ch_assembly)
        ch_vs2 = VIRSORTER2.out.virus.join(VIRSORTER2.out.scores)
    } else {
        ch_vs2 = ch_assembly.map { meta, _assembly -> [meta, [], []] }
    }

    // ---- 3. Candidate pool (one FASTA + one table for all samples) ----
    SELECT_CANDIDATES(ch_assembly.join(ch_genomad).join(ch_vs2))

    ch_candidates_fasta = SELECT_CANDIDATES.out.fasta
        .map { _meta, fasta -> fasta }
        .collectFile(name: 'candidates.fna', storeDir: "${params.outdir}/03_candidates", sort: true)
    ch_candidates_table = SELECT_CANDIDATES.out.table
        .map { _meta, table -> table }
        .collectFile(name: 'candidates.tsv', storeDir: "${params.outdir}/03_candidates", keepHeader: true, skip: 1, sort: true)

    // ---- 4-5. CheckV on the pool, then filtering ----
    CHECKV(ch_candidates_fasta, file(params.checkv_db, checkIfExists: true))
    FILTER_VIRAL_CONTIGS(
        CHECKV.out.quality_summary,
        CHECKV.out.viruses,
        CHECKV.out.proviruses,
        ch_candidates_table
    )

    // ---- 6. vOTU clustering ----
    VOTU_CLUSTERING(FILTER_VIRAL_CONTIGS.out.fasta)

    // ---- 7-8. Mapping + vOTU abundance ----
    ch_coverm_input = ch_input
        .map { meta, _assembly, r1, r2 -> [meta, r1, r2] }
        .combine(VOTU_CLUSTERING.out.votus)
    COVERM_VOTU(ch_coverm_input)
    VOTU_ABUNDANCE_TABLES(COVERM_VOTU.out.abundance.map { _meta, tsv -> tsv }.collect())

    // ---- 9. Host prediction ----
    // With --host_mags, the default DB plus your MAGs becomes a new full
    // iPHoP DB, and predictions use that DB (it already contains every
    // default host, so a second run against the default DB adds nothing).
    if (params.host_mags) {
        ch_host_mags = channel
            .fromPath("${params.host_mags}/*.${params.host_mags_ext}", checkIfExists: true)
            .collect()
        GTDBTK_DE_NOVO(ch_host_mags)
        IPHOP_ADD_TO_DB(ch_host_mags, GTDBTK_DE_NOVO.out.gtdbtk, file(params.iphop_db, checkIfExists: true))
        ch_iphop_db = IPHOP_ADD_TO_DB.out.db
    } else {
        ch_iphop_db = channel.value(file(params.iphop_db, checkIfExists: true))
    }

    SPLIT_VOTUS(VOTU_CLUSTERING.out.votus)
    IPHOP_PREDICT(SPLIT_VOTUS.out.chunks.flatten(), ch_iphop_db)

    IPHOP_PREDICT.out.genus
        .collectFile(name: 'Host_prediction_to_genus_m90.csv', storeDir: "${params.outdir}/09_iphop", keepHeader: true, skip: 1, sort: true)
    IPHOP_PREDICT.out.genome
        .collectFile(name: 'Host_prediction_to_genome_m90.csv', storeDir: "${params.outdir}/09_iphop", keepHeader: true, skip: 1, sort: true)
}
