#!/usr/bin/env nextflow

// ============================================================
// VIRION 3 - Modular Plug-and-Play Pipeline
//
// Architecture: every stage has a paired LOADER sub-workflow
// that reads existing outputs from disk — same channel shapes
// as if the upstream processor had just run. Entry workflows
// pick the right loader for their starting point and wire
// only the stages they need.
//
// Entry workflows (each fully independent):
//
//   Full pipeline (start from raw assemblies):
//     VIRUS_DETECTION  →  detection + contig extraction
//     CHECKV           →  + CheckV quality filtering
//     CLUSTERING        →  + per-sample ANI clustering
//     MAPPING          →  + read mapping → sorted BAM
//     ABUNDANCE        →  + CoverM vOTU table
//     ENHANCED         →  + cross-sample enhanced virome
//     MICRODIVERSITY   →  + MetaPop microdiversity (full run)
//     FULL_PIPELINE    →  alias for MICRODIVERSITY
//
//   Resume from existing outputs (skip upstream steps):
//     FROM_CHECKV       →  load CheckV FASTAs   → cluster → map → abundance → enhanced
//     FROM_CLUSTERING   →  load cluster FASTAs  → map → abundance → enhanced
//     FROM_MAPPING      →  load BAM files       → abundance → enhanced
//     FROM_ABUNDANCE    →  load BAM + cluster   → enhanced → microdiversity
//     FROM_ENHANCED     →  load enhanced BAM/ref → microdiversity only
//
// Parameters for loaders (all default null — only set what you need):
//   --checkv_dir      directory with <sid>_final_viral_contigs.fasta files
//   --cluster_dir     directory with <sid>.self-blastn.clusters.fna files
//   --bam_dir         directory with <sid>_sorted.bam files
//   --enhanced_ref    path to the single enhanced clusters.fna file
//   --enhanced_bam    path to the single enhanced sorted BAM file
// ============================================================

// ------ Default Parameters ------
params.threads            = 40
params.detect             = "genomad"              // 'genomad' | 'both'
params.input_csv          = "${launchDir}/data/samples.csv"
params.outdir             = "${launchDir}/results"
params.virus              = "${params.outdir}/virus"
params.vs2                = "${params.virus}/vs2"
params.genomad            = "${params.virus}/genomad"
params.genomad_db         = "/fs/project/PAS1117/ricardo/genomad_db"
params.virus_length       = 5000
params.virus_score        = 0.7
params.virus_completeness = 0
params.combined           = "${params.virus}/combined"
params.checkv             = "${params.virus}/checkv"
params.cluster            = "${params.virus}/cluster"
params.mapping            = "${params.virus}/mapping"
params.abundance          = "${params.virus}/abundance"
params.microdiversity     = "${params.virus}/microdiversity"
params.skip_microdiversity = false                 // set to true to skip MetaPop
params.length             = 5000                   // min contig length for FILTER_CHECKV
params.short_reads        = null                   // assembly filename used as short-read representative
params.hybrid_reads       = null                   // assembly filename for hybrid method (optional)
params.long_reads         = null                   // assembly filename for long-read method (optional)

// ------ Loader input params (null = not loading from disk) ------
params.checkv_dir         = null   // dir with *_final_viral_contigs.fasta
params.cluster_dir        = null   // dir with *.self-blastn.clusters.fna
params.cluster_suffix     = "_final_viral_contigs.fasta.self-blastn.clusters.fna"  // suffix stripped from filename to recover sample_id
params.bam_dir            = null   // dir with *_sorted.bam
params.enhanced_ref       = null   // single path: enhanced clusters.fna
params.enhanced_bam       = null   // single path: enhanced sorted BAM

// ------ Module Imports ------
include { VIRSORTER }                              from './modules/virus_detection.nf'
include { GENOMAD_END_TO_END }                     from './modules/virus_detection.nf'
include { COMBINE_VS2_GENOMAD }                    from './modules/virus_detection.nf'
include { SELECT_CONTIGS }                         from './modules/virus_detection.nf'
include { SELECT_CONTIGS_GENOMAD }                 from './modules/virus_detection.nf'
include { EXTRACT_CONTIGS }                        from './modules/virus_detection.nf'
include { CHECKV_END_TO_END }                      from './modules/checkv.nf'
include { FILTER_CHECKV }                          from './modules/checkv.nf'
include { CLUSTERING }                             from './modules/clustering.nf'
include { CLUSTERING as CLUSTERING_ENH }           from './modules/clustering.nf'
include { CONCAT_FASTA }                           from './modules/clustering.nf'
include { MAPPING_LONG_READS }                     from './modules/clustering.nf'
include { MAPPING_SHORT_READS }                    from './modules/clustering.nf'
include { MAPPING_SHORT_READS as MAPPING_SHORT_READS_ENH } from './modules/clustering.nf'
include { SAMTOOLS }                               from './modules/clustering.nf'
include { SAMTOOLS as SAMTOOLS_ENH }               from './modules/clustering.nf'
include { COVERM_ABUNDANCE }                       from './modules/clustering.nf'
include { COVERM_ABUNDANCE as COVERM_ABUNDANCE_ENH } from './modules/clustering.nf'
include { CREATE_ABUND_TB }                        from './modules/clustering.nf'
include { METAPOP_MICRODIVERSITY }                 from './modules/microdiversity.nf'


// ============================================================
//  LOADERS
//  Each loader reads existing outputs from disk and emits the
//  exact same channel shapes as the paired processor sub-workflow.
//  This lets any entry workflow skip upstream steps cleanly.
// ============================================================

// ------------------------------------------------------------
// LOADER: samples CSV → structured channel
// Always required (reads are needed for mapping even when
// upstream steps are skipped).
// Emits: tuple(sample_id, read_path, read_type, raw_read1, raw_read2)
// ------------------------------------------------------------
workflow LOAD_SAMPLES {
    main:
        samples_ch = Channel
            .fromPath(params.input_csv)
            .splitCsv(header: true)
            .map { row ->
                // Classify as paired-end if raw_read2 is present and non-empty
                def read_type = (row.raw_read2 && row.raw_read2.trim() != "") ? "paired" : "single"
                tuple(
                    row.sample_id,
                    file(row.read_path),
                    read_type,
                    row.raw_read1               ? file(row.raw_read1) : null,
                    row.raw_read2?.trim() != "" ? file(row.raw_read2) : null
                )
            }
    emit:
        samples = samples_ch   // tuple(sample_id, read_path, read_type, raw_read1, raw_read2)
}

// ------------------------------------------------------------
// LOADER: existing CheckV filtered FASTAs from --checkv_dir
// Replaces: DETECT_VIRUSES → RUN_CHECKV
// File pattern: <sample_id>_final_viral_contigs.fasta
// Emits: tuple(sample_id, fasta)  — same as RUN_CHECKV.out.filtered_fasta
// ------------------------------------------------------------
workflow LOAD_CHECKV_FASTAS {
    main:
        if (!params.checkv_dir) { error "--checkv_dir is required for this entry workflow" }

        filtered_fasta_ch = Channel
            .fromPath("${params.checkv_dir}/*_final_viral_contigs.fasta")
            .map { fasta ->
                // Strip the suffix to recover sample_id
                def sid = fasta.name.replaceAll(/_final_viral_contigs\.fasta$/, "")
                tuple(sid, fasta)
            }
    emit:
        filtered_fasta = filtered_fasta_ch   // tuple(sample_id, fasta)
}

// ------------------------------------------------------------
// LOADER: existing cluster FASTAs from --cluster_dir
// Replaces: DETECT_VIRUSES → RUN_CHECKV → RUN_CLUSTERING
// File pattern: <sample_id>.self-blastn.clusters.fna
// Emits: cluster_fasta -> tuple(sample_id, fasta)
//        fasta_list    -> list<path>   (for RUN_ENHANCED)
// ------------------------------------------------------------
workflow LOAD_CLUSTER_FASTAS {
    take:
        samples_ch   // tuple(sample_id, ...) — used to broadcast the shared reference to every sample
    main:
        if (!params.cluster_dir) { error "--cluster_dir is required for this entry workflow" }

        // Exactly one shared cluster fasta per --cluster_dir (e.g. concatenated_all...clusters.fna)
        shared_fasta_ch = Channel
            .fromPath("${params.cluster_dir}/**/*${params.cluster_suffix}")
            .first()

        // Broadcast the shared reference to every sample_id present in samples_ch
        cluster_fasta_ch = samples_ch
            .map { sid, rp, rt, r1, r2 -> sid }
            .combine(shared_fasta_ch)

        // Single-element list for RUN_ENHANCED's CONCAT_FASTA input
        fasta_list_ch = shared_fasta_ch.map { f -> [f] }
    emit:
        cluster_fasta = cluster_fasta_ch
        fasta_list    = fasta_list_ch
}

// ------------------------------------------------------------
// LOADER: existing sorted BAM files from --bam_dir
// Replaces: everything up through RUN_MAPPING
// File pattern: <sample_id>_sorted.bam
// Emits: bam_files -> tuple(sample_id, bam)  — same as RUN_MAPPING.out.bam_files
// ------------------------------------------------------------
workflow LOAD_BAM_FILES {
    main:
        if (!params.bam_dir) { error "--bam_dir is required for this entry workflow" }

        bam_files_ch = Channel
            .fromPath("${params.bam_dir}/*_sorted.bam")
            .map { bam ->
                // Strip the suffix to recover sample_id
                def sid = bam.name.replaceAll(/_sorted\.bam$/, "")
                tuple(sid, bam)
            }
    emit:
        bam_files = bam_files_ch   // tuple(sample_id, bam)
}

// ------------------------------------------------------------
// LOADER: existing enhanced virome outputs from --enhanced_ref / --enhanced_bam
// Replaces: RUN_ENHANCED
// Emits: enhanced_ref -> path
//        enhanced_bam -> tuple('enhanced', bam)
// ------------------------------------------------------------
workflow LOAD_ENHANCED {
    main:
        if (!params.enhanced_ref) { error "--enhanced_ref is required for this entry workflow" }
        if (!params.enhanced_bam) { error "--enhanced_bam is required for this entry workflow" }

        // Wrap single files into the channel shapes RUN_MICRODIVERSITY expects
        enhanced_ref_ch = Channel.value(file(params.enhanced_ref))
        enhanced_bam_ch = Channel.value(tuple('enhanced', file(params.enhanced_bam)))

    emit:
        enhanced_ref = enhanced_ref_ch   // path
        enhanced_bam = enhanced_bam_ch   // tuple('enhanced', bam)
}


// ============================================================
//  PROCESSORS
//  Each processor runs one pipeline stage and emits its outputs.
//  Processors never load files from disk — that is the loader's job.
// ============================================================

// ------------------------------------------------------------
// PROCESSOR: Virus Detection
// Runs VirSorter2 and/or geNomad, selects viral contigs
// Emits: filt_fasta -> tuple(sample_id, fasta)
// ------------------------------------------------------------
workflow RUN_DETECTION {
    take:
        samples_ch   // tuple(sample_id, read_path, read_type, raw_read1, raw_read2)

    main:
        // Detection tools only need the assembly FASTA
        initial_ch    = samples_ch.map { sid, rp, rt, r1, r2 -> tuple(sid, rp) }
        genomad_db_ch = Channel.value(params.genomad_db)

        if (params.detect == 'both') {
            virsorter_ch      = VIRSORTER(initial_ch)
            genomad_ch        = GENOMAD_END_TO_END(initial_ch, genomad_db_ch)
            genomad_virus_fna = genomad_ch.genomad_virus

            // Pair VS2 and geNomad scores by sample_id then combine
            paired_ch = virsorter_ch.viral_score
                .join(genomad_ch.virus_summary, by: 0)
                .map { sid, vs2, gmd -> tuple(sid, vs2, gmd) }

            combined_tbl_ch = COMBINE_VS2_GENOMAD(paired_ch)

            filter_input_ch = initial_ch
                .combine(combined_tbl_ch.combined_df.flatten())
                .map { sid, rp, combined -> tuple(sid, combined, rp) }

            selected_ch = SELECT_CONTIGS(filter_input_ch)

        } else if (params.detect == 'genomad') {
            genomad_ch        = GENOMAD_END_TO_END(initial_ch, genomad_db_ch)
            genomad_virus_fna = genomad_ch.genomad_virus

            filter_input_ch = initial_ch
                .join(genomad_ch.virus_summary, by: 0)
                .map { sid, rp, summary -> tuple(sid, summary, rp) }

            selected_ch = SELECT_CONTIGS_GENOMAD(filter_input_ch)

        } else {
            error "Invalid --detect value: '${params.detect}'. Use 'genomad' or 'both'."
        }

        // Extract contig sequences using the selected contig ID list
        extract_in_ch = selected_ch
            .join(initial_ch, by: 0)
            .join(genomad_virus_fna, by: 0)
            .map { sid, contigs, rp, gvirus -> tuple(sid, contigs, rp, gvirus) }

        filtered_ch = EXTRACT_CONTIGS(extract_in_ch)

    emit:
        filt_fasta = filtered_ch.filt_fasta   // tuple(sample_id, fasta)
}

// ------------------------------------------------------------
// PROCESSOR: CheckV Quality Assessment
// Emits: filtered_fasta -> tuple(sample_id, fasta)
// ------------------------------------------------------------
workflow RUN_CHECKV {
    take:
        filt_fasta_ch   // tuple(sample_id, fasta) from RUN_DETECTION or LOAD_CHECKV_FASTAS

    main:
        checkv_ch = CHECKV_END_TO_END(filt_fasta_ch)

        // Pair quality summary with its combined FASTA for the filter step
        checkv_pair_ch = checkv_ch.quality_summary
            .join(checkv_ch.checkv_fasta, by: 0)
            .map { sid, qsum, fasta -> tuple(sid, qsum, fasta) }

        filter_check_ch = FILTER_CHECKV(checkv_pair_ch)

    emit:
        filtered_fasta = filter_check_ch.filtered_checkv_fasta   // tuple(sample_id, fasta)
}

// ------------------------------------------------------------
// PROCESSOR: Per-sample ANI Clustering
// Emits: cluster_fasta -> tuple(sample_id, fasta)
//        fasta_list    -> list<path>
// ------------------------------------------------------------
workflow RUN_CLUSTERING {
    take:
        filtered_fasta_ch   // tuple(sample_id, fasta) from RUN_CHECKV or LOAD_CHECKV_FASTAS

    main:
        clustering_ch = CLUSTERING(filtered_fasta_ch)

        // Collect all cluster FASTAs into a flat list for RUN_ENHANCED
        fasta_list_ch = clustering_ch.cluster_fasta.map { sid, f -> f }.collect()

    emit:
        cluster_fasta = clustering_ch.cluster_fasta   // tuple(sample_id, fasta)
        fasta_list    = fasta_list_ch                  // list<path>
}

// ------------------------------------------------------------
// PROCESSOR: Read Mapping → Sorted BAM
// Routes long / short reads to the appropriate mapper
// Emits: bam_files -> tuple(sample_id, bam)
// ------------------------------------------------------------
workflow RUN_MAPPING {
    take:
        samples_ch        // tuple(sample_id, read_path, read_type, raw_read1, raw_read2)
        cluster_fasta_ch  // tuple(sample_id, fasta) — per-sample cluster reference

    main:
        // Build mapping input by joining sample metadata with its cluster reference
        mapping_input_ch = samples_ch
            .join(cluster_fasta_ch, by: 0)
            .map { sid, rp, rt, r1, r2, cfasta ->
                if (rt == "paired" && r1 && r2)   { tuple(sid, cfasta, r1, r2, "short") }
                else if (r1)                       { tuple(sid, cfasta, r1,     "long")  }
                else                               { tuple(sid, cfasta, rp,     "long")  }
            }

        // Route by read type
        long_reads_ch  = mapping_input_ch
            .filter { it.last() == "long"  }
            .map    { sid, cf, ref, _t -> tuple(sid, cf, ref) }

        short_reads_ch = mapping_input_ch
            .filter { it.last() == "short" }
            .map    { sid, cf, r1, r2, _t -> tuple(sid, cf, r1, r2) }

        // Map → SAM → sorted BAM
        all_sam_ch  = MAPPING_LONG_READS(long_reads_ch).sam_file
            .mix(MAPPING_SHORT_READS(short_reads_ch).sam_file)
        samtools_ch = SAMTOOLS(all_sam_ch)

    emit:
        bam_files = samtools_ch.bam_file   // tuple(sample_id, bam)
}

// ------------------------------------------------------------
// PROCESSOR: CoverM Abundance + vOTU Table
// Emits: coverm_output -> tuple(sample_id, tsv)
// ------------------------------------------------------------
workflow RUN_ABUNDANCE {
    take:
        bam_files_ch   // tuple(sample_id, bam) from RUN_MAPPING or LOAD_BAM_FILES

    main:
        coverm_ch    = COVERM_ABUNDANCE(bam_files_ch)
        coverm_paths = coverm_ch.coverm_output.map { it[1] }.collect()
        CREATE_ABUND_TB(coverm_paths)

    emit:
        coverm_output = coverm_ch.coverm_output   // tuple(sample_id, tsv)
}

// ------------------------------------------------------------
// PROCESSOR: Enhanced (cross-sample) Virome
// Pools all cluster FASTAs, re-clusters, maps short reads
// Emits: enhanced_ref -> path
//        enhanced_bam -> tuple('enhanced', bam)
// ------------------------------------------------------------
workflow RUN_ENHANCED {
    take:
        samples_ch     // tuple(sample_id, read_path, read_type, raw_read1, raw_read2)
        fasta_list_ch  // list<path> — all per-sample cluster FASTAs

    main:
        // Concatenate all per-sample FASTAs and re-cluster cross-sample
        concat_ch         = CONCAT_FASTA(fasta_list_ch)
        enhanced_input_ch = concat_ch.concat_fasta.map { ref -> tuple('enhanced', ref) }
        clustering_all_ch = CLUSTERING_ENH(enhanced_input_ch)

        // Strip the sample_id wrapper — there is only one enhanced reference
        enhanced_ref_ch = clustering_all_ch.cluster_fasta
            .filter { sid, ref -> sid == 'enhanced' }
            .map    { sid, ref -> ref }

        // Use the first available paired-end sample as the short-read representative
        first_paired_sid_ch = samples_ch
            .filter { sid, rp, rt, r1, r2 -> rt == 'paired' && r1 && r2 }
            .map    { sid, rp, rt, r1, r2 -> sid }
            .first()

        // Build mapping input: (enhanced, ref, r1, r2)
        enhanced_map_in_ch = first_paired_sid_ch
            .combine(samples_ch)
            .filter { pid, sid, rp, rt, r1, r2 -> sid == pid && rt == 'paired' && r1 && r2 }
            .map    { pid, sid, rp, rt, r1, r2 -> tuple('enhanced', r1, r2) }
            .combine(enhanced_ref_ch)
            .map    { sid, r1, r2, ref -> tuple(sid, ref, r1, r2) }

        // Map → SAM → sorted BAM → CoverM
        enh_sam_ch      = MAPPING_SHORT_READS_ENH(enhanced_map_in_ch).sam_file
        samtools_enh_ch = SAMTOOLS_ENH(enh_sam_ch)
        COVERM_ABUNDANCE_ENH(samtools_enh_ch.bam_file)

    emit:
        enhanced_ref = enhanced_ref_ch              // path
        enhanced_bam = samtools_enh_ch.bam_file     // tuple('enhanced', bam)
}

// ------------------------------------------------------------
// PROCESSOR: MetaPop Microdiversity
// Runs MetaPop for short, enhanced, and optionally hybrid/long methods
// ------------------------------------------------------------
workflow RUN_MICRODIVERSITY {
    take:
        samples_ch       // tuple(sample_id, read_path, read_type, raw_read1, raw_read2)
        bam_files_ch     // tuple(sample_id, bam)  — per-sample BAMs
        cluster_fasta_ch // tuple(sample_id, fasta) — per-sample cluster references
        enhanced_ref_ch  // path — enhanced cluster reference
        enhanced_bam_ch  // tuple('enhanced', bam)

    main:
        // --- Enhanced method pack (always included) ---
        enhanced_bams_ch = enhanced_bam_ch.map { sid, bam -> bam }.collect()
        enhanced_pack    = enhanced_bams_ch
            .combine(enhanced_ref_ch)
            .map { bams, ref -> tuple('enhanced', bams, ref) }

        // --- Short method: resolve sample_id from --short_reads filename ---
        short_reads_name = file(params.short_reads).name

        short_sid_ch = samples_ch
            .filter { sid, rp, rt, r1, r2 -> rp.name == short_reads_name }
            .map    { sid, rp, rt, r1, r2 -> sid }
            .first()

        short_bams_ch = short_sid_ch
            .combine(bam_files_ch)
            .filter { ssid, sid, bam -> sid == ssid }
            .map    { ssid, sid, bam -> bam }
            .collect()

        short_ref_ch = short_sid_ch
            .combine(cluster_fasta_ch)
            .filter { ssid, sid, ref -> sid == ssid }
            .map    { ssid, sid, ref -> ref }

        short_pack = short_bams_ch
            .combine(short_ref_ch)
            .map { bams, ref -> tuple('short', bams, ref) }

        methods_ch = short_pack.mix(enhanced_pack)

        // --- Hybrid method (optional) ---
        if (params.hybrid_reads) {
            hybrid_name = file(params.hybrid_reads).name

            hybrid_sid_ch = samples_ch
                .filter { sid, rp, rt, r1, r2 -> rp.name == hybrid_name }
                .map    { sid, rp, rt, r1, r2 -> sid }
                .first()

            hybrid_bams_ch = hybrid_sid_ch
                .combine(bam_files_ch)
                .filter { hsid, sid, bam -> sid == hsid }
                .map    { hsid, sid, bam -> bam }
                .collect()

            hybrid_ref_ch = hybrid_sid_ch
                .combine(cluster_fasta_ch)
                .filter { hsid, sid, ref -> sid == hsid }
                .map    { hsid, sid, ref -> ref }

            methods_ch = methods_ch.mix(
                hybrid_bams_ch.combine(hybrid_ref_ch)
                    .map { bams, ref -> tuple('hybrid', bams, ref) }
            )
        }

        // --- Long-read method (optional) ---
        if (params.long_reads) {
            long_name = file(params.long_reads).name

            long_sid_ch = samples_ch
                .filter { sid, rp, rt, r1, r2 -> rp.name == long_name }
                .map    { sid, rp, rt, r1, r2 -> sid }
                .first()

            long_bams_ch = long_sid_ch
                .combine(bam_files_ch)
                .filter { lsid, sid, bam -> sid == lsid }
                .map    { lsid, sid, bam -> bam }
                .collect()

            long_ref_ch = long_sid_ch
                .combine(cluster_fasta_ch)
                .filter { lsid, sid, ref -> sid == lsid }
                .map    { lsid, sid, ref -> ref }

            methods_ch = methods_ch.mix(
                long_bams_ch.combine(long_ref_ch)
                    .map { bams, ref -> tuple('long', bams, ref) }
            )
        }

        // Run MetaPop for each method
        METAPOP_MICRODIVERSITY(methods_ch)
}


// ============================================================
//  ENTRY WORKFLOWS — FULL CHAIN (start from raw assemblies)
//  Each stops at the named stage; all prior stages also run.
// ============================================================

workflow VIRUS_DETECTION {
    log.info "=== VIRUS_DETECTION ==="
    LOAD_SAMPLES()
    RUN_DETECTION(LOAD_SAMPLES.out.samples)
}

workflow CHECKV {
    log.info "=== CHECKV ==="
    LOAD_SAMPLES()
    RUN_DETECTION(LOAD_SAMPLES.out.samples)
    RUN_CHECKV(RUN_DETECTION.out.filt_fasta)
}

workflow STEP_CLUSTERING {
    log.info "=== STEP_CLUSTERING ==="
    LOAD_SAMPLES()
    RUN_DETECTION(LOAD_SAMPLES.out.samples)
    RUN_CHECKV(RUN_DETECTION.out.filt_fasta)
    RUN_CLUSTERING(RUN_CHECKV.out.filtered_fasta)
}


workflow MAPPING {
    log.info "=== MAPPING ==="
    LOAD_SAMPLES()
    RUN_DETECTION(LOAD_SAMPLES.out.samples)
    RUN_CHECKV(RUN_DETECTION.out.filt_fasta)
    RUN_CLUSTERING(RUN_CHECKV.out.filtered_fasta)
    RUN_MAPPING(LOAD_SAMPLES.out.samples, RUN_CLUSTERING.out.cluster_fasta)
}

workflow ABUNDANCE {
    log.info "=== ABUNDANCE ==="
    LOAD_SAMPLES()
    RUN_DETECTION(LOAD_SAMPLES.out.samples)
    RUN_CHECKV(RUN_DETECTION.out.filt_fasta)
    RUN_CLUSTERING(RUN_CHECKV.out.filtered_fasta)
    RUN_MAPPING(LOAD_SAMPLES.out.samples, RUN_CLUSTERING.out.cluster_fasta)
    RUN_ABUNDANCE(RUN_MAPPING.out.bam_files)
}

workflow ENHANCED {
    log.info "=== ENHANCED ==="
    LOAD_SAMPLES()
    RUN_DETECTION(LOAD_SAMPLES.out.samples)
    RUN_CHECKV(RUN_DETECTION.out.filt_fasta)
    RUN_CLUSTERING(RUN_CHECKV.out.filtered_fasta)
    RUN_MAPPING(LOAD_SAMPLES.out.samples, RUN_CLUSTERING.out.cluster_fasta)
    RUN_ABUNDANCE(RUN_MAPPING.out.bam_files)
    RUN_ENHANCED(LOAD_SAMPLES.out.samples, RUN_CLUSTERING.out.fasta_list)
}

// Entry: full pipeline — detection → CheckV → clustering → mapping → abundance → enhanced → (MetaPop)
// Use --skip_microdiversity true to stop before MetaPop
workflow FULL_PIPELINE {
    log.info "=== FULL PIPELINE ==="
    LOAD_SAMPLES()
    RUN_DETECTION(LOAD_SAMPLES.out.samples)
    RUN_CHECKV(RUN_DETECTION.out.filt_fasta)
    RUN_CLUSTERING(RUN_CHECKV.out.filtered_fasta)
    RUN_MAPPING(LOAD_SAMPLES.out.samples, RUN_CLUSTERING.out.cluster_fasta)
    RUN_ABUNDANCE(RUN_MAPPING.out.bam_files)
    RUN_ENHANCED(LOAD_SAMPLES.out.samples, RUN_CLUSTERING.out.fasta_list)
    if (!params.skip_microdiversity) {
        log.info "=== FULL PIPELINE + Microdiversity ==="
        RUN_MICRODIVERSITY(
            LOAD_SAMPLES.out.samples,
            RUN_MAPPING.out.bam_files,
            RUN_CLUSTERING.out.cluster_fasta,
            RUN_ENHANCED.out.enhanced_ref,
            RUN_ENHANCED.out.enhanced_bam
        )
    } else {
        log.info "=== FULL PIPELINE - skipping Microdiversity ==="
    }
}


// ============================================================
//  ENTRY WORKFLOWS — RESUME FROM EXISTING OUTPUTS
//  Use loaders to inject existing files; run only what is needed.
// ============================================================

// Resume from existing CheckV filtered FASTAs
// Requires: --checkv_dir
// nextflow run main.nf -entry FROM_CHECKV --checkv_dir results/virus/checkv
workflow FROM_CHECKV {
    log.info "=== FROM_CHECKV → clustering → mapping → abundance → enhanced ==="
    LOAD_SAMPLES()
    LOAD_CHECKV_FASTAS()
    RUN_CLUSTERING(LOAD_CHECKV_FASTAS.out.filtered_fasta)
    RUN_MAPPING(LOAD_SAMPLES.out.samples, RUN_CLUSTERING.out.cluster_fasta)
    RUN_ABUNDANCE(RUN_MAPPING.out.bam_files)
    RUN_ENHANCED(LOAD_SAMPLES.out.samples, RUN_CLUSTERING.out.fasta_list)
}

// Resume from existing per-sample cluster FASTAs
// Requires: --cluster_dir
// nextflow run main.nf -entry FROM_CLUSTERING --cluster_dir results/virus/cluster
workflow FROM_CLUSTERING {
    log.info "=== FROM_CLUSTERING → mapping → abundance → enhanced ==="
    LOAD_SAMPLES()
    LOAD_CLUSTER_FASTAS(LOAD_SAMPLES.out.samples)
    RUN_MAPPING(LOAD_SAMPLES.out.samples, LOAD_CLUSTER_FASTAS.out.cluster_fasta)
    RUN_ABUNDANCE(RUN_MAPPING.out.bam_files)
    RUN_ENHANCED(LOAD_SAMPLES.out.samples, LOAD_CLUSTER_FASTAS.out.fasta_list)
}

// Resume from existing sorted BAM files
// Requires: --bam_dir  --cluster_dir
// nextflow run main.nf -entry FROM_MAPPING --bam_dir results/virus/mapping --cluster_dir results/virus/cluster
workflow FROM_MAPPING {
    log.info "=== FROM_MAPPING → abundance → enhanced ==="
    LOAD_SAMPLES()
    LOAD_BAM_FILES()
    LOAD_CLUSTER_FASTAS(LOAD_SAMPLES.out.samples)
    RUN_ABUNDANCE(LOAD_BAM_FILES.out.bam_files)
    RUN_ENHANCED(LOAD_SAMPLES.out.samples, LOAD_CLUSTER_FASTAS.out.fasta_list)
}

// Resume from existing BAMs + cluster FASTAs; run enhanced then microdiversity
// Requires: --bam_dir  --cluster_dir
// nextflow run main.nf -entry FROM_ABUNDANCE --bam_dir ... --cluster_dir ...
workflow FROM_ABUNDANCE {
    log.info "=== FROM_ABUNDANCE → enhanced → microdiversity ==="
    LOAD_SAMPLES()
    LOAD_BAM_FILES()
    LOAD_CLUSTER_FASTAS(LOAD_SAMPLES.out.samples)
    RUN_ENHANCED(LOAD_SAMPLES.out.samples, LOAD_CLUSTER_FASTAS.out.fasta_list)
    // Skip MetaPop if --skip_microdiversity is set
    if (!params.skip_microdiversity) {
        RUN_MICRODIVERSITY(
            LOAD_SAMPLES.out.samples,
            LOAD_BAM_FILES.out.bam_files,
            LOAD_CLUSTER_FASTAS.out.cluster_fasta,
            RUN_ENHANCED.out.enhanced_ref,
            RUN_ENHANCED.out.enhanced_bam
        )
    }
}

// Resume from existing enhanced virome outputs; run microdiversity only
// Requires: --bam_dir  --cluster_dir  --enhanced_ref  --enhanced_bam
// nextflow run main.nf -entry FROM_ENHANCED \
//     --bam_dir ...  --cluster_dir ...  --enhanced_ref ...  --enhanced_bam ...
workflow FROM_ENHANCED {
    log.info "=== FROM_ENHANCED → microdiversity only ==="
    LOAD_SAMPLES()
    LOAD_BAM_FILES()
    LOAD_CLUSTER_FASTAS(LOAD_SAMPLES.out.samples)
    LOAD_ENHANCED()
    // Skip MetaPop if --skip_microdiversity is set
    if (!params.skip_microdiversity) {
        RUN_MICRODIVERSITY(
            LOAD_SAMPLES.out.samples,
            LOAD_BAM_FILES.out.bam_files,
            LOAD_CLUSTER_FASTAS.out.cluster_fasta,
            LOAD_ENHANCED.out.enhanced_ref,
            LOAD_ENHANCED.out.enhanced_bam
        )
    }
}

// ============================================================
// Default entry (no -entry flag)
// ============================================================
workflow {
    FULL_PIPELINE()
}
