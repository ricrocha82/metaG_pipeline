#!/usr/bin/env nextflow
/*
 * metaG_pipeline: umbrella pipeline chaining all three stages in one run.
 *
 *   pipeline1-qc-assembly   -> QC, trimming, host/PhiX removal, taxonomy,
 *                              assembly, evaluation, gene prediction/annotation
 *   pipeline2-binning-mags  -> binning, refinement, cross-sample dereplication,
 *                              MAG taxonomy, MAG abundance
 *   pipeline3-viromics      -> virus detection, CheckV, clustering,
 *                              abundance, microdiversity
 *
 * pipeline2 and pipeline3 both branch off pipeline1's assemblies + clean
 * reads and run in parallel -- pipeline3 doesn't depend on pipeline2's
 * bins, and pipeline2 doesn't depend on pipeline3's viral contigs. Their
 * outputs are wired directly in-memory (no intermediate CSV round-trip)
 * via each pipeline's QC_ASSEMBLY/BINNING_MAGS/VIROMICS subworkflow.
 *
 * Each stage remains independently runnable via its own main.nf (see each
 * subdirectory's README) -- e.g. `nextflow run pipeline2-binning-mags/main.nf
 * --binning_input ...` still works standalone, unchanged. Params were
 * renamed to be pipeline-scoped (qc_*, binning_*, viromics_*) specifically
 * so this umbrella can run all three in one Nextflow session without one
 * pipeline's config silently clobbering another's (see root README).
 *
 * Use --skip_binning / --skip_viromics to run only a subset of stages.
 */

nextflow.enable.dsl = 2

include { QC_ASSEMBLY  } from './pipeline1-qc-assembly/main.nf'
include { BINNING_MAGS } from './pipeline2-binning-mags/main.nf'
include { VIROMICS     } from './pipeline3-viromics/main.nf'

workflow {

    // ---- pipeline1: QC, trimming, assembly, annotation (always runs) ----
    QC_ASSEMBLY(null)

    // ---- pipeline2: binning + MAGs, fed directly from pipeline1's output ----
    if (!params.skip_binning) {
        BINNING_MAGS(QC_ASSEMBLY.out.assemblies)
    } else {
        log.info "=== Skipping pipeline2-binning-mags (--skip_binning true) ==="
    }

    // ---- pipeline3: viromics, fed directly from pipeline1's output ----
    if (!params.skip_viromics) {
        // Adapt pipeline1's [meta, assembly, r1, r2] shape into pipeline3's
        // expected samples_ch shape: tuple(sample_id, read_path, read_type, raw_read1, raw_read2)
        ch_viromics_samples = QC_ASSEMBLY.out.assemblies
            .map { meta, asm, r1, r2 -> tuple(meta.id, asm, 'paired', r1, r2) }

        VIROMICS(ch_viromics_samples)
    } else {
        log.info "=== Skipping pipeline3-viromics (--skip_viromics true) ==="
    }
}
