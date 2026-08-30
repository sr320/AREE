#!/usr/bin/env nextflow
/*
 * AREE proteomics reanalysis / harmonization workflow
 *
 * STATUS: processed_results_harmonization executes end to end and is covered
 * by current-version CI. The raw abundance-matrix path remains unexecuted. Real Channel factories, real
 * `if (params.mode == ...)` branching, and real per-process pandas/R logic
 * (see modules/proteomics/*.nf headers for per-step honesty notes). This
 * workflow has NOT been executed end-to-end against real data in this build
 * -- no raw proteomics data or compute budget was available in this
 * environment. See README.md "What has and has not been run".
 *
 * Two modes (params.mode), matching docs/design.md section 7:
 *
 *   raw_reanalysis:
 *     params.proteomics.raw_abundance_matrix (wide protein x sample matrix)
 *     params.proteomics.sample_sheet (sample_id -> group mapping)
 *       -> HARMONIZE_INPUT -> NORMALIZE -> MISSINGNESS_REPORT
 *       -> DIFFERENTIAL_ABUNDANCE -> ID_TRANSLATION -> STANDARDIZE_OUTPUT
 *       -> EMIT_MANIFEST -> RENDER_REPORT
 *
 *   processed_results_harmonization:
 *     params.proteomics.processed_results (already-shaped protein-level
 *     differential abundance table)
 *       -> STANDARDIZE_OUTPUT (validate/reshape; compute missingness_percent
 *          from source file if present, else leave null)
 *       -> EMIT_MANIFEST (flags that raw QC/normalization were not
 *          independently verified)
 *       -> RENDER_REPORT
 *
 * IMPORTANT simplification (see README.md): "raw" in this workflow means a
 * raw/unfiltered wide protein (or peptide) abundance matrix, e.g. exported
 * from MaxQuant/Skyline/DIA-NN, NOT raw mass-spectrometry spectra. Spectral
 * search/reprocessing (MaxQuant/Skyline-scale infrastructure) is explicitly
 * out of scope for this scaffold.
 */

nextflow.enable.dsl = 2

include { HARMONIZE_INPUT }        from '../../modules/proteomics/harmonize_input.nf'
include { NORMALIZE }              from '../../modules/proteomics/normalize.nf'
include { MISSINGNESS_REPORT }     from '../../modules/proteomics/missingness_report.nf'
include { DIFFERENTIAL_ABUNDANCE } from '../../modules/proteomics/differential_abundance.nf'
include { ID_TRANSLATION }         from '../../modules/proteomics/id_translation.nf'
include { STANDARDIZE_OUTPUT }     from '../../modules/proteomics/standardize_output.nf'
include { EMIT_MANIFEST }          from '../../modules/proteomics/emit_manifest.nf'
include { RENDER_REPORT }          from '../../modules/proteomics/render_report.nf'

// Function, not a bare assignment: strict DSL2 forbids script-level statements.
def workflowVersion() { '0.1.0-scaffold' }

def requireParam(String name, value) {
    if (!value) {
        exit 1, "Missing required parameter: ${name} (mode=${params.mode})"
    }
    return value
}

workflow {

    if (!params.proteomics?.study_id) {
        exit 1, "params.proteomics.study_id is required"
    }
    if (!params.proteomics?.comparison_id) {
        exit 1, "params.proteomics.comparison_id is required"
    }

    study_id_ch      = Channel.value(params.proteomics.study_id)
    comparison_id_ch = Channel.value(params.proteomics.comparison_id)
    mode_ch          = Channel.value(params.mode)
    start_time_ch    = Channel.value(workflow.start.toString())

    report_template = file("${projectDir}/assets/report_template.qmd")

    if (params.mode == 'raw_reanalysis') {

        raw_matrix_path   = requireParam('proteomics.raw_abundance_matrix', params.proteomics.raw_abundance_matrix)
        sample_sheet_path = requireParam('proteomics.sample_sheet', params.proteomics.sample_sheet)
        id_mapping_path   = requireParam('proteomics.id_mapping_table', params.proteomics.id_mapping_table)

        raw_matrix_ch   = Channel.fromPath(raw_matrix_path, checkIfExists: true)
        sample_sheet_ch = Channel.fromPath(sample_sheet_path, checkIfExists: true)
        id_mapping_ch   = Channel.fromPath(id_mapping_path, checkIfExists: true)

        HARMONIZE_INPUT(study_id_ch, comparison_id_ch, raw_matrix_ch, sample_sheet_ch)

        NORMALIZE(HARMONIZE_INPUT.out.long_table)

        MISSINGNESS_REPORT(NORMALIZE.out.normalized_long)

        DIFFERENTIAL_ABUNDANCE(
            NORMALIZE.out.normalized_long,
            MISSINGNESS_REPORT.out.per_protein
        )

        ID_TRANSLATION(DIFFERENTIAL_ABUNDANCE.out.results, id_mapping_ch)

        STANDARDIZE_OUTPUT(
            ID_TRANSLATION.out.translated,
            mode_ch,
            MISSINGNESS_REPORT.out.per_protein.map { sid, cid, f -> f }
        )

        manifest_qc_input = MISSINGNESS_REPORT.out.summary_json

    } else if (params.mode == 'processed_results_harmonization') {

        processed_results_path = requireParam('proteomics.processed_results', params.proteomics.processed_results)
        processed_results_ch   = Channel.fromPath(processed_results_path, checkIfExists: true)

        processed_tuple_ch = study_id_ch
            .combine(comparison_id_ch)
            .combine(processed_results_ch)
            .map { sid, cid, f -> tuple(sid, cid, f) }

        // No independently-computed missingness file exists in this mode;
        // pass a static, checked-in placeholder (assets/NO_FILE) so the
        // process's optional input channel is satisfied without claiming a
        // real missingness report. Modules recognize this literal filename
        // and treat it as "not supplied" (see standardize_output.nf /
        // emit_manifest.nf).
        no_missingness_file = file("${projectDir}/assets/NO_FILE", checkIfExists: true)

        STANDARDIZE_OUTPUT(
            processed_tuple_ch,
            mode_ch,
            Channel.value(no_missingness_file)
        )

        manifest_qc_input = Channel.value(no_missingness_file)

    } else {
        exit 1, "Unknown params.mode '${params.mode}': expected 'raw_reanalysis' or 'processed_results_harmonization'"
    }

    EMIT_MANIFEST(
        study_id_ch,
        comparison_id_ch,
        mode_ch,
        STANDARDIZE_OUTPUT.out.standardized.map { sid, cid, f -> f },
        manifest_qc_input,
        Channel.value(workflowVersion()),
        start_time_ch
    )

    RENDER_REPORT(
        study_id_ch,
        comparison_id_ch,
        STANDARDIZE_OUTPUT.out.standardized.map { sid, cid, f -> f },
        EMIT_MANIFEST.out.manifest,
        Channel.value(report_template)
    )

    // Neither `params` nor `workflow` resolves inside a closure nested in the
    // workflow body, so snapshot both into locals the handler can capture.
    def _params = params
    def _wf = workflow
    // Moved inside the entry workflow: strict DSL2 (Nextflow 25.10+) rejects
    // `workflow.onComplete` as a script-level statement.
    workflow.onComplete {
        log.info """
        AREE proteomics workflow scaffold complete.
        mode        : ${_params.mode}
        study_id    : ${_params.proteomics?.study_id}
        comparison  : ${_params.proteomics?.comparison_id}
        outdir      : ${_params.outdir}
        status      : ${_wf.success ? 'OK' : 'FAILED'}
        NOTE: this scaffold has not been validated end-to-end against real data.
        """.stripIndent()
    }
}
