#!/usr/bin/env nextflow
/*
 * AREE metabolomics workflow scaffold
 * ------------------------------------
 * STATUS: unexecuted structural scaffold. See README.md "What has and has
 * not been run" before treating any output of a real `nextflow run` of this
 * file as validated. Every process below has real, syntactically checked
 * script logic (pandas / limma / Quarto) but none has been run against a
 * real container runtime in this build.
 *
 * Implements CLAUDE.md Layer 2.D (Metabolomics):
 *   metabolite feature table intake -> annotation confidence tracking ->
 *   normalization and QC -> differential abundance -> pathway/metabolite-
 *   class mapping -> standardized feature-level evidence output.
 *
 * Two modes (params.mode), matching docs/design.md section 7:
 *
 *   raw_reanalysis:
 *     params.raw_feature_table (features x samples intensity TSV/CSV, e.g.
 *     XCMS/MZmine aligned feature table) + params.sample_sheet
 *       -> FEATURE_TABLE_INTAKE
 *       -> ANNOTATION_CONFIDENCE
 *       -> NORMALIZE_QC
 *       -> DIFFERENTIAL_ABUNDANCE
 *       -> PATHWAY_MAPPING
 *       -> STANDARDIZE_OUTPUT
 *       -> EMIT_MANIFEST
 *       -> RENDER_REPORT
 *
 *   processed_results_harmonization (default):
 *     params.processed_results (already-shaped differential feature table)
 *       -> STANDARDIZE_OUTPUT (validate/reshape only)
 *       -> EMIT_MANIFEST (flags that raw QC/normalization were not
 *          independently verified)
 *       -> RENDER_REPORT
 *
 * The output of STANDARDIZE_OUTPUT in either mode is the exact TSV schema
 * consumed by:
 *   aree harmonize --study STUDY_ID --input <standardized_tsv>
 * (see src/harmonize/metabolomics.py and README.md of this workflow).
 */

nextflow.enable.dsl = 2

include { FEATURE_TABLE_INTAKE  } from '../../modules/metabolomics/feature_table_intake.nf'
include { ANNOTATION_CONFIDENCE } from '../../modules/metabolomics/annotation_confidence.nf'
include { NORMALIZE_QC          } from '../../modules/metabolomics/normalize_qc.nf'
include { DIFFERENTIAL_ABUNDANCE } from '../../modules/metabolomics/differential_abundance.nf'
include { PATHWAY_MAPPING       } from '../../modules/metabolomics/pathway_mapping.nf'
include { STANDARDIZE_OUTPUT    } from '../../modules/metabolomics/standardize_output.nf'
include { EMIT_MANIFEST         } from '../../modules/metabolomics/emit_manifest.nf'
include { RENDER_REPORT         } from '../../modules/metabolomics/render_report.nf'

workflow {

    log.info """
    AREE metabolomics workflow (scaffold)
    --------------------------------------
    mode          : ${params.mode}
    study_id      : ${params.metabolomics.study_id}
    comparison_id : ${params.metabolomics.comparison_id}
    outdir        : ${params.outdir}
    """.stripIndent()

    if (!(params.mode in ['raw_reanalysis', 'processed_results_harmonization'])) {
        error "params.mode must be 'raw_reanalysis' or 'processed_results_harmonization', got: ${params.mode}"
    }

    study_id_ch      = Channel.value(params.metabolomics.study_id)
    comparison_id_ch = Channel.value(params.metabolomics.comparison_id)
    mode_ch          = Channel.value(params.mode)
    report_template  = Channel.fromPath("${projectDir}/assets/report_template.qmd", checkIfExists: true)

    if (params.mode == 'raw_reanalysis') {

        if (!params.metabolomics.raw_feature_table) {
            error "params.metabolomics.raw_feature_table is required when params.mode == 'raw_reanalysis'"
        }
        if (!params.metabolomics.sample_sheet) {
            error "params.metabolomics.sample_sheet is required when params.mode == 'raw_reanalysis'"
        }

        raw_feature_table_ch = Channel.fromPath(params.metabolomics.raw_feature_table, checkIfExists: true)
        sample_sheet_ch      = Channel.fromPath(params.metabolomics.sample_sheet, checkIfExists: true)

        annotation_map_ch = params.metabolomics.metabolite_annotation_map
            ? Channel.fromPath(params.metabolomics.metabolite_annotation_map, checkIfExists: true)
            : Channel.fromPath("${projectDir}/assets/NO_FILE")

        class_map_ch = params.metabolomics.metabolite_class_map
            ? Channel.fromPath(params.metabolomics.metabolite_class_map, checkIfExists: true)
            : Channel.fromPath("${projectDir}/assets/NO_FILE")

        FEATURE_TABLE_INTAKE(
            study_id_ch,
            comparison_id_ch,
            raw_feature_table_ch,
            sample_sheet_ch
        )

        ANNOTATION_CONFIDENCE(
            study_id_ch,
            comparison_id_ch,
            FEATURE_TABLE_INTAKE.out.validated_table.map { it[2] },
            annotation_map_ch
        )

        NORMALIZE_QC(
            study_id_ch,
            comparison_id_ch,
            ANNOTATION_CONFIDENCE.out.annotated_table.map { it[2] },
            sample_sheet_ch,
            Channel.value(params.metabolomics.qc_cv_threshold)
        )

        DIFFERENTIAL_ABUNDANCE(
            study_id_ch,
            comparison_id_ch,
            NORMALIZE_QC.out.normalized_table.map { it[2] },
            sample_sheet_ch
        )

        PATHWAY_MAPPING(
            study_id_ch,
            comparison_id_ch,
            DIFFERENTIAL_ABUNDANCE.out.results.map { it[2] },
            class_map_ch
        )

        raw_results_for_standardize = PATHWAY_MAPPING.out.classified_table.map { it[2] }

        qc_json_for_manifest = FEATURE_TABLE_INTAKE.out.intake_qc
            .mix(ANNOTATION_CONFIDENCE.out.annotation_qc)
            .mix(NORMALIZE_QC.out.normalization_qc)
            .mix(PATHWAY_MAPPING.out.mapping_qc)
            .collect()

        input_files_for_checksum = raw_feature_table_ch
            .mix(sample_sheet_ch)
            .collect()

        versions_ch = FEATURE_TABLE_INTAKE.out.versions
            .mix(ANNOTATION_CONFIDENCE.out.versions)
            .mix(NORMALIZE_QC.out.versions)
            .mix(DIFFERENTIAL_ABUNDANCE.out.versions)
            .mix(PATHWAY_MAPPING.out.versions)

    } else {

        if (!params.metabolomics.processed_results) {
            error "params.metabolomics.processed_results is required when params.mode == 'processed_results_harmonization'"
        }

        processed_results_ch = Channel.fromPath(params.metabolomics.processed_results, checkIfExists: true)

        raw_results_for_standardize = processed_results_ch

        qc_json_for_manifest = Channel.empty().collect()
        input_files_for_checksum = processed_results_ch.collect()
        versions_ch = Channel.empty()
    }

    STANDARDIZE_OUTPUT(
        study_id_ch,
        comparison_id_ch,
        raw_results_for_standardize
    )

    EMIT_MANIFEST(
        study_id_ch,
        comparison_id_ch,
        mode_ch,
        STANDARDIZE_OUTPUT.out.standardized_tsv.map { it[2] },
        input_files_for_checksum,
        qc_json_for_manifest
    )

    RENDER_REPORT(
        study_id_ch,
        comparison_id_ch,
        mode_ch,
        STANDARDIZE_OUTPUT.out.standardized_tsv.map { it[2] },
        EMIT_MANIFEST.out.manifest,
        report_template
    )

    versions_ch
        .mix(STANDARDIZE_OUTPUT.out.versions)
        .mix(EMIT_MANIFEST.out.versions)
        .collectFile(name: 'versions_combined.yml', storeDir: "${params.outdir}/metabolomics/pipeline_info")
}
