#!/usr/bin/env nextflow
/*
 * AREE methylation / WGBS / EM-seq workflow
 *
 * STATUS: unexecuted structural scaffold. See README.md in this directory
 * for a full, honest description of what has and has not been run. In
 * short: this is real DSL2 wiring with real tool CLI invocations in the
 * per-process modules, but no process in this workflow has actually been
 * executed against FASTQ, BAM, or a real reference genome in this build.
 * The `processed_results_harmonization` branch *could* be run today against
 * the shipped demo TSV (data/demo/methylation/*_demo.tsv) since it has no
 * heavyweight bioinformatics dependency beyond Python's standard library and
 * Quarto, but that execution has not been performed as part of this task
 * either.
 *
 * Branches on params.mode:
 *   - raw_reanalysis:
 *       FASTQC -> TRIM_GALORE -> BISMARK_GENOME_PREPARATION + BISMARK_ALIGN
 *       -> BISMARK_DEDUPLICATE -> BISMARK_METHYLATION_EXTRACTOR
 *       -> COVERAGE_FILTER -> DMR_METHYLKIT -> ANNOTATE_REGIONS
 *       -> STANDARDIZE_OUTPUT -> EMIT_MANIFEST -> RENDER_REPORT
 *   - processed_results_harmonization:
 *       STANDARDIZE_OUTPUT (validate/reshape processed_results TSV)
 *       -> EMIT_MANIFEST (with raw-QC-not-verified warning)
 *       -> RENDER_REPORT
 */

nextflow.enable.dsl = 2

include { FASTQC }                          from '../../modules/methylation/fastqc.nf'
include { TRIM_GALORE }                     from '../../modules/methylation/trim_galore.nf'
include { BISMARK_GENOME_PREPARATION }      from '../../modules/methylation/bismark_align.nf'
include { BISMARK_ALIGN }                   from '../../modules/methylation/bismark_align.nf'
include { BISMARK_DEDUPLICATE }             from '../../modules/methylation/bismark_align.nf'
include { BISMARK_METHYLATION_EXTRACTOR }   from '../../modules/methylation/bismark_methylation_extractor.nf'
include { COVERAGE_FILTER }                 from '../../modules/methylation/coverage_filter.nf'
include { DMR_METHYLKIT }                   from '../../modules/methylation/dmr_methylkit.nf'
include { ANNOTATE_REGIONS }                from '../../modules/methylation/annotate_regions.nf'
include { STANDARDIZE_OUTPUT }              from '../../modules/methylation/standardize_output.nf'
include { EMIT_MANIFEST }                   from '../../modules/methylation/emit_manifest.nf'
include { RENDER_REPORT }                   from '../../modules/methylation/render_report.nf'

def paramsAsJson() {
    // Explicit, auditable parameter dump for the manifest -- not "defaults",
    // per docs/design.md Sec. 4 provenance model.
    def dump = [
        mode                 : params.mode,
        study_id             : params.methylation?.study_id,
        comparison_id        : params.methylation?.comparison_id,
        min_coverage         : params.methylation?.min_coverage,
        qvalue_cutoff        : params.methylation?.qvalue_cutoff,
        meth_diff_cutoff     : params.methylation?.meth_diff_cutoff,
        dmr_mode             : params.methylation?.dmr_mode,
        tile_size            : params.methylation?.tile_size,
        tile_step            : params.methylation?.tile_step,
        promoter_upstream_bp : params.methylation?.promoter_upstream_bp,
        genome_assembly      : params.genome_assembly,
    ]
    return groovy.json.JsonOutput.toJson(dump)
}

workflow {

    def study_id      = params.methylation.study_id
    def comparison_id = params.methylation.comparison_id
    def run_start_iso = workflow.start.toString()
    def report_template = file("${projectDir}/assets/report_template.qmd")

    if (params.mode == 'raw_reanalysis') {

        // ---- Raw-data reanalysis branch ------------------------------------
        // Expects paired-end FASTQ, a genome FASTA, and a GTF annotation.
        reads_ch = Channel
            .fromFilePairs(params.reads, checkIfExists: true)

        genome_fasta_ch = Channel.fromPath(params.genome_fasta, checkIfExists: true)
        annotation_gtf_ch = Channel.fromPath(params.annotation_gtf, checkIfExists: true)

        FASTQC(reads_ch)
        TRIM_GALORE(reads_ch)

        BISMARK_GENOME_PREPARATION(genome_fasta_ch)

        aligned_ch = BISMARK_ALIGN(
            TRIM_GALORE.out.trimmed_reads,
            BISMARK_GENOME_PREPARATION.out.genome_dir.first()
        )

        deduped_ch = BISMARK_DEDUPLICATE(aligned_ch.bam)

        extracted_ch = BISMARK_METHYLATION_EXTRACTOR(
            deduped_ch.dedup_bam,
            BISMARK_GENOME_PREPARATION.out.genome_dir.first()
        )

        filtered_ch = COVERAGE_FILTER(
            extracted_ch.cx_report,
            params.methylation.min_coverage
        )

        // Collect all per-sample filtered CX reports into ordered lists so
        // DMR_METHYLKIT can build methylKit's sample.id/treatment vectors.
        // sample_sheet maps sample_id -> treatment (0=control, 1=treatment),
        // supplied via params.methylation.sample_sheet (a simple CSV:
        // sample_id,treatment).
        sample_sheet_ch = Channel
            .fromPath(params.methylation.sample_sheet, checkIfExists: true)
            .splitCsv(header: true)
            .map { row -> tuple(row.sample_id, row.treatment as Integer) }

        dmr_inputs_ch = filtered_ch.filtered_report
            .join(sample_sheet_ch)
            .toSortedList { a, b -> a[0] <=> b[0] }
            .map { rows ->
                def sample_ids = rows.collect { it[0] }
                def cx_reports = rows.collect { it[1] }
                def treatments = rows.collect { it[2] }
                tuple(sample_ids, cx_reports, treatments)
            }

        dmr_ch = dmr_inputs_ch.map { sample_ids, cx_reports, treatments ->
            tuple(study_id, comparison_id, cx_reports, sample_ids, treatments)
        }

        dmr_result_ch = DMR_METHYLKIT(
            dmr_ch.map { it[0] },
            dmr_ch.map { it[1] },
            dmr_ch.map { it[2] },
            dmr_ch.map { it[3] },
            dmr_ch.map { it[4] },
            params.methylation.min_coverage,
            params.methylation.qvalue_cutoff,
            params.methylation.meth_diff_cutoff,
            params.methylation.dmr_mode,
            params.methylation.tile_size,
            params.methylation.tile_step
        )

        annotated_ch = ANNOTATE_REGIONS(
            study_id,
            comparison_id,
            dmr_result_ch.dmr_table,
            annotation_gtf_ch,
            params.methylation.promoter_upstream_bp
        )

        standardized_ch = STANDARDIZE_OUTPUT(
            study_id,
            comparison_id,
            annotated_ch.annotated_table,
            params.mode
        )

        qc_files_ch = dmr_result_ch.qc_metrics
            .concat(annotated_ch.qc_metrics)
            .concat(standardized_ch.qc_metrics)
            .collect()

        warnings_ch = Channel.value([])

        manifest_ch = EMIT_MANIFEST(
            study_id,
            comparison_id,
            params.mode,
            standardized_ch.standardized_table,
            qc_files_ch,
            paramsAsJson(),
            warnings_ch,
            workflow.manifest.version ?: '0.1.0-scaffold',
            run_start_iso
        )

        RENDER_REPORT(
            study_id,
            comparison_id,
            standardized_ch.standardized_table,
            manifest_ch.manifest,
            report_template,
            params.mode
        )

    } else if (params.mode == 'processed_results_harmonization') {

        // ---- Processed-results harmonization branch ------------------------
        // Expects a single DMR-shaped TSV; raw QC/coverage are not
        // independently verifiable in this mode (flagged explicitly below).
        processed_ch = Channel.fromPath(params.methylation.processed_results, checkIfExists: true)

        standardized_ch = STANDARDIZE_OUTPUT(
            study_id,
            comparison_id,
            processed_ch,
            params.mode
        )

        qc_files_ch = standardized_ch.qc_metrics.collect()

        warnings_ch = Channel.value([
            'raw_data_not_available: this study was processed in ' +
            'processed_results_harmonization mode; read QC, alignment, ' +
            'coverage, and DML/DMR-calling steps were performed by the ' +
            'original study authors and were not re-run or independently ' +
            'verified by AREE.'
        ])

        manifest_ch = EMIT_MANIFEST(
            study_id,
            comparison_id,
            params.mode,
            standardized_ch.standardized_table,
            qc_files_ch,
            paramsAsJson(),
            warnings_ch,
            workflow.manifest.version ?: '0.1.0-scaffold',
            run_start_iso
        )

        RENDER_REPORT(
            study_id,
            comparison_id,
            standardized_ch.standardized_table,
            manifest_ch.manifest,
            report_template,
            params.mode
        )

    } else {
        error "Unknown params.mode '${params.mode}': expected 'raw_reanalysis' or 'processed_results_harmonization'"
    }
}

workflow.onComplete {
    log.info """
    AREE methylation workflow (${params.mode}) finished.
    Study:      ${params.methylation?.study_id}
    Comparison: ${params.methylation?.comparison_id}
    Output dir: ${params.outdir}/methylation
    Success:    ${workflow.success}
    NOTE: this is a scaffold workflow; see workflows/methylation/README.md
    for what has and has not actually been executed/verified.
    """
}
