#!/usr/bin/env nextflow
/*
 * AREE :: RNA-seq reanalysis / harmonization workflow
 *
 * STATUS: structurally complete DSL2 scaffold. NOT executed against real
 * FASTQ or real processed tables in this build — no compute budget or real
 * public data was available in this environment. See README.md "What has
 * and has not been run" before treating any output as validated.
 *
 * Implements CLAUDE.md Layer 2 section A (RNA-seq):
 *   FASTQ QC -> adapter trimming -> pseudoalignment/quant -> sample QC ->
 *   differential expression -> effect-size table -> standardized output ->
 *   manifest + report.
 *
 * Two modes (docs/design.md section 7), selected by params.mode:
 *   - raw_reanalysis:                 FASTQ in, full pipeline runs.
 *   - processed_results_harmonization: an existing DE-like TSV is validated
 *     and reshaped directly; a warning is recorded that raw QC was not
 *     independently verified.
 *
 * Output in either mode: a standardized TSV with columns
 *   gene_id  baseMean  log2FoldChange  lfcSE  stat  pvalue  padj
 * which is exactly what `aree harmonize --study STUDY_ID --input <tsv>`
 * (src/harmonize/rnaseq.py::harmonize_rnaseq) expects to read.
 */

nextflow.enable.dsl = 2

include { FASTQC }                              from '../../modules/rnaseq/fastqc.nf'
include { TRIM_FASTP }                          from '../../modules/rnaseq/trim_fastp.nf'
include { SALMON_INDEX; SALMON_QUANT }          from '../../modules/rnaseq/align_quant_salmon.nf'
include { SAMPLE_QC }                           from '../../modules/rnaseq/sample_qc.nf'
include { DIFFERENTIAL_EXPRESSION_DESEQ2 }      from '../../modules/rnaseq/differential_expression_deseq2.nf'
include { STANDARDIZE_OUTPUT }                  from '../../modules/rnaseq/standardize_output.nf'
include { EMIT_MANIFEST }                       from '../../modules/rnaseq/emit_manifest.nf'
include { RENDER_REPORT }                       from '../../modules/rnaseq/render_report.nf'

def VALID_MODES = ['raw_reanalysis', 'processed_results_harmonization']

workflow {

    if (!(params.mode in VALID_MODES)) {
        error "params.mode must be one of ${VALID_MODES}, got: ${params.mode}"
    }

    log.info """
    AREE RNA-seq workflow
    ----------------------
    mode           : ${params.mode}
    study_id       : ${params.rnaseq?.study_id}
    comparison_id  : ${params.rnaseq?.comparison_id}
    outdir         : ${params.outdir}
    NOTE: unexecuted scaffold in this build (see README.md).
    """.stripIndent()

    study_id      = params.rnaseq?.study_id ?: params.study_id
    comparison_id = params.rnaseq?.comparison_id ?: params.comparison_id
    workflow_start_iso = workflow.start.toString()

    report_template = file("${projectDir}/assets/report_template.qmd")

    if (params.mode == 'raw_reanalysis') {

        // ---- raw_reanalysis branch --------------------------------------
        // Expect params.reads as a glob resolving paired-end FASTQ, e.g.:
        //   params.reads = "data/raw/GIGAS_HEAT01/*_R{1,2}.fastq.gz"
        if (!params.reads) {
            error "raw_reanalysis mode requires --reads (glob for fromFilePairs), e.g. 'reads/*_R{1,2}.fastq.gz'"
        }

        read_pairs_ch = Channel.fromFilePairs(params.reads, checkIfExists: true)

        FASTQC(read_pairs_ch)
        TRIM_FASTP(read_pairs_ch)

        // Build or reuse a Salmon index.
        if (params.salmon_index) {
            salmon_index_ch = Channel.fromPath(params.salmon_index, checkIfExists: true)
        } else {
            if (!params.transcript_fasta || !params.gtf) {
                error "raw_reanalysis mode without --salmon_index requires --transcript_fasta and --gtf"
            }
            SALMON_INDEX(
                Channel.fromPath(params.transcript_fasta, checkIfExists: true),
                Channel.fromPath(params.gtf, checkIfExists: true)
            )
            salmon_index_ch = SALMON_INDEX.out.index
        }

        gtf_ch = Channel.fromPath(params.gtf, checkIfExists: true)

        SALMON_QUANT(
            TRIM_FASTP.out.trimmed_reads,
            salmon_index_ch.first(),
            gtf_ch.first()
        )

        // Aggregate QC logs (FastQC zips + fastp jsons + salmon quant dirs)
        // into one MultiQC run.
        qc_inputs_ch = FASTQC.out.zip.map { sid, f -> f }
            .mix(TRIM_FASTP.out.json.map { sid, f -> f })
            .mix(SALMON_QUANT.out.quant_dir.map { sid, d -> d })
            .collect()

        SAMPLE_QC(qc_inputs_ch)

        if (!params.tx2gene) {
            error "raw_reanalysis mode requires --tx2gene (transcript-to-gene map TSV) for DESeq2 gene-level summarization"
        }
        if (!params.sample_sheet) {
            error "raw_reanalysis mode requires --sample_sheet (TSV with sample_id, condition, quant_subdir columns)"
        }

        quant_dirs_ch = SALMON_QUANT.out.quant_dir.map { sid, d -> d }.collect()

        DIFFERENTIAL_EXPRESSION_DESEQ2(
            study_id,
            comparison_id,
            quant_dirs_ch,
            Channel.fromPath(params.tx2gene, checkIfExists: true),
            Channel.fromPath(params.sample_sheet, checkIfExists: true)
        )

        raw_de_ch = DIFFERENTIAL_EXPRESSION_DESEQ2.out.results.map { sid, cid, tsv -> tsv }
        workflow_warnings = []

        STANDARDIZE_OUTPUT(study_id, comparison_id, raw_de_ch)

    } else {

        // ---- processed_results_harmonization branch ---------------------
        // Expect params.rnaseq.processed_results (see config/demo.config)
        // or params.processed_results as a fallback.
        processed_path = params.rnaseq?.processed_results ?: params.processed_results
        if (!processed_path) {
            error "processed_results_harmonization mode requires params.rnaseq.processed_results (or --processed_results) pointing at an existing DE-like TSV"
        }

        processed_ch = Channel.fromPath(processed_path, checkIfExists: true)

        workflow_warnings = [
            "raw_qc_not_independently_verified: this study was harmonized from " +
            "an existing processed results table; FastQC/fastp/Salmon QC metrics " +
            "were not generated or verified by this workflow run (see docs/design.md section 7)."
        ]

        STANDARDIZE_OUTPUT(study_id, comparison_id, processed_ch)
    }

    // ---- shared tail: manifest + report (both modes) --------------------
    standardized_ch = STANDARDIZE_OUTPUT.out.standardized_tsv.map { sid, cid, tsv -> tsv }

    // Groovy list -> JSON array literal string, safely embeddable into the
    // Python source generated inside EMIT_MANIFEST's script block.
    warnings_json = groovy.json.JsonOutput.toJson(workflow_warnings)

    EMIT_MANIFEST(
        study_id,
        comparison_id,
        params.mode,
        standardized_ch,
        workflow_start_iso,
        warnings_json
    )

    RENDER_REPORT(
        study_id,
        comparison_id,
        params.mode,
        standardized_ch,
        EMIT_MANIFEST.out.manifest,
        report_template
    )
}

workflow.onComplete {
    log.info """
    AREE RNA-seq workflow finished
    -------------------------------
    mode      : ${params.mode}
    status    : ${workflow.success ? 'OK' : 'FAILED'}
    outdir    : ${params.outdir}
    duration  : ${workflow.duration}
    NOTE: this run reflects DSL2 wiring/validation only in this build unless
    executed against real containers and real inputs outside this environment.
    """.stripIndent()
}
