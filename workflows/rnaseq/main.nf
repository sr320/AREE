#!/usr/bin/env nextflow
/*
 * AREE :: RNA-seq reanalysis / harmonization workflow
 *
 * STATUS: executed end to end against subsampled real FASTQ and against the
 * bundled processed-results demo. CI also exercises the raw DAG with tiny
 * paired FASTQ fixtures under `-stub-run`. Containers and full-depth biology
 * remain unverified; see README.md "What has and has not been run".
 *
 * Implements CLAUDE.md Layer 2 section A (RNA-seq):
 *   FASTQ QC -> adapter trimming -> pseudoalignment/quant -> sample QC ->
 *   differential expression -> effect-size table -> standardized output ->
 *   manifest + report. Optionally (--viral_reference), per-sample viral read
 *   QC, and curated sample exclusion before DE (--exclude_samples).
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
include { VIRAL_INDEX; VIRAL_QUANT; VIRAL_SUMMARY } from '../../modules/rnaseq/viral_load.nf'

// Nextflow's strict DSL2 syntax (25.10+) forbids bare statements at script
// level, so this is a function declaration rather than a top-level assignment.
def validModes() {
    ['raw_reanalysis', 'processed_results_harmonization']
}

// Nextflow saves each finished task's cache entry with Kryo, which needs
// java.base/java.util opened to it. The `nextflow` launcher script passes
// --add-opens for this; `java -jar ~/.nextflow/framework/<ver>/nextflow-<ver>-one.jar`
// does not. Without it every cache write fails silently inside Nextflow's
// writer thread: the run succeeds, but -resume never reuses a task. This is
// what made the 2026-10-12 DUPERRET2025_STARVE_RNASEQ comparisons each
// recompute every library (see README.md, "Several comparisons of one study").
def cacheEntriesSerializable() {
    Collections.getModule().isOpen('java.util', nextflow.Session.getModule())
}

workflow {

    if (!(params.mode in validModes())) {
        error "params.mode must be one of ${validModes()}, got: ${params.mode}"
    }

    if (!cacheEntriesSerializable()) {
        def why = "this JVM does not open java.base/java.util to Nextflow, so no task " +
            "will be written to the resume cache. Launch Nextflow with the `nextflow` " +
            "launcher script (e.g. from `conda activate aree`), not with " +
            "`java -jar nextflow-<ver>-one.jar`; see workflows/rnaseq/README.md, " +
            "\"Several comparisons of one study\"."
        if (workflow.resume) {
            error "-resume cannot reuse any task: ${why}"
        }
        log.warn "A later -resume will not reuse this run: ${why}"
    }

    log.info """
    AREE RNA-seq workflow
    ----------------------
    mode           : ${params.mode}
    study_id       : ${params.rnaseq?.study_id ?: params.study_id}
    comparison_id  : ${params.rnaseq?.comparison_id ?: params.comparison_id}
    outdir         : ${params.outdir}
    """.stripIndent()

    study_id      = params.rnaseq?.study_id ?: params.study_id
    comparison_id = params.rnaseq?.comparison_id ?: params.comparison_id
    workflow_start_iso = workflow.start.toString()

    report_template = file("${projectDir}/assets/report_template.qmd")
    no_viral_qc     = file("${projectDir}/assets/NO_VIRAL_QC")
    no_deseq2_design = file("${projectDir}/assets/NO_DESEQ2_DESIGN")
    collect_versions_py = file("${projectDir}/../../modules/shared/collect_versions.py", checkIfExists: true)

    if (params.mode == 'raw_reanalysis') {

        // ---- raw_reanalysis branch --------------------------------------
        // Expect params.reads as a glob resolving paired-end FASTQ, e.g.:
        //   params.reads = "data/raw/GIGAS_HEAT01/*_R{1,2}.fastq.gz"
        // or, with --single_end, one FASTQ per library, e.g.:
        //   params.reads = "data/raw/DUPERRET2025_STARVE_RNASEQ/ERR*.fastq.gz"
        // Single-end sample ids are the file name without .fastq/.fq(.gz), so
        // ENA downloads are keyed by run accession, which the DESeq2 step
        // matches against the sample sheet's run_accession column.
        if (!params.reads) {
            error "raw_reanalysis mode requires --reads (glob for fromFilePairs, or one file per library with --single_end), e.g. 'reads/*_R{1,2}.fastq.gz'"
        }

        if (params.single_end) {
            read_pairs_ch = Channel.fromPath(params.reads, checkIfExists: true)
                .map { f -> [ f.name.replaceFirst(/\.(fastq|fq)(\.gz)?$/, ''), [ f ] ] }
        } else {
            read_pairs_ch = Channel.fromFilePairs(params.reads, checkIfExists: true)
        }

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
        // Each process's versions.yml, merged by EMIT_MANIFEST into the
        // manifest's software_versions. A process added to this DAG must add
        // its .out.versions here too (tests/test_workflow_versions.py checks).
        versions_ch = FASTQC.out.versions
            .mix(TRIM_FASTP.out.versions)
        if (!params.salmon_index) {
            versions_ch = versions_ch.mix(SALMON_INDEX.out.versions)
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
        versions_ch = versions_ch
            .mix(SALMON_QUANT.out.versions)
            .mix(SAMPLE_QC.out.versions)

        if (!params.tx2gene) {
            error "raw_reanalysis mode requires --tx2gene (transcript-to-gene map TSV) for DESeq2 gene-level summarization"
        }
        if (!params.sample_sheet) {
            error "raw_reanalysis mode requires --sample_sheet (TSV with sample_id and condition columns; quant_subdir optional)"
        }
        // The two condition levels are named explicitly rather than assumed to
        // be "control"/"treatment": a real sample sheet carries the study's own
        // group labels, and a whole-BioProject sheet may describe groups this
        // contrast does not use.
        if (!params.control_level || !params.treatment_level) {
            error "raw_reanalysis mode requires --control_level and --treatment_level naming the two `condition` values to contrast"
        }

        quant_dirs_ch = SALMON_QUANT.out.quant_dir.map { sid, d -> d }.collect()

        DIFFERENTIAL_EXPRESSION_DESEQ2(
            study_id,
            comparison_id,
            params.control_level,
            params.treatment_level,
            params.exclude_samples ?: '',
            params.replicate_unit ?: '',
            params.covariates ?: '',
            params.tag_seq as boolean,
            quant_dirs_ch,
            Channel.fromPath(params.tx2gene, checkIfExists: true),
            Channel.fromPath(params.sample_sheet, checkIfExists: true)
        )

        raw_de_ch = DIFFERENTIAL_EXPRESSION_DESEQ2.out.results.map { sid, cid, tsv -> tsv }
        deseq2_design_ch = DIFFERENTIAL_EXPRESSION_DESEQ2.out.design
        versions_ch = versions_ch.mix(DIFFERENTIAL_EXPRESSION_DESEQ2.out.versions)
        workflow_warnings = []

        // Optional viral read QC: did every challenged animal carry the
        // pathogen, and is every control free of it?
        if (params.viral_reference) {
            viral_fasta_ch = Channel.fromPath(params.viral_reference, checkIfExists: true).first()
            VIRAL_INDEX(viral_fasta_ch)
            VIRAL_QUANT(TRIM_FASTP.out.trimmed_reads, VIRAL_INDEX.out.index.first())
            VIRAL_SUMMARY(
                study_id,
                comparison_id,
                params.control_level,
                params.treatment_level,
                params.exclude_samples ?: '',
                params.viral_detection_per_million,
                viral_fasta_ch,
                Channel.fromPath(params.sample_sheet, checkIfExists: true).first(),
                VIRAL_QUANT.out.counts.collect()
            )
            viral_summary_ch = VIRAL_SUMMARY.out.json
            versions_ch = versions_ch
                .mix(VIRAL_INDEX.out.versions)
                .mix(VIRAL_QUANT.out.versions)
                .mix(VIRAL_SUMMARY.out.versions)
        } else {
            viral_summary_ch = Channel.value(no_viral_qc)
        }

        STANDARDIZE_OUTPUT(study_id, comparison_id, raw_de_ch)

        // The reference files this run read, checksummed into the manifest. A
        // bundle fetched with `aree fetch-reference` carries
        // reference_provenance.json beside them; it is staged too, so the
        // manifest can confirm the files still match what was fetched.
        def ref_paths = [params.transcript_fasta, params.gtf, params.tx2gene]
            .findAll { it }
            .collect { file(it, checkIfExists: true) }
        if (ref_paths) {
            def bundle_provenance = ref_paths[0].parent.resolve('reference_provenance.json')
            if (bundle_provenance.exists()) {
                ref_paths << bundle_provenance
            }
        }
        reference_files_ch = Channel.value(ref_paths)

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
        viral_summary_ch = Channel.value(no_viral_qc)
        deseq2_design_ch = Channel.value(no_deseq2_design)
        versions_ch = Channel.empty()
        reference_files_ch = Channel.value([])
    }

    // ---- shared tail: manifest + report (both modes) --------------------
    standardized_ch = STANDARDIZE_OUTPUT.out.standardized_tsv.map { sid, cid, tsv -> tsv }

    // Groovy list -> JSON array literal string, safely embeddable into the
    // Python source generated inside EMIT_MANIFEST's script block.
    warnings_json = groovy.json.JsonOutput.toJson(workflow_warnings)

    versions_ch = versions_ch.mix(STANDARDIZE_OUTPUT.out.versions)

    EMIT_MANIFEST(
        study_id,
        comparison_id,
        params.mode,
        standardized_ch,
        workflow_start_iso,
        warnings_json,
        viral_summary_ch,
        deseq2_design_ch,
        reference_files_ch,
        versions_ch.collect(),
        collect_versions_py
    )

    RENDER_REPORT(
        study_id,
        comparison_id,
        params.mode,
        standardized_ch,
        EMIT_MANIFEST.out.manifest,
        report_template
    )

    // Neither `params` nor `workflow` resolves inside a closure nested in the
    // workflow body, so snapshot both into locals the handler can capture.
    def _params = params
    def _wf = workflow
    // Moved inside the entry workflow: strict DSL2 (Nextflow 25.10+) rejects
    // `workflow.onComplete` as a script-level statement.
    workflow.onComplete {
        log.info """
        AREE RNA-seq workflow finished
        -------------------------------
        mode      : ${_params.mode}
        status    : ${_wf.success ? 'OK' : 'FAILED'}
        outdir    : ${_params.outdir}
        duration  : ${_wf.duration}
        NOTE: interpret biological output only when the run used full-depth,
        quality-controlled inputs; stub runs validate wiring only.
        """.stripIndent()
    }
}
