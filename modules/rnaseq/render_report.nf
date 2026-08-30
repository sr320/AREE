// Render the human-readable HTML QC/DE report via Quarto, per CLAUDE.md
// Layer 2 requirement ("a human-readable HTML or Quarto report").
//
// STATUS: structurally complete DSL2 process with a real `quarto render`
// invocation against a real parameterized .qmd template
// (workflows/rnaseq/assets/report_template.qmd). Not executed in this
// build — no Quarto/R runtime invoked here.

process RENDER_REPORT {
    tag "${study_id}:${comparison_id}"
    label 'process_low'
    container 'ghcr.io/quarto-dev/quarto:1.5.57'
    publishDir "${params.outdir}/rnaseq/reports", mode: params.publish_mode

    input:
    val study_id
    val comparison_id
    val mode
    path standardized_tsv
    path manifest_json
    path report_template

    output:
    path "${study_id}_${comparison_id}_rnaseq_report.html", emit: html
    path "versions.yml", emit: versions

    script:
    // Real Quarto CLI invocation with parameter injection. Not executed here.
    """
    quarto render ${report_template} \\
        --to html \\
        --output ${study_id}_${comparison_id}_rnaseq_report.html \\
        -P study_id:${study_id} \\
        -P comparison_id:${comparison_id} \\
        -P mode:${mode} \\
        -P standardized_tsv:${standardized_tsv} \\
        -P manifest_json:${manifest_json}

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        quarto: \$(quarto --version)
    END_VERSIONS
    """

    stub:
    """
    printf '<html><body>AREE RNA-seq stub report</body></html>\n' > ${study_id}_${comparison_id}_rnaseq_report.html
    echo '${task.process}:' > versions.yml
    """
}
