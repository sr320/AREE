// Render a human-readable HTML report from the standardized DMR table and
// run manifest, via Quarto.
//
// STATUS: structurally complete DSL2 process invoking a real `quarto render`
// command line against a real .qmd template
// (workflows/methylation/assets/report_template.qmd). NOT executed in this
// build — no Quarto installation was exercised as part of this task. The
// template itself is syntactically valid R/Quarto and would render given the
// standardized TSV + manifest JSON this workflow produces.

process RENDER_REPORT {
    tag "${study_id}:${comparison_id}"
    label 'process_low'
    container 'ghcr.io/quarto-dev/quarto:1.5.57'
    publishDir "${params.outdir}/methylation/report", mode: params.publish_mode

    input:
    val study_id
    val comparison_id
    path standardized_table
    path manifest_json
    path report_template
    val mode

    output:
    path "${study_id}_${comparison_id}_report.html", emit: html
    path "versions.yml", emit: versions

    script:
    // Quarto parameterized render: -P passes params through to the .qmd's
    // YAML `params:` block. --output pins a predictable, study-scoped
    // filename.
    """
    quarto render ${report_template} \\
        --to html \\
        --output ${study_id}_${comparison_id}_report.html \\
        -P study_id:${study_id} \\
        -P comparison_id:${comparison_id} \\
        -P standardized_table:${standardized_table} \\
        -P manifest_json:${manifest_json} \\
        -P mode:${mode}

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        quarto: \$(quarto --version)
    END_VERSIONS
    """
}
