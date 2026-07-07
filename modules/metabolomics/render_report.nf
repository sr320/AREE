// Render the human-readable HTML report from the standardized metabolite
// feature table + provenance manifest, using the Quarto document at
// workflows/metabolomics/assets/report_template.qmd.
//
// STATUS: real Quarto/R Markdown template with real ggplot2/dplyr/knitr
// logic (would render correctly given a real Quarto + R + package
// environment and real input files). NOT executed in this build — no
// container runtime invoked here.

process RENDER_REPORT {
    tag "${study_id}:${comparison_id}"
    label 'process_low'
    container 'ghcr.io/quarto-dev/quarto:1.5.57'
    publishDir "${params.outdir}/metabolomics/report", mode: params.publish_mode

    input:
    val study_id
    val comparison_id
    val mode
    path standardized_tsv
    path manifest_json
    path report_template

    output:
    path "${study_id}_${comparison_id}_metabolomics_report.html", emit: html

    script:
    """
    cp ${report_template} report.qmd

    quarto render report.qmd \\
        --to html \\
        --output ${study_id}_${comparison_id}_metabolomics_report.html \\
        -P study_id:"${study_id}" \\
        -P comparison_id:"${comparison_id}" \\
        -P mode:"${mode}" \\
        -P standardized_tsv:"${standardized_tsv}" \\
        -P manifest_json:"${manifest_json}"
    """

    stub:
    """
    touch ${study_id}_${comparison_id}_metabolomics_report.html
    """
}
