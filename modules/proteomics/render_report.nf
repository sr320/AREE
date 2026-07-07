// Render the human-readable HTML QC/results report from the standardized
// output table + manifest, via Quarto.
//
// STATUS: structurally complete DSL2 process invoking a real Quarto render
// command against a real .qmd template (workflows/proteomics/assets/report_template.qmd,
// which contains executable R/ggplot2 chunks: missingness histogram, volcano
// plot, top-hits table). NOT executed in this build -- no Quarto/R runtime
// was invoked here; see workflows/proteomics/README.md.

process RENDER_REPORT {
    tag "${study_id}:${comparison_id}"
    label 'process_low'
    container 'ghcr.io/quarto-dev/quarto:1.5.57'
    publishDir "${params.outdir}/proteomics/report", mode: params.publish_mode

    input:
    val study_id
    val comparison_id
    path standardized_table
    path manifest_json
    path report_template

    output:
    path "${study_id}_${comparison_id}_proteomics_report.html", emit: html

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail

    cp "${report_template}" report_template.qmd

    quarto render report_template.qmd \\
      --to html \\
      --output "${study_id}_${comparison_id}_proteomics_report.html" \\
      -P standardized_table:"${standardized_table}" \\
      -P manifest_json:"${manifest_json}" \\
      -P study_id:"${study_id}" \\
      -P comparison_id:"${comparison_id}"
    """

    stub:
    """
    touch ${study_id}_${comparison_id}_proteomics_report.html
    """
}
