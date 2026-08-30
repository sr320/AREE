// Sample-level QC aggregation across FastQC/fastp/Salmon outputs via MultiQC
// (raw_reanalysis mode only).
//
// STATUS: structurally complete DSL2 process with a real MultiQC command
// line. NOT executed in this build — see workflows/rnaseq/README.md.

process SAMPLE_QC {
    tag "cohort"
    label 'process_low'
    container 'biocontainers/multiqc:1.21--pyhdfd78af_0'
    publishDir "${params.outdir}/rnaseq/multiqc", mode: params.publish_mode

    input:
    path qc_files, stageAs: 'inputs/*'

    output:
    path "multiqc_report.html", emit: report
    // MultiQC names its data directory after the report file, so `-n
    // multiqc_report.html` produces multiqc_report_data/. Declaring
    // "multiqc_data" made the task fail as a missing output even though
    // MultiQC had exited 0 and written everything.
    path "multiqc_report_data", emit: data
    path "versions.yml", emit: versions

    script:
    // Real MultiQC aggregation over whatever FastQC/fastp/Salmon logs were
    // staged into ./inputs. Not executed here.
    """
    multiqc inputs -n multiqc_report.html

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        multiqc: \$(multiqc --version | sed -e 's/multiqc, version //')
    END_VERSIONS
    """

    stub:
    """
    touch multiqc_report.html
    mkdir multiqc_report_data
    touch multiqc_report_data/multiqc_data.json
    echo '${task.process}:' > versions.yml
    """
}
