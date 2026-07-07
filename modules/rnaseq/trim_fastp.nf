// Adapter trimming / read QC filtering with fastp (raw_reanalysis mode only).
//
// STATUS: structurally complete DSL2 process with a real fastp command line.
// NOT executed in this build — see workflows/rnaseq/README.md.

process TRIM_FASTP {
    tag "${sample_id}"
    label 'process_medium'
    container 'biocontainers/fastp:0.23.4--h5f740d0_0'
    publishDir "${params.outdir}/rnaseq/trimmed", mode: params.publish_mode

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_trimmed_R{1,2}.fastq.gz"), emit: trimmed_reads
    tuple val(sample_id), path("${sample_id}.fastp.json"), emit: json
    tuple val(sample_id), path("${sample_id}.fastp.html"), emit: html
    path "versions.yml", emit: versions

    script:
    // Real paired-end fastp invocation. Not executed here.
    """
    fastp \\
        -i ${reads[0]} -I ${reads[1]} \\
        -o ${sample_id}_trimmed_R1.fastq.gz -O ${sample_id}_trimmed_R2.fastq.gz \\
        --detect_adapter_for_pe \\
        --thread ${task.cpus} \\
        --json ${sample_id}.fastp.json \\
        --html ${sample_id}.fastp.html \\
        --report_title "${sample_id} fastp report"

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        fastp: \$(fastp --version 2>&1 | sed -e 's/fastp //g')
    END_VERSIONS
    """
}
