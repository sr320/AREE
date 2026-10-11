// Adapter trimming / read QC filtering with fastp (raw_reanalysis mode only).
//
// STATUS: structurally complete DSL2 process with a real fastp command line.
// NOT executed in this build — see workflows/rnaseq/README.md.

process TRIM_FASTP {
    tag "${sample_id}"
    label 'process_medium'
    container 'quay.io/biocontainers/fastp:0.23.4--h5f740d0_0'
    // Trimmed FASTQ are tens of GB per study and are regenerable from the raw
    // reads, so they are only ever symlinked into outdir; `--publish_mode copy`
    // applies to the small fastp reports alone. Copying them once filled a
    // workstation disk.
    publishDir "${params.outdir}/rnaseq/trimmed", mode: 'symlink', pattern: '*.fastq.gz'
    publishDir "${params.outdir}/rnaseq/trimmed", mode: params.publish_mode, pattern: '*.fastp.{json,html}'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_trimmed_R{1,2}.fastq.gz"), emit: trimmed_reads
    tuple val(sample_id), path("${sample_id}.fastp.json"), emit: json
    tuple val(sample_id), path("${sample_id}.fastp.html"), emit: html
    path "versions.yml", emit: versions

    script:
    // Paired-end by default. --single_end takes one FASTQ per sample;
    // --tag_seq (3' QuantSeq-style libraries) also trims poly-A/poly-X tails,
    // which those libraries carry at the 3' end of most reads.
    def tag_seq_opts = params.tag_seq ? '--trim_poly_x' : ''
    if (params.single_end) {
    """
    fastp \\
        -i ${reads[0]} \\
        -o ${sample_id}_trimmed_R1.fastq.gz \\
        ${tag_seq_opts} \\
        --thread ${task.cpus} \\
        --json ${sample_id}.fastp.json \\
        --html ${sample_id}.fastp.html \\
        --report_title "${sample_id} fastp report"

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        fastp: \$(fastp --version 2>&1 | sed -e 's/fastp //g')
    END_VERSIONS
    """
    } else {
    """
    fastp \\
        -i ${reads[0]} -I ${reads[1]} \\
        -o ${sample_id}_trimmed_R1.fastq.gz -O ${sample_id}_trimmed_R2.fastq.gz \\
        --detect_adapter_for_pe \\
        ${tag_seq_opts} \\
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

    stub:
    """
    touch ${sample_id}_trimmed_R1.fastq.gz ${params.single_end ? '' : "${sample_id}_trimmed_R2.fastq.gz"}
    echo '{}' > ${sample_id}.fastp.json
    touch ${sample_id}.fastp.html
    echo '${task.process}:' > versions.yml
    """
}
