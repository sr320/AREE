// FASTQ quality control (raw_reanalysis mode only).
//
// STATUS: structurally complete DSL2 process with a real FastQC command line.
// NOT executed in this build — no real/synthetic FASTQ fixtures ship in this
// repository (see workflows/rnaseq/README.md and config/demo.config header).

process FASTQC {
    tag "${sample_id}"
    label 'process_low'
    container 'quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0'
    publishDir "${params.outdir}/rnaseq/fastqc", mode: params.publish_mode

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("*_fastqc.zip"), emit: zip
    tuple val(sample_id), path("*_fastqc.html"), emit: html
    path "versions.yml", emit: versions

    script:
    // Real invocation a bioinformatician would run against paired-end reads.
    // Not executed here — see README "What has and has not been run".
    """
    fastqc -o . -t ${task.cpus} ${reads}

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        fastqc: \$(fastqc --version | sed 's/FastQC v//')
    END_VERSIONS
    """

    stub:
    """
    touch ${sample_id}_R1_fastqc.zip ${sample_id}_R2_fastqc.zip
    touch ${sample_id}_R1_fastqc.html ${sample_id}_R2_fastqc.html
    echo '${task.process}:' > versions.yml
    """
}
