// FASTQ quality control for bisulfite/EM-seq reads (raw_reanalysis mode only).
//
// STATUS: structurally complete DSL2 process with a real FastQC command line.
// NOT executed in this build — no real/synthetic bisulfite FASTQ fixtures
// ship in this repository (see workflows/methylation/README.md and
// config/demo.config header). FastQC does not need to know reads are
// bisulfite-converted; the same invocation used for RNA-seq applies here.

process FASTQC {
    tag "${sample_id}"
    label 'process_low'
    container 'biocontainers/fastqc:0.12.1--hdfd78af_0'
    publishDir "${params.outdir}/methylation/fastqc", mode: params.publish_mode

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("*_fastqc.zip"), emit: zip
    tuple val(sample_id), path("*_fastqc.html"), emit: html
    path "versions.yml", emit: versions

    script:
    // Real invocation for paired-end bisulfite/EM-seq reads. Not executed
    // here — see README "What has and has not been run".
    """
    fastqc -o . -t ${task.cpus} ${reads}

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        fastqc: \$(fastqc --version | sed 's/FastQC v//')
    END_VERSIONS
    """
}
