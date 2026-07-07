// Adapter/quality trimming with Trim Galore in bisulfite-aware mode
// (raw_reanalysis mode only).
//
// STATUS: structurally complete DSL2 process with a real Trim Galore command
// line for paired-end bisulfite/EM-seq data. NOT executed in this build —
// see workflows/methylation/README.md.

process TRIM_GALORE {
    tag "${sample_id}"
    label 'process_medium'
    container 'biocontainers/trim-galore:0.6.10--hdfd78af_0'
    publishDir "${params.outdir}/methylation/trimmed", mode: params.publish_mode

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_val_1.fq.gz"), path("${sample_id}_val_2.fq.gz"), emit: trimmed_reads
    tuple val(sample_id), path("*trimming_report.txt"), emit: report
    path "versions.yml", emit: versions

    script:
    // Real paired-end Trim Galore invocation. --paired triggers Cutadapt
    // paired-end trimming; no bisulfite-specific flag is required at this
    // step (bisulfite awareness is applied downstream by Bismark alignment).
    // Not executed here.
    """
    trim_galore \\
        --paired \\
        --cores ${task.cpus} \\
        --basename ${sample_id} \\
        --output_dir . \\
        ${reads[0]} ${reads[1]}

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        trim_galore: \$(trim_galore --version | grep version | sed 's/.*version //')
    END_VERSIONS
    """
}
