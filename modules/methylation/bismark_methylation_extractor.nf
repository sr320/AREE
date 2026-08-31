// Per-cytosine methylation calling from deduplicated Bismark alignments
// (raw_reanalysis mode only).
//
// STATUS: structurally complete DSL2 process with a real
// bismark_methylation_extractor command line. NOT executed in this build —
// see workflows/methylation/README.md. --cytosine_report requires the
// genome directory used for alignment so it can report every cytosine in
// the reference, not just covered ones.

process BISMARK_METHYLATION_EXTRACTOR {
    tag "${sample_id}"
    label 'process_high'
    container 'quay.io/biocontainers/bismark:0.24.2--hdfd78af_0'
    publishDir "${params.outdir}/methylation/methylation_calls", mode: params.publish_mode

    input:
    tuple val(sample_id), path(dedup_bam)
    path genome_dir

    output:
    tuple val(sample_id), path("*.CX_report.txt.gz"), emit: cx_report
    tuple val(sample_id), path("*splitting_report.txt"), emit: splitting_report
    tuple val(sample_id), path("*.bedGraph.gz"), emit: bedgraph
    path "versions.yml", emit: versions

    script:
    // Real methylation extraction with genome-wide cytosine report
    // (--CX reports CpG/CHG/CHH contexts; relevant for non-CpG methylation
    // signal sometimes reported in invertebrate WGBS studies). Not executed
    // here.
    """
    bismark_methylation_extractor \\
        --comprehensive \\
        --paired-end \\
        --cytosine_report \\
        --CX \\
        --genome_folder ${genome_dir} \\
        --bedGraph \\
        --gzip \\
        --multicore ${Math.max((task.cpus / 3) as int, 1)} \\
        ${dedup_bam}

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        bismark: \$(bismark_methylation_extractor --version | grep 'Bismark Version' | sed 's/.*Bismark Version: v//')
    END_VERSIONS
    """
}
