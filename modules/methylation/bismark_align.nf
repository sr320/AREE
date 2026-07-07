// Bismark genome preparation, alignment, and PCR-duplicate removal
// (raw_reanalysis mode only).
//
// STATUS: structurally complete DSL2 processes with real Bismark CLI
// invocations. NOT executed in this build — no reference genome build or
// bisulfite FASTQ fixtures are staged in this repository (see
// workflows/methylation/README.md). Genome preparation is split into its own
// process so it runs once per reference and is reused across all samples via
// Nextflow's implicit value-channel broadcast.

process BISMARK_GENOME_PREPARATION {
    label 'process_high'
    container 'biocontainers/bismark:0.24.2--hdfd78af_1'
    publishDir "${params.outdir}/methylation/bismark_genome", mode: params.publish_mode

    input:
    path genome_fasta

    output:
    path "bismark_genome", emit: genome_dir
    path "versions.yml", emit: versions

    script:
    // Real bismark_genome_preparation invocation. Bowtie2 is the default
    // aligner backend. Not executed here.
    """
    mkdir -p bismark_genome
    cp ${genome_fasta} bismark_genome/
    bismark_genome_preparation --bowtie2 bismark_genome

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        bismark: \$(bismark --version | grep 'Bismark Version' | sed 's/.*Bismark Version: v//')
    END_VERSIONS
    """
}

process BISMARK_ALIGN {
    tag "${sample_id}"
    label 'process_high'
    container 'biocontainers/bismark:0.24.2--hdfd78af_1'
    publishDir "${params.outdir}/methylation/bismark_align", mode: params.publish_mode

    input:
    tuple val(sample_id), path(read1), path(read2)
    path genome_dir

    output:
    tuple val(sample_id), path("${sample_id}_bismark_bt2_pe.bam"), emit: bam
    tuple val(sample_id), path("*_PE_report.txt"), emit: report
    path "versions.yml", emit: versions

    script:
    // Real paired-end Bismark/Bowtie2 alignment against a prepared genome
    // directory. Bismark names outputs from the first read file by default;
    // --basename pins a predictable sample-scoped output name. Not executed
    // here.
    """
    bismark \\
        --genome ${genome_dir} \\
        -1 ${read1} \\
        -2 ${read2} \\
        --bowtie2 \\
        --multicore ${Math.max((task.cpus / 3) as int, 1)} \\
        --basename ${sample_id}_bismark_bt2 \\
        -o .

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        bismark: \$(bismark --version | grep 'Bismark Version' | sed 's/.*Bismark Version: v//')
    END_VERSIONS
    """
}

process BISMARK_DEDUPLICATE {
    tag "${sample_id}"
    label 'process_medium'
    container 'biocontainers/bismark:0.24.2--hdfd78af_1'
    publishDir "${params.outdir}/methylation/bismark_dedup", mode: params.publish_mode

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}*.deduplicated.bam"), emit: dedup_bam
    tuple val(sample_id), path("*.deduplication_report.txt"), emit: report
    path "versions.yml", emit: versions

    script:
    // Real deduplicate_bismark invocation for paired-end BAMs. Not executed
    // here.
    """
    deduplicate_bismark --bam --paired ${bam}

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        bismark: \$(deduplicate_bismark --version | grep 'Bismark Version' | sed 's/.*Bismark Version: v//')
    END_VERSIONS
    """
}
