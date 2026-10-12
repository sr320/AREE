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
    // applies to the small fastp reports and versions.yml alone. Copying them
    // once filled a workstation disk.
    publishDir "${params.outdir}/rnaseq/trimmed", mode: 'symlink', pattern: '*.fastq.gz'
    publishDir "${params.outdir}/rnaseq/trimmed", mode: params.publish_mode, pattern: '{*.fastp.json,*.fastp.html,versions.yml}'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_trimmed_R{1,2}.fastq.gz"), emit: trimmed_reads
    tuple val(sample_id), path("${sample_id}.fastp.json"), emit: json
    tuple val(sample_id), path("${sample_id}.fastp.html"), emit: html
    path "versions.yml", emit: versions

    script:
    // Paired-end or single-end, decided by how many files the sample carries
    // (main.nf sends one with --single_end). A single file stages as a Path,
    // not a list, so normalize before indexing. Single-end output is R1 only.
    //
    // Single-end adapters are given explicitly (params.se_adapter), not
    // auto-detected. fastp's single-end detection guesses from overrepresented
    // sequence: on DUPERRET2025_STARVE_RNASEQ it picked an abundant oyster
    // transcript as the "adapter" and trimmed it from about 1% of reads, and
    // each library gets its own guess. Paired-end detection uses read overlap
    // and is reliable, so it is kept.
    def fq = reads instanceof Collection ? reads as List : [ reads ]
    def se_adapter = params.se_adapter == 'auto' ? '' : "--adapter_sequence ${params.se_adapter}"
    // 3' Tag-Seq reads run into the poly(A) tail; trim homopolymer tails so
    // they neither fail to map nor map to A-rich sequence.
    def poly_x = params.library_protocol == 'three_prime_tag' ? '--trim_poly_x' : ''
    def io = fq.size() == 1 \
        ? "-i ${fq[0]} -o ${sample_id}_trimmed_R1.fastq.gz ${se_adapter}" \
        : "-i ${fq[0]} -I ${fq[1]} -o ${sample_id}_trimmed_R1.fastq.gz -O ${sample_id}_trimmed_R2.fastq.gz --detect_adapter_for_pe"
    """
    fastp \\
        ${io} \\
        ${poly_x} \\
        --thread ${task.cpus} \\
        --json ${sample_id}.fastp.json \\
        --html ${sample_id}.fastp.html \\
        --report_title "${sample_id} fastp report"

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        fastp: \$(fastp --version 2>&1 | sed -e 's/fastp //g')
    END_VERSIONS
    """

    stub:
    def n = reads instanceof Collection ? (reads as List).size() : 1
    def trimmed = n == 1 ? "${sample_id}_trimmed_R1.fastq.gz" : "${sample_id}_trimmed_R1.fastq.gz ${sample_id}_trimmed_R2.fastq.gz"
    """
    touch ${trimmed}
    echo '{}' > ${sample_id}.fastp.json
    touch ${sample_id}.fastp.html
    echo '${task.process}:' > versions.yml
    """
}
