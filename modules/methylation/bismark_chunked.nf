// Chunked Bismark alignment (raw_reanalysis mode, opt-in through
// params.methylation.align_chunk_pairs).
//
// Deep libraries take a day or more to align as one task. On a preemptible
// scheduler partition (e.g. Hyak's ckpt, where a CPU job may run at most a
// few hours and can be killed and restarted at any time) such a task never
// finishes. Splitting each trimmed library into fixed-size chunks of read
// pairs makes every alignment task short, lets a preempted chunk be retried
// on its own, and lets the chunks of one library align in parallel.
//
//   SPLIT_FASTQ_PAIRS     trimmed R1/R2 -> N chunk pairs (seqkit split2)
//   BISMARK_ALIGN_CHUNK   one Bismark run per chunk pair (same command as
//                         BISMARK_ALIGN)
//   BISMARK_MERGE_CHUNKS  samtools cat of the chunk BAMs in chunk order, so
//                         mates stay adjacent for deduplicate_bismark; the
//                         chunk reports are summed into one PE report
//
// Deduplication runs on the merged BAM, so duplicates that fall in different
// chunks are still removed. The merged BAM is equivalent to an unchunked
// run's BAM up to record order, because Bismark aligns each pair
// independently.

process SPLIT_FASTQ_PAIRS {
    tag "${sample_id}"
    label 'process_medium'
    container 'quay.io/biocontainers/seqkit:2.8.2--h9ee0642_0'

    input:
    tuple val(sample_id), path(read1), path(read2)
    val chunk_pairs

    output:
    tuple val(sample_id), path("chunks/*_val_1.part_*.fq.gz"), path("chunks/*_val_2.part_*.fq.gz"), emit: chunks
    path "versions.yml", emit: versions

    script:
    """
    seqkit split2 \\
        -1 ${read1} \\
        -2 ${read2} \\
        --by-size ${chunk_pairs} \\
        --threads ${task.cpus} \\
        --out-dir chunks

    n1=\$(ls chunks/*_val_1.part_*.fq.gz | wc -l)
    n2=\$(ls chunks/*_val_2.part_*.fq.gz | wc -l)
    [ "\$n1" -eq "\$n2" ] || { echo "R1/R2 chunk counts differ: \$n1 vs \$n2" >&2; exit 1; }

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        seqkit: \$(seqkit version | sed 's/seqkit v//')
    END_VERSIONS
    """
}

process BISMARK_ALIGN_CHUNK {
    tag "${sample_id}:${chunk_id}"
    label 'process_high'
    container 'quay.io/biocontainers/bismark:0.24.2--hdfd78af_0'

    input:
    tuple val(sample_id), val(chunk_id), path(read1), path(read2)
    path genome_dir

    output:
    tuple val(sample_id), val(chunk_id), path("${sample_id}.${chunk_id}_bismark_bt2_pe.bam"), path("${sample_id}.${chunk_id}_bismark_bt2_PE_report.txt"), emit: chunk_bam
    path "versions.yml", emit: versions

    script:
    """
    bismark \\
        --genome ${genome_dir} \\
        -1 ${read1} \\
        -2 ${read2} \\
        --bowtie2 \\
        --multicore ${Math.max((task.cpus / 3) as int, 1)} \\
        ${params.methylation?.bismark_align_args ?: ''} \\
        -o .

    mv *_bismark_bt2_pe.bam ${sample_id}.${chunk_id}_bismark_bt2_pe.bam
    mv *_bismark_bt2_PE_report.txt ${sample_id}.${chunk_id}_bismark_bt2_PE_report.txt

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        bismark: \$(bismark --version | grep 'Bismark Version' | sed 's/.*Bismark Version: v//')
    END_VERSIONS
    """
}

process BISMARK_MERGE_CHUNKS {
    tag "${sample_id}"
    label 'process_low'
    container 'quay.io/biocontainers/bismark:0.24.2--hdfd78af_0'
    publishDir "${params.outdir}/methylation/bismark_align", mode: params.publish_mode

    input:
    tuple val(sample_id), val(chunk_ids), path(bams), path(reports)

    output:
    tuple val(sample_id), path("${sample_id}_bismark_bt2_pe.bam"), emit: bam
    tuple val(sample_id), path("${sample_id}_bismark_bt2_PE_report.txt"), emit: report
    path "${sample_id}_chunk_reports", emit: chunk_reports
    path "versions.yml", emit: versions

    script:
    // Chunk order is restored from the zero-padded chunk ids, not from the
    // order the tasks finished in.
    def order = [chunk_ids, bams, reports].transpose().sort { it[0] }
    def bam_list = order.collect { it[1] }.join(' ')
    """
    samtools cat -o ${sample_id}_bismark_bt2_pe.bam ${bam_list}

    mkdir ${sample_id}_chunk_reports
    cp ${reports.join(' ')} ${sample_id}_chunk_reports/
    sum_bismark_pe_reports.sh ${sample_id} ${sample_id}_chunk_reports/*_PE_report.txt \\
        > ${sample_id}_bismark_bt2_PE_report.txt

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        samtools: \$(samtools --version | head -n1 | sed 's/samtools //')
    END_VERSIONS
    """
}
