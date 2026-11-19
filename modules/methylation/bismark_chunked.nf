// Chunked trimming, alignment and deduplication (raw_reanalysis mode,
// opt-in through params.methylation.align_chunk_pairs).
//
// Deep libraries take a day or more to align as one task. On a preemptible
// scheduler partition (e.g. Hyak's ckpt, where a CPU job may run at most a
// few hours and can be killed and restarted at any time) such a task never
// finishes. Here each library is aligned as fixed-size chunks of read pairs,
// so every task is short, a preempted chunk is retried on its own, and the
// chunks of one library align in parallel.
//
// The second design goal is a small footprint on shared storage. Every
// intermediate about the size of the reads (chunk FASTQs, trimmed reads, the
// merged BAM, Bismark's temporary files) lives only in the task's node-local
// scratch directory (`scratch true`); only the chunk BAMs and the
// deduplicated BAM are written back. The first chunked version persisted
// trimmed reads, chunk FASTQs and the merged BAM: ~4.4x the raw reads on
// shared disk, against ~1.8x here.
//
//   COUNT_READ_PAIRS      raw R1 -> number of read pairs
//   TRIM_ALIGN_CHUNK      raw R1/R2 records [first, last] -> Trim Galore ->
//                         Bismark, on node-local disk -> chunk BAM
//   MERGE_DEDUP_CHUNKS    samtools cat of the chunk BAMs in chunk order (mates
//                         stay adjacent) -> deduplicate_bismark, on node-local
//                         disk -> deduplicated BAM; chunk reports summed into
//                         one per-library Bismark PE report
//
// TRIM_ALIGN_CHUNK needs Trim Galore and Bismark in one environment. No
// public image has both, so it declares a pinned conda environment, which
// -profile klone enables. Natively (-profile local) the tools on PATH are
// used, as for the other steps.

process COUNT_READ_PAIRS {
    tag "${sample_id}"
    label 'process_low'
    container 'quay.io/biocontainers/bismark:0.24.2--hdfd78af_0'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path(reads), stdout, emit: counted

    script:
    """
    gzip -dc ${reads[0]} | awk 'END { printf "%d", NR / 4 }'
    """
}

process TRIM_ALIGN_CHUNK {
    tag "${sample_id}:${chunk_id}"
    label 'process_high'
    conda 'bioconda::bismark=0.24.2 bioconda::trim-galore=0.6.10 bioconda::samtools=1.18 conda-forge::pigz=2.8'
    scratch true

    input:
    tuple val(sample_id), val(chunk_id), val(first_pair), val(last_pair), path(reads)
    path genome_dir

    output:
    tuple val(sample_id), val(chunk_id), path("${sample_id}.${chunk_id}_bismark_bt2_pe.bam"), path("${sample_id}.${chunk_id}_bismark_bt2_PE_report.txt"), path("${sample_id}.${chunk_id}*trimming_report.txt"), emit: chunk_bam
    path "versions.yml", emit: versions

    script:
    // Records [first_pair, last_pair] (1-based) of both mates; the two files
    // list pairs in the same order, so the mates stay in step. awk stops
    // reading after the last record, which ends gzip with SIGPIPE: the
    // subshell's `|| true` keeps that from failing the task under pipefail.
    // Trim Galore and Bismark then run as in TRIM_GALORE and BISMARK_ALIGN.
    def first_line = ((first_pair as long) - 1) * 4 + 1
    def last_line  = (last_pair as long) * 4
    def n_pairs    = (last_pair as long) - (first_pair as long) + 1
    def p          = "${sample_id}.${chunk_id}"
    """
    extract() {
        ( gzip -dc "\$1" || true ) | awk -v s=${first_line} -v e=${last_line} 'NR > e { exit } NR >= s' \\
            | pigz -1 -p 2 > "\$2"
    }
    extract ${reads[0]} ${p}_R1.fq.gz
    extract ${reads[1]} ${p}_R2.fq.gz
    n1=\$(( \$(gzip -dc ${p}_R1.fq.gz | wc -l) / 4 ))
    n2=\$(( \$(gzip -dc ${p}_R2.fq.gz | wc -l) / 4 ))
    [ "\$n1" -eq ${n_pairs} ] && [ "\$n2" -eq ${n_pairs} ] \\
        || { echo "chunk has \$n1 / \$n2 pairs, expected ${n_pairs}" >&2; exit 1; }

    trim_galore \\
        --paired \\
        --cores ${Math.max((task.cpus / 3) as int, 1)} \\
        --basename ${p} \\
        --output_dir . \\
        ${p}_R1.fq.gz ${p}_R2.fq.gz
    rm ${p}_R1.fq.gz ${p}_R2.fq.gz

    bismark \\
        --genome ${genome_dir} \\
        -1 ${p}_val_1.fq.gz \\
        -2 ${p}_val_2.fq.gz \\
        --bowtie2 \\
        --multicore ${Math.max((task.cpus / 3) as int, 1)} \\
        ${params.methylation?.bismark_align_args ?: ''} \\
        -o .
    rm ${p}_val_1.fq.gz ${p}_val_2.fq.gz

    mv *_bismark_bt2_pe.bam ${p}_bismark_bt2_pe.bam
    mv *_bismark_bt2_PE_report.txt ${p}_bismark_bt2_PE_report.txt

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        trim_galore: \$(trim_galore --version | grep version | sed 's/.*version //')
        cutadapt: \$(cutadapt --version)
        bismark: \$(bismark --version | grep 'Bismark Version' | sed 's/.*Bismark Version: v//')
    END_VERSIONS
    """
}

process MERGE_DEDUP_CHUNKS {
    tag "${sample_id}"
    label 'process_medium'
    container 'quay.io/biocontainers/bismark:0.24.2--hdfd78af_0'
    scratch true
    publishDir "${params.outdir}/methylation/bismark_dedup", mode: params.publish_mode, pattern: "*{deduplicated.bam,deduplication_report.txt}"
    publishDir "${params.outdir}/methylation/bismark_align", mode: params.publish_mode, pattern: "*{_PE_report.txt,_chunk_reports}"

    input:
    tuple val(sample_id), val(chunk_ids), path(bams), path(reports), path(trim_reports)

    output:
    tuple val(sample_id), path("${sample_id}_bismark_bt2_pe.deduplicated.bam"), emit: dedup_bam
    tuple val(sample_id), path("${sample_id}_bismark_bt2_pe.deduplication_report.txt"), emit: dedup_report
    tuple val(sample_id), path("${sample_id}_bismark_bt2_PE_report.txt"), emit: report
    path "${sample_id}_chunk_reports", emit: chunk_reports
    path "versions.yml", emit: versions

    script:
    // Chunk order is restored from the zero-padded chunk ids, not from the
    // order the tasks finished in. The merged BAM exists only in scratch.
    def order = [chunk_ids, bams].transpose().sort { it[0] }
    def bam_list = order.collect { it[1] }.join(' ')
    """
    samtools cat -o ${sample_id}_bismark_bt2_pe.bam ${bam_list}
    deduplicate_bismark --bam --paired ${sample_id}_bismark_bt2_pe.bam
    rm ${sample_id}_bismark_bt2_pe.bam

    mkdir ${sample_id}_chunk_reports
    cp ${reports.join(' ')} ${trim_reports.flatten().join(' ')} ${sample_id}_chunk_reports/
    sum_bismark_pe_reports.sh ${sample_id} ${sample_id}_chunk_reports/*_PE_report.txt \\
        > ${sample_id}_bismark_bt2_PE_report.txt

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        samtools: \$(samtools --version | head -n1 | sed 's/samtools //')
        bismark: \$(deduplicate_bismark --version | grep -oE 'Version: v[0-9.]+' | head -1 | sed 's/Version: v//')
    END_VERSIONS
    """
}
