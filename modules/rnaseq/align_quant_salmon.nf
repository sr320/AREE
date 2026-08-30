// Pseudoalignment + transcript/gene-level quantification with Salmon
// (raw_reanalysis mode only).
//
// STATUS: structurally complete DSL2 process with real Salmon command lines
// (index build is optional/conditional on params.salmon_index being unset).
// NOT executed in this build — see workflows/rnaseq/README.md.

process SALMON_INDEX {
    tag "${genome_fasta.baseName}"
    label 'process_high'
    container 'combinelab/salmon:1.10.3'

    input:
    path genome_fasta
    path gtf

    output:
    path "salmon_index", emit: index
    path "versions.yml", emit: versions

    script:
    // Real decoy-aware-style index build would normally also fold in the
    // genome as a decoy; kept minimal here since this is unexecuted scaffold.
    """
    salmon index \\
        -t ${genome_fasta} \\
        -i salmon_index \\
        -k 31 \\
        -p ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        salmon: \$(salmon --version | sed 's/salmon //')
    END_VERSIONS
    """

    stub:
    """
    mkdir salmon_index
    touch salmon_index/versionInfo.json
    echo '${task.process}:' > versions.yml
    """
}

process SALMON_QUANT {
    tag "${sample_id}"
    label 'process_high'
    container 'combinelab/salmon:1.10.3'
    publishDir "${params.outdir}/rnaseq/salmon", mode: params.publish_mode

    input:
    tuple val(sample_id), path(reads)
    path salmon_index
    path gtf

    output:
    tuple val(sample_id), path("${sample_id}"), emit: quant_dir
    tuple val(sample_id), path("${sample_id}/quant.sf"), emit: quant_sf
    path "versions.yml", emit: versions

    script:
    // Real paired-end Salmon quant invocation with GTF-based gene-level
    // summarization via --geneMap. Not executed here.
    """
    salmon quant \\
        -i ${salmon_index} \\
        -l A \\
        -1 ${reads[0]} -2 ${reads[1]} \\
        -p ${task.cpus} \\
        --geneMap ${gtf} \\
        --validateMappings \\
        --gcBias \\
        -o ${sample_id}

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        salmon: \$(salmon --version | sed 's/salmon //')
    END_VERSIONS
    """

    stub:
    """
    mkdir ${sample_id}
    printf 'Name\tLength\tEffectiveLength\tTPM\tNumReads\nTX1\t100\t80\t1\t1\n' > ${sample_id}/quant.sf
    echo '${task.process}:' > versions.yml
    """
}
