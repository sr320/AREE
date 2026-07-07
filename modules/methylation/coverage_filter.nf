// Coverage filtering of per-cytosine reports prior to DML/DMR calling
// (raw_reanalysis mode only).
//
// STATUS: structurally complete DSL2 process. The filtering logic itself is
// a real, minimal awk implementation of a coverage threshold over a Bismark
// CX/coverage report (columns: chrom, position, strand, count_methylated,
// count_unmethylated, context, trinucleotide). It is intentionally simple —
// production use would typically hand this off to methylKit's own
// filterByCoverage() (see dmr_methylkit.nf, which re-applies coverage
// filtering there as the authoritative step). This process exists mainly to
// (a) emit an early per-sample coverage QC metric before pooling samples for
// DMR calling, and (b) demonstrate that coverage filtering is a distinct,
// auditable pipeline stage per docs/design.md Sec. 4 (parameters must be
// explicit, not buried inside a downstream black box). NOT executed in this
// build.

process COVERAGE_FILTER {
    tag "${sample_id}"
    label 'process_medium'
    container 'biocontainers/bismark:0.24.2--hdfd78af_1'
    publishDir "${params.outdir}/methylation/coverage_filtered", mode: params.publish_mode

    input:
    tuple val(sample_id), path(cx_report)
    val min_coverage

    output:
    tuple val(sample_id), path("${sample_id}.cov_filtered.CX_report.txt.gz"), emit: filtered_report
    tuple val(sample_id), path("${sample_id}.coverage_qc.tsv"), emit: coverage_qc
    path "versions.yml", emit: versions

    script:
    // Bismark CX report columns: chrom pos strand count_meth count_unmeth
    // context trinucleotide. Coverage = count_meth + count_unmeth. Cytosines
    // below min_coverage are dropped from the filtered report; summary
    // counts are written to a small per-sample QC table consumed later by
    // emit_manifest.nf / render_report.nf.
    """
    zcat ${cx_report} | awk -v min_cov=${min_coverage} '
        BEGIN { OFS="\\t"; total=0; pass=0; cov_sum=0 }
        {
            cov = \$4 + \$5
            total++
            cov_sum += cov
            if (cov >= min_cov) {
                pass++
                print
            }
        }
        END {
            mean_cov = (total > 0) ? cov_sum/total : 0
            pct_pass = (total > 0) ? 100.0*pass/total : 0
            printf "sample_id\\tmin_coverage_threshold\\tn_cytosines_total\\tn_cytosines_passing\\tpct_passing\\tmean_coverage_unfiltered\\n" > "/dev/stderr"
            printf "${sample_id}\\t%d\\t%d\\t%d\\t%.2f\\t%.2f\\n", min_cov, total, pass, pct_pass, mean_cov > "/dev/stderr"
        }
    ' > ${sample_id}.cov_filtered.CX_report.txt 2> ${sample_id}.coverage_qc.tsv

    gzip -f ${sample_id}.cov_filtered.CX_report.txt

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        awk: \$(awk --version 2>&1 | head -n1 | sed 's/,.*//')
    END_VERSIONS
    """
}
