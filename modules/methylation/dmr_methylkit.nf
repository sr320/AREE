// DML/DMR calling with methylKit (raw_reanalysis mode only).
//
// STATUS: structurally complete DSL2 process wrapping a real, syntactically
// correct methylKit R script (methRead -> filterByCoverage -> unite ->
// calculateDiffMeth -> getMethylDiff). This is the core scientific logic of
// the raw-mode methylation workflow. NOT executed in this build — no real
// per-sample CX reports exist to feed it (see workflows/methylation/README.md).
// Tile-based DMR calling (tileMethylCounts) is used rather than single-base
// DML calling by default, controlled by params.methylation.dmr_mode, because
// most published oyster WGBS resilience studies report region-level, not
// single-cytosine, differential methylation.

process DMR_METHYLKIT {
    tag "${study_id}:${comparison_id}"
    label 'process_high'
    container 'bioconductor/bioconductor_docker:RELEASE_3_18'
    publishDir "${params.outdir}/methylation/dmr", mode: params.publish_mode

    input:
    val study_id
    val comparison_id
    path cx_reports          // one filtered CX report per sample, all samples in one channel list
    val sample_ids            // matching list of sample_id strings
    val treatment_labels      // matching list of 0/1 treatment/control flags (methylKit "treatment" vector)
    val min_coverage
    val qvalue_cutoff
    val meth_diff_cutoff
    val dmr_mode              // "tile" or "base"
    val tile_size
    val tile_step

    output:
    path "${study_id}_${comparison_id}_dmr_raw.tsv", emit: dmr_table
    path "${study_id}_${comparison_id}_methylkit_qc.tsv", emit: qc_metrics
    path "versions.yml", emit: versions

    script:
    // sample_ids / treatment_labels / cx_reports are positionally matched by
    // the caller (see main.nf assembly of this process's inputs). The R
    // script below is a real methylKit DMR-calling script, simplified only
    // in that it assumes CpG context and a single pairwise comparison
    // (treatment vs control), consistent with the demo comparison design in
    // registry/studies/.
    """
    cat <<-'EOF_R' > run_methylkit.R
    library(methylKit)

    sample_ids  <- strsplit("${sample_ids.join(',')}", ",")[[1]]
    treatments  <- as.integer(strsplit("${treatment_labels.join(',')}", ",")[[1]])
    cx_files    <- strsplit("${cx_reports.join(',')}", ",")[[1]]
    min_cov     <- as.integer(${min_coverage})
    qvalue_cut  <- as.numeric(${qvalue_cutoff})
    meth_diff_cut <- as.numeric(${meth_diff_cutoff})
    dmr_mode    <- "${dmr_mode}"
    tile_size   <- as.integer(${tile_size})
    tile_step   <- as.integer(${tile_step})

    stopifnot(length(sample_ids) == length(treatments))
    stopifnot(length(sample_ids) == length(cx_files))

    # methRead expects a list of per-sample file paths and matching sample IDs
    # / treatment vector (1 = treatment/exposed, 0 = control/sham).
    obj <- methRead(
        as.list(cx_files),
        sample.id   = as.list(sample_ids),
        assembly    = "${params.genome_assembly ?: 'unspecified'}",
        treatment   = treatments,
        context     = "CpG",
        pipeline    = "bismarkCytosineReport",
        mincov      = min_cov
    )

    # Coverage filtering is re-applied here as the authoritative step
    # (percentile-based high-coverage outlier removal + min-coverage floor),
    # independent of the earlier per-sample COVERAGE_FILTER QC pass.
    filtered <- filterByCoverage(
        obj,
        lo.count = min_cov,
        lo.perc  = NULL,
        hi.count = NULL,
        hi.perc  = 99.9
    )

    normalized <- normalizeCoverage(filtered)

    if (dmr_mode == "tile") {
        tiled  <- tileMethylCounts(normalized, win.size = tile_size, step.size = tile_step, cov.bases = min_cov)
        united <- methylKit::unite(tiled, destrand = FALSE)
    } else {
        united <- methylKit::unite(normalized, destrand = FALSE)
    }

    diff <- calculateDiffMeth(united, mc.cores = ${task.cpus})

    # getMethylDiff with difference/qvalue thresholds applied explicitly and
    # recorded in the manifest (see emit_manifest.nf) rather than left as
    # implicit defaults.
    dmr_sig <- getMethylDiff(diff, difference = meth_diff_cut, qvalue = qvalue_cut)
    dmr_all <- getData(diff)

    region_ids <- paste0(
        "DMR_", seq_len(nrow(dmr_all))
    )

    out <- data.frame(
        region_id          = region_ids,
        chrom              = dmr_all\$chr,
        start              = dmr_all\$start,
        end                = dmr_all\$end,
        meth_diff_percent  = round(dmr_all\$meth.diff, 3),
        qvalue             = signif(dmr_all\$qvalue, 4),
        n_samples          = length(sample_ids),
        passes_significance_filter = (dmr_all\$qvalue <= qvalue_cut) & (abs(dmr_all\$meth.diff) >= meth_diff_cut),
        stringsAsFactors = FALSE
    )

    write.table(out, file = "${study_id}_${comparison_id}_dmr_raw.tsv",
                sep = "\\t", quote = FALSE, row.names = FALSE)

    qc <- data.frame(
        study_id            = "${study_id}",
        comparison_id       = "${comparison_id}",
        n_samples           = length(sample_ids),
        n_candidate_regions = nrow(dmr_all),
        n_significant_regions = nrow(dmr_sig),
        dmr_mode            = dmr_mode,
        tile_size           = ifelse(dmr_mode == "tile", tile_size, NA),
        min_coverage        = min_cov,
        qvalue_cutoff       = qvalue_cut,
        meth_diff_cutoff    = meth_diff_cut
    )
    write.table(qc, file = "${study_id}_${comparison_id}_methylkit_qc.tsv",
                sep = "\\t", quote = FALSE, row.names = FALSE)
    EOF_R

    Rscript run_methylkit.R

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        r-base: \$(R --version | head -n1 | sed 's/R version //; s/ .*//')
        bioconductor-methylkit: \$(Rscript -e 'cat(as.character(packageVersion("methylKit")))')
    END_VERSIONS
    """
}
