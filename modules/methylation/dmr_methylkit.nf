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
    val replicate_unit        // design-sheet column naming the unit of replication, or ''
    path design_sheet         // TSV mapping run_accession/sample_id -> unit (placeholder when unused)

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
    unit_col    <- trimws("${replicate_unit}")
    design_path <- "${design_sheet}"

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

    # Unit-level replication. When the phenotype or treatment is assigned to a
    # group of libraries (a family, a tank), the libraries are not independent
    # replicates of it. Sum each unit's methylated / unmethylated counts per
    # region into one pseudo-library, then test with the McCullagh-Nelder
    # overdispersion correction and an F test, so that the error comes from
    # unit-to-unit variation. Without it, summed read counts would make the
    # logistic test treat every read as independent.
    n_units <- length(sample_ids)
    overdispersion <- "none"
    test_used <- "Chisq"
    if (nzchar(unit_col)) {
        design <- read.delim(design_path, colClasses = "character", check.names = FALSE)
        if (!unit_col %in% names(design)) {
            stop("replicate_unit column '", unit_col, "' is not in the design sheet")
        }
        key <- if ("run_accession" %in% names(design)) "run_accession" else "sample_id"
        units <- design[[unit_col]][match(sample_ids, design[[key]])]
        if (any(is.na(units) | !nzchar(units))) {
            stop("no '", unit_col, "' for: ", paste(sample_ids[is.na(units) | !nzchar(units)], collapse = ", "))
        }
        arms_per_unit <- tapply(treatments, units, function(x) length(unique(x)))
        if (any(arms_per_unit > 1)) {
            stop("replicate_unit(s) spanning both arms: ",
                 paste(names(arms_per_unit)[arms_per_unit > 1], collapse = ", "))
        }
        unit_levels <- sort(unique(units))
        unit_trt <- as.integer(tapply(treatments, units, function(x) x[1])[unit_levels])
        if (any(table(factor(unit_trt, levels = c(0, 1))) < 2)) {
            stop("each arm needs at least 2 replicate units; got ",
                 paste(sprintf("%s=%d", unit_levels, unit_trt), collapse = ", "))
        }
        d <- getData(united)
        agg <- d[, c("chr", "start", "end", "strand")]
        for (i in seq_along(unit_levels)) {
            cols <- which(units == unit_levels[i])
            agg[[paste0("coverage", i)]] <- rowSums(d[, united@coverage.index[cols], drop = FALSE])
            agg[[paste0("numCs", i)]]    <- rowSums(d[, united@numCs.index[cols], drop = FALSE])
            agg[[paste0("numTs", i)]]    <- rowSums(d[, united@numTs.index[cols], drop = FALSE])
        }
        message("Pooled by '", unit_col, "': ", length(unit_levels), " units from ",
                length(sample_ids), " libraries (",
                paste(sprintf("%s:%d=%d", unit_levels, unit_trt, as.integer(table(units)[unit_levels])),
                      collapse = ", "), ")")
        united <- new("methylBase", agg,
                      sample.ids     = as.character(unit_levels),
                      assembly       = united@assembly,
                      context        = united@context,
                      treatment      = unit_trt,
                      coverage.index = 5 + 3 * (seq_along(unit_levels) - 1),
                      numCs.index    = 6 + 3 * (seq_along(unit_levels) - 1),
                      numTs.index    = 7 + 3 * (seq_along(unit_levels) - 1),
                      destranded     = united@destranded,
                      resolution     = united@resolution)
        n_units <- length(unit_levels)
        overdispersion <- "MN"
        test_used <- "F"
        diff <- calculateDiffMeth(united, overdispersion = "MN", test = "F", mc.cores = ${task.cpus})
    } else {
        diff <- calculateDiffMeth(united, mc.cores = ${task.cpus})
    }

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
        n_replicate_units  = n_units,
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
        meth_diff_cutoff    = meth_diff_cut,
        replicate_unit      = ifelse(nzchar(unit_col), unit_col, NA),
        n_replicate_units   = n_units,
        overdispersion      = overdispersion,
        test                = test_used
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
