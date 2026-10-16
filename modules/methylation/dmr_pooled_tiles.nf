// Tile-level differential methylation on pooled replicate units
// (raw_reanalysis mode only), selected by params.methylation.dmr_engine =
// 'limma' (recommended) or 'dss'.
//
// Written for designs with very few replicate units per arm, where
// methylKit's per-region McCullagh-Nelder overdispersion is unstable: with 2
// families per arm, a region whose families happen to agree gets an estimated
// dispersion below 1, methylKit resets it to 1, and the test falls back to a
// binomial test on summed read counts (CIBNOR2021_HEAT_RRSS, first pass:
// 296 of 397 tiles at q < 0.05 differed by 5 points or less).
//
// Two tests on the same pooled tiles:
//   limma (default): moderated t on per-unit M-values (logit of methylation),
//     variance shrunk across tiles by empirical Bayes, with small-sample df.
//   dss: DSS::DMLtest, beta-binomial with shrunken dispersion, Wald test.
// On synthetic data with 2 vs 2 families, family-level variation and 20 true
// regions among 400 tiles: limma found 4/20 with 1/380 false positives, and
// 5.0% of null tiles had p < 0.05 (calibrated); DSS found 17/20 but with
// 57/380 false positives, because its Wald test has no small-sample
// correction. Prefer limma unless there are several units per arm.
//
// Steps, per library CpG report (Bismark cytosine report, coverage-filtered
// by COVERAGE_FILTER at params.methylation.min_coverage):
//   1. destrand: a minus-strand CpG at position p is the same CpG as the
//      plus-strand C at p - 1; their counts are summed.
//   2. pool by replicate unit (e.g. family), when replicate_unit is set;
//      otherwise each library is its own unit.
//   3. keep CpGs covered by >= unit_min_coverage reads in EVERY unit.
//   4. sum kept CpGs into non-overlapping tiles of tile_size bp; keep tiles
//      with >= tile_min_cpgs CpGs.
//   5. the chosen test on the tiles, treatment units vs control units.
//      meth_diff_percent is treatment minus control mean methylation (%);
//      qvalue is the Benjamini-Hochberg FDR (limma) or DSS's FDR.

process DMR_POOLED_TILES {
    tag "${study_id}:${comparison_id}"
    label 'process_high'
    container 'bioconductor/bioconductor_docker:RELEASE_3_18'
    publishDir "${params.outdir}/methylation/dmr", mode: params.publish_mode

    input:
    val study_id
    val comparison_id
    path cx_reports          // one coverage-filtered cytosine report per library
    val sample_ids
    val treatment_labels     // 0 = control, 1 = treatment, matched to sample_ids
    val unit_min_coverage
    val qvalue_cutoff
    val meth_diff_cutoff
    val tile_size
    val tile_min_cpgs
    val replicate_unit
    path design_sheet
    val dmr_test              // 'limma' or 'dss'

    output:
    path "${study_id}_${comparison_id}_dmr_raw.tsv", emit: dmr_table
    path "${study_id}_${comparison_id}_methylkit_qc.tsv", emit: qc_metrics
    path "versions.yml", emit: versions

    script:
    """
    cat <<-'EOF_R' > run_pooled_tiles.R
    suppressPackageStartupMessages({
        library(data.table)
    })
    dmr_test      <- "${dmr_test}"
    stopifnot(dmr_test %in% c("limma", "dss"))

    sample_ids    <- strsplit("${sample_ids.join(',')}", ",")[[1]]
    treatments    <- as.integer(strsplit("${treatment_labels.join(',')}", ",")[[1]])
    cx_files      <- strsplit("${cx_reports.join(',')}", ",")[[1]]
    unit_min_cov  <- as.integer(${unit_min_coverage})
    qvalue_cut    <- as.numeric(${qvalue_cutoff})
    meth_diff_cut <- as.numeric(${meth_diff_cutoff})
    tile_size     <- as.integer(${tile_size})
    tile_min_cpgs <- as.integer(${tile_min_cpgs})
    unit_col      <- trimws("${replicate_unit}")
    design_path   <- "${design_sheet}"

    stopifnot(length(sample_ids) == length(treatments), length(sample_ids) == length(cx_files))
    # Nextflow stages files under their own names; match each report to its
    # sample by the sample id it starts with.
    cx_files <- vapply(sample_ids, function(s) {
        hit <- cx_files[startsWith(basename(cx_files), paste0(s, "."))]
        if (length(hit) != 1) stop("expected one report for ", s, ", found ", length(hit))
        hit
    }, character(1))

    units <- sample_ids
    if (nzchar(unit_col)) {
        design <- read.delim(design_path, colClasses = "character", check.names = FALSE)
        if (!unit_col %in% names(design)) stop("replicate_unit column '", unit_col, "' is not in the design sheet")
        key <- if ("run_accession" %in% names(design)) "run_accession" else "sample_id"
        units <- design[[unit_col]][match(sample_ids, design[[key]])]
        if (any(is.na(units) | !nzchar(units))) {
            stop("no '", unit_col, "' for: ", paste(sample_ids[is.na(units) | !nzchar(units)], collapse = ", "))
        }
    }
    arms_per_unit <- tapply(treatments, units, function(x) length(unique(x)))
    if (any(arms_per_unit > 1)) {
        stop("replicate_unit(s) spanning both arms: ", paste(names(arms_per_unit)[arms_per_unit > 1], collapse = ", "))
    }
    unit_levels <- sort(unique(units))
    unit_trt <- as.integer(tapply(treatments, units, function(x) x[1])[unit_levels])
    if (any(table(factor(unit_trt, levels = c(0, 1))) < 2)) {
        stop("each arm needs at least 2 replicate units; got ",
             paste(sprintf("%s=%d", unit_levels, unit_trt), collapse = ", "))
    }

    read_destranded <- function(f) {
        d <- fread(cmd = paste("gzip -dc", shQuote(f)), header = FALSE, select = 1:5,
                   col.names = c("chr", "pos", "strand", "m", "u"))
        d[strand == "-", pos := pos - 1L]
        d[, .(m = sum(m), u = sum(u)), by = .(chr, pos)]
    }

    pooled <- list()
    for (u in unit_levels) {
        libs <- which(units == u)
        d <- rbindlist(lapply(cx_files[libs], read_destranded))
        d <- d[, .(m = sum(m), u = sum(u)), by = .(chr, pos)]
        d <- d[(m + u) >= unit_min_cov]
        setnames(d, c("m", "u"), paste0(c("m_", "u_"), u))
        pooled[[u]] <- d
        message(sprintf("unit %s: %d libraries, %d destranded CpGs at >= %dx", u, length(libs), nrow(d), unit_min_cov))
    }
    wide <- Reduce(function(a, b) merge(a, b, by = c("chr", "pos")), pooled)
    message(sprintf("CpGs at >= %dx in all %d units: %d", unit_min_cov, length(unit_levels), nrow(wide)))

    wide[, tile := (pos - 1L) %/% tile_size]
    mcols <- paste0("m_", unit_levels); ucols <- paste0("u_", unit_levels)
    tiles <- wide[, c(list(n_cpg = .N), lapply(.SD, sum)), by = .(chr, tile), .SDcols = c(mcols, ucols)]
    tiles <- tiles[n_cpg >= tile_min_cpgs]
    message(sprintf("tiles of %d bp with >= %d CpGs: %d", tile_size, tile_min_cpgs, nrow(tiles)))
    if (nrow(tiles) == 0) stop("no tiles pass the coverage and CpG-count filters")

    M   <- as.matrix(tiles[, ..mcols])
    Cov <- M + as.matrix(tiles[, ..ucols])
    colnames(M) <- colnames(Cov) <- unit_levels
    treat_units   <- unit_levels[unit_trt == 1]
    control_units <- unit_levels[unit_trt == 0]
    tiles[, start := as.integer(tile * tile_size + 1L)]
    df_prior <- NA_real_

    if (dmr_test == "limma") {
        suppressPackageStartupMessages(library(limma))
        beta <- (M + 0.5) / (Cov + 1)
        mval <- log2(beta / (1 - beta))
        X <- model.matrix(~ factor(unit_trt, levels = c(0, 1)))
        fit <- eBayes(lmFit(mval, X))
        tt <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
        df_prior <- fit\$df.prior
        diff <- rowMeans(beta[, treat_units, drop = FALSE]) - rowMeans(beta[, control_units, drop = FALSE])
        res <- data.table(
            chr  = tiles\$chr, pos = tiles\$start,
            # treatment minus control, as for methylKit's meth.diff
            diff = diff,
            # SE on the methylation-difference scale implied by limma's
            # moderated t, so that diff / se reproduces the test statistic
            # when the effect is pooled across studies.
            diff.se = ifelse(tt\$t != 0, abs(diff / tt\$t), NA_real_),
            pval = tt\$P.Value,
            fdr  = tt\$adj.P.Val
        )
    } else {
        suppressPackageStartupMessages({ library(DSS); library(bsseq) })
        bs <- BSseq(chr = tiles\$chr, pos = tiles\$start, M = M, Cov = Cov, sampleNames = unit_levels)
        # group1 = treatment, so diff = treatment - control.
        dml <- as.data.table(DMLtest(bs, group1 = treat_units, group2 = control_units, smoothing = FALSE))
        res <- dml[tiles[, .(chr, start)], on = .(chr, pos = start), nomatch = 0L]
    }

    out <- data.frame(
        region_id         = paste0("DMR_", seq_len(nrow(res))),
        chrom             = res\$chr,
        start             = res\$pos,
        end               = res\$pos + tile_size - 1L,
        meth_diff_percent = round(100 * res\$diff, 3),
        qvalue            = signif(res\$fdr, 4),
        pvalue            = signif(res\$pval, 4),
        se_meth_diff_percent = signif(100 * res\$diff.se, 4),
        analysis_method   = ifelse(dmr_test == "limma", "limma_pooled_tiles", "DSS_pooled_tiles"),
        n_samples         = length(sample_ids),
        n_replicate_units = length(unit_levels),
        passes_significance_filter = !is.na(res\$fdr) & res\$fdr <= qvalue_cut &
                                     abs(100 * res\$diff) >= meth_diff_cut,
        stringsAsFactors  = FALSE
    )
    write.table(out, file = "${study_id}_${comparison_id}_dmr_raw.tsv",
                sep = "\\t", quote = FALSE, row.names = FALSE)

    qc <- data.frame(
        study_id              = "${study_id}",
        comparison_id         = "${comparison_id}",
        n_samples             = length(sample_ids),
        n_candidate_regions   = nrow(out),
        n_significant_regions = sum(out\$passes_significance_filter),
        dmr_mode              = "tile",
        tile_size             = tile_size,
        min_coverage          = unit_min_cov,
        qvalue_cutoff         = qvalue_cut,
        meth_diff_cutoff      = meth_diff_cut,
        replicate_unit        = ifelse(nzchar(unit_col), unit_col, NA),
        n_replicate_units     = length(unit_levels),
        overdispersion        = ifelse(dmr_test == "limma", "empirical-Bayes moderated variance", "DSS shrunken beta-binomial"),
        test                  = ifelse(dmr_test == "limma", "limma moderated t on unit M-values, BH", "DSS::DMLtest (Wald)"),
        engine                = dmr_test,
        df_prior              = df_prior,
        destranded            = TRUE,
        tile_min_cpgs         = tile_min_cpgs,
        n_cpgs_all_units      = nrow(wide)
    )
    write.table(qc, file = "${study_id}_${comparison_id}_methylkit_qc.tsv",
                sep = "\\t", quote = FALSE, row.names = FALSE)
    EOF_R

    Rscript run_pooled_tiles.R

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        r-base: \$(R --version | head -n1 | sed 's/R version //; s/ .*//')
        bioconductor-limma: \$(Rscript -e 'cat(as.character(packageVersion("limma")))')
        bioconductor-dss: \$(Rscript -e 'cat(as.character(packageVersion("DSS")))')
        bioconductor-bsseq: \$(Rscript -e 'cat(as.character(packageVersion("bsseq")))')
        r-data.table: \$(Rscript -e 'cat(as.character(packageVersion("data.table")))')
    END_VERSIONS
    """
}
