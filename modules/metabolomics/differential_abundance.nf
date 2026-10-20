// Differential abundance testing on normalized (log2) metabolite feature
// intensities (raw_reanalysis mode only), using limma's moderated t-test
// (empirical Bayes shrinkage of variance) — the standard approach for
// small-n omics intensity matrices, widely used for metabolomics feature
// tables as well as microarray/proteomics data.
//
// With --metabolomics.replicate_unit (a sample-sheet column, e.g. tank), log2
// values are averaged within each unit before the fit, so the units are the
// replicates; a unit may not sit in both arms. Missing values (below
// detection) are left out of each unit's mean. A feature is tested only if it
// is detected in at least --metabolomics.min_detected_per_arm replicates in
// each arm; every feature, tested or not, is listed in the missingness report.
//
// lfcSE is limma's posterior standard error of log2FC
// (sqrt(s2.post) * stdev.unscaled), so the evidence record carries one.
//
// STATUS: executed on the real targeted panel of DUPERRET2025_POMS_STARVE
// (2026-10-16) and by the raw-mode CI smoke test.
process DIFFERENTIAL_ABUNDANCE {
    tag "${study_id}:${comparison_id}"
    label 'process_medium'
    container 'bioconductor/bioconductor_docker:RELEASE_3_18'
    publishDir "${params.outdir}/metabolomics/differential_abundance", mode: params.publish_mode

    input:
    val study_id
    val comparison_id
    path normalized_table
    path sample_sheet
    val replicate_unit
    val min_detected_per_arm
    val acquisition

    output:
    tuple val(study_id), val(comparison_id), path("${study_id}_${comparison_id}_diffabundance_raw.tsv"), emit: results
    path "${study_id}_${comparison_id}_missingness.tsv", emit: missingness
    path "${study_id}_${comparison_id}_diffabundance_qc.json", emit: da_qc
    path "versions.yml", emit: versions

    script:
    """
    cat <<-'EOF' > run_limma.R
    #!/usr/bin/env Rscript
    # AREE metabolomics differential abundance (limma moderated t-test)
    suppressPackageStartupMessages({
        library(limma)
        library(jsonlite)
    })

    args <- list(
        normalized_table = "${normalized_table}",
        sample_sheet     = "${sample_sheet}",
        study_id         = "${study_id}",
        comparison_id    = "${comparison_id}",
        replicate_unit   = "${replicate_unit ?: ''}",
        min_detected     = as.integer("${min_detected_per_arm}"),
        acquisition      = "${acquisition}",
        out_tsv          = "${study_id}_${comparison_id}_diffabundance_raw.tsv",
        out_missing      = "${study_id}_${comparison_id}_missingness.tsv",
        out_qc           = "${study_id}_${comparison_id}_diffabundance_qc.json"
    )
    if (!(args\$acquisition %in% c("untargeted", "targeted"))) {
        stop("acquisition must be 'untargeted' or 'targeted', got '", args\$acquisition, "'")
    }

    features <- read.delim(args\$normalized_table, stringsAsFactors = FALSE, check.names = FALSE)
    samples  <- read.delim(args\$sample_sheet, stringsAsFactors = FALSE, check.names = FALSE)
    samples\$condition <- factor(samples\$condition, levels = c("control", "treatment"))
    if (any(is.na(samples\$condition))) {
        stop("sample sheet condition must be 'control' or 'treatment' for every sample")
    }

    meta_cols <- intersect(
        c("feature_id", "putative_metabolite_name", "annotation_confidence_level", "platform"),
        colnames(features)
    )
    sample_ids <- intersect(samples\$sample_id, colnames(features))
    stopifnot(length(sample_ids) >= 2)
    samples <- samples[match(sample_ids, samples\$sample_id), ]
    expr <- as.matrix(features[, sample_ids, drop = FALSE])
    rownames(expr) <- features\$feature_id

    unit_col <- trimws(args\$replicate_unit)
    if (nzchar(unit_col)) {
        if (!(unit_col %in% colnames(samples))) {
            stop("replicate_unit column '", unit_col, "' is not in the sample sheet")
        }
        units <- as.character(samples[[unit_col]])
        if (any(is.na(units) | units == "")) {
            stop("replicate_unit column '", unit_col, "' is empty for some samples")
        }
        arms_per_unit <- tapply(as.character(samples\$condition), units, function(x) length(unique(x)))
        if (any(arms_per_unit > 1)) {
            stop("replicate_unit(s) spanning both arms: ",
                 paste(names(arms_per_unit)[arms_per_unit > 1], collapse = ", "))
        }
        unit_ids <- unique(units)
        expr <- sapply(unit_ids, function(u) {
            m <- rowMeans(expr[, units == u, drop = FALSE], na.rm = TRUE)
            m[is.nan(m)] <- NA
            m
        })
        expr <- matrix(expr, nrow = nrow(features), dimnames = list(features\$feature_id, unit_ids))
        condition <- factor(as.character(samples\$condition[match(unit_ids, units)]),
                            levels = c("control", "treatment"))
    } else {
        condition <- samples\$condition
    }

    detected <- !is.na(expr)
    n_control   <- rowSums(detected[, condition == "control", drop = FALSE])
    n_treatment <- rowSums(detected[, condition == "treatment", drop = FALSE])
    tested <- n_control >= args\$min_detected & n_treatment >= args\$min_detected
    reason <- ifelse(tested, "tested",
              ifelse(n_control == 0 & n_treatment == 0, "not_detected",
              ifelse(n_control < args\$min_detected & n_treatment < args\$min_detected,
                     "detected_too_rarely_in_both_arms",
              ifelse(n_control < args\$min_detected, "detected_too_rarely_in_control",
                     "detected_too_rarely_in_treatment"))))
    missingness <- data.frame(features[, meta_cols, drop = FALSE],
                              n_replicates_control = sum(condition == "control"),
                              n_detected_control = n_control,
                              n_replicates_treatment = sum(condition == "treatment"),
                              n_detected_treatment = n_treatment,
                              status = reason, check.names = FALSE, stringsAsFactors = FALSE)
    write.table(missingness, args\$out_missing, sep = "\\t", quote = FALSE, row.names = FALSE, na = "NA")
    if (sum(tested) < 2) stop("fewer than 2 features pass min_detected_per_arm = ", args\$min_detected)

    design <- model.matrix(~condition)
    fit <- eBayes(lmFit(expr[tested, , drop = FALSE], design))
    coef <- "conditiontreatment"
    tt <- topTable(fit, coef = coef, number = Inf, sort.by = "none")
    out <- data.frame(
        feature_id      = rownames(tt),
        log2FC          = tt\$logFC,
        lfcSE           = sqrt(fit\$s2.post) * fit\$stdev.unscaled[, coef],
        pvalue          = tt\$P.Value,
        padj            = tt\$adj.P.Val,
        analysis_method = paste0(args\$acquisition, "_MS_limma_moderated_t"),
        stringsAsFactors = FALSE
    )
    # Carry through annotation metadata already attached upstream so
    # downstream steps do not need to re-join it.
    for (col in setdiff(meta_cols, "feature_id")) {
        out[[col]] <- features[[col]][match(out\$feature_id, features\$feature_id)]
    }
    out <- out[order(out\$padj, na.last = TRUE), ]
    write.table(out, args\$out_tsv, sep = "\\t", quote = FALSE, row.names = FALSE, na = "NA")

    qc <- list(
        study_id = args\$study_id,
        comparison_id = args\$comparison_id,
        replicate_unit = if (nzchar(unit_col)) unit_col else NULL,
        n_replicates_control = sum(condition == "control"),
        n_replicates_treatment = sum(condition == "treatment"),
        min_detected_per_arm = args\$min_detected,
        n_features_in = nrow(features),
        n_features_tested = sum(tested),
        n_features_not_tested = as.list(table(reason[!tested])),
        n_padj_below_0.05 = sum(out\$padj < 0.05, na.rm = TRUE),
        limma_df_prior = unname(fit\$df.prior[1]),
        analysis_method = out\$analysis_method[1]
    )
    writeLines(toJSON(qc, auto_unbox = TRUE, pretty = TRUE, null = "null"), args\$out_qc)
    cat(sprintf(
        "AREE limma differential abundance complete for %s:%s: %d of %d features tested on %d vs %d replicates\\n",
        args\$study_id, args\$comparison_id, sum(tested), nrow(features),
        sum(condition == "control"), sum(condition == "treatment")
    ))
    EOF

    Rscript run_limma.R

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        r-base: \$(Rscript -e 'cat(as.character(getRversion()))')
        bioconductor-limma: \$(Rscript -e 'cat(as.character(packageVersion("limma")))')
        r-jsonlite: \$(Rscript -e 'cat(as.character(packageVersion("jsonlite")))')
    END_VERSIONS
    """

    stub:
    """
    touch ${study_id}_${comparison_id}_diffabundance_raw.tsv
    touch ${study_id}_${comparison_id}_missingness.tsv
    touch ${study_id}_${comparison_id}_diffabundance_qc.json
    touch versions.yml
    """
}
