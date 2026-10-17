// Differential protein abundance with limma's moderated t-test on normalized
// log2 abundances, raw_reanalysis mode only.
//
// A protein is tested only if it is detected (non-missing) in at least
// --proteomics.min_detected_per_arm samples in each arm; every protein, tested
// or not, is listed with its status in the *_detection_status.tsv table.
// lfcSE is limma's posterior standard error of log2FC
// (sqrt(s2.post) * stdev.unscaled), so evidence records carry one and can be
// pooled. --proteomics.quantification (intensity | spectral_count) is written
// into analysis_method so spectral-count and intensity evidence stay
// distinguishable downstream.
//
// STATUS: executed by the raw-mode CI smoke test (tests/fixtures/proteomics_raw).
process DIFFERENTIAL_ABUNDANCE {
    tag "${study_id}:${comparison_id}"
    label 'process_medium'
    container 'bioconductor/bioconductor_docker:RELEASE_3_18'
    publishDir "${params.outdir}/proteomics/differential_abundance", mode: params.publish_mode

    input:
    tuple val(study_id), val(comparison_id), path(normalized_long_table)
    val min_detected_per_arm
    val quantification

    output:
    tuple val(study_id), val(comparison_id), path("${study_id}_${comparison_id}_limma_results.tsv"), emit: results
    path "${study_id}_${comparison_id}_detection_status.tsv", emit: detection_status
    path "${study_id}_${comparison_id}_diffabundance_qc.json", emit: da_qc
    path "versions.yml", emit: versions

    script:
    """
    cat <<-'EOF' > run_limma.R
    #!/usr/bin/env Rscript
    suppressPackageStartupMessages({
        library(limma)
        library(jsonlite)
    })

    args <- list(
        table          = "${normalized_long_table}",
        study_id       = "${study_id}",
        comparison_id  = "${comparison_id}",
        min_detected   = as.integer("${min_detected_per_arm}"),
        quantification = "${quantification}",
        out_tsv        = "${study_id}_${comparison_id}_limma_results.tsv",
        out_status     = "${study_id}_${comparison_id}_detection_status.tsv",
        out_qc         = "${study_id}_${comparison_id}_diffabundance_qc.json"
    )
    if (!(args\$quantification %in% c("intensity", "spectral_count"))) {
        stop("proteomics.quantification must be 'intensity' or 'spectral_count', got '", args\$quantification, "'")
    }

    long_df <- read.delim(args\$table, stringsAsFactors = FALSE, check.names = FALSE)
    stopifnot(all(c("protein_accession", "sample_id", "group", "log2_abundance_normalized") %in% colnames(long_df)))

    samples <- unique(long_df[, c("sample_id", "group")])
    if (any(duplicated(samples\$sample_id))) stop("a sample_id has more than one group")
    condition <- factor(samples\$group, levels = c("control", "treatment"))
    if (any(is.na(condition))) stop("group must be 'control' or 'treatment' for every sample")

    proteins <- unique(long_df\$protein_accession)
    expr <- matrix(NA_real_, nrow = length(proteins), ncol = nrow(samples),
                   dimnames = list(proteins, samples\$sample_id))
    expr[cbind(match(long_df\$protein_accession, proteins), match(long_df\$sample_id, samples\$sample_id))] <-
        long_df\$log2_abundance_normalized

    detected <- !is.na(expr)
    n_control   <- rowSums(detected[, condition == "control", drop = FALSE])
    n_treatment <- rowSums(detected[, condition == "treatment", drop = FALSE])
    tested <- n_control >= args\$min_detected & n_treatment >= args\$min_detected
    status <- ifelse(tested, "tested",
              ifelse(n_control == 0 & n_treatment == 0, "not_detected",
              ifelse(n_control < args\$min_detected & n_treatment < args\$min_detected,
                     "detected_too_rarely_in_both_arms",
              ifelse(n_control < args\$min_detected, "detected_too_rarely_in_control",
                     "detected_too_rarely_in_treatment"))))
    write.table(data.frame(protein_accession = proteins,
                           n_replicates_control = sum(condition == "control"),
                           n_detected_control = n_control,
                           n_replicates_treatment = sum(condition == "treatment"),
                           n_detected_treatment = n_treatment,
                           status = status, stringsAsFactors = FALSE),
                args\$out_status, sep = "\\t", quote = FALSE, row.names = FALSE)
    if (sum(tested) < 2) stop("fewer than 2 proteins pass min_detected_per_arm = ", args\$min_detected)

    design <- model.matrix(~condition)
    fit <- eBayes(lmFit(expr[tested, , drop = FALSE], design))
    coef <- "conditiontreatment"
    tt <- topTable(fit, coef = coef, number = Inf, sort.by = "none")
    out <- data.frame(
        protein_accession    = rownames(tt),
        log2FC               = tt\$logFC,
        lfcSE                = sqrt(fit\$s2.post) * fit\$stdev.unscaled[, coef],
        pvalue               = tt\$P.Value,
        padj                 = tt\$adj.P.Val,
        ave_expr             = tt\$AveExpr,
        t_stat               = tt\$t,
        n_detected_control   = n_control[rownames(tt)],
        n_detected_treatment = n_treatment[rownames(tt)],
        analysis_method      = paste0(args\$quantification, "_limma_moderated_t"),
        stringsAsFactors = FALSE
    )
    write.table(out, args\$out_tsv, sep = "\\t", quote = FALSE, row.names = FALSE, na = "NA")

    qc <- list(
        study_id = args\$study_id,
        comparison_id = args\$comparison_id,
        n_replicates_control = sum(condition == "control"),
        n_replicates_treatment = sum(condition == "treatment"),
        min_detected_per_arm = args\$min_detected,
        n_proteins_in = length(proteins),
        n_proteins_tested = sum(tested),
        n_proteins_not_tested = as.list(table(status[!tested])),
        n_padj_below_0.05 = sum(out\$padj < 0.05, na.rm = TRUE),
        n_padj_below_0.1 = sum(out\$padj < 0.1, na.rm = TRUE),
        limma_df_prior = unname(fit\$df.prior[1]),
        analysis_method = out\$analysis_method[1]
    )
    writeLines(toJSON(qc, auto_unbox = TRUE, pretty = TRUE, null = "null"), args\$out_qc)
    cat(sprintf("AREE proteomics limma: %d of %d proteins tested on %d vs %d samples\\n",
                sum(tested), length(proteins), sum(condition == "control"), sum(condition == "treatment")))
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
    touch ${study_id}_${comparison_id}_limma_results.tsv
    touch ${study_id}_${comparison_id}_detection_status.tsv
    touch ${study_id}_${comparison_id}_diffabundance_qc.json
    touch versions.yml
    """
}
