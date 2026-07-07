// Differential abundance testing on normalized (log2) metabolite feature
// intensities (raw_reanalysis mode only), using limma's moderated t-test
// (empirical Bayes shrinkage of variance) — the standard approach for
// small-n omics intensity matrices, widely used for metabolomics feature
// tables as well as microarray/proteomics data.
//
// STATUS: structurally complete DSL2 process. The R script below is real,
// syntactically valid limma usage (lmFit -> eBayes -> topTable) that would
// run correctly given a real normalized intensity matrix and sample sheet.
// NOT executed in this build — no compute / no real normalized-matrix input
// exists here (see workflows/metabolomics/README.md).

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

    output:
    tuple val(study_id), val(comparison_id), path("${study_id}_${comparison_id}_diffabundance_raw.tsv"), emit: results
    path "versions.yml", emit: versions

    script:
    // Real limma moderated-t-test skeleton. Output columns are written to
    // already match the STANDARDIZE_OUTPUT target schema (feature_id,
    // log2FC, pvalue, padj) plus metadata columns carried through so that
    // step is largely a pass-through/validation in raw mode.
    """
    cat <<-'EOF' > run_limma.R
    #!/usr/bin/env Rscript
    # AREE metabolomics differential abundance (limma moderated t-test)
    #
    # NOTE: structurally complete, unexecuted in this build. Written to run
    # correctly against a real log2-normalized feature x sample intensity
    # matrix and a real two-group sample sheet.
    suppressPackageStartupMessages({
        library(limma)
    })

    args <- list(
        normalized_table = "${normalized_table}",
        sample_sheet     = "${sample_sheet}",
        study_id         = "${study_id}",
        comparison_id    = "${comparison_id}",
        out_tsv          = "${study_id}_${comparison_id}_diffabundance_raw.tsv"
    )

    features <- read.delim(args\$normalized_table, stringsAsFactors = FALSE, check.names = FALSE)
    samples  <- read.delim(args\$sample_sheet, stringsAsFactors = FALSE)
    samples\$condition <- factor(samples\$condition, levels = c("control", "treatment"))

    meta_cols <- intersect(
        c("feature_id", "putative_metabolite_name", "annotation_confidence_level"),
        colnames(features)
    )
    sample_ids <- intersect(samples\$sample_id, colnames(features))
    stopifnot(length(sample_ids) >= 2)

    samples <- samples[match(sample_ids, samples\$sample_id), ]

    expr <- as.matrix(features[, sample_ids, drop = FALSE])
    rownames(expr) <- features\$feature_id

    design <- model.matrix(~condition, data = samples)

    fit <- lmFit(expr, design)
    fit <- eBayes(fit)

    tt <- topTable(fit, coef = "conditiontreatment", number = Inf, sort.by = "none")

    out <- data.frame(
        feature_id = rownames(tt),
        log2FC     = tt\$logFC,
        pvalue     = tt\$P.Value,
        padj       = tt\$adj.P.Val,
        stringsAsFactors = FALSE
    )

    # Carry through annotation metadata already attached upstream so
    # downstream steps do not need to re-join it.
    if ("putative_metabolite_name" %in% meta_cols) {
        out\$putative_metabolite_name <- features\$putative_metabolite_name[match(out\$feature_id, features\$feature_id)]
    }
    if ("annotation_confidence_level" %in% meta_cols) {
        out\$annotation_confidence_level <- features\$annotation_confidence_level[match(out\$feature_id, features\$feature_id)]
    }

    out <- out[order(out\$padj, na.last = TRUE), ]

    write.table(out, args\$out_tsv, sep = "\\t", quote = FALSE, row.names = FALSE, na = "NA")

    cat(sprintf(
        "AREE limma differential abundance complete for %s:%s — %d features tested\\n",
        args\$study_id, args\$comparison_id, nrow(out)
    ))
    EOF

    Rscript run_limma.R

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        r-base: \$(Rscript -e 'cat(as.character(getRversion()))')
        bioconductor-limma: \$(Rscript -e 'cat(as.character(packageVersion("limma")))')
    END_VERSIONS
    """

    stub:
    """
    touch ${study_id}_${comparison_id}_diffabundance_raw.tsv
    touch versions.yml
    """
}
