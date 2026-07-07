// Differential protein abundance via limma, appropriate for proteomics
// abundance matrices (as opposed to count-based RNA-seq DE tools) because it
// handles continuous, approximately-normal (post log2/median-normalization)
// data and tolerates missing values in individual features without requiring
// a full count model.
//
// STATUS: structurally complete DSL2 process wrapping a real, syntactically
// valid R/limma script (pivot to wide protein x sample matrix -> lmFit ->
// eBayes -> topTable). NOT executed against real data in this build — no
// compute budget / real abundance matrix available (see
// workflows/proteomics/README.md). raw_reanalysis mode only.
//
// NA handling note: limma's lmFit() tolerates NAs in individual cells of the
// data matrix (features with too few non-NA observations in a group simply
// get NA/less-precise statistics rather than failing the whole fit), which is
// why limma is preferred here over tools that require a complete matrix. We
// additionally drop proteins missing in >50% of samples before fitting
// (`max_missing_frac` below) since a linear model on a mostly-imputed/absent
// feature is not meaningful; this threshold is a documented, adjustable
// simplification, not a hidden default.

process DIFFERENTIAL_ABUNDANCE {
    tag "${study_id}:${comparison_id}"
    label 'process_medium'
    container 'bioconductor/bioconductor_docker:RELEASE_3_18'
    publishDir "${params.outdir}/proteomics/differential_abundance", mode: params.publish_mode

    input:
    tuple val(study_id), val(comparison_id), path(normalized_long_table)
    tuple val(study_id2), val(comparison_id2), path(missingness_per_protein)

    output:
    tuple val(study_id), val(comparison_id), path("${study_id}_${comparison_id}_limma_results.tsv"), emit: results
    path "${study_id}_${comparison_id}_limma.RData", emit: rdata
    path "versions.yml", emit: versions

    script:
    """
    #!/usr/bin/env Rscript

    suppressMessages({
      library(limma)
      library(tidyr)
      library(dplyr)
    })

    study_id <- "${study_id}"
    comparison_id <- "${comparison_id}"
    max_missing_frac <- 0.5  # drop proteins missing in >50% of samples before fitting

    long_df <- read.delim("${normalized_long_table}", sep = "\\t", stringsAsFactors = FALSE)
    stopifnot(all(c("protein_accession", "sample_id", "group", "log2_abundance_normalized") %in% colnames(long_df)))

    miss_df <- read.delim("${missingness_per_protein}", sep = "\\t", stringsAsFactors = FALSE)

    keep_proteins <- miss_df\$protein_accession[miss_df\$missingness_percent <= (max_missing_frac * 100)]
    cat(sprintf(
      "Retaining %d / %d proteins with <= %.0f%% missingness for differential abundance\\n",
      length(keep_proteins), nrow(miss_df), max_missing_frac * 100
    ))
    long_df <- long_df[long_df\$protein_accession %in% keep_proteins, ]

    # Pivot long -> wide: rows = proteins, columns = samples. limma expects a
    # numeric matrix with proteins as rows and samples as columns; NAs are
    # tolerated in individual cells by lmFit().
    wide_df <- long_df %>%
      select(protein_accession, sample_id, log2_abundance_normalized) %>%
      pivot_wider(names_from = sample_id, values_from = log2_abundance_normalized)

    protein_ids <- wide_df\$protein_accession
    expr_matrix <- as.matrix(wide_df[, -1, drop = FALSE])
    rownames(expr_matrix) <- protein_ids

    sample_groups <- long_df %>%
      distinct(sample_id, group) %>%
      arrange(match(sample_id, colnames(expr_matrix)))
    stopifnot(identical(sample_groups\$sample_id, colnames(expr_matrix)))

    group_levels <- unique(sample_groups\$group)
    if (length(group_levels) != 2) {
      stop(sprintf(
        "DIFFERENTIAL_ABUNDANCE expects exactly 2 groups (treatment vs control), found: %s",
        paste(group_levels, collapse = ", ")
      ))
    }

    # Convention: alphabetically first non-'control' label (or simply the
    # second level) is treated as treatment; 'control' as reference. Study
    # authors should confirm this against their sample sheet's `group` coding.
    control_label <- if ("control" %in% group_levels) "control" else group_levels[1]
    treatment_label <- setdiff(group_levels, control_label)[1]

    group_factor <- factor(sample_groups\$group, levels = c(control_label, treatment_label))
    design <- model.matrix(~group_factor)
    colnames(design) <- c("Intercept", "treatment_vs_control")

    fit <- lmFit(expr_matrix, design)
    fit <- eBayes(fit)

    tt <- topTable(fit, coef = "treatment_vs_control", number = Inf, sort.by = "none")
    tt\$protein_accession <- rownames(tt)

    out <- tt %>%
      transmute(
        protein_accession = protein_accession,
        log2FC = logFC,
        pvalue = P.Value,
        padj = adj.P.Val,
        ave_expr = AveExpr,
        t_stat = t
      )

    out_path <- sprintf("%s_%s_limma_results.tsv", study_id, comparison_id)
    write.table(out, out_path, sep = "\\t", quote = FALSE, row.names = FALSE)

    save(fit, tt, file = sprintf("%s_%s_limma.RData", study_id, comparison_id))

    r_version <- paste(R.version\$major, R.version\$minor, sep = ".")
    limma_version <- as.character(packageVersion("limma"))
    versions_lines <- c(
      "DIFFERENTIAL_ABUNDANCE:",
      sprintf("    r_base: \\"%s\\"", r_version),
      sprintf("    limma: \\"%s\\"", limma_version)
    )
    writeLines(versions_lines, "versions.yml")

    cat(sprintf("limma differential abundance complete: %d proteins tested\\n", nrow(out)))
    """

    stub:
    """
    touch ${study_id}_${comparison_id}_limma_results.tsv
    touch ${study_id}_${comparison_id}_limma.RData
    touch versions.yml
    """
}
