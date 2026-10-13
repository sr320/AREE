// Differential expression via DESeq2, from Salmon per-sample quant.sf files
// (raw_reanalysis mode only).
//
// STATUS: structurally complete DSL2 process. The R script embedded below is
// a real, syntactically valid DESeq2 analysis (tximport -> DESeqDataSet ->
// results()) that would run correctly given real quant.sf files, a real
// tx2gene map, and a real sample sheet. It has NOT been executed in this
// build (no compute / no real quant.sf inputs exist here) — see
// workflows/rnaseq/README.md "What has and has not been run".

process DIFFERENTIAL_EXPRESSION_DESEQ2 {
    tag "${study_id}:${comparison_id}"
    label 'process_medium'
    container 'bioconductor/bioconductor_docker:RELEASE_3_18'
    publishDir "${params.outdir}/rnaseq/deseq2", mode: params.publish_mode

    input:
    val study_id
    val comparison_id
    val control_level
    val treatment_level
    val exclude_samples
    val replicate_unit
    val covariates
    val tag_seq
    path quant_dirs, stageAs: 'quants/*'
    path tx2gene
    path sample_sheet

    output:
    tuple val(study_id), val(comparison_id), path("${study_id}_${comparison_id}_deseq2_raw.tsv"), emit: results
    path "${study_id}_${comparison_id}_deseq2.RData", emit: rdata
    path "${study_id}_${comparison_id}_deseq2_design.json", emit: design
    path "versions.yml", emit: versions

    script:
    // Real DESeq2 skeleton (tximport -> DESeqDataSetFromTximport -> results()).
    // Column names in the output (gene_id, baseMean, log2FoldChange, lfcSE,
    // stat, pvalue, padj) are deliberately written to already match the
    // STANDARDIZE_OUTPUT target schema so that downstream step is a
    // pass-through/validation in raw mode. NOT executed in this build.
    """
    cat <<-'EOF' > run_deseq2.R
    #!/usr/bin/env Rscript
    # AREE RNA-seq differential expression (DESeq2)
    #
    # NOTE: structurally complete, unexecuted in this build. Written to run
    # correctly against real Salmon quant.sf outputs and a real sample sheet.
    suppressPackageStartupMessages({
        library(tximport)
        library(DESeq2)
    })

    args <- list(
        sample_sheet   = "${sample_sheet}",
        tx2gene        = "${tx2gene}",
        quant_dir      = "quants",
        study_id       = "${study_id}",
        comparison_id  = "${comparison_id}",
        control_level   = "${control_level}",
        treatment_level = "${treatment_level}",
        exclude_samples = "${exclude_samples}",
        replicate_unit  = "${replicate_unit}",
        covariates      = "${covariates}",
        tag_seq         = ${tag_seq ? 'TRUE' : 'FALSE'},
        out_tsv        = "${study_id}_${comparison_id}_deseq2_raw.tsv",
        out_rdata      = "${study_id}_${comparison_id}_deseq2.RData",
        out_design     = "${study_id}_${comparison_id}_deseq2_design.json"
    )

    # sample_sheet.tsv is expected to have columns: sample_id, condition,
    # quant_subdir (name of the per-sample Salmon output directory), and
    # optionally covariate columns named by --covariates (family, batch, ...).
    samples <- read.delim(args\$sample_sheet, stringsAsFactors = FALSE)

    # The two condition levels are named by the caller rather than assumed to be
    # literally "control"/"treatment". A real sample sheet carries the study's
    # own group labels (e.g. Midori_Control / Midori_France), and a study whose
    # groups are not named control/treatment must not silently produce an
    # all-NA factor and a misleading DESeq2 error.
    if (!"condition" %in% names(samples)) {
        stop("sample sheet has no 'condition' column: ", args\$sample_sheet)
    }
    # Each arm is one `condition` level or several, comma-separated. Several
    # levels pool into one arm, e.g. three single-tank pH levels on each side
    # of a tipping point, so that no single tank decides the contrast.
    split_levels <- function(x) {
        lv <- trimws(strsplit(x, ",")[[1]])
        lv[nzchar(lv)]
    }
    control_levels   <- split_levels(args\$control_level)
    treatment_levels <- split_levels(args\$treatment_level)
    if (length(control_levels) == 0 || length(treatment_levels) == 0) {
        stop("control_level and treatment_level must each name at least one condition level")
    }
    shared <- intersect(control_levels, treatment_levels)
    if (length(shared) > 0) {
        stop("condition level(s) in both arms: ", paste(shared, collapse = ", "))
    }
    # A one-level arm keeps its own name, so single-level runs are unchanged.
    # A pooled arm gets a syntactic name joining its levels.
    arm_label <- function(lv) if (length(lv) == 1) lv else paste(lv, collapse = "_or_")
    control_arm   <- arm_label(control_levels)
    treatment_arm <- arm_label(treatment_levels)

    present <- unique(samples\$condition)
    for (lvl in c(control_levels, treatment_levels)) {
        if (!(lvl %in% present)) {
            stop(sprintf(
                "condition level '%s' is not present in the sample sheet. Levels found: %s",
                lvl, paste(present, collapse = ", ")))
        }
    }

    # Samples excluded by curation (e.g. a control that carries the pathogen),
    # named by sample_id. Each must exist in the sheet, so a typo cannot
    # silently exclude nothing. The reason is recorded in the study registry.
    excluded <- trimws(strsplit(args\$exclude_samples, ",")[[1]])
    excluded <- excluded[nzchar(excluded)]
    unknown <- setdiff(excluded, samples\$sample_id)
    if (length(unknown) > 0) {
        stop("exclude_samples names sample_id(s) not in the sample sheet: ",
             paste(unknown, collapse = ", "))
    }
    if (length(excluded) > 0) {
        message("Excluding by curation: ", paste(excluded, collapse = ", "))
        samples <- samples[!(samples\$sample_id %in% excluded), ]
    }

    # A sample sheet may describe a whole BioProject; this contrast uses only
    # the two groups named for it.
    samples <- samples[samples\$condition %in% c(control_levels, treatment_levels), ]
    samples\$source_condition <- samples\$condition
    samples\$condition <- factor(
        ifelse(samples\$condition %in% control_levels, control_arm, treatment_arm),
        levels = c(control_arm, treatment_arm))
    if (length(control_levels) > 1 || length(treatment_levels) > 1) {
        per_level <- table(samples\$source_condition)
        message("Pooled arms: ", control_arm, " vs ", treatment_arm, " (",
                paste(sprintf("%s=%d", names(per_level), per_level), collapse = ", "), ")")
    }

    # Pseudo-bulk. When treatment is applied per tank (or family, pen...),
    # animals within a unit are not independent replicates of it. Naming that
    # column sums each unit's counts into one library, so DESeq2 estimates
    # dispersion and lfcSE from unit-to-unit variation, which is the right
    # error for the treatment. Each unit must sit wholly inside one arm.
    unit_col <- trimws(args\$replicate_unit)
    pseudo_bulk <- nzchar(unit_col)
    if (pseudo_bulk) {
        if (!unit_col %in% names(samples)) {
            stop("replicate_unit column '", unit_col, "' is not in the sample sheet")
        }
        units <- as.character(samples[[unit_col]])
        if (any(is.na(units) | !nzchar(units))) {
            stop("replicate_unit column '", unit_col, "' is empty for some samples")
        }
        arms_per_unit <- tapply(as.character(samples\$condition), units, function(x) length(unique(x)))
        if (any(arms_per_unit > 1)) {
            stop("replicate_unit(s) spanning both arms: ",
                 paste(names(arms_per_unit)[arms_per_unit > 1], collapse = ", "))
        }
        samples\$replicate_unit <- units
    }

    # Covariates: sample-sheet columns fitted as factors ahead of condition,
    # e.g. family or sequencing batch (~ family + condition). The contrast is
    # still condition only. A covariate adjusts for a known source of variation;
    # it does not make correlated libraries independent. If one replicate unit
    # (tank) holds libraries at several covariate levels (both families), the
    # libraries are still clustered by unit. Use --replicate_unit, or one run
    # per covariate level, for that design (workflows/rnaseq/README.md).
    covariates_requested <- split_levels(args\$covariates)
    reserved <- c("condition", "sample_id", "source_condition", "replicate_unit", "quant_subdir")
    if (pseudo_bulk) reserved <- c(reserved, unit_col)
    clash <- intersect(covariates_requested, reserved)
    if (length(clash) > 0) {
        stop("covariates cannot include ", paste(clash, collapse = ", "),
             " (the contrast, an identifier, or the replicate unit)")
    }
    unknown_cov <- setdiff(covariates_requested, names(samples))
    if (length(unknown_cov) > 0) {
        stop("covariate column(s) not in the sample sheet: ", paste(unknown_cov, collapse = ", "))
    }
    for (cv in covariates_requested) {
        v <- as.character(samples[[cv]])
        if (any(is.na(v) | !nzchar(v))) {
            stop("covariate column '", cv, "' is empty for some samples in this contrast")
        }
        samples[[cv]] <- v
        if (pseudo_bulk) {
            levels_per_unit <- tapply(v, samples\$replicate_unit, function(x) length(unique(x)))
            if (any(levels_per_unit > 1)) {
                stop("covariate '", cv, "' varies within replicate_unit '", unit_col, "' (",
                     paste(names(levels_per_unit)[levels_per_unit > 1], collapse = ", "),
                     "); summed units cannot carry it. Run once per level instead.")
            }
        }
    }
    # A covariate with one level here (e.g. family in a single-family run)
    # carries no information; it is dropped, and the drop is recorded.
    n_levels <- vapply(covariates_requested, function(cv) length(unique(samples[[cv]])), integer(1))
    covariates_dropped <- covariates_requested[n_levels < 2]
    covariates_used    <- covariates_requested[n_levels >= 2]
    if (length(covariates_dropped) > 0) {
        message("Covariate(s) constant in this contrast, dropped: ", paste(covariates_dropped, collapse = ", "))
    }
    design_formula <- as.formula(paste("~", paste(c(covariates_used, "condition"), collapse = " + ")))
    message("DESeq2 design: ", deparse(design_formula))

    n_per_group <- if (pseudo_bulk) {
        table(unique(samples[, c("replicate_unit", "condition")])\$condition)
    } else {
        table(samples\$condition)
    }
    if (any(n_per_group < 2)) {
        stop("each group needs at least 2 replicates; got ",
             paste(sprintf("%s=%d", names(n_per_group), n_per_group), collapse = ", "))
    }
    if (any(n_per_group < 3)) {
        warning("a group has fewer than 3 replicates; dispersion estimates will be unreliable")
    }

    # quant_subdir is optional. Salmon output directories are named for the
    # read-file key, which is the run accession for archive downloads and the
    # sample_id for hand-named files, so use whichever directory exists rather
    # than requiring the curator to add the column.
    if (!"quant_subdir" %in% names(samples)) {
        samples\$quant_subdir <- samples\$sample_id
        if ("run_accession" %in% names(samples)) {
            by_run <- file.exists(file.path(args\$quant_dir, samples\$run_accession, "quant.sf"))
            samples\$quant_subdir[by_run] <- samples\$run_accession[by_run]
        }
    }

    quant_files <- file.path(args\$quant_dir, samples\$quant_subdir, "quant.sf")
    names(quant_files) <- samples\$sample_id
    missing_quant <- quant_files[!file.exists(quant_files)]
    if (length(missing_quant) > 0) {
        stop("missing Salmon quant.sf for: ", paste(names(missing_quant), collapse = ", "))
    }

    # A headerless tx2gene would otherwise lose its first transcript silently.
    tx2gene <- read.delim(args\$tx2gene, header = TRUE, stringsAsFactors = FALSE)
    if (!all(c("transcript_id", "gene_id") %in% names(tx2gene))) {
        tx2gene <- read.delim(args\$tx2gene, header = FALSE, stringsAsFactors = FALSE)
        if (ncol(tx2gene) < 2) stop("tx2gene must have at least 2 columns: ", args\$tx2gene)
        tx2gene <- tx2gene[, 1:2]
    }
    names(tx2gene)[1:2] <- c("transcript_id", "gene_id")

    # tximport(ignoreTxVersion = TRUE) strips the version suffix from the ids
    # read out of quant.sf but does NOT touch tx2gene, so a map built from a
    # versioned annotation (every NCBI GTF) silently stops matching. Strip
    # both sides so the two are compared on the same footing.
    tx2gene\$transcript_id <- sub("\\\\..*", "", tx2gene\$transcript_id)
    tx2gene <- unique(tx2gene)

    # 3' tag libraries (QuantSeq-style) give about one read per transcript, so
    # a gene's count does not scale with transcript length. Their counts are
    # used as they are, with no tximport length offset or length scaling.
    if (args\$tag_seq) message("tag_seq: DESeq2 on counts without transcript-length offsets")
    # Which length model was fitted, for the manifest: the one property that
    # decides whether these fold changes can be compared with other studies'.
    length_model <- if (args\$tag_seq) {
        "none: counts without length offsets (tag_seq)"
    } else if (pseudo_bulk) {
        "lengthScaledTPM counts summed per replicate unit"
    } else {
        "tximport average transcript length offsets"
    }

    if (pseudo_bulk) {
        # Length-scaled counts can be summed across animals; tximport's
        # per-sample length offsets cannot, so they are folded in here.
        txi <- tximport(quant_files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE,
                        countsFromAbundance = if (args\$tag_seq) "no" else "lengthScaledTPM")
        unit_counts <- t(rowsum(t(txi\$counts), group = samples\$replicate_unit))
        unit_info <- unique(samples[, c("replicate_unit", "condition", covariates_used)])
        unit_info\$n_libraries <- as.integer(table(samples\$replicate_unit)[unit_info\$replicate_unit])
        rownames(unit_info) <- unit_info\$replicate_unit
        unit_info <- unit_info[colnames(unit_counts), ]
        message("Pseudo-bulk by '", unit_col, "': ", ncol(unit_counts), " units from ",
                nrow(samples), " libraries (",
                paste(sprintf("%s:%s=%d", unit_info\$condition, unit_info\$replicate_unit,
                              unit_info\$n_libraries), collapse = ", "), ")")
        col_data <- unit_info
    } else {
        txi <- tximport(quant_files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)
        col_data <- samples
    }

    # Covariates enter as factors. Check the model before DESeq2 sees it, so a
    # covariate confounded with condition (e.g. one family per arm) fails with
    # its name rather than DESeq2's generic full-rank error.
    for (cv in covariates_used) col_data[[cv]] <- factor(col_data[[cv]])
    mm <- model.matrix(design_formula, data = as.data.frame(col_data))
    if (qr(mm)\$rank < ncol(mm)) {
        stop("design ", deparse(design_formula), " is not full rank: covariate(s) ",
             paste(covariates_used, collapse = ", "), " are confounded with condition or each other")
    }
    residual_df <- nrow(mm) - ncol(mm)
    if (residual_df < 1) {
        stop("design ", deparse(design_formula), " leaves no residual degrees of freedom (",
             nrow(mm), " libraries, ", ncol(mm), " coefficients)")
    }

    dds <- if (pseudo_bulk) {
        DESeqDataSetFromMatrix(countData = round(unit_counts), colData = col_data, design = design_formula)
    } else if (args\$tag_seq) {
        rownames(col_data) <- col_data\$sample_id
        DESeqDataSetFromMatrix(countData = round(txi\$counts), colData = col_data[colnames(txi\$counts), , drop = FALSE],
                               design = design_formula)
    } else {
        DESeqDataSetFromTximport(txi, colData = col_data, design = design_formula)
    }

    # Minimal prefiltering: drop genes with essentially no signal anywhere.
    dds <- dds[rowSums(counts(dds)) >= 10, ]

    dds <- DESeq(dds)

    res <- results(
        dds,
        contrast = c("condition", treatment_arm, control_arm),
        alpha    = 0.05
    )

    res_df <- as.data.frame(res)
    res_df\$gene_id <- rownames(res_df)

    out <- res_df[, c("gene_id", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
    out <- out[order(out\$padj, na.last = TRUE), ]

    write.table(out, args\$out_tsv, sep = "\\t", quote = FALSE, row.names = FALSE, na = "NA")
    save(dds, res, file = args\$out_rdata)

    # What was actually fitted, for the manifest. Written without jsonlite so
    # the step needs nothing beyond DESeq2 and tximport.
    json_str <- function(x) paste0('"', gsub('"', '\\\\"', x), '"')
    json_arr <- function(x) paste0("[", paste(vapply(x, json_str, ""), collapse = ", "), "]")
    writeLines(c(
        "{",
        paste0('  "design": ', json_str(paste(deparse(design_formula), collapse = "")), ","),
        paste0('  "covariates_requested": ', json_arr(covariates_requested), ","),
        paste0('  "covariates_used": ', json_arr(covariates_used), ","),
        paste0('  "covariates_dropped_constant": ', json_arr(covariates_dropped), ","),
        paste0('  "tag_seq": ', if (args\$tag_seq) "true" else "false", ","),
        paste0('  "length_model": ', json_str(length_model), ","),
        paste0('  "n_libraries_or_units": ', nrow(mm), ","),
        paste0('  "residual_df": ', residual_df),
        "}"
    ), args\$out_design)

    cat(sprintf(
        "AREE DESeq2 step complete for %s:%s — %d genes tested\\n",
        args\$study_id, args\$comparison_id, nrow(out)
    ))
    EOF

    Rscript run_deseq2.R

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        r-base: \$(Rscript -e 'cat(as.character(getRversion()))')
        bioconductor-deseq2: \$(Rscript -e 'cat(as.character(packageVersion("DESeq2")))')
        bioconductor-tximport: \$(Rscript -e 'cat(as.character(packageVersion("tximport")))')
    END_VERSIONS
    """

    stub:
    """
    printf 'gene_id\tbaseMean\tlog2FoldChange\tlfcSE\tstat\tpvalue\tpadj\nGENE1\t10\t1\t0.2\t5\t0.001\t0.01\n' > ${study_id}_${comparison_id}_deseq2_raw.tsv
    touch ${study_id}_${comparison_id}_deseq2.RData
    printf '{"design": "~condition", "covariates_requested": [], "covariates_used": [], "covariates_dropped_constant": [], "n_libraries_or_units": 0, "residual_df": 0, "stub": true}\n' > ${study_id}_${comparison_id}_deseq2_design.json
    echo '${task.process}:' > versions.yml
    """

}
