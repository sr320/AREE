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
    path quant_dirs, stageAs: 'quants/*'
    path tx2gene
    path sample_sheet

    output:
    tuple val(study_id), val(comparison_id), path("${study_id}_${comparison_id}_deseq2_raw.tsv"), emit: results
    path "${study_id}_${comparison_id}_deseq2.RData", emit: rdata
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
        out_tsv        = "${study_id}_${comparison_id}_deseq2_raw.tsv",
        out_rdata      = "${study_id}_${comparison_id}_deseq2.RData"
    )

    # sample_sheet.tsv is expected to have columns: sample_id, condition,
    # quant_subdir (name of the per-sample Salmon output directory), and
    # optionally covariates (batch, family_line, etc.).
    samples <- read.delim(args\$sample_sheet, stringsAsFactors = FALSE)

    # The two condition levels are named by the caller rather than assumed to be
    # literally "control"/"treatment". A real sample sheet carries the study's
    # own group labels (e.g. Midori_Control / Midori_France), and a study whose
    # groups are not named control/treatment must not silently produce an
    # all-NA factor and a misleading DESeq2 error.
    if (!"condition" %in% names(samples)) {
        stop("sample sheet has no 'condition' column: ", args\$sample_sheet)
    }
    present <- unique(samples\$condition)
    for (lvl in c(args\$control_level, args\$treatment_level)) {
        if (!(lvl %in% present)) {
            stop(sprintf(
                "condition level '%s' is not present in the sample sheet. Levels found: %s",
                lvl, paste(present, collapse = ", ")))
        }
    }

    # A sample sheet may describe a whole BioProject; this contrast uses only
    # the two groups named for it.
    samples <- samples[samples\$condition %in% c(args\$control_level, args\$treatment_level), ]
    samples\$condition <- factor(samples\$condition,
                                levels = c(args\$control_level, args\$treatment_level))

    n_per_group <- table(samples\$condition)
    if (any(n_per_group < 2)) {
        stop("each group needs at least 2 replicates; got ",
             paste(sprintf("%s=%d", names(n_per_group), n_per_group), collapse = ", "))
    }
    if (any(n_per_group < 3)) {
        warning("a group has fewer than 3 replicates; dispersion estimates will be unreliable")
    }

    # quant_subdir is optional: Salmon output directories are named for the
    # sample by default, so fall back to sample_id rather than requiring the
    # curator to duplicate the column.
    if (!"quant_subdir" %in% names(samples)) {
        samples\$quant_subdir <- samples\$sample_id
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

    txi <- tximport(quant_files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)

    dds <- DESeqDataSetFromTximport(
        txi,
        colData = samples,
        design  = ~condition
    )

    # Minimal prefiltering: drop genes with essentially no signal anywhere.
    dds <- dds[rowSums(counts(dds)) >= 10, ]

    dds <- DESeq(dds)

    res <- results(
        dds,
        contrast = c("condition", args\$treatment_level, args\$control_level),
        alpha    = 0.05
    )

    res_df <- as.data.frame(res)
    res_df\$gene_id <- rownames(res_df)

    out <- res_df[, c("gene_id", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
    out <- out[order(out\$padj, na.last = TRUE), ]

    write.table(out, args\$out_tsv, sep = "\\t", quote = FALSE, row.names = FALSE, na = "NA")
    save(dds, res, file = args\$out_rdata)

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

}
