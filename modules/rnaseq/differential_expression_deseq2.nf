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
        out_tsv        = "${study_id}_${comparison_id}_deseq2_raw.tsv",
        out_rdata      = "${study_id}_${comparison_id}_deseq2.RData"
    )

    # sample_sheet.tsv is expected to have columns: sample_id, condition,
    # quant_subdir (name of the per-sample Salmon output directory), and
    # optionally covariates (batch, family_line, etc.).
    samples <- read.delim(args\$sample_sheet, stringsAsFactors = FALSE)
    samples\$condition <- factor(samples\$condition, levels = c("control", "treatment"))

    quant_files <- file.path(args\$quant_dir, samples\$quant_subdir, "quant.sf")
    names(quant_files) <- samples\$sample_id
    stopifnot(all(file.exists(quant_files)))

    tx2gene <- read.delim(args\$tx2gene, header = TRUE, stringsAsFactors = FALSE)

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
        contrast = c("condition", "treatment", "control"),
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
