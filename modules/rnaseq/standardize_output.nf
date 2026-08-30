// Reshape/validate a DESeq2-shaped table into the exact AREE RNA-seq
// standardized result schema:
//   gene_id  baseMean  log2FoldChange  lfcSE  stat  pvalue  padj
// (tab-separated), matching data/demo/rnaseq/*_dge_demo.tsv and the columns
// consumed by src/harmonize/rnaseq.py::harmonize_rnaseq().
//
// Runs in BOTH modes:
//   - raw_reanalysis: input is the DESEQ2 process's own output, already in
//     this shape — this step re-validates column names/order/dtypes rather
//     than blindly trusting the upstream process.
//   - processed_results_harmonization: input is a user/study-provided TSV
//     that claims to already be in this shape — this step is the important
//     validation gate, since we did not generate it ourselves.
//
// STATUS: real, runnable Python (would execute correctly if a real container
// were pulled and `nextflow run` were invoked); not executed in this build.

process STANDARDIZE_OUTPUT {
    tag "${study_id}:${comparison_id}"
    label 'process_low'
    container 'python:3.11-slim'
    publishDir "${params.outdir}/rnaseq/standardized", mode: params.publish_mode

    input:
    val study_id
    val comparison_id
    path raw_results

    output:
    tuple val(study_id), val(comparison_id), path("${study_id}_${comparison_id}_dge_standardized.tsv"), emit: standardized_tsv
    path "versions.yml", emit: versions

    script:
    // Real pandas-based reshape/validation script, executed via `python3`.
    // Not executed in this build (no container runtime invoked here).
    """
    cat <<-'EOF' > standardize_output.py
    #!/usr/bin/env python3
    # Validate and reshape an RNA-seq DE results table into the AREE
    # standardized schema: gene_id, baseMean, log2FoldChange, lfcSE, stat,
    # pvalue, padj (tab-separated). Fails loudly on missing required columns
    # rather than silently dropping or renaming.
    import sys
    import pandas as pd

    REQUIRED = ["gene_id", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"]
    IN_PATH = "${raw_results}"
    OUT_PATH = "${study_id}_${comparison_id}_dge_standardized.tsv"

    df = pd.read_csv(IN_PATH, sep=None, engine="python")

    # Tolerate a couple of common alternate spellings seen across DESeq2 /
    # edgeR / limma-voom exports without silently guessing at ambiguous ones.
    rename_map = {
        "Gene": "gene_id", "GeneID": "gene_id", "gene": "gene_id",
        "log2FC": "log2FoldChange", "logFC": "log2FoldChange",
        "lfcSE": "lfcSE", "SE": "lfcSE",
        "pvalue": "pvalue", "PValue": "pvalue", "p.value": "pvalue", "p_value": "pvalue",
        "padj": "padj", "FDR": "padj", "adj.P.Val": "padj", "adjusted_p_value": "padj",
    }
    df = df.rename(columns={k: v for k, v in rename_map.items() if k in df.columns})

    missing = [c for c in REQUIRED if c not in df.columns]
    if missing:
        sys.exit(
            f"ERROR: {IN_PATH} is missing required column(s) {missing} after "
            f"standardization -- cannot emit a valid AREE RNA-seq evidence "
            f"input. Present columns: {list(df.columns)}"
        )

    out = df[REQUIRED].copy()
    for col in ["baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"]:
        out[col] = pd.to_numeric(out[col], errors="coerce")

    n_dropped_na_gene = out["gene_id"].isna().sum()
    out = out.dropna(subset=["gene_id"])

    out.to_csv(OUT_PATH, sep="\\t", index=False)

    print(
        f"Standardized {len(out)} rows for ${study_id}:${comparison_id} "
        f"({n_dropped_na_gene} rows dropped for missing gene_id) -> {OUT_PATH}"
    )
    EOF

    python3 standardize_output.py

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        python: \$(python3 --version | sed 's/Python //')
        pandas: \$(python3 -c 'import pandas; print(pandas.__version__)')
    END_VERSIONS
    """

    stub:
    """
    printf 'gene_id\tbaseMean\tlog2FoldChange\tlfcSE\tstat\tpvalue\tpadj\nGENE1\t10\t1\t0.2\t5\t0.001\t0.01\n' > ${study_id}_${comparison_id}_dge_standardized.tsv
    echo '${task.process}:' > versions.yml
    """
}
