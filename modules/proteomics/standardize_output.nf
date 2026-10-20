// Reshape/validate results into the exact standardized protein-level evidence
// schema consumed by `aree harmonize` (src/harmonize/proteomics.py) and
// matching data/demo/proteomics/*_protein_abundance_demo.tsv:
//
//     protein_accession  gene_symbol  log2FC  lfcSE  pvalue  padj
//     missingness_percent  analysis_method
//
// raw_reanalysis passes the id-translated limma table plus the missingness
// report; lfcSE and analysis_method come from limma. The
// processed_results_harmonization branch passes an already-shaped external
// table through the same validation; lfcSE, missingness_percent and
// analysis_method are kept when the source reports them and left null when it
// does not (never imputed).
//
// STATUS: processed mode runs in the demo CI; raw mode in the raw-mode CI
// smoke test (tests/fixtures/proteomics_raw).

process STANDARDIZE_OUTPUT {
    tag "${study_id}:${comparison_id}"
    label 'process_low'
    container 'python:3.11-slim'
    publishDir "${params.outdir}/proteomics/standardize_output", mode: params.publish_mode

    input:
    tuple val(study_id), val(comparison_id), path(input_table)
    val mode
    path missingness_per_protein_optional

    output:
    tuple val(study_id), val(comparison_id), path("${study_id}_${comparison_id}_protein_abundance_standardized.tsv"), emit: standardized
    path "${study_id}_${comparison_id}_standardize_output.log", emit: log
    path "versions.yml", emit: versions

    script:
    """
    #!/usr/bin/env python3
    import sys
    import pandas as pd

    study_id = "${study_id}"
    comparison_id = "${comparison_id}"
    mode = "${mode}"
    missingness_path = "${missingness_per_protein_optional}"

    target_cols = ["protein_accession", "gene_symbol", "log2FC", "lfcSE", "pvalue", "padj",
                   "missingness_percent", "analysis_method"]
    log_lines = []

    def log(msg):
        print(msg)
        log_lines.append(str(msg))

    log(f"AREE proteomics standardize_output: study={study_id} comparison={comparison_id} mode={mode}")

    df = pd.read_csv("${input_table}", sep="\\t")
    log(f"Loaded input table with columns: {list(df.columns)} ({len(df)} rows)")

    if mode == "raw_reanalysis":
        # Expect id_translation output: protein_accession, gene_symbol,
        # log2FC, pvalue, padj, ave_expr, t_stat, id_mapping_confidence, ...
        required = {"protein_accession", "gene_symbol", "log2FC", "lfcSE", "pvalue", "padj", "analysis_method"}
        missing = required - set(df.columns)
        if missing:
            log(f"ERROR: raw_reanalysis standardize_output missing columns {sorted(missing)}")
            sys.exit(1)

        if missingness_path and missingness_path not in ("", "NO_FILE"):
            miss_df = pd.read_csv(missingness_path, sep="\\t")
            df = df.merge(
                miss_df[["protein_accession", "missingness_percent"]],
                on="protein_accession",
                how="left",
            )
            n_missing_join = int(df["missingness_percent"].isna().sum())
            if n_missing_join:
                log(f"WARNING: {n_missing_join} proteins had no missingness_report match; leaving missingness_percent null")
        else:
            log("WARNING: no missingness_report file supplied; missingness_percent set to null")
            df["missingness_percent"] = pd.NA

    elif mode == "processed_results_harmonization":
        # Already-shaped external result table. Validate required columns
        # exist; if the source table happens to include raw per-sample
        # abundance columns we cannot assume that here (schema varies study to
        # study), so missingness_percent is either already present in the
        # source file or explicitly left null (never fabricated/imputed) —
        # this mirrors docs/design.md section 7's requirement that processed
        # mode not claim independently-verified QC.
        required = {"protein_accession", "log2FC", "pvalue", "padj"}
        missing = required - set(df.columns)
        if missing:
            log(f"ERROR: processed_results_harmonization input missing required columns {sorted(missing)}")
            sys.exit(1)
        if "gene_symbol" not in df.columns:
            log("WARNING: input has no gene_symbol column; falling back to protein_accession as placeholder")
            df["gene_symbol"] = df["protein_accession"]
        for col in ("lfcSE", "missingness_percent", "analysis_method"):
            if col not in df.columns:
                log(f"WARNING: input has no {col} column; source study did not report it -> set to null")
                df[col] = pd.NA
    else:
        log(f"ERROR: unknown mode '{mode}'")
        sys.exit(1)

    out = df[target_cols].copy()

    # Type coercion / sanity checks, without silently dropping rows.
    for col in ["log2FC", "lfcSE", "pvalue", "padj", "missingness_percent"]:
        out[col] = pd.to_numeric(out[col], errors="coerce")

    n_bad_pvalue = int(((out["pvalue"] < 0) | (out["pvalue"] > 1)).sum(skipna=True))
    if n_bad_pvalue:
        log(f"WARNING: {n_bad_pvalue} rows have pvalue outside [0,1] after coercion")

    out_path = f"{study_id}_{comparison_id}_protein_abundance_standardized.tsv"
    out.to_csv(out_path, sep="\\t", index=False)
    log(f"Wrote standardized output: {out_path} ({len(out)} rows, columns={target_cols})")
    log("This file is the intended --input for: aree harmonize --study "
        f"{study_id} --input {out_path}")

    with open(f"{study_id}_{comparison_id}_standardize_output.log", "w") as fh:
        fh.write("\\n".join(log_lines) + "\\n")

    with open("versions.yml", "w") as fh:
        fh.write("STANDARDIZE_OUTPUT:\\n")
        fh.write(f"    python: \\"{sys.version.split()[0]}\\"\\n")
        fh.write(f"    pandas: \\"{pd.__version__}\\"\\n")
    """

    stub:
    """
    touch ${study_id}_${comparison_id}_protein_abundance_standardized.tsv
    touch ${study_id}_${comparison_id}_standardize_output.log
    touch versions.yml
    """
}
