// Reshape/validate a metabolomics results table into the exact AREE
// metabolomics standardized result schema:
//   feature_id  putative_metabolite_name  annotation_confidence_level
//   metabolite_class  log2FC  pvalue  padj
// (tab-separated), matching
// data/demo/metabolomics/GIGAS_GROW06_hypoxia_vs_control_metabolite_features_demo.tsv
// and the columns consumed by src/harmonize/metabolomics.py::harmonize_metabolomics().
//
// Runs in BOTH modes:
//   - raw_reanalysis: input is PATHWAY_MAPPING's own output, already close to
//     this shape — this step re-validates column names/order/dtypes and
//     range-checks annotation_confidence_level rather than blindly trusting
//     the upstream process.
//   - processed_results_harmonization: input is a user/study-provided TSV
//     that claims to already be in this shape — this step is the important
//     validation gate, since we did not generate it ourselves. This is
//     exactly the shape of data/demo/metabolomics/*_metabolite_features_demo.tsv.
//
// STATUS: real, runnable pandas logic (would execute correctly given a real
// container). NOT executed in this build.

process STANDARDIZE_OUTPUT {
    tag "${study_id}:${comparison_id}"
    label 'process_low'
    container 'python:3.11-slim'
    publishDir "${params.outdir}/metabolomics/standardized", mode: params.publish_mode

    input:
    val study_id
    val comparison_id
    path raw_results

    output:
    tuple val(study_id), val(comparison_id), path("${study_id}_${comparison_id}_metabolite_features_standardized.tsv"), emit: standardized_tsv
    path "versions.yml", emit: versions

    script:
    """
    cat <<-'EOF' > standardize_output.py
    #!/usr/bin/env python3
    ${'"""'}
    Validate and reshape a metabolomics differential-abundance table into the
    AREE standardized schema: feature_id, putative_metabolite_name,
    annotation_confidence_level, metabolite_class, log2FC, pvalue, padj
    (tab-separated). Fails loudly on missing required columns rather than
    silently dropping or renaming. This is the single gate both raw_reanalysis
    and processed_results_harmonization modes pass through, so downstream
    src/harmonize/metabolomics.py never needs to know which mode produced the
    file.
    ${'"""'}
    import sys

    import pandas as pd

    REQUIRED = [
        "feature_id",
        "putative_metabolite_name",
        "annotation_confidence_level",
        "metabolite_class",
        "log2FC",
        "pvalue",
        "padj",
    ]
    IN_PATH = "${raw_results}"
    OUT_PATH = "${study_id}_${comparison_id}_metabolite_features_standardized.tsv"

    df = pd.read_csv(IN_PATH, sep=None, engine="python")

    # Tolerate a couple of common alternate spellings without silently
    # guessing at ambiguous ones.
    rename_map = {
        "Feature": "feature_id", "feature": "feature_id", "FeatureID": "feature_id",
        "metabolite_name": "putative_metabolite_name", "name": "putative_metabolite_name",
        "confidence_level": "annotation_confidence_level", "msi_level": "annotation_confidence_level",
        "class": "metabolite_class", "chemical_class": "metabolite_class",
        "log2FoldChange": "log2FC", "logFC": "log2FC",
        "p.value": "pvalue", "PValue": "pvalue", "p_value": "pvalue",
        "FDR": "padj", "adj.P.Val": "padj", "adjusted_p_value": "padj",
    }
    df = df.rename(columns={k: v for k, v in rename_map.items() if k in df.columns})

    missing = [c for c in REQUIRED if c not in df.columns]
    if missing:
        sys.exit(
            f"ERROR: {IN_PATH} is missing required column(s) {missing} after "
            f"standardization — cannot emit a valid AREE metabolomics evidence "
            f"input. Present columns: {list(df.columns)}"
        )

    out = df[REQUIRED].copy()

    out["annotation_confidence_level"] = pd.to_numeric(
        out["annotation_confidence_level"], errors="coerce"
    )
    n_bad_level = int(
        (~out["annotation_confidence_level"].isin([1, 2, 3, 4])).sum()
    )
    if n_bad_level:
        print(
            f"WARNING: {n_bad_level} row(s) have an annotation_confidence_level "
            f"outside {{1,2,3,4}} or missing; defaulting them to 4 (unknown).",
            file=sys.stderr,
        )
        bad_mask = ~out["annotation_confidence_level"].isin([1, 2, 3, 4])
        out.loc[bad_mask, "annotation_confidence_level"] = 4
    out["annotation_confidence_level"] = out["annotation_confidence_level"].astype(int)

    for col in ["log2FC", "pvalue", "padj"]:
        out[col] = pd.to_numeric(out[col], errors="coerce")

    out["metabolite_class"] = out["metabolite_class"].fillna("unknown")

    n_dropped_na_id = out["feature_id"].isna().sum()
    out = out.dropna(subset=["feature_id"])

    out.to_csv(OUT_PATH, sep="\\t", index=False)

    print(
        f"Standardized {len(out)} metabolite feature rows for "
        f"${study_id}:${comparison_id} ({n_dropped_na_id} rows dropped for "
        f"missing feature_id) -> {OUT_PATH}"
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
    touch ${study_id}_${comparison_id}_metabolite_features_standardized.tsv
    touch versions.yml
    """
}
