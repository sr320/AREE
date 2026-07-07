// Compute per-protein and per-sample missingness from the normalized
// long-format abundance table. Per-protein missingness feeds directly into
// the `missingness_percent` column of the standardized evidence-ready TSV
// (see src/harmonize/proteomics.py, which treats missingness_percent >= 10
// as a trigger for the `processed_only` quality flag).
//
// STATUS: structurally complete DSL2 process with a real pandas script
// (isna()/groupby-based missingness calculation). NOT executed against real
// data in this build. raw_reanalysis mode only.

process MISSINGNESS_REPORT {
    tag "${study_id}:${comparison_id}"
    label 'process_low'
    container 'python:3.11-slim'
    publishDir "${params.outdir}/proteomics/missingness_report", mode: params.publish_mode

    input:
    tuple val(study_id), val(comparison_id), path(normalized_long_table)

    output:
    tuple val(study_id), val(comparison_id), path("${study_id}_${comparison_id}_missingness_per_protein.tsv"), emit: per_protein
    path "${study_id}_${comparison_id}_missingness_per_sample.tsv", emit: per_sample
    path "${study_id}_${comparison_id}_missingness_summary.json", emit: summary_json
    path "versions.yml", emit: versions

    script:
    """
    #!/usr/bin/env python3
    import json
    import sys
    import pandas as pd

    study_id = "${study_id}"
    comparison_id = "${comparison_id}"

    df = pd.read_csv("${normalized_long_table}", sep="\\t")
    required = {"protein_accession", "sample_id", "log2_abundance_normalized"}
    missing_cols = required - set(df.columns)
    if missing_cols:
        print(f"ERROR: missingness_report input missing columns: {sorted(missing_cols)}", file=sys.stderr)
        sys.exit(1)

    n_samples_total = df["sample_id"].nunique()

    # Per-protein missingness: fraction of samples where this protein's
    # normalized abundance is NA, expressed as a percent (matches the demo
    # schema's `missingness_percent` column, e.g. 2.0 == 2%).
    per_protein = (
        df.groupby("protein_accession")["log2_abundance_normalized"]
        .apply(lambda s: 100.0 * s.isna().mean())
        .reset_index(name="missingness_percent")
    )
    per_protein["n_samples_total"] = n_samples_total

    # Per-sample missingness: fraction of proteins missing in a given sample
    # (useful for detecting failed/low-quality runs).
    per_sample = (
        df.groupby("sample_id")["log2_abundance_normalized"]
        .apply(lambda s: 100.0 * s.isna().mean())
        .reset_index(name="missingness_percent")
    )
    n_proteins_total = df["protein_accession"].nunique()
    per_sample["n_proteins_total"] = n_proteins_total

    per_protein_path = f"{study_id}_{comparison_id}_missingness_per_protein.tsv"
    per_sample_path = f"{study_id}_{comparison_id}_missingness_per_sample.tsv"
    per_protein.to_csv(per_protein_path, sep="\\t", index=False)
    per_sample.to_csv(per_sample_path, sep="\\t", index=False)

    summary = {
        "study_id": study_id,
        "comparison_id": comparison_id,
        "n_proteins_total": int(n_proteins_total),
        "n_samples_total": int(n_samples_total),
        "mean_missingness_percent_per_protein": float(per_protein["missingness_percent"].mean()),
        "max_missingness_percent_per_protein": float(per_protein["missingness_percent"].max()),
        "n_proteins_over_10pct_missing": int((per_protein["missingness_percent"] >= 10).sum()),
    }
    with open(f"{study_id}_{comparison_id}_missingness_summary.json", "w") as fh:
        json.dump(summary, fh, indent=2)

    print(f"Missingness summary: {summary}")

    with open("versions.yml", "w") as fh:
        fh.write("MISSINGNESS_REPORT:\\n")
        fh.write(f"    python: \\"{sys.version.split()[0]}\\"\\n")
        fh.write(f"    pandas: \\"{pd.__version__}\\"\\n")
    """

    stub:
    """
    touch ${study_id}_${comparison_id}_missingness_per_protein.tsv
    touch ${study_id}_${comparison_id}_missingness_per_sample.tsv
    touch ${study_id}_${comparison_id}_missingness_summary.json
    touch versions.yml
    """
}
