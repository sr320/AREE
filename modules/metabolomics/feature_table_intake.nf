// Intake and structural validation of a raw/unfiltered metabolite feature
// intensity table (raw_reanalysis mode only) — e.g. an XCMS `featureTable`
// export or an MZmine "aligned feature list" CSV/TSV: rows = features,
// columns = samples (intensities), plus feature metadata (m/z, retention
// time, etc.).
//
// IMPORTANT SIMPLIFICATION (see workflows/metabolomics/README.md): "raw" in
// this workflow means this raw/unfiltered feature-by-sample intensity table,
// NOT raw mzML spectral data. Reprocessing mzML (peak picking, alignment,
// retention-time correction) is out of scope for AREE and is not implemented
// anywhere in this repository — a lab would run XCMS/MZmine/etc. upstream and
// feed the resulting feature table in here.
//
// STATUS: real, runnable pandas logic (would execute correctly given a real
// python:3.11-slim container and a real feature table). NOT executed in this
// build — no compute / no real feature-table fixture exists here; the demo
// config only exercises processed_results_harmonization mode (see
// config/demo.config header).

process FEATURE_TABLE_INTAKE {
    tag "${study_id}:${comparison_id}"
    label 'process_low'
    container 'python:3.11-slim'
    publishDir "${params.outdir}/metabolomics/intake", mode: params.publish_mode

    input:
    val study_id
    val comparison_id
    path raw_feature_table
    path sample_sheet

    output:
    tuple val(study_id), val(comparison_id), path("${study_id}_${comparison_id}_features_validated.tsv"), emit: validated_table
    path "${study_id}_${comparison_id}_intake_qc.json", emit: intake_qc
    path "versions.yml", emit: versions

    script:
    // Real pandas structural validation:
    //  - reads the feature table (tab or comma separated, auto-detected)
    //  - reads the sample sheet (sample_id, condition[, qc_pool] columns)
    //  - checks every sample_id in the sample sheet has a matching intensity
    //    column in the feature table
    //  - checks feature_id is present and unique
    //  - coerces intensity columns to numeric, reporting how many values
    //    were non-numeric/coerced to NaN
    //  - fails loudly (non-zero exit) rather than silently dropping columns
    """
    cat <<-'EOF' > feature_table_intake.py
    #!/usr/bin/env python3
    ${'"""'}
    AREE metabolomics feature-table intake.

    Validates the structural contract of a raw/unfiltered feature-by-sample
    intensity table before any normalization or statistics are applied.
    ${'"""'}
    import json
    import sys

    import pandas as pd

    STUDY_ID = "${study_id}"
    COMPARISON_ID = "${comparison_id}"
    FEATURE_TABLE_PATH = "${raw_feature_table}"
    SAMPLE_SHEET_PATH = "${sample_sheet}"
    OUT_TABLE = "${study_id}_${comparison_id}_features_validated.tsv"
    OUT_QC = "${study_id}_${comparison_id}_intake_qc.json"

    def sniff_sep(path):
        with open(path) as fh:
            header = fh.readline()
        return "\\t" if header.count("\\t") >= header.count(",") else ","

    features = pd.read_csv(FEATURE_TABLE_PATH, sep=sniff_sep(FEATURE_TABLE_PATH))
    samples = pd.read_csv(SAMPLE_SHEET_PATH, sep=sniff_sep(SAMPLE_SHEET_PATH))

    warnings = []

    if "feature_id" not in features.columns:
        sys.exit(
            f"ERROR: {FEATURE_TABLE_PATH} has no 'feature_id' column — cannot "
            f"identify metabolite features. Present columns: {list(features.columns)}"
        )

    n_before = len(features)
    dup_mask = features["feature_id"].duplicated(keep=False)
    if dup_mask.any():
        warnings.append(
            f"{dup_mask.sum()} duplicate feature_id value(s) found; keeping first "
            f"occurrence of each"
        )
        features = features[~features["feature_id"].duplicated(keep='first')]

    required_sample_cols = {"sample_id", "condition"}
    missing_sheet_cols = required_sample_cols - set(samples.columns)
    if missing_sheet_cols:
        sys.exit(
            f"ERROR: sample sheet {SAMPLE_SHEET_PATH} is missing required "
            f"column(s) {sorted(missing_sheet_cols)}"
        )

    sample_ids = samples["sample_id"].astype(str).tolist()
    missing_intensity_cols = [s for s in sample_ids if s not in features.columns]
    if missing_intensity_cols:
        sys.exit(
            f"ERROR: sample sheet lists sample_id(s) {missing_intensity_cols} "
            f"with no matching intensity column in {FEATURE_TABLE_PATH}. "
            f"Present columns: {list(features.columns)}"
        )

    n_coerced = 0
    for col in sample_ids:
        before_na = features[col].isna().sum()
        features[col] = pd.to_numeric(features[col], errors="coerce")
        after_na = features[col].isna().sum()
        n_coerced += max(after_na - before_na, 0)

    if n_coerced:
        warnings.append(f"{n_coerced} non-numeric intensity value(s) coerced to NaN")

    is_qc_pool = samples.get("qc_pool")
    n_qc_pools = int(is_qc_pool.fillna(False).astype(bool).sum()) if is_qc_pool is not None else 0

    features.to_csv(OUT_TABLE, sep="\\t", index=False)

    qc = {
        "study_id": STUDY_ID,
        "comparison_id": COMPARISON_ID,
        "n_features_in": n_before,
        "n_features_out": len(features),
        "n_samples": len(sample_ids),
        "n_qc_pool_samples": n_qc_pools,
        "n_intensity_values_coerced_to_nan": int(n_coerced),
        "warnings": warnings,
    }
    with open(OUT_QC, "w") as fh:
        json.dump(qc, fh, indent=2)

    print(
        f"Intake validated {len(features)} features x {len(sample_ids)} samples "
        f"for {STUDY_ID}:{COMPARISON_ID} ({len(warnings)} warning(s))"
    )
    EOF

    python3 feature_table_intake.py

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        python: \$(python3 --version | sed 's/Python //')
        pandas: \$(python3 -c 'import pandas; print(pandas.__version__)')
    END_VERSIONS
    """

    stub:
    """
    touch ${study_id}_${comparison_id}_features_validated.tsv
    touch ${study_id}_${comparison_id}_intake_qc.json
    touch versions.yml
    """
}
