// Normalize feature intensities and compute QC metrics (raw_reanalysis mode
// only).
//
// Normalization: total-ion-current (TIC) normalization per sample (each
// sample's intensities are scaled so the sample's total intensity equals the
// across-sample median total intensity), followed by log2(x + 1) transform.
// This is a standard, simple untargeted-metabolomics normalization choice —
// not the only defensible one (median-fold-change / PQN are common
// alternatives) — documented here so a curator can swap it out per study.
//
// QC: if the sample sheet marks any samples as QC pooled samples
// (`qc_pool` column == true/1), computes the coefficient of variation (CV)
// of each feature across those pooled injections and reports the fraction of
// features passing a configurable CV threshold (--qc_cv_threshold, default
// 0.30) — a standard untargeted-LCMS QC metric. If no QC pool samples are
// present, this metric is explicitly reported as not_applicable rather than
// silently omitted.
//
// STATUS: real, runnable pandas logic. NOT executed in this build.

process NORMALIZE_QC {
    tag "${study_id}:${comparison_id}"
    label 'process_medium'
    container 'python:3.11-slim'
    publishDir "${params.outdir}/metabolomics/normalized", mode: params.publish_mode

    input:
    val study_id
    val comparison_id
    path annotated_table
    path sample_sheet
    val qc_cv_threshold

    output:
    tuple val(study_id), val(comparison_id), path("${study_id}_${comparison_id}_features_normalized.tsv"), emit: normalized_table
    path "${study_id}_${comparison_id}_normalization_qc.json", emit: normalization_qc
    path "versions.yml", emit: versions

    script:
    """
    cat <<-'EOF' > normalize_qc.py
    #!/usr/bin/env python3
    ${'"""'}
    TIC normalization + log2 transform of metabolite feature intensities,
    plus pooled-QC-sample CV reporting.
    ${'"""'}
    import json

    import numpy as np
    import pandas as pd

    STUDY_ID = "${study_id}"
    COMPARISON_ID = "${comparison_id}"
    IN_PATH = "${annotated_table}"
    SAMPLE_SHEET_PATH = "${sample_sheet}"
    CV_THRESHOLD = float("${qc_cv_threshold}")
    OUT_TABLE = "${study_id}_${comparison_id}_features_normalized.tsv"
    OUT_QC = "${study_id}_${comparison_id}_normalization_qc.json"

    df = pd.read_csv(IN_PATH, sep="\\t")
    samples = pd.read_csv(SAMPLE_SHEET_PATH, sep="\\t" if "\\t" in open(SAMPLE_SHEET_PATH).readline() else ",")

    meta_cols = [c for c in ["feature_id", "putative_metabolite_name", "annotation_confidence_level"] if c in df.columns]
    sample_ids = samples["sample_id"].astype(str).tolist()
    intensity_cols = [c for c in sample_ids if c in df.columns]

    intensities = df[intensity_cols].astype(float)

    # Total-ion-current normalization: scale each sample column so its sum
    # equals the median per-sample sum across all samples.
    sample_totals = intensities.sum(axis=0)
    median_total = sample_totals.median()
    scaling_factors = median_total / sample_totals.replace(0, np.nan)
    normalized = intensities.mul(scaling_factors, axis=1)

    log2_normalized = np.log2(normalized.clip(lower=0) + 1)

    out = pd.concat([df[meta_cols], log2_normalized], axis=1)
    out.to_csv(OUT_TABLE, sep="\\t", index=False)

    qc_pool_col = samples.get("qc_pool")
    if qc_pool_col is not None:
        qc_pool_col = qc_pool_col.fillna(False).astype(bool)
        qc_pool_samples = samples.loc[qc_pool_col, "sample_id"].astype(str).tolist()
        qc_pool_samples = [s for s in qc_pool_samples if s in log2_normalized.columns]
    else:
        qc_pool_samples = []

    if len(qc_pool_samples) >= 2:
        pool_data = normalized[qc_pool_samples]  # CV computed on linear scale, not log
        pool_mean = pool_data.mean(axis=1)
        pool_sd = pool_data.std(axis=1, ddof=1)
        cv = (pool_sd / pool_mean.replace(0, np.nan)).abs()
        n_pass = int((cv <= CV_THRESHOLD).sum())
        pct_pass = round(100.0 * n_pass / len(cv), 2) if len(cv) else None
        qc_pool_metric = {
            "n_qc_pool_samples": len(qc_pool_samples),
            "cv_threshold": CV_THRESHOLD,
            "n_features_passing_cv_threshold": n_pass,
            "pct_features_passing_cv_threshold": pct_pass,
        }
    else:
        qc_pool_metric = {
            "n_qc_pool_samples": len(qc_pool_samples),
            "cv_threshold": CV_THRESHOLD,
            "n_features_passing_cv_threshold": "not_applicable",
            "pct_features_passing_cv_threshold": "not_applicable",
            "note": "fewer than 2 QC pooled samples available; pooled-sample CV QC not computed",
        }

    qc = {
        "study_id": STUDY_ID,
        "comparison_id": COMPARISON_ID,
        "normalization_method": "total_ion_current_then_log2",
        "n_features": len(df),
        "n_samples_normalized": len(intensity_cols),
        "sample_scaling_factors": {k: (None if pd.isna(v) else round(float(v), 4)) for k, v in scaling_factors.items()},
        "qc_pool_metric": qc_pool_metric,
    }
    with open(OUT_QC, "w") as fh:
        json.dump(qc, fh, indent=2)

    print(
        f"Normalized {len(intensity_cols)} samples x {len(df)} features for "
        f"{STUDY_ID}:{COMPARISON_ID}; QC pool metric: {qc_pool_metric}"
    )
    EOF

    python3 normalize_qc.py

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        python: \$(python3 --version | sed 's/Python //')
        pandas: \$(python3 -c 'import pandas; print(pandas.__version__)')
        numpy: \$(python3 -c 'import numpy; print(numpy.__version__)')
    END_VERSIONS
    """

    stub:
    """
    touch ${study_id}_${comparison_id}_features_normalized.tsv
    touch ${study_id}_${comparison_id}_normalization_qc.json
    touch versions.yml
    """
}
