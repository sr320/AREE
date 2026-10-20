// Normalize feature intensities and compute QC metrics (raw_reanalysis mode
// only).
//
// Normalization (--metabolomics.normalization):
//   tic (default)      each sample's intensities are scaled so its total equals
//                      the across-sample median total, then log2(x + 1). A
//                      standard, simple untargeted-metabolomics choice.
//   median             log2, then each sample is centred on its median
//                      detected feature (plus the grand median, to keep the
//                      scale). Suits targeted panels with no usable standard.
//   internal_standard  log2, then each sample's mean log2 internal-standard
//                      area is subtracted (plus its across-sample mean). The
//                      usual choice for targeted MRM panels that spike
//                      standards into every sample.
// Every method runs separately within each value of the optional `platform`
// column, so panels from different instruments are never scaled against each
// other. Internal-standard rows (`is_internal_standard` true) are dropped from
// the output: they are not metabolites to test.
//
// Detection limit (--metabolomics.detection_limit): values at or below it are
// below detection and become missing (NaN), not a tiny intensity. Without it,
// 'tic' keeps its original behaviour; 'median' and 'internal_standard' treat
// non-positive values as missing.
//
// QC: if the sample sheet marks any samples as QC pooled samples
// (`qc_pool` column == true/1), computes the coefficient of variation (CV)
// of each feature across those pooled injections and reports the fraction of
// features passing a configurable CV threshold (--qc_cv_threshold, default
// 0.30) — a standard untargeted-LCMS QC metric. If no QC pool samples are
// present, this metric is explicitly reported as not_applicable rather than
// silently omitted.
//
// STATUS: 'tic' runs in no real study yet; 'internal_standard' with a
// detection limit is what DUPERRET2025_POMS_STARVE's targeted panel uses.

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
    val normalization
    val detection_limit

    output:
    tuple val(study_id), val(comparison_id), path("${study_id}_${comparison_id}_features_normalized.tsv"), emit: normalized_table
    path "${study_id}_${comparison_id}_normalization_qc.json", emit: normalization_qc
    path "versions.yml", emit: versions

    script:
    """
    cat <<-'EOF' > normalize_qc.py
    #!/usr/bin/env python3
    ${'"""'}
    Per-platform normalization + log2 transform of metabolite feature
    intensities, plus pooled-QC-sample CV reporting.
    ${'"""'}
    import json
    import sys

    import numpy as np
    import pandas as pd

    STUDY_ID = "${study_id}"
    COMPARISON_ID = "${comparison_id}"
    IN_PATH = "${annotated_table}"
    SAMPLE_SHEET_PATH = "${sample_sheet}"
    CV_THRESHOLD = float("${qc_cv_threshold}")
    METHOD = "${normalization}"
    DETECTION_LIMIT = None if "${detection_limit}" in ("", "null") else float("${detection_limit}")
    OUT_TABLE = "${study_id}_${comparison_id}_features_normalized.tsv"
    OUT_QC = "${study_id}_${comparison_id}_normalization_qc.json"

    if METHOD not in ("tic", "median", "internal_standard"):
        sys.exit(f"ERROR: normalization must be 'tic', 'median' or 'internal_standard', got {METHOD!r}")

    df = pd.read_csv(IN_PATH, sep="\t")
    samples = pd.read_csv(SAMPLE_SHEET_PATH, sep="\t" if "\t" in open(SAMPLE_SHEET_PATH).readline() else ",")

    meta_cols = [c for c in ["feature_id", "putative_metabolite_name", "annotation_confidence_level", "platform"]
                 if c in df.columns]
    sample_ids = samples["sample_id"].astype(str).tolist()
    intensity_cols = [c for c in sample_ids if c in df.columns]

    intensities = df[intensity_cols].astype(float)
    is_std = (df["is_internal_standard"].astype(str).str.lower().isin(["true", "1", "yes"])
              if "is_internal_standard" in df.columns else pd.Series(False, index=df.index))
    platform = df["platform"].astype(str) if "platform" in df.columns else pd.Series("all", index=df.index)

    n_below_detection = 0
    if DETECTION_LIMIT is not None:
        below = intensities <= DETECTION_LIMIT
        n_below_detection = int(below[~is_std].sum().sum())
        intensities = intensities.mask(below)

    normalized_log2 = pd.DataFrame(index=df.index, columns=intensity_cols, dtype=float)
    normalized = pd.DataFrame(index=df.index, columns=intensity_cols, dtype=float)  # linear, for the pooled-QC CV
    scaling = {}
    for plat in sorted(platform.unique()):
        rows = platform == plat
        feats = rows & ~is_std
        block = intensities.loc[rows]
        if METHOD == "tic":
            totals = intensities.loc[feats].sum(axis=0)
            factors = totals.median() / totals.replace(0, np.nan)
            scaled = block.mul(factors, axis=1)
            normalized.loc[rows] = scaled
            if DETECTION_LIMIT is None:
                normalized_log2.loc[rows] = np.log2(scaled.clip(lower=0) + 1)
            else:
                normalized_log2.loc[rows] = np.log2(scaled)
            scaling[plat] = {k: (None if pd.isna(v) else round(float(v), 6)) for k, v in factors.items()}
            continue
        logged = np.log2(block.where(block > 0))
        if METHOD == "median":
            ref = logged.loc[feats[rows]].median(axis=0)
        else:
            std_rows = is_std[rows]
            if not std_rows.any():
                sys.exit(f"ERROR: normalization 'internal_standard' but platform {plat!r} has no "
                         "row with is_internal_standard true")
            ref = logged.loc[std_rows].mean(axis=0)
            if ref.isna().any():
                sys.exit(f"ERROR: platform {plat!r}: no detected internal standard in sample(s) "
                         f"{list(ref[ref.isna()].index)}")
        offset = ref - ref.mean()
        normalized_log2.loc[rows] = logged.sub(offset, axis=1)
        normalized.loc[rows] = np.exp2(normalized_log2.loc[rows])
        scaling[plat] = {k: round(float(v), 6) for k, v in offset.items()}

    keep = ~is_std
    log2_normalized = normalized_log2.loc[keep]
    out = pd.concat([df.loc[keep, meta_cols], log2_normalized], axis=1)
    out.to_csv(OUT_TABLE, sep="\t", index=False, na_rep="NA")

    qc_pool_col = samples.get("qc_pool")
    if qc_pool_col is not None:
        qc_pool_col = qc_pool_col.fillna(False).astype(bool)
        qc_pool_samples = samples.loc[qc_pool_col, "sample_id"].astype(str).tolist()
        qc_pool_samples = [s for s in qc_pool_samples if s in log2_normalized.columns]
    else:
        qc_pool_samples = []

    if len(qc_pool_samples) >= 2:
        pool_data = normalized.loc[keep, qc_pool_samples]  # CV computed on linear scale, not log
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
        "normalization_method": {"tic": "total_ion_current_then_log2", "median": "log2_then_sample_median",
                                 "internal_standard": "log2_then_internal_standard"}[METHOD]
                                + ("_per_platform" if "platform" in df.columns else ""),
        "detection_limit": DETECTION_LIMIT,
        "n_values_below_detection": n_below_detection,
        "n_features": int(keep.sum()),
        "n_internal_standards_excluded": int(is_std.sum()),
        "internal_standards": df.loc[is_std, "feature_id"].astype(str).tolist(),
        "n_samples_normalized": len(intensity_cols),
        # tic: linear scaling factor per sample; median/internal_standard: the
        # log2 offset subtracted from each sample. Keyed by platform.
        "sample_scaling": scaling,
        "qc_pool_metric": qc_pool_metric,
    }
    with open(OUT_QC, "w") as fh:
        json.dump(qc, fh, indent=2)

    print(
        f"Normalized ({METHOD}) {len(intensity_cols)} samples x {int(keep.sum())} features for "
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
