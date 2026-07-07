// Assign or validate a per-feature annotation_confidence_level using a
// simplified Metabolomics Standards Initiative (MSI)-style 1-4 scale:
//   1 = confirmed by authentic chemical standard (RT + MS/MS match)
//   2 = putatively annotated (spectral library match, no standard run)
//   3 = putative class / substructure only
//   4 = unknown / unannotated feature
//
// If the input feature table already carries an `annotation_confidence_level`
// column (as produced by many XCMS/MZmine + annotation-tool pipelines, e.g.
// CAMERA, MS-DIAL, GNPS), that column is validated and passed through.
// Otherwise, a lookup table (--metabolite_annotation_map, TSV with columns
// feature_id, putative_metabolite_name, annotation_confidence_level) is
// joined in, and any feature still unmatched is conservatively assigned
// level 4 ("unknown") — never silently guessed at a stronger level.
//
// STATUS: real, runnable pandas logic. NOT executed in this build.

process ANNOTATION_CONFIDENCE {
    tag "${study_id}:${comparison_id}"
    label 'process_low'
    container 'python:3.11-slim'
    publishDir "${params.outdir}/metabolomics/annotation", mode: params.publish_mode

    input:
    val study_id
    val comparison_id
    path validated_table
    path annotation_map

    output:
    tuple val(study_id), val(comparison_id), path("${study_id}_${comparison_id}_features_annotated.tsv"), emit: annotated_table
    path "${study_id}_${comparison_id}_annotation_qc.json", emit: annotation_qc
    path "versions.yml", emit: versions

    script:
    // annotation_map may be a real TSV path or the literal sentinel
    // "NO_FILE" (used when the input table already has the two required
    // columns and no external lookup is supplied).
    """
    cat <<-'EOF' > annotation_confidence.py
    #!/usr/bin/env python3
    ${'"""'}
    Assign/validate MSI-style annotation_confidence_level (1-4) per feature.

    Precedence:
      1. Values already present in the input table win (pass-through + validate).
      2. Missing values are filled from --metabolite_annotation_map by feature_id.
      3. Anything still missing after both is conservatively set to level 4
         (unknown) with putative_metabolite_name left as the original feature_id.
    ${'"""'}
    import json

    import pandas as pd

    STUDY_ID = "${study_id}"
    COMPARISON_ID = "${comparison_id}"
    IN_PATH = "${validated_table}"
    MAP_PATH = "${annotation_map}"
    OUT_TABLE = "${study_id}_${comparison_id}_features_annotated.tsv"
    OUT_QC = "${study_id}_${comparison_id}_annotation_qc.json"
    VALID_LEVELS = {1, 2, 3, 4}

    df = pd.read_csv(IN_PATH, sep="\\t")

    if "putative_metabolite_name" not in df.columns:
        df["putative_metabolite_name"] = pd.NA
    if "annotation_confidence_level" not in df.columns:
        df["annotation_confidence_level"] = pd.NA

    n_from_input = df["annotation_confidence_level"].notna().sum()

    if MAP_PATH not in ("NO_FILE", "", None):
        amap = pd.read_csv(MAP_PATH, sep="\\t")
        required = {"feature_id", "putative_metabolite_name", "annotation_confidence_level"}
        missing_cols = required - set(amap.columns)
        if missing_cols:
            raise SystemExit(
                f"ERROR: annotation map {MAP_PATH} missing required column(s) "
                f"{sorted(missing_cols)}"
            )
        amap = amap.set_index("feature_id")
        need_fill = df["annotation_confidence_level"].isna()
        for idx in df.index[need_fill]:
            fid = df.at[idx, "feature_id"]
            if fid in amap.index:
                df.at[idx, "putative_metabolite_name"] = amap.at[fid, "putative_metabolite_name"]
                df.at[idx, "annotation_confidence_level"] = amap.at[fid, "annotation_confidence_level"]

    n_from_map = df["annotation_confidence_level"].notna().sum() - n_from_input

    still_missing = df["annotation_confidence_level"].isna()
    n_defaulted_unknown = int(still_missing.sum())
    df.loc[still_missing, "annotation_confidence_level"] = 4
    df.loc[still_missing & df["putative_metabolite_name"].isna(), "putative_metabolite_name"] = (
        "unknown_" + df.loc[still_missing & df["putative_metabolite_name"].isna(), "feature_id"].astype(str)
    )

    df["annotation_confidence_level"] = df["annotation_confidence_level"].astype(int)
    invalid_levels = ~df["annotation_confidence_level"].isin(VALID_LEVELS)
    if invalid_levels.any():
        raise SystemExit(
            f"ERROR: {int(invalid_levels.sum())} feature(s) have an "
            f"annotation_confidence_level outside {sorted(VALID_LEVELS)}"
        )

    df.to_csv(OUT_TABLE, sep="\\t", index=False)

    level_counts = df["annotation_confidence_level"].value_counts().sort_index().to_dict()
    qc = {
        "study_id": STUDY_ID,
        "comparison_id": COMPARISON_ID,
        "n_features": len(df),
        "n_confidence_from_input_table": int(n_from_input),
        "n_confidence_from_annotation_map": int(n_from_map),
        "n_defaulted_to_unknown_level_4": n_defaulted_unknown,
        "level_counts": {str(k): int(v) for k, v in level_counts.items()},
    }
    with open(OUT_QC, "w") as fh:
        json.dump(qc, fh, indent=2)

    print(
        f"Annotation confidence resolved for {len(df)} features "
        f"({STUDY_ID}:{COMPARISON_ID}); level counts: {level_counts}"
    )
    EOF

    python3 annotation_confidence.py

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        python: \$(python3 --version | sed 's/Python //')
        pandas: \$(python3 -c 'import pandas; print(pandas.__version__)')
    END_VERSIONS
    """

    stub:
    """
    touch ${study_id}_${comparison_id}_features_annotated.tsv
    touch ${study_id}_${comparison_id}_annotation_qc.json
    touch versions.yml
    """
}
