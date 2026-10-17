// Put protein abundances on a log2 scale and (optionally) median-normalize
// samples, raw_reanalysis mode only.
//
// --proteomics.input_scale
//     auto   : the original heuristic. Max value < 100 is taken as already
//              log2; otherwise log2(x + 1). Wrong for spectral-count ratios
//              such as NSAF (all values < 1, but linear), so set it explicitly
//              for those.
//     linear : log2(x). Zeros and values at or below detection_limit become
//              missing (not detected), never log2 of a pseudocount.
//     log2   : used as is.
// --proteomics.detection_limit
//     Values at or below it are missing. Defaults to 0 for linear input, so a
//     zero spectral count means "not detected", not "very low".
// --proteomics.normalization
//     median : subtract each sample's median log2 abundance (over detected
//              proteins), re-centred on the grand median.
//     none   : leave as is, e.g. for NSAF, which is already normalized within
//              each sample.
//
// STATUS: executed by the raw-mode CI smoke test (tests/fixtures/proteomics_raw).

process NORMALIZE {
    tag "${study_id}:${comparison_id}"
    label 'process_low'
    container 'python:3.11-slim'
    publishDir "${params.outdir}/proteomics/normalize", mode: params.publish_mode

    input:
    tuple val(study_id), val(comparison_id), path(long_abundance_table)
    val input_scale
    val detection_limit
    val normalization

    output:
    tuple val(study_id), val(comparison_id), path("${study_id}_${comparison_id}_normalized_long.tsv"), emit: normalized_long
    path "${study_id}_${comparison_id}_normalization_qc.tsv", emit: qc
    path "versions.yml", emit: versions

    script:
    """
    #!/usr/bin/env python3
    import sys
    import numpy as np
    import pandas as pd

    study_id = "${study_id}"
    comparison_id = "${comparison_id}"
    input_scale = "${input_scale}"
    detection_limit = "${detection_limit}"
    normalization = "${normalization}"

    if input_scale not in ("auto", "linear", "log2"):
        sys.exit(f"ERROR: proteomics.input_scale must be auto, linear or log2, got {input_scale!r}")
    if normalization not in ("median", "none"):
        sys.exit(f"ERROR: proteomics.normalization must be median or none, got {normalization!r}")

    df = pd.read_csv("${long_abundance_table}", sep="\\t")
    required = {"protein_accession", "sample_id", "abundance", "group"}
    missing = required - set(df.columns)
    if missing:
        sys.exit(f"ERROR: normalize input missing columns: {sorted(missing)}")
    df["abundance"] = pd.to_numeric(df["abundance"], errors="coerce")

    if input_scale == "auto":
        max_val = df["abundance"].max(skipna=True)
        input_scale = "log2" if pd.notna(max_val) and max_val < 100 else "auto_linear_pseudocount"

    limit = None if detection_limit in ("", "null", "None") else float(detection_limit)
    if input_scale == "linear":
        limit = 0.0 if limit is None else limit
        detected = df["abundance"] > limit
        df["log2_abundance"] = np.log2(df["abundance"].where(detected))
        transform_applied = f"log2(x); x <= {limit} treated as not detected (missing)"
    elif input_scale == "log2":
        df["log2_abundance"] = df["abundance"]
        if limit is not None:
            df.loc[df["abundance"] <= limit, "log2_abundance"] = np.nan
        transform_applied = "none (input declared or detected as log2)"
    else:
        df["log2_abundance"] = np.log2(df["abundance"].clip(lower=0) + 1.0)
        transform_applied = "log2(x + 1) (auto: input looked linear)"

    sample_medians = df.groupby("sample_id")["log2_abundance"].median()
    grand_median = sample_medians.median()
    if normalization == "median":
        offset = sample_medians - grand_median
    else:
        offset = sample_medians * 0.0
    df["log2_abundance_normalized"] = df["log2_abundance"] - df["sample_id"].map(offset)

    qc = pd.DataFrame({
        "sample_id": sample_medians.index,
        "n_detected": df.groupby("sample_id")["log2_abundance"].count().reindex(sample_medians.index).values,
        "pre_norm_median_log2": sample_medians.values,
        "offset_applied": offset.values,
    })
    qc["grand_median_log2"] = grand_median
    qc["transform_applied"] = transform_applied
    qc["normalization"] = normalization

    df.to_csv(f"{study_id}_{comparison_id}_normalized_long.tsv", sep="\\t", index=False)
    qc.to_csv(f"{study_id}_{comparison_id}_normalization_qc.tsv", sep="\\t", index=False)
    print(f"Normalization complete: transform={transform_applied}, normalization={normalization}, n_samples={len(sample_medians)}")

    with open("versions.yml", "w") as fh:
        fh.write("NORMALIZE:\\n")
        fh.write(f"    python: \\"{sys.version.split()[0]}\\"\\n")
        fh.write(f"    pandas: \\"{pd.__version__}\\"\\n")
        fh.write(f"    numpy: \\"{np.__version__}\\"\\n")
    """

    stub:
    """
    touch ${study_id}_${comparison_id}_normalized_long.tsv
    touch ${study_id}_${comparison_id}_normalization_qc.tsv
    touch versions.yml
    """
}
