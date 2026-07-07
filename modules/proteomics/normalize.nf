// Normalize long-format protein abundance values: log2 transform (if not
// already log-space) + median normalization across samples, a standard
// proteomics normalization approach that is robust to the large dynamic
// range and non-normality of raw intensities.
//
// STATUS: structurally complete DSL2 process with a real pandas script
// (median-of-medians centering per sample). NOT executed against real data
// in this build. raw_reanalysis mode only.

process NORMALIZE {
    tag "${study_id}:${comparison_id}"
    label 'process_low'
    container 'python:3.11-slim'
    publishDir "${params.outdir}/proteomics/normalize", mode: params.publish_mode

    input:
    tuple val(study_id), val(comparison_id), path(long_abundance_table)

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

    df = pd.read_csv("${long_abundance_table}", sep="\\t")
    required = {"protein_accession", "sample_id", "abundance", "group"}
    missing = required - set(df.columns)
    if missing:
        print(f"ERROR: normalize input missing columns: {sorted(missing)}", file=sys.stderr)
        sys.exit(1)

    df["abundance"] = pd.to_numeric(df["abundance"], errors="coerce")

    # Heuristic: if the max observed abundance is large (>100), assume linear
    # intensity space and log2-transform (adding a small pseudocount to
    # tolerate zeros); otherwise assume already log-space.
    max_val = df["abundance"].max(skipna=True)
    already_log = pd.notna(max_val) and max_val < 100
    if already_log:
        df["log2_abundance"] = df["abundance"]
        transform_applied = "none (input already appeared to be log-space)"
    else:
        pseudocount = 1.0
        df["log2_abundance"] = np.log2(df["abundance"].clip(lower=0) + pseudocount)
        transform_applied = f"log2(x + {pseudocount})"

    # Median normalization: subtract each sample's median log2 abundance,
    # then re-center on the grand median so overall scale is preserved.
    sample_medians = df.groupby("sample_id")["log2_abundance"].median()
    grand_median = sample_medians.median()
    median_offset = df["sample_id"].map(sample_medians) - grand_median
    df["log2_abundance_normalized"] = df["log2_abundance"] - median_offset

    qc = pd.DataFrame({
        "sample_id": sample_medians.index,
        "pre_norm_median_log2": sample_medians.values,
        "median_offset_applied": (sample_medians - grand_median).values,
    })
    qc["grand_median_log2"] = grand_median
    qc["transform_applied"] = transform_applied

    out_path = f"{study_id}_{comparison_id}_normalized_long.tsv"
    df.to_csv(out_path, sep="\\t", index=False)

    qc_path = f"{study_id}_{comparison_id}_normalization_qc.tsv"
    qc.to_csv(qc_path, sep="\\t", index=False)

    print(f"Normalization complete: transform={transform_applied}, n_samples={len(sample_medians)}")

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
