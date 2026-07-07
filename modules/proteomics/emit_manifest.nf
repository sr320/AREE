// Emit a machine-readable provenance manifest (JSON) and a versions.yml for
// the proteomics workflow run, per docs/design.md section 4 (provenance
// model) and section 7 (raw vs. processed modes).
//
// STATUS: structurally complete DSL2 process with a real Python script
// (hashlib sha256 checksums, parameter dump, qc metrics read from upstream
// JSON when available). NOT executed against real data in this build.
//
// Provenance note (docs/design.md section 4): this manifest documents
// *infra-level* provenance (what ran, with what parameters/versions/inputs,
// when). It is distinct from *evidence-level* provenance recorded per-row in
// the harmonized evidence table by src/harmonize/proteomics.py
// (source_file, workflow_version, date_generated, generated_by). The two are
// linked by workflow_version + study_id/comparison_id + input checksums, not
// merged into one record, so an evidence row can always be traced back to
// the run that produced its source file.

process EMIT_MANIFEST {
    tag "${study_id}:${comparison_id}"
    label 'process_low'
    container 'python:3.11-slim'
    publishDir "${params.outdir}/proteomics/manifest", mode: params.publish_mode

    input:
    val study_id
    val comparison_id
    val mode
    path standardized_table
    path missingness_summary_json_optional
    val workflow_version
    val start_time_iso

    output:
    path "${study_id}_${comparison_id}_manifest.json", emit: manifest
    path "versions.yml", emit: versions

    script:
    """
    #!/usr/bin/env python3
    import hashlib
    import json
    import sys
    import pandas as pd

    study_id = "${study_id}"
    comparison_id = "${comparison_id}"
    mode = "${mode}"
    workflow_version = "${workflow_version}"
    start_time_iso = "${start_time_iso}"
    missingness_summary_path = "${missingness_summary_json_optional}"

    def sha256_of(path):
        h = hashlib.sha256()
        with open(path, "rb") as fh:
            for chunk in iter(lambda: fh.read(8192), b""):
                h.update(chunk)
        return h.hexdigest()

    standardized_path = "${standardized_table}"
    checksum = sha256_of(standardized_path)

    df = pd.read_csv(standardized_path, sep="\\t")
    n_proteins = int(len(df))
    mean_missingness = None
    if "missingness_percent" in df.columns and df["missingness_percent"].notna().any():
        mean_missingness = float(df["missingness_percent"].mean())

    qc_metrics = {
        "n_proteins": n_proteins,
        "mean_missingness_percent": mean_missingness,
    }

    warnings = []
    if mean_missingness is None:
        warnings.append(
            "missingness_percent unavailable or entirely null for this comparison; "
            "source study likely did not report per-protein missingness."
        )
    if mode == "processed_results_harmonization":
        warnings.append(
            "processed_results_harmonization mode: raw QC/normalization/differential-abundance "
            "steps were NOT independently re-run or verified by this workflow; results are "
            "taken from the externally-supplied processed table as-is (see docs/design.md section 7)."
        )
    if mean_missingness is not None and mean_missingness >= 10:
        warnings.append(
            f"mean per-protein missingness ({mean_missingness:.1f}%) is >= 10%; "
            "downstream harmonize step will apply a 'processed_only'-style quality flag "
            "per src/harmonize/proteomics.py thresholding logic."
        )

    missingness_summary = None
    if missingness_summary_path and missingness_summary_path not in ("", "NO_FILE"):
        try:
            with open(missingness_summary_path) as fh:
                missingness_summary = json.load(fh)
        except FileNotFoundError:
            missingness_summary = None

    manifest = {
        "workflow_name": "aree/workflows/proteomics",
        "workflow_version": workflow_version,
        "mode": mode,
        "study_id": study_id,
        "comparison_id": comparison_id,
        "date_generated": start_time_iso,
        "generated_by": f"automated:aree-proteomics-workflow@{workflow_version}",
        "provenance_level": "infra-level (see docs/design.md section 4; distinct from per-row evidence-level provenance)",
        "inputs": {
            "standardized_table": {
                "path": standardized_path,
                "sha256": checksum,
            },
        },
        "parameters": {
            "mode": mode,
        },
        "software_versions": {
            "python": "3.11 (python:3.11-slim container)",
            "pandas": pd.__version__,
            "limma": "bioconductor RELEASE_3_18 (see containers/README.md; exact limma version pinned at container build time, not independently re-verified in this scaffold)",
        },
        "qc_metrics": qc_metrics,
        "missingness_summary": missingness_summary,
        "warnings": warnings,
    }

    manifest_path = f"{study_id}_{comparison_id}_manifest.json"
    with open(manifest_path, "w") as fh:
        json.dump(manifest, fh, indent=2)

    print(f"Wrote manifest: {manifest_path}")
    print(json.dumps(manifest, indent=2))

    with open("versions.yml", "w") as fh:
        fh.write("EMIT_MANIFEST:\\n")
        fh.write(f"    python: \\"{sys.version.split()[0]}\\"\\n")
        fh.write(f"    pandas: \\"{pd.__version__}\\"\\n")
    """

    stub:
    """
    touch ${study_id}_${comparison_id}_manifest.json
    touch versions.yml
    """
}
