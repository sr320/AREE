// Emit a machine-readable provenance manifest (JSON) and a top-level
// versions.yml for a metabolomics workflow run, per the provenance model in
// docs/design.md section 4.
//
// Distinguishes two levels of provenance, per docs/design.md section 7:
//   - "infra-level" provenance: which files/params/containers/tool-versions
//     this Nextflow run itself used (this manifest's job).
//   - "evidence-level" provenance: the per-row source_file/workflow_version/
//     date_generated/generated_by fields stamped onto each evidence record
//     by src/harmonize/metabolomics.py when the standardized TSV is later
//     harmonized via `aree harmonize`. This manifest is a required INPUT to
//     that step's audit trail, not a replacement for it.
//
// In processed_results_harmonization mode, this manifest explicitly flags
// that upstream normalization/QC (TIC normalization, pooled-QC-sample CV)
// could NOT be independently verified, since no raw feature table passed
// through this workflow's own NORMALIZE_QC step.
//
// STATUS: real, runnable Python logic (sha256 checksums, JSON manifest).
// NOT executed in this build.

process EMIT_MANIFEST {
    tag "${study_id}:${comparison_id}"
    label 'process_low'
    container 'python:3.11-slim'
    publishDir "${params.outdir}/metabolomics/manifest", mode: params.publish_mode

    input:
    val study_id
    val comparison_id
    val mode
    path standardized_tsv
    path input_files_for_checksum
    path qc_json_files

    output:
    path "${study_id}_${comparison_id}_manifest.json", emit: manifest
    path "versions.yml", emit: versions

    script:
    // qc_json_files may be an empty list (processed_results_harmonization
    // mode never produces intake/annotation/normalization QC JSON since
    // those steps do not run) — the Python script tolerates zero, one, or
    // many QC files.
    """
    cat <<-'EOF' > emit_manifest.py
    #!/usr/bin/env python3
    ${'"""'}
    Emit AREE metabolomics workflow provenance manifest + versions.yml.
    ${'"""'}
    import glob
    import hashlib
    import json
    import os
    from datetime import datetime, timezone

    STUDY_ID = "${study_id}"
    COMPARISON_ID = "${comparison_id}"
    MODE = "${mode}"
    STANDARDIZED_TSV = "${standardized_tsv}"
    OUT_MANIFEST = f"{STUDY_ID}_{COMPARISON_ID}_manifest.json"

    def sha256_file(path):
        h = hashlib.sha256()
        with open(path, "rb") as fh:
            for chunk in iter(lambda: fh.read(8192), b""):
                h.update(chunk)
        return h.hexdigest()

    input_files = sorted(set(glob.glob("*")) - {OUT_MANIFEST, "emit_manifest.py", "versions.yml"})
    # Only checksum real, readable files staged by Nextflow for this task
    # (excludes directories and the script itself).
    checksums = {
        f: sha256_file(f)
        for f in input_files
        if os.path.isfile(f) and not f.endswith(".json") and not f.endswith(".py")
    }

    qc_metrics = {}
    warnings = []

    for qc_path in sorted(glob.glob("*_qc.json")) + sorted(glob.glob("*_qc_metrics.json")):
        try:
            with open(qc_path) as fh:
                qc_metrics[qc_path] = json.load(fh)
        except (OSError, json.JSONDecodeError) as exc:
            warnings.append(f"could not parse QC file {qc_path}: {exc}")

    n_features = None
    n_conf_1_or_2 = None
    if os.path.isfile(STANDARDIZED_TSV):
        with open(STANDARDIZED_TSV) as fh:
            header = fh.readline().rstrip("\\n").split("\\t")
            rows = [line.rstrip("\\n").split("\\t") for line in fh]
        n_features = len(rows)
        if "annotation_confidence_level" in header:
            idx = header.index("annotation_confidence_level")
            try:
                n_conf_1_or_2 = sum(1 for r in rows if len(r) > idx and r[idx] in ("1", "2"))
            except (ValueError, IndexError):
                n_conf_1_or_2 = None

    pct_conf_1_or_2 = (
        round(100.0 * n_conf_1_or_2 / n_features, 2)
        if n_features and n_conf_1_or_2 is not None
        else None
    )

    if MODE == "processed_results_harmonization":
        warnings.append(
            "processed_results_harmonization mode: upstream normalization "
            "(TIC + log2) and pooled-QC-sample CV metrics were NOT "
            "independently verified by this workflow run — they were computed "
            "(if at all) by the source study's own pipeline. See "
            "docs/design.md section 7 and quality_flags on downstream "
            "evidence records (identifier_mapping_uncertain / raw QC "
            "unverified)."
        )
    if n_features is None:
        warnings.append(f"could not read standardized output table {STANDARDIZED_TSV} to compute qc_metrics")

    manifest = {
        "workflow_name": "aree-metabolomics",
        "workflow_version": "0.1.0-scaffold",
        "mode": MODE,
        "study_id": STUDY_ID,
        "comparison_id": COMPARISON_ID,
        "date_generated": datetime.now(timezone.utc).isoformat(),
        "generated_by": "automated:aree-metabolomics@0.1.0-scaffold",
        "provenance_level": "infra-level (see docs/design.md section 4); "
                             "per-row evidence-level provenance is stamped "
                             "separately by src/harmonize/metabolomics.py",
        "parameters": {
            "mode": MODE,
        },
        "input_checksums_sha256": checksums,
        "qc_metrics": {
            "n_features": n_features,
            "n_annotation_confidence_level_1_or_2": n_conf_1_or_2,
            "pct_annotation_confidence_level_1_or_2": pct_conf_1_or_2,
            "upstream_step_qc": qc_metrics,
        },
        "warnings": warnings,
    }

    with open(OUT_MANIFEST, "w") as fh:
        json.dump(manifest, fh, indent=2, sort_keys=True, default=str)

    print(f"Wrote manifest {OUT_MANIFEST} for {STUDY_ID}:{COMPARISON_ID} (mode={MODE})")
    EOF

    python3 emit_manifest.py

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        python: \$(python3 --version | sed 's/Python //')
        workflow_mode: "${mode}"
        pipeline_software_versions:
            python_utility_layer: "python:3.11-slim"
            r_differential_abundance: "bioconductor/bioconductor_docker:RELEASE_3_18"
            report_render: "ghcr.io/quarto-dev/quarto:1.5.57"
    END_VERSIONS
    """

    stub:
    """
    touch ${study_id}_${comparison_id}_manifest.json
    touch versions.yml
    """
}
