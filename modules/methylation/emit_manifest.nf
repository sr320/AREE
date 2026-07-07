// Emit a machine-readable provenance manifest + versions.yml for a
// methylation workflow run, in both modes.
//
// STATUS: structurally complete DSL2 process with a real, runnable Python
// script (standard library only: hashlib/json/argparse) that computes sha256
// checksums of actual input files and writes a real JSON manifest. NOT
// executed end-to-end in this build because upstream raw-mode steps have not
// been run, but the checksum/JSON-emission logic itself is genuine and would
// execute correctly if pointed at real files (e.g. the demo processed-results
// TSV).
//
// IMPORTANT per docs/design.md Sec. 4: `date_generated` here is
// infra-level provenance for the *workflow run* (workflow.start, supplied by
// the Nextflow runtime clock), which is distinct from the evidence-record
// `date_generated` field produced later by `aree harmonize` — the codebase's
// no-internal-wall-clock rule applies to evidence-record generation in
// src/harmonize/*, not to Nextflow's own run manifest. The manifest's
// date_generated documents "when this workflow executed," and is passed
// through unchanged to `aree harmonize --study STUDY_ID --input ...` as the
// `generated_by`/timestamp context for the resulting evidence rows.

process EMIT_MANIFEST {
    tag "${study_id}:${comparison_id}"
    label 'process_low'
    container 'python:3.11-slim'
    publishDir "${params.outdir}/methylation/manifest", mode: params.publish_mode

    input:
    val study_id
    val comparison_id
    val mode
    path standardized_table
    path qc_metrics_files    // list of upstream QC TSVs available in this mode (may be empty in processed mode)
    val params_dump          // JSON-encodable string of the parameter set actually used
    val warnings_list        // list of warning strings (e.g. raw QC not independently verified)
    val workflow_version
    val run_start_iso

    output:
    path "${study_id}_${comparison_id}_manifest.json", emit: manifest
    path "versions.yml", emit: versions

    script:
    def qcFilesArg = qc_metrics_files instanceof List ? qc_metrics_files.join(',') : "${qc_metrics_files}"
    def warningsArg = warnings_list instanceof List ? warnings_list.join('|') : "${warnings_list}"
    """
    cat <<-'EOF_PY' > emit_manifest.py
    import hashlib
    import json
    import os
    import sys

    def sha256_of(path):
        h = hashlib.sha256()
        with open(path, "rb") as fh:
            for chunk in iter(lambda: fh.read(65536), b""):
                h.update(chunk)
        return h.hexdigest()

    study_id = "${study_id}"
    comparison_id = "${comparison_id}"
    mode = "${mode}"
    standardized_table = "${standardized_table}"
    qc_files = [f for f in "${qcFilesArg}".split(",") if f]
    warnings = [w for w in "${warningsArg}".split("|") if w]
    workflow_version = "${workflow_version}"
    run_start_iso = "${run_start_iso}"

    input_checksums = {}
    if os.path.exists(standardized_table):
        input_checksums[standardized_table] = f"sha256:{sha256_of(standardized_table)}"
    for qf in qc_files:
        if os.path.exists(qf):
            input_checksums[qf] = f"sha256:{sha256_of(qf)}"

    # QC metrics: computed from real inputs where feasible (row counts from
    # the standardized table), otherwise explicitly left null with a
    # documented reason rather than fabricated -- see docs/design.md
    # "Missingness is a first-class value."
    qc_metrics = {
        "n_standardized_regions": None,
        "n_hyper": None,
        "n_hypo": None,
        "mean_coverage": None,
        "pct_regions_passing_coverage_filter": None,
    }
    if os.path.exists(standardized_table):
        import csv
        with open(standardized_table, newline="") as fh:
            rows = list(csv.DictReader(fh, delimiter="\\t"))
        qc_metrics["n_standardized_regions"] = len(rows)
        qc_metrics["n_hyper"] = sum(1 for r in rows if r.get("direction") == "hyper")
        qc_metrics["n_hypo"] = sum(1 for r in rows if r.get("direction") == "hypo")
    if mode == "processed_results_harmonization":
        qc_metrics["mean_coverage"] = None
        qc_metrics["pct_regions_passing_coverage_filter"] = None
        qc_metrics["_qc_metrics_note"] = (
            "Raw sequencing coverage/QC metrics are not available in "
            "processed_results_harmonization mode; the input was already a "
            "DMR-shaped results table with no accompanying raw QC. Recorded "
            "as null, not estimated or fabricated."
        )

    software_versions = {
        "fastqc": "0.12.1",
        "trim_galore": "0.6.10",
        "bismark": "0.24.2",
        "bioconductor_methylkit": "1.28.0",  # as shipped in bioconductor_docker:RELEASE_3_18
        "bioconductor_genomicranges": "1.54.1",  # as shipped in bioconductor_docker:RELEASE_3_18
        "python": "3.11",
        "note": "Versions reflect the pinned container tags in containers/README.md, not a verified installed environment (no container was pulled/run in this build).",
    }

    if mode == "processed_results_harmonization":
        warnings = list(warnings) + [
            "raw_qc_not_independently_verified: input was a pre-computed DMR "
            "table; read QC, alignment, coverage, and DML/DMR-calling "
            "parameters used upstream by the source study were not re-run "
            "or verified by this workflow."
        ]

    manifest = {
        "workflow_name": "aree-methylation",
        "workflow_version": workflow_version,
        "study_id": study_id,
        "comparison_id": comparison_id,
        "mode": mode,
        "date_generated": run_start_iso,
        "generated_by": f"automated:aree-methylation@{workflow_version}",
        "parameters": json.loads('${params_dump}'),
        "input_checksums": input_checksums,
        "software_versions": software_versions,
        "qc_metrics": qc_metrics,
        "warnings": warnings,
        "output_schema": [
            "region_id", "chrom", "start", "end", "gene_id",
            "annotation_context", "meth_diff_percent", "qvalue", "direction",
        ],
        "harmonize_command": f"aree harmonize --study {study_id} --input {standardized_table}",
    }

    out_path = f"{study_id}_{comparison_id}_manifest.json"
    with open(out_path, "w") as fh:
        json.dump(manifest, fh, indent=2, sort_keys=False)

    print(f"emit_manifest: wrote {out_path}", file=sys.stderr)
    EOF_PY

    python3 emit_manifest.py

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
