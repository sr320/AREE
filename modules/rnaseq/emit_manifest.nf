// Emit a machine-readable provenance manifest (JSON) and a versions.yml file
// for the RNA-seq workflow run, per CLAUDE.md Layer 2 requirements ("Every
// workflow must output: a machine-readable manifest, parameters used,
// software versions, input checksums, QC metrics, warnings or failure
// states, ...").
//
// NOTE on provenance timestamps: docs/design.md section 4 states the AREE
// *Python* codebase deliberately avoids internal wall-clock calls so
// evidence-record timestamps are always caller-supplied and reproducible.
// This process is infrastructure-level (workflow-run) provenance, not an
// evidence-record field, and Nextflow's own `workflow.start` is the
// appropriate, already-reproducible source for "when did this run start" —
// it is captured once by the Nextflow engine itself, not invented here.
//
// STATUS: real, runnable Python (sha256 checksums via hashlib, real JSON
// manifest); not executed in this build.

process EMIT_MANIFEST {
    tag "${study_id}:${comparison_id}"
    label 'process_low'
    container 'python:3.11-slim'
    publishDir "${params.outdir}/rnaseq/manifests", mode: params.publish_mode

    input:
    val study_id
    val comparison_id
    val mode
    path standardized_tsv
    val workflow_start_iso
    val extra_warnings
    path viral_summary   // VIRAL_SUMMARY JSON, or assets/NO_VIRAL_QC when viral QC was not run
    // Every upstream process's versions.yml (mixed + collected in main.nf),
    // renamed on staging because they all share one file name.
    path process_versions, stageAs: 'process_versions/versions_*.yml'
    path collect_versions_py  // modules/shared/collect_versions.py

    output:
    path "${study_id}_${comparison_id}_manifest.json", emit: manifest
    path "versions.yml", emit: versions

    script:
    // software_versions is built from the versions.yml each upstream process
    // wrote, i.e. from the tools that actually ran, native or containerized.
    // It used to be a hard-coded copy of the container tags, which misreported
    // every native (-profile local) run.
    //
    // extra_warnings is a Groovy list already rendered to a JSON array string
    // by the calling workflow (see main.nf), so it can be dropped straight
    // into the Python literal below without re-escaping.
    """
    python3 ${collect_versions_py} software_versions.json process_versions/*.yml

    cat <<-'EOF' > emit_manifest.py
    #!/usr/bin/env python3
    # Build the AREE RNA-seq workflow provenance manifest.
    import hashlib
    import json
    import platform

    def sha256_of(path):
        h = hashlib.sha256()
        with open(path, "rb") as fh:
            for chunk in iter(lambda: fh.read(1 << 20), b""):
                h.update(chunk)
        return h.hexdigest()

    standardized_path = "${standardized_tsv}"

    # Curated sample exclusions, recorded as run parameters; the reasons live in
    # the study registry (comparisons[].excluded_samples).
    excluded_samples = [s.strip() for s in "${params.exclude_samples ?: ''}".split(",") if s.strip()]

    # The contrast as run: each arm may pool several condition levels, and
    # replicate_unit says whether DESeq2 was fitted on animals or on units
    # (e.g. tanks). Without these a reader cannot tell what an lfcSE means.
    def _levels(value):
        return [s.strip() for s in value.split(",") if s.strip()]
    control_levels = _levels("${params.control_level ?: ''}")
    treatment_levels = _levels("${params.treatment_level ?: ''}")
    replicate_unit = "${params.replicate_unit ?: ''}".strip() or None

    viral_path = "${viral_summary}"
    viral = None
    warnings = ${extra_warnings}
    if viral_path != "NO_VIRAL_QC":
        with open(viral_path) as fh:
            viral = json.load(fh)
        warnings = warnings + viral.get("warnings", [])

    with open("software_versions.json") as fh:
        collected = json.load(fh)
    software_versions = collected["software_versions"]
    # This process's own interpreter; its versions.yml is written after the
    # manifest, so it cannot be among the collected files.
    software_versions["EMIT_MANIFEST"] = {"python": platform.python_version()}
    warnings = warnings + collected["warnings"]
    if excluded_samples:
        warnings.append(
            "samples_excluded_by_curation: " + ", ".join(excluded_samples)
            + " (see the study registry for the reason)")

    manifest = {
        "workflow_name": "aree-rnaseq",
        "workflow_version": "${workflow.manifest.version ?: '0.1.0-scaffold'}",
        "mode": "${mode}",
        "study_id": "${study_id}",
        "comparison_id": "${comparison_id}",
        "input_files": [
            {
                "name": standardized_path,
                "sha256": sha256_of(standardized_path),
            }
        ],
        "parameters": {
            "mode": "${mode}",
            "study_id": "${study_id}",
            "comparison_id": "${comparison_id}",
            "outdir": "${params.outdir}",
            "exclude_samples": excluded_samples,
            "control_levels": control_levels,
            "treatment_levels": treatment_levels,
            "replicate_unit": replicate_unit,
            "single_end": ${params.single_end ? 'True' : 'False'},
            "se_adapter": "${params.se_adapter}" if ${params.single_end ? 'True' : 'False'} else None,
            "viral_reference": "${params.viral_reference ?: ''}" or None,
        },
        # Tools as reported by each process that ran, keyed by process name.
        "software_versions": software_versions,
        "software_versions_not_captured": {
            "RENDER_REPORT": "renders this manifest, so it runs afterwards; its versions.yml is published to rnaseq/reports/",
        },
        # Null for a native run (-profile local). The images a container
        # profile would pull are declared on each process in modules/rnaseq/;
        # inside a container, software_versions reports that image's tools.
        "container_engine": "${workflow.containerEngine ?: ''}" or None,
        "reference": {
            "genome_assembly": "${params.rnaseq?.genome_assembly ?: 'not_specified'}",
            "annotation_version": "${params.rnaseq?.annotation_version ?: 'not_specified'}",
        },
        "workflow_start": "${workflow_start_iso}",
        "qc_metrics": {
            # Placeholder structure: in raw_reanalysis mode this would be
            # populated from parsed FastQC/fastp/MultiQC/Salmon logs. Left as
            # an explicit null/placeholder rather than a fabricated number.
            "mean_reads_per_sample": None,
            "mean_mapping_rate_percent": None,
            "n_samples": None,
            "note": "QC metrics are not populated in processed_results_harmonization mode or in this unexecuted scaffold run.",
        },
        # Per-sample viral read counts when --viral_reference was given; null
        # means the check was not run, not that no virus was found.
        "viral_load": viral,
        "warnings": warnings,
    }

    with open(f"${study_id}_${comparison_id}_manifest.json", "w") as fh:
        json.dump(manifest, fh, indent=2, sort_keys=False)

    print(f"Wrote manifest for ${study_id}:${comparison_id}")
    EOF

    python3 emit_manifest.py

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    // The stub still merges the collected versions.yml files so a -stub-run
    // checks that every process's versions channel reaches the manifest.
    """
    python3 ${collect_versions_py} software_versions.json process_versions/*.yml
    python3 -c 'import json; v = json.load(open("software_versions.json")); json.dump({"study_id": "${study_id}", "comparison_id": "${comparison_id}", "mode": "${mode}", "stub": True, "software_versions": dict(v["software_versions"], EMIT_MANIFEST={}), "software_versions_not_captured": {"RENDER_REPORT": "runs after the manifest"}, "warnings": v["warnings"]}, open("${study_id}_${comparison_id}_manifest.json", "w"), indent=2)'
    echo '${task.process}:' > versions.yml
    """
}
