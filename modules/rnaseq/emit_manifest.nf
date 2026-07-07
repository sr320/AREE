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

    output:
    path "${study_id}_${comparison_id}_manifest.json", emit: manifest
    path "versions.yml", emit: versions

    script:
    // extra_warnings is a Groovy list already rendered to a JSON array string
    // by the calling workflow (see main.nf), so it can be dropped straight
    // into the Python literal below without re-escaping.
    """
    cat <<-'EOF' > emit_manifest.py
    #!/usr/bin/env python3
    # Build the AREE RNA-seq workflow provenance manifest. Real hashlib-based
    # sha256 checksum + real JSON emission. Not executed in this build.
    import hashlib
    import json

    def sha256_of(path):
        h = hashlib.sha256()
        with open(path, "rb") as fh:
            for chunk in iter(lambda: fh.read(1 << 20), b""):
                h.update(chunk)
        return h.hexdigest()

    standardized_path = "${standardized_tsv}"

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
        },
        "software_versions": {
            # Static, representative versions matching containers/README.md
            # image tags. These are NOT introspected from a running
            # container in this build (no container runtime invoked here) --
            # they are declared to match the pinned image tags.
            "fastqc": "0.12.1",
            "fastp": "0.23.4",
            "salmon": "1.10.3",
            "star": "2.7.11b",
            "multiqc": "1.21",
            "bioconductor_docker": "RELEASE_3_18",
            "deseq2": "1.42.0",
            "tximport": "1.30.0",
            "quarto": "1.5.57",
            "python": "3.11",
        },
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
        "warnings": ${extra_warnings},
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
}
