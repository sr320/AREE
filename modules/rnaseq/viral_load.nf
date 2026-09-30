// Per-sample viral read QC for pathogen-challenge RNA-seq (raw_reanalysis mode,
// optional: runs only when params.viral_reference is set).
//
// Counts, for every library, the read pairs that Salmon maps to a viral genome.
// It answers two questions a host-only reanalysis cannot: did each challenged
// animal actually carry the virus, and is every control free of it? In
// CALLA2026_OSHV this is how a control carrying OsHV-1 at infected-animal
// levels (DO_CT_4) and three uninfected "challenged" animals were found.
//
// The count is a QC measure, not a titre: it depends on library depth and
// composition and on how close the reference is to the challenge strain, so a
// divergent strain lowers counts in infected animals. It cannot put reads into
// an uninfected one.
//
// STATUS: executed natively (-profile local, Homebrew Salmon) on the full-depth
// CALLA2026_OSHV Midori comparisons; containers not verified.

process VIRAL_INDEX {
    tag "${viral_fasta.baseName}"
    label 'process_low'
    container 'combinelab/salmon:1.10.3'

    input:
    path viral_fasta

    output:
    path "viral_index", emit: index
    path "versions.yml", emit: versions

    script:
    """
    salmon index -t ${viral_fasta} -i viral_index -k 31 -p ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        salmon: \$(salmon --version | sed 's/salmon //')
    END_VERSIONS
    """

    stub:
    """
    mkdir viral_index
    touch viral_index/versionInfo.json
    echo '${task.process}:' > versions.yml
    """
}

process VIRAL_QUANT {
    tag "${sample_id}"
    label 'process_medium'
    container 'combinelab/salmon:1.10.3'

    input:
    tuple val(sample_id), path(reads)
    path viral_index

    output:
    path "${sample_id}.viral.tsv", emit: counts
    path "versions.yml", emit: versions

    script:
    // Every read pair is used (no subsampling), so the fraction is not biased
    // by the flowcell position of the first reads in the file.
    """
    salmon quant -i ${viral_index} -l A \\
        -1 ${reads[0]} -2 ${reads[1]} \\
        -p ${task.cpus} -o viral_quant

    processed=\$(sed -n 's/.*"num_processed": *\\([0-9]*\\).*/\\1/p' viral_quant/aux_info/meta_info.json)
    viral=\$(awk -F'\\t' 'NR>1 {s += \$5} END {printf "%.0f", s}' viral_quant/quant.sf)
    printf 'library_id\\tpairs_processed\\tviral_pairs\\n%s\\t%s\\t%s\\n' "${sample_id}" "\$processed" "\$viral" > ${sample_id}.viral.tsv

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        salmon: \$(salmon --version | sed 's/salmon //')
    END_VERSIONS
    """

    stub:
    """
    printf 'library_id\\tpairs_processed\\tviral_pairs\\n%s\\t1000\\t0\\n' "${sample_id}" > ${sample_id}.viral.tsv
    echo '${task.process}:' > versions.yml
    """
}

process VIRAL_SUMMARY {
    tag "${study_id}:${comparison_id}"
    label 'process_low'
    container 'python:3.11-slim'
    publishDir "${params.outdir}/rnaseq/viral_load", mode: params.publish_mode

    input:
    val study_id
    val comparison_id
    val control_level
    val treatment_level
    val exclude_samples
    val threshold_per_million
    path viral_fasta
    path sample_sheet
    path counts

    output:
    path "${study_id}_${comparison_id}_viral_load.tsv", emit: tsv
    path "${study_id}_${comparison_id}_viral_load.json", emit: json
    path "versions.yml", emit: versions

    script:
    """
    cat <<-'EOF' > viral_summary.py
    #!/usr/bin/env python3
    # Join per-library viral counts to the sample sheet and flag controls that
    # carry virus and challenged samples that do not.
    import csv
    import glob
    import hashlib
    import json

    control, treatment = "${control_level}", "${treatment_level}"
    threshold = float("${threshold_per_million}")
    excluded = {s.strip() for s in "${exclude_samples}".split(",") if s.strip()}

    with open("${sample_sheet}", newline="") as fh:
        sheet = list(csv.DictReader(fh, delimiter="\\t"))
    # Reads are keyed by file name (usually the run accession); the sample sheet
    # may name that key in quant_subdir, run_accession, or sample_id.
    by_key = {}
    for row in sheet:
        for col in ("quant_subdir", "run_accession", "sample_id"):
            if row.get(col):
                by_key.setdefault(row[col], row)

    rows, warnings = [], []
    for path in sorted(glob.glob("*.viral.tsv")):
        with open(path, newline="") as fh:
            rec = next(csv.DictReader(fh, delimiter="\\t"))
        meta = by_key.get(rec["library_id"], {})
        processed, viral = int(rec["pairs_processed"]), int(rec["viral_pairs"])
        per_million = viral / processed * 1e6 if processed else float("nan")
        sample_id = meta.get("sample_id", rec["library_id"])
        condition = meta.get("condition", "")
        detected = per_million >= threshold
        rows.append({
            "sample_id": sample_id, "library_id": rec["library_id"], "condition": condition,
            "pairs_processed": processed, "viral_pairs": viral,
            "viral_per_million": round(per_million, 2), "virus_detected": detected,
            "excluded_from_de": sample_id in excluded,
        })
        note = " (excluded from DE)" if sample_id in excluded else ""
        if condition == control and detected:
            warnings.append(
                f"viral_reads_in_control: {sample_id} ({condition}) has {per_million:,.1f} "
                f"viral pairs per million, at or above the {threshold:g} threshold{note}")
        if condition == treatment and not detected:
            warnings.append(
                f"no_viral_reads_in_challenged_sample: {sample_id} ({condition}) has "
                f"{per_million:,.1f} viral pairs per million, below the {threshold:g} threshold{note}")
        if not condition:
            warnings.append(f"viral_qc_unmatched_library: {rec['library_id']} is not in the sample sheet")

    cols = ["sample_id", "library_id", "condition", "pairs_processed", "viral_pairs",
            "viral_per_million", "virus_detected", "excluded_from_de"]
    with open("${study_id}_${comparison_id}_viral_load.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\\t", lineterminator="\\n")
        w.writeheader()
        w.writerows(rows)

    with open("${viral_fasta}", "rb") as fh:
        ref_sha = hashlib.sha256(fh.read()).hexdigest()
    with open("${viral_fasta}") as fh:
        ref_header = fh.readline().strip().lstrip(">")
    summary = {
        "reference": {"file": "${viral_fasta}", "sha256": ref_sha, "header": ref_header},
        "threshold_viral_pairs_per_million": threshold,
        "samples": rows,
        "warnings": warnings,
    }
    with open("${study_id}_${comparison_id}_viral_load.json", "w") as fh:
        json.dump(summary, fh, indent=2)
    print(f"{len(rows)} libraries, {len(warnings)} viral QC warning(s)")
    EOF

    python3 viral_summary.py

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    printf 'sample_id\\tlibrary_id\\tcondition\\tpairs_processed\\tviral_pairs\\tviral_per_million\\tvirus_detected\\texcluded_from_de\\n' > ${study_id}_${comparison_id}_viral_load.tsv
    printf '{"reference": {}, "threshold_viral_pairs_per_million": ${threshold_per_million}, "samples": [], "warnings": [], "stub": true}\\n' > ${study_id}_${comparison_id}_viral_load.json
    echo '${task.process}:' > versions.yml
    """
}
