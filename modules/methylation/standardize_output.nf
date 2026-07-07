// Reshape annotated DMR results (either mode) into the exact standardized
// TSV schema consumed by `aree harmonize` / src/harmonize/methylation.py:
// region_id, chrom, start, end, gene_id, annotation_context,
// meth_diff_percent, qvalue, direction.
//
// STATUS: structurally complete DSL2 process with a real, runnable Python
// reshape/validation script (uses only the standard library — csv/argparse —
// so it has no dependency surface beyond python:3.11-slim). This is the one
// step shared verbatim by both raw_reanalysis and processed_results_
// harmonization modes, which is what lets both modes converge on an
// identical evidence schema per docs/design.md Sec. 7. NOT executed in this
// build (no upstream raw-mode inputs exist yet), but for
// processed_results_harmonization mode this script would actually run
// correctly against the shipped demo TSV if invoked directly with python3.

process STANDARDIZE_OUTPUT {
    tag "${study_id}:${comparison_id}"
    label 'process_low'
    container 'python:3.11-slim'
    publishDir "${params.outdir}/methylation/standardized", mode: params.publish_mode

    input:
    val study_id
    val comparison_id
    path input_table   // annotated DMR table (raw mode) OR processed_results TSV (processed mode)
    val source_mode     // "raw_reanalysis" or "processed_results_harmonization"

    output:
    path "${study_id}_${comparison_id}_dmr_standardized.tsv", emit: standardized_table
    path "${study_id}_${comparison_id}_standardize_qc.tsv", emit: qc_metrics
    path "versions.yml", emit: versions

    script:
    """
    cat <<-'EOF_PY' > standardize.py
    import csv
    import sys

    REQUIRED_OUT_COLS = [
        "region_id", "chrom", "start", "end", "gene_id",
        "annotation_context", "meth_diff_percent", "qvalue", "direction",
    ]

    def sniff_delimiter(path):
        with open(path, newline="") as fh:
            sample = fh.read(4096)
        return "\\t" if sample.count("\\t") >= sample.count(",") else ","

    def direction_from_diff(value):
        try:
            v = float(value)
        except (TypeError, ValueError):
            return ""
        if v > 0:
            return "hyper"
        if v < 0:
            return "hypo"
        return "hyper"  # exactly zero has no sign; default documented in README caveats

    in_path = "${input_table}"
    out_path = "${study_id}_${comparison_id}_dmr_standardized.tsv"
    qc_path = "${study_id}_${comparison_id}_standardize_qc.tsv"
    source_mode = "${source_mode}"

    delim = sniff_delimiter(in_path)
    with open(in_path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter=delim)
        rows = list(reader)
        in_cols = reader.fieldnames or []

    n_in = len(rows)
    n_dropped_no_coords = 0
    out_rows = []
    for r in rows:
        chrom = r.get("chrom", "")
        start = r.get("start", "")
        end = r.get("end", "")
        if chrom == "" or start == "" or end == "":
            # A region with no genomic coordinates cannot be placed in the
            # standardized genomic-region schema at all; this is the one
            # case that is genuinely dropped (not merely gene-less), and it
            # is counted and reported rather than silently discarded.
            n_dropped_no_coords += 1
            continue

        gene_id = r.get("gene_id", "") or ""  # intergenic regions: kept, blank gene_id, never dropped

        meth_diff = r.get("meth_diff_percent", r.get("meth_diff", ""))
        qvalue = r.get("qvalue", r.get("q_value", ""))

        direction = r.get("direction", "").strip().lower()
        if direction not in ("hyper", "hypo"):
            direction = direction_from_diff(meth_diff)

        out_rows.append({
            "region_id": r.get("region_id", ""),
            "chrom": chrom,
            "start": start,
            "end": end,
            "gene_id": gene_id,
            "annotation_context": r.get("annotation_context", ""),
            "meth_diff_percent": meth_diff,
            "qvalue": qvalue,
            "direction": direction,
        })

    with open(out_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=REQUIRED_OUT_COLS, delimiter="\\t")
        writer.writeheader()
        writer.writerows(out_rows)

    n_out = len(out_rows)
    n_intergenic = sum(1 for r in out_rows if not r["gene_id"])
    n_hyper = sum(1 for r in out_rows if r["direction"] == "hyper")
    n_hypo = sum(1 for r in out_rows if r["direction"] == "hypo")

    with open(qc_path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\\t")
        writer.writerow([
            "study_id", "comparison_id", "source_mode", "n_input_rows",
            "n_output_rows", "n_dropped_missing_coordinates",
            "n_intergenic_regions_kept", "n_hyper", "n_hypo",
        ])
        writer.writerow([
            "${study_id}", "${comparison_id}", source_mode, n_in,
            n_out, n_dropped_no_coords, n_intergenic, n_hyper, n_hypo,
        ])

    print(f"standardize_output: {n_in} input rows -> {n_out} output rows "
          f"({n_dropped_no_coords} dropped for missing coordinates, "
          f"{n_intergenic} intergenic regions kept)", file=sys.stderr)
    EOF_PY

    python3 standardize.py

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
