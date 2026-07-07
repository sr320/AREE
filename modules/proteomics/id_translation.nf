// Translate protein accessions (e.g. UniProt-style accessions) to gene
// symbols (and, when available, other standardized identifiers) via a join
// against an identifier crosswalk table.
//
// STATUS: structurally complete DSL2 process with a real pandas merge/join
// script. NOT executed against real data in this build. The default
// `--id_mapping_table` points at the repository's illustrative demo
// crosswalk (data/mappings/gene_id_crosswalk.tsv), which is explicitly
// SYNTHETIC (see header of that file) and self-consistent only with
// data/demo/* accessions — it is not a real UniProt idmapping resource. A
// production run must supply a real crosswalk (e.g. exported from UniProt's
// idmapping service or the oyster reference annotation) via
// `--id_mapping_table`.

process ID_TRANSLATION {
    tag "${study_id}:${comparison_id}"
    label 'process_low'
    container 'python:3.11-slim'
    publishDir "${params.outdir}/proteomics/id_translation", mode: params.publish_mode

    input:
    tuple val(study_id), val(comparison_id), path(diffabund_results)
    path id_mapping_table

    output:
    tuple val(study_id), val(comparison_id), path("${study_id}_${comparison_id}_id_translated.tsv"), emit: translated
    path "${study_id}_${comparison_id}_id_translation_report.json", emit: report_json
    path "versions.yml", emit: versions

    script:
    """
    #!/usr/bin/env python3
    import json
    import sys
    import pandas as pd

    study_id = "${study_id}"
    comparison_id = "${comparison_id}"

    results = pd.read_csv("${diffabund_results}", sep="\\t")
    if "protein_accession" not in results.columns:
        print("ERROR: differential abundance results missing 'protein_accession' column", file=sys.stderr)
        sys.exit(1)

    mapping = pd.read_csv("${id_mapping_table}", sep="\\t", comment="#")
    required_map_cols = {"uniprot_accession", "gene_symbol"}
    missing_map_cols = required_map_cols - set(mapping.columns)
    if missing_map_cols:
        print(
            f"ERROR: id_mapping_table missing required columns {sorted(missing_map_cols)}; "
            f"found columns: {list(mapping.columns)}",
            file=sys.stderr,
        )
        sys.exit(1)

    mapping_small = mapping[["uniprot_accession", "gene_symbol"]].drop_duplicates()

    # Flag any accession that maps to >1 distinct gene_symbol (many-to-one from
    # the accession's perspective is fine; one-to-many from accession -> symbol
    # is a mapping-confidence concern) before the join, so it is visible in the
    # report even though the merge below will just take all matching rows.
    dup_counts = mapping_small.groupby("uniprot_accession")["gene_symbol"].nunique()
    ambiguous_accessions = sorted(dup_counts[dup_counts > 1].index.tolist())

    merged = results.merge(
        mapping_small,
        left_on="protein_accession",
        right_on="uniprot_accession",
        how="left",
    )

    unresolved_mask = merged["gene_symbol"].isna()
    n_unresolved = int(unresolved_mask.sum())
    if n_unresolved > 0:
        print(f"WARNING: {n_unresolved} / {len(merged)} protein_accession values had no gene_symbol match")
        # Preserve the row; downstream schema still needs a value, so fall back
        # to the original accession as a placeholder gene_symbol and mark it.
        merged.loc[unresolved_mask, "gene_symbol"] = merged.loc[unresolved_mask, "protein_accession"]

    merged["id_mapping_confidence"] = "exact"
    merged.loc[unresolved_mask, "id_mapping_confidence"] = "unresolved"
    merged.loc[merged["protein_accession"].isin(ambiguous_accessions), "id_mapping_confidence"] = "one_to_many_ortholog"

    merged = merged.drop(columns=["uniprot_accession"])

    out_path = f"{study_id}_{comparison_id}_id_translated.tsv"
    merged.to_csv(out_path, sep="\\t", index=False)

    report = {
        "study_id": study_id,
        "comparison_id": comparison_id,
        "id_mapping_table": "${id_mapping_table}",
        "n_input_rows": int(len(results)),
        "n_unresolved": n_unresolved,
        "n_ambiguous_accessions": len(ambiguous_accessions),
        "ambiguous_accessions": ambiguous_accessions,
    }
    with open(f"{study_id}_{comparison_id}_id_translation_report.json", "w") as fh:
        json.dump(report, fh, indent=2)

    print(f"ID translation complete: {report}")

    with open("versions.yml", "w") as fh:
        fh.write("ID_TRANSLATION:\\n")
        fh.write(f"    python: \\"{sys.version.split()[0]}\\"\\n")
        fh.write(f"    pandas: \\"{pd.__version__}\\"\\n")
    """

    stub:
    """
    touch ${study_id}_${comparison_id}_id_translated.tsv
    touch ${study_id}_${comparison_id}_id_translation_report.json
    touch versions.yml
    """
}
