// Map each feature's putative_metabolite_name to a metabolite_class (and,
// where the lookup table provides one, a pathway) using a lookup/crosswalk
// TSV supplied via --metabolite_class_map.
//
// NOTE ON data/mappings/: this repository's existing crosswalk
// (data/mappings/gene_id_crosswalk.tsv) is a GENE identifier crosswalk
// (NCBI/Ensembl/UniProt/locus/symbol/orthogroup) for RNA-seq/proteomics —
// it has no metabolite entries and is not applicable here. No
// metabolite-class or pathway crosswalk currently ships in
// data/mappings/, so this module requires an explicit
// --metabolite_class_map TSV (columns: putative_metabolite_name,
// metabolite_class[, pathway]). A small illustrative crosswalk is expected
// to be added under data/mappings/ (e.g. metabolite_class_map.tsv) as a
// follow-up — see workflows/metabolomics/README.md. Features with no match
// are assigned metabolite_class "unknown" rather than dropped or guessed.
//
// STATUS: real, runnable pandas logic. NOT executed in this build.

process PATHWAY_MAPPING {
    tag "${study_id}:${comparison_id}"
    label 'process_low'
    container 'python:3.11-slim'
    publishDir "${params.outdir}/metabolomics/pathway_mapping", mode: params.publish_mode

    input:
    val study_id
    val comparison_id
    path diffabundance_table
    path metabolite_class_map

    output:
    tuple val(study_id), val(comparison_id), path("${study_id}_${comparison_id}_diffabundance_classified.tsv"), emit: classified_table
    path "${study_id}_${comparison_id}_pathway_mapping_qc.json", emit: mapping_qc
    path "versions.yml", emit: versions

    script:
    // metabolite_class_map may be a real TSV path or the literal sentinel
    // "NO_FILE" (used when no crosswalk is supplied; everything falls back
    // to metabolite_class "unknown").
    """
    cat <<-'EOF' > pathway_mapping.py
    #!/usr/bin/env python3
    ${'"""'}
    Map putative_metabolite_name -> metabolite_class (and optional pathway)
    via a lookup TSV. Unmatched features get metabolite_class "unknown" —
    never silently dropped from the table.
    ${'"""'}
    import json

    import pandas as pd

    STUDY_ID = "${study_id}"
    COMPARISON_ID = "${comparison_id}"
    IN_PATH = "${diffabundance_table}"
    MAP_PATH = "${metabolite_class_map}"
    OUT_TABLE = "${study_id}_${comparison_id}_diffabundance_classified.tsv"
    OUT_QC = "${study_id}_${comparison_id}_pathway_mapping_qc.json"

    df = pd.read_csv(IN_PATH, sep="\\t")

    if "putative_metabolite_name" not in df.columns:
        raise SystemExit(
            f"ERROR: {IN_PATH} has no 'putative_metabolite_name' column — "
            f"cannot map to metabolite class/pathway."
        )

    if MAP_PATH not in ("NO_FILE", "", None):
        class_map = pd.read_csv(MAP_PATH, sep="\\t")
        if "putative_metabolite_name" not in class_map.columns or "metabolite_class" not in class_map.columns:
            raise SystemExit(
                f"ERROR: metabolite class map {MAP_PATH} must have columns "
                f"'putative_metabolite_name' and 'metabolite_class' "
                f"(optionally 'pathway'). Present columns: {list(class_map.columns)}"
            )
        merge_cols = ["putative_metabolite_name", "metabolite_class"] + (
            ["pathway"] if "pathway" in class_map.columns else []
        )
        df = df.merge(class_map[merge_cols].drop_duplicates("putative_metabolite_name"), on="putative_metabolite_name", how="left")
        n_from_map = df["metabolite_class"].notna().sum()
    else:
        df["metabolite_class"] = pd.NA
        n_from_map = 0

    if "pathway" not in df.columns:
        df["pathway"] = pd.NA

    n_unmatched = int(df["metabolite_class"].isna().sum())
    df["metabolite_class"] = df["metabolite_class"].fillna("unknown")

    df.to_csv(OUT_TABLE, sep="\\t", index=False)

    class_counts = df["metabolite_class"].value_counts().to_dict()
    qc = {
        "study_id": STUDY_ID,
        "comparison_id": COMPARISON_ID,
        "metabolite_class_map_used": MAP_PATH if MAP_PATH not in ("NO_FILE", "", None) else None,
        "n_features": len(df),
        "n_matched_via_map": int(n_from_map),
        "n_unmatched_defaulted_to_unknown": n_unmatched,
        "metabolite_class_counts": {str(k): int(v) for k, v in class_counts.items()},
    }
    with open(OUT_QC, "w") as fh:
        json.dump(qc, fh, indent=2)

    print(
        f"Pathway/class mapping complete for {STUDY_ID}:{COMPARISON_ID}: "
        f"{n_from_map} matched, {n_unmatched} defaulted to 'unknown'"
    )
    EOF

    python3 pathway_mapping.py

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        python: \$(python3 --version | sed 's/Python //')
        pandas: \$(python3 -c 'import pandas; print(pandas.__version__)')
    END_VERSIONS
    """

    stub:
    """
    touch ${study_id}_${comparison_id}_diffabundance_classified.tsv
    touch ${study_id}_${comparison_id}_pathway_mapping_qc.json
    touch versions.yml
    """
}
