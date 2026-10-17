// Attach a gene symbol to each protein accession for the report and the
// standardized table, by a join against an identifier crosswalk.
//
// Two routes, whichever the crosswalk supports: UniProt accessions
// (uniprot_accession) and, for crosswalks built with
// src/mappings/build_crosswalk.py >= 1.1.0, legacy gene-model IDs such as the
// oyster_v9 CGI_ models (legacy_locus_tags). A ProteinProphet protein group
// ("CGI_1,CGI_2") is matched only if every member names the same gene.
//
// This is display annotation only. The identifier resolution that evidence
// records carry (feature_id_standardized, mapping_confidence) is done later by
// `aree harmonize` (src/harmonize/proteomics.py) against the study's
// registered crosswalk. The default --proteomics.id_mapping_table is the
// SYNTHETIC demo crosswalk; a real run passes the species crosswalk.
//
// STATUS: executed by the raw-mode CI smoke test (tests/fixtures/proteomics_raw).

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
        sys.exit("ERROR: differential abundance results missing 'protein_accession' column")

    mapping = pd.read_csv("${id_mapping_table}", sep="\\t", comment="#", dtype=str).fillna("")
    if "gene_symbol" not in mapping.columns:
        sys.exit(f"ERROR: id_mapping_table has no gene_symbol column; found {list(mapping.columns)}")

    # key -> set of symbols, from every route the crosswalk carries.
    routes = {"uniprot_accession": "uniprot", "uniprot_accessions_all": "uniprot", "legacy_locus_tags": "legacy_locus_tag"}
    symbol_of, route_of = {}, {}
    for col, route in routes.items():
        if col not in mapping.columns:
            continue
        for keys, sym in zip(mapping[col], mapping["gene_symbol"]):
            for key in (k for k in keys.split(";") if k):
                symbol_of.setdefault(key, set()).add(sym)
                route_of[key] = route

    def translate(acc):
        members = [m.strip() for m in str(acc).replace(";", ",").split(",") if m.strip()]
        found = [symbol_of.get(m) for m in members]
        if not members or any(f is None for f in found):
            return acc, "unresolved", None
        symbols = set().union(*found)
        if len(symbols) != 1:
            return acc, "ambiguous", None
        return symbols.pop(), "resolved", route_of[members[0]]

    translated = results["protein_accession"].map(translate)
    results["gene_symbol"] = [t[0] for t in translated]
    results["id_translation_status"] = [t[1] for t in translated]
    results["id_translation_route"] = [t[2] for t in translated]

    out_path = f"{study_id}_{comparison_id}_id_translated.tsv"
    results.to_csv(out_path, sep="\\t", index=False)
    report = {
        "study_id": study_id,
        "comparison_id": comparison_id,
        "id_mapping_table": "${id_mapping_table}",
        "routes_available": sorted({r for c, r in routes.items() if c in mapping.columns}),
        "n_input_rows": int(len(results)),
        "n_by_status": {k: int(v) for k, v in results["id_translation_status"].value_counts().items()},
        "note": "display annotation only; evidence-level identifier resolution is done by aree harmonize",
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
