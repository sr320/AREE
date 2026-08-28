"""Harmonize a protein abundance / differential-abundance table into evidence records."""
from __future__ import annotations

import pandas as pd

from .identifiers import resolve_identifier
from .schema import EVIDENCE_COLUMNS, compute_quality_flags, make_evidence_id, molecular_direction_from_effect, source_file_ref


def harmonize_proteomics(
    study: dict,
    comparison: dict,
    results_path,
    *,
    workflow_version: str,
    date_generated: str,
    generated_by: str,
    analysis_method: str = "MSstats",
) -> pd.DataFrame:
    df = pd.read_csv(results_path, sep="\t")
    source_ref = source_file_ref(results_path)
    rows = []

    for _, r in df.iterrows():
        raw_id = r["protein_accession"]
        resolved = resolve_identifier(raw_id, "uniprot_accession")
        effect_size = float(r["log2FC"]) if pd.notna(r.get("log2FC")) else None
        direction = molecular_direction_from_effect(effect_size)

        quality_flags = compute_quality_flags(study, comparison, resolved.mapping_confidence)
        missingness = r.get("missingness_percent")
        if pd.notna(missingness) and float(missingness) >= 10 and "processed_only" not in quality_flags:
            quality_flags = sorted(set(quality_flags) | {"processed_only"})

        rows.append({
            "evidence_id": make_evidence_id(study["study_id"], comparison["comparison_id"], raw_id, analysis_method),
            "study_id": study["study_id"],
            "comparison_id": comparison["comparison_id"],
            "simulated": bool(study.get("simulated", False)),
            "feature_id_original": raw_id,
            "feature_id_standardized": resolved.feature_id_standardized,
            "feature_type": "protein",
            "orthogroup_id": resolved.orthogroup_id,
            "species": study["species"],
            "genome_assembly": study["genome_assembly"],
            "annotation_version": study.get("annotation_version"),
            "annotation_context": None,
            "molecular_direction": direction,
            "effect_size": effect_size,
            "effect_size_type": "log2FoldChange",
            "standard_error": None,
            "ci_lower": None,
            "ci_upper": None,
            "p_value": float(r["pvalue"]) if pd.notna(r.get("pvalue")) else None,
            "adjusted_p_value": float(r["padj"]) if pd.notna(r.get("padj")) else None,
            "sample_size": comparison["sample_size"],
            "tissue": comparison["tissue"],
            "life_stage": comparison["life_stage"],
            "stressor": comparison["stressor_standardized"],
            "phenotype": comparison["phenotype"],
            "phenotype_direction": comparison["phenotype_direction"],
            "analysis_method": analysis_method,
            "mapping_confidence": resolved.mapping_confidence,
            "quality_flags": quality_flags,
            "source_file": source_ref,
            "workflow_version": workflow_version,
            "date_generated": date_generated,
            "generated_by": generated_by,
        })

    return pd.DataFrame(rows, columns=EVIDENCE_COLUMNS)
