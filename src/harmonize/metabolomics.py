"""Harmonize a metabolomics feature table into evidence records.

Metabolite identifiers are not part of the gene/protein crosswalk hierarchy;
instead, mapping_confidence is derived from the metabolomics annotation
confidence level (a Metabolomics Standards Initiative-style 1-4 scale, where
1 = confirmed against an authentic standard and 3-4 = tentative/unknown).
This mapping is intentionally conservative: only level 1-2 features get a
standardized identifier at all.
"""
from __future__ import annotations

import pandas as pd

from .schema import EVIDENCE_COLUMNS, compute_quality_flags, make_evidence_id, molecular_direction_from_effect, source_file_ref

LEVEL_TO_CONFIDENCE = {1: "exact", 2: "inferred", 3: "unresolved", 4: "unresolved"}


def harmonize_metabolomics(
    study: dict,
    comparison: dict,
    results_path,
    *,
    workflow_version: str,
    date_generated: str,
    generated_by: str,
    analysis_method: str = "untargeted_LCMS",
) -> pd.DataFrame:
    df = pd.read_csv(results_path, sep="\t")
    source_ref = source_file_ref(results_path)
    rows = []

    for _, r in df.iterrows():
        raw_id = r["feature_id"]
        level = int(r["annotation_confidence_level"]) if pd.notna(r.get("annotation_confidence_level")) else 4
        confidence = LEVEL_TO_CONFIDENCE.get(level, "unresolved")
        standardized = r.get("putative_metabolite_name") if confidence != "unresolved" else None

        effect_size = float(r["log2FC"]) if pd.notna(r.get("log2FC")) else None
        direction = molecular_direction_from_effect(effect_size)

        quality_flags = compute_quality_flags(study, comparison, confidence)
        if level >= 3 and "identifier_mapping_uncertain" not in quality_flags:
            quality_flags = sorted(set(quality_flags) | {"identifier_mapping_uncertain"})

        rows.append({
            "evidence_id": make_evidence_id(study["study_id"], comparison["comparison_id"], raw_id, analysis_method),
            "study_id": study["study_id"],
            "comparison_id": comparison["comparison_id"],
            "simulated": bool(study.get("simulated", False)),
            "feature_id_original": raw_id,
            "feature_id_standardized": standardized,
            "feature_type": "metabolite_feature",
            "orthogroup_id": None,
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
            "mapping_confidence": confidence,
            "quality_flags": quality_flags,
            "source_file": source_ref,
            "workflow_version": workflow_version,
            "date_generated": date_generated,
            "generated_by": generated_by,
        })

    return pd.DataFrame(rows, columns=EVIDENCE_COLUMNS)
