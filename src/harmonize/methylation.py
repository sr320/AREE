"""Harmonize a DML/DMR table into evidence records.

Regions annotated to a gene are resolved via the gene identifier crosswalk;
intergenic regions (no gene_id) are kept as genomic_region features with
mapping_confidence 'unresolved' rather than dropped, per the no-silent-loss
provenance requirement.
"""
from __future__ import annotations

import pandas as pd

from .identifiers import ResolvedIdentifier, resolve_identifier
from .schema import (
    EVIDENCE_COLUMNS,
    compute_quality_flags,
    make_evidence_id,
    source_file_ref,
    study_reference_fields,
)

DIRECTION_MAP = {"hyper": "hyper", "hypo": "hypo"}


def harmonize_methylation(
    study: dict,
    comparison: dict,
    results_path,
    *,
    workflow_version: str,
    date_generated: str,
    generated_by: str,
    analysis_method: str = "methylKit",
) -> pd.DataFrame:
    df = pd.read_csv(results_path, sep="\t")
    reference_fields = study_reference_fields(study)
    source_ref = source_file_ref(results_path)
    rows = []

    for _, r in df.iterrows():
        gene_id = r.get("gene_id")
        has_gene = pd.notna(gene_id) and str(gene_id).strip() != ""
        if has_gene:
            resolved = resolve_identifier(gene_id, "ncbi_gene_id")
            feature_id_original = gene_id
            feature_type = "methylation_region"
        else:
            resolved = ResolvedIdentifier(None, "unresolved", None)
            feature_id_original = r["region_id"]
            feature_type = "genomic_region"

        meth_diff = float(r["meth_diff_percent"]) if pd.notna(r.get("meth_diff_percent")) else None

        rows.append({
            "evidence_id": make_evidence_id(study["study_id"], comparison["comparison_id"], feature_id_original, analysis_method),
            "study_id": study["study_id"],
            "comparison_id": comparison["comparison_id"],
            "simulated": bool(study.get("simulated", False)),
            "feature_id_original": feature_id_original,
            "feature_id_standardized": resolved.feature_id_standardized,
            "feature_type": feature_type,
            "orthogroup_id": resolved.orthogroup_id,
            **reference_fields,
            "annotation_context": r.get("annotation_context"),
            "molecular_direction": DIRECTION_MAP.get(r.get("direction"), "ambiguous"),
            "effect_size": meth_diff,
            "effect_size_type": "methylation_diff_percent",
            "standard_error": None,
            "ci_lower": None,
            "ci_upper": None,
            "p_value": None,
            "adjusted_p_value": float(r["qvalue"]) if pd.notna(r.get("qvalue")) else None,
            "sample_size": comparison["sample_size"],
            "tissue": comparison["tissue"],
            "life_stage": comparison["life_stage"],
            "stressor": comparison["stressor_standardized"],
            "phenotype": comparison["phenotype"],
            "phenotype_direction": comparison["phenotype_direction"],
            "analysis_method": analysis_method,
            "mapping_confidence": resolved.mapping_confidence,
            "quality_flags": compute_quality_flags(study, comparison, resolved.mapping_confidence),
            "source_file": source_ref,
            "workflow_version": workflow_version,
            "date_generated": date_generated,
            "generated_by": generated_by,
        })

    return pd.DataFrame(rows, columns=EVIDENCE_COLUMNS)
