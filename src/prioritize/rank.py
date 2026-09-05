"""Candidate tiering: hard gates first, score second.

A candidate's numeric score never promotes it past these gates — this is the
mechanism that keeps "significant in one study" from reading as "validated"
(see docs/design.md section 6).

Ranking is a single grouped pass over the evidence table. Each candidate's
evidence subset is looked up from a precomputed index instead of re-filtering
the whole table per candidate, which is what makes a genome-wide pool of
30,000+ candidates rank in seconds rather than tens of minutes.
"""
from __future__ import annotations

import ast
from functools import lru_cache

import pandas as pd

from common import load_vocab

from .scoring import (
    candidate_score,
    compute_components,
    phenotype_relevance,
    significance_p_value,
)

HIGH_PRIORITY_MIN_STUDIES = 2
HIGH_PRIORITY_MIN_DIRECTION_CONSISTENCY = 0.7
HIGH_PRIORITY_MIN_QUALITY_SCORE = 0.4
# Two genome-wide studies agree in sign for roughly half of all null genes, so
# study count and direction consistency alone cannot define the top tier. The
# pooled effect must also survive family-wise false-discovery control.
HIGH_PRIORITY_MAX_ADJUSTED_P = 0.05
# Multi-omics convergence counts molecular *layers* (transcriptomics, DNA
# methylation, proteomics, ...) that each carry a significant record for the
# same feature under the SAME phenotype. A gene differentially expressed under
# heat and differentially methylated under pathogen challenge is two
# observations about two questions, not convergent evidence for either.
MULTI_OMICS_MIN_LAYERS = 2
LAYER_SUPPORT_MAX_ADJUSTED_P = 0.05

# Order in which tiers are presented: strongest evidence first.
TIER_ORDER = ["high_priority_cross_study", "multi_omics_convergence", "emerging"]

# What KIND of evidence supports a candidate, independent of how strong it is.
# The tier answers "how well replicated?"; this answers "replicated evidence of
# what?" — and they are not the same question. Two studies agreeing that a gene
# responds to OsHV-1 infection is well-replicated *disease-response* evidence;
# it is not evidence that the gene marks a resilient animal, because neither
# study measured survival. Keeping these orthogonal is what stops a strong tier
# from being read as a resilience claim (see docs/resilience_vs_exposure.md).
EVIDENCE_CLASS_BY_RELEVANCE = {
    "resilience": "resilience_associated",
    "disease": "disease_associated",
    "stress_response": "stress_response",
    "exposure_only": "exposure_only",
}
DEFAULT_EVIDENCE_CLASS = "exposure_only"

PARTITION_COLUMNS = ["_simulated", "_species_taxid"]
CANDIDATE_GROUP_COLUMNS = ["feature_id_standardized", "phenotype", "feature_type", *PARTITION_COLUMNS]
FEATURE_GROUP_COLUMNS = ["feature_id_standardized", *PARTITION_COLUMNS]


def _parse_flags(value) -> list:
    if isinstance(value, list):
        return value
    if isinstance(value, str) and value.startswith("["):
        return ast.literal_eval(value)
    return []


def _partition_value(value) -> str:
    if isinstance(value, bool):
        return str(value).lower()
    return str(value).strip().lower()


def add_partition_columns(evidence_df: pd.DataFrame) -> pd.DataFrame:
    """Return the evidence with normalized `_simulated` / `_species_taxid` columns.

    `simulated` arrives as a bool in memory but as the text "True"/"False" from a
    TSV, and `species_taxid` as int or text; normalizing once here is what lets
    every downstream grouping key be compared by plain equality. Idempotent.
    """
    if all(col in evidence_df.columns for col in PARTITION_COLUMNS):
        return evidence_df
    return evidence_df.assign(
        _simulated=evidence_df["simulated"].map(_partition_value),
        _species_taxid=evidence_df["species_taxid"].map(_partition_value),
    )


@lru_cache(maxsize=1)
def molecular_layers() -> dict:
    """feature_type -> molecular_layer, from the feature_types vocabulary."""
    return {
        term_id: term.get("molecular_layer", term_id)
        for term_id, term in load_vocab("feature_types").items()
    }


def layer_of(feature_type: str) -> str:
    return molecular_layers().get(feature_type, str(feature_type))


def evidence_class(meta_row: dict) -> str:
    """What kind of evidence this candidate rests on, from its phenotype."""
    return EVIDENCE_CLASS_BY_RELEVANCE.get(
        phenotype_relevance(meta_row["phenotype"]), DEFAULT_EVIDENCE_CLASS
    )


def pooled_stressors(meta_row: dict) -> list:
    """The standardized stressors behind the pooled estimate.

    `distinct_stressors` holds pipe-separated names, not a count, unlike its
    two similarly-named siblings — see docs/interpreting_meta_analysis.md.
    """
    return sorted({s for s in str(meta_row.get("distinct_stressors") or "").split("|") if s})


def context_replication(meta_row: dict) -> str:
    """Whether independent studies replicate a result across biological contexts.

    `k_studies` counts studies; this counts *contexts*. Two studies of the same
    stressor in the same tissue at the same life stage replicate the assay, not
    the biology's generality — a distinction that disappears if only the study
    count is reported.
    """
    if meta_row["k_studies"] < HIGH_PRIORITY_MIN_STUDIES:
        return "single_study"
    varies = (
        len(pooled_stressors(meta_row)) > 1
        or meta_row["distinct_tissues"] > 1
        or meta_row["distinct_life_stages"] > 1
    )
    return "multi_context" if varies else "single_context"


def candidate_key(meta_row: dict) -> tuple:
    return (
        str(meta_row["feature_id_standardized"]),
        meta_row["phenotype"],
        meta_row["feature_type"],
        _partition_value(meta_row["simulated"]),
        _partition_value(meta_row["species_taxid"]),
    )


class EvidenceIndex:
    """Positional indices into a normalized evidence frame, built once per run."""

    def __init__(self, evidence_df: pd.DataFrame):
        evidence = add_partition_columns(evidence_df)
        mapped = evidence[evidence["feature_id_standardized"].notna()].copy()
        mapped["feature_id_standardized"] = mapped["feature_id_standardized"].astype(str)
        mapped["_layer"] = mapped["feature_type"].map(layer_of)
        adjusted = (
            pd.to_numeric(mapped["adjusted_p_value"], errors="coerce")
            if "adjusted_p_value" in mapped.columns
            else pd.Series(float("nan"), index=mapped.index)
        )
        mapped["_significant"] = adjusted.notna() & (adjusted <= LAYER_SUPPORT_MAX_ADJUSTED_P)
        self.evidence = mapped
        self._by_candidate = mapped.groupby(CANDIDATE_GROUP_COLUMNS, sort=False).indices
        self._by_feature = mapped.groupby(FEATURE_GROUP_COLUMNS, sort=False).indices
        # Layers with at least one significant record, per feature + phenotype + partition.
        significant = mapped[mapped["_significant"]]
        self._supported_layers = (
            significant.groupby(["feature_id_standardized", "phenotype", *PARTITION_COLUMNS])["_layer"]
            .agg(set).to_dict()
        )

    def _slice(self, index, key) -> pd.DataFrame:
        positions = index.get(key)
        if positions is None:
            return self.evidence.iloc[0:0]
        return self.evidence.iloc[positions]

    def for_candidate(self, meta_row: dict) -> pd.DataFrame:
        """All records for this feature / phenotype / feature type / partition."""
        return self._slice(self._by_candidate, candidate_key(meta_row))

    def contributing(self, meta_row: dict) -> pd.DataFrame:
        """The records that actually entered the pooled estimate."""
        subset = self.for_candidate(meta_row)
        ids = meta_row.get("contributing_evidence_ids")
        if not ids or ids != ids:
            return subset
        wanted = set(str(ids).split("|"))
        return subset[subset["evidence_id"].astype(str).isin(wanted)]

    def other_layers(self, meta_row: dict) -> pd.DataFrame:
        """Records for the same feature and partition from OTHER molecular layers,
        under any phenotype. Callers use `phenotype` and `_significant` on the
        result to separate convergent support from merely adjacent evidence."""
        key = candidate_key(meta_row)
        same_feature = self._slice(self._by_feature, (key[0], key[3], key[4]))
        return same_feature[same_feature["_layer"] != layer_of(meta_row["feature_type"])]

    def supporting_layers(self, meta_row: dict) -> set:
        """Molecular layers with a significant record for this feature under
        this candidate's own phenotype and partition. The candidate's own layer
        also counts when its pooled estimate is significant, even if no single
        contributing record is."""
        key = candidate_key(meta_row)
        layers = set(self._supported_layers.get((key[0], key[1], key[3], key[4]), set()))
        if significance_p_value(meta_row) <= LAYER_SUPPORT_MAX_ADJUSTED_P:
            layers.add(layer_of(meta_row["feature_type"]))
        return layers


def build_candidates(meta_df: pd.DataFrame, evidence_df: pd.DataFrame) -> pd.DataFrame:
    if len(meta_df) == 0:
        return pd.DataFrame()

    index = EvidenceIndex(evidence_df)

    rows = []
    for meta_dict in meta_df.to_dict("records"):
        subset = index.contributing(meta_dict)
        mapping_confidences = sorted(subset["mapping_confidence"].unique())
        quality_flags = sorted({f for flags in subset["quality_flags"] for f in _parse_flags(flags)})
        own_layer = layer_of(meta_dict["feature_type"])
        supporting = index.supporting_layers(meta_dict)
        n_supporting_layers = len(supporting)

        meta_dict["n_supporting_layers"] = n_supporting_layers
        components = compute_components(meta_dict, mapping_confidences, quality_flags)
        score = candidate_score(components)

        is_high_priority = (
            meta_dict["k_studies"] >= HIGH_PRIORITY_MIN_STUDIES
            and components["phenotype_relevance_score"] > 0.1
            and meta_dict["direction_consistency"] >= HIGH_PRIORITY_MIN_DIRECTION_CONSISTENCY
            and components["quality_score"] >= HIGH_PRIORITY_MIN_QUALITY_SCORE
            and significance_p_value(meta_dict) <= HIGH_PRIORITY_MAX_ADJUSTED_P
        )
        # The candidate's own layer must be among the supporters: two other
        # layers converging on a gene this candidate itself shows no signal for
        # is not evidence for this candidate.
        is_multi_omics = own_layer in supporting and n_supporting_layers >= MULTI_OMICS_MIN_LAYERS

        if is_high_priority:
            tier = "high_priority_cross_study"
        elif is_multi_omics:
            tier = "multi_omics_convergence"
        else:
            tier = "emerging"

        rows.append({
            **meta_dict,
            "score": score,
            "tier": tier,
            "evidence_class": evidence_class(meta_dict),
            "context_replication": context_replication(meta_dict),
            "pooled_stressors": "|".join(pooled_stressors(meta_dict)),
            "is_high_priority": is_high_priority,
            "is_multi_omics_convergence": is_multi_omics,
            "molecular_layer": own_layer,
            "n_supporting_layers": n_supporting_layers,
            "supporting_layers": "|".join(sorted(supporting)),
            "mapping_confidences": "|".join(mapping_confidences),
            "quality_flags_union": "|".join(quality_flags),
            **{f"component_{k}": v for k, v in components.items()},
        })

    out = pd.DataFrame(rows)
    return sort_candidates(out)


def sort_candidates(candidates: pd.DataFrame) -> pd.DataFrame:
    """Strongest tier first, then score descending within tier."""
    if len(candidates) == 0:
        return candidates
    tier_rank = candidates["tier"].map({t: i for i, t in enumerate(TIER_ORDER)}).fillna(len(TIER_ORDER))
    return (
        candidates.assign(_tier_rank=tier_rank)
        .sort_values(["_tier_rank", "score"], ascending=[True, False])
        .drop(columns="_tier_rank")
    )
