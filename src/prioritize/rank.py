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
from mappings.metabolite_links import links_for_partition

from .scoring import (
    candidate_score,
    compute_components,
    conservative_relevance,
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

# A metabolite cannot share a gene identifier, so metabolomics is linked to the
# gene-level layers biochemically instead (metabolite -> KEGG reaction -> KO ->
# C. gigas gene; src/mappings/metabolite_links.py). That is annotation, not a
# second measurement of the same molecule, so it is a separate, weaker tier:
# never merged into `multi_omics_convergence`, never counted in
# n_supporting_layers, and never part of the numeric score
# (docs/metabolite_gene_links.md).
METABOLITE_LAYER = "metabolomics"

# Order in which tiers are presented: strongest evidence first.
TIER_ORDER = [
    "high_priority_cross_study",
    "multi_omics_convergence",
    "metabolite_linked_convergence",
    "emerging",
]

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

# What a comparison's two arms differ in (schemas/study.schema.json). It is part
# of the candidate key, as of the meta-analysis grouping: a gene that differs
# between resilient and susceptible animals at baseline and a gene that responds
# to heat are two different candidates, even under one phenotype.
CONTRAST_TYPES = [
    "baseline_phenotype_difference",
    "phenotype_difference_under_exposure",
    "exposure_response",
    "phenotype_by_exposure_interaction",
]
BASELINE_CONTRAST = "baseline_phenotype_difference"
# For rows and frames built by hand (tests). Harmonized evidence always carries
# a contrast type, and the meta-analysis refuses evidence that does not.
DEFAULT_CONTRAST_TYPE = "exposure_response"

PARTITION_COLUMNS = ["_simulated", "_species_taxid"]
CANDIDATE_GROUP_COLUMNS = ["feature_id_standardized", "phenotype", "feature_type", *PARTITION_COLUMNS,
                           "contrast_type"]
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


def contrast_type_of(row: dict) -> str:
    value = row.get("contrast_type")
    return value if isinstance(value, str) and value else DEFAULT_CONTRAST_TYPE


def with_contrast_type(evidence_df: pd.DataFrame) -> pd.DataFrame:
    """The evidence with a contrast_type column, defaulted only where it is absent
    altogether (frames built by hand); never overwrites a recorded value."""
    if "contrast_type" in evidence_df.columns:
        return evidence_df
    return evidence_df.assign(contrast_type=DEFAULT_CONTRAST_TYPE)


def contrast_suffix(row: dict) -> str:
    """Card-name suffix that keeps two contrast types of one feature apart.

    Empty for exposure_response, so the names of existing cards do not change.
    """
    contrast = contrast_type_of(row)
    return "" if contrast == DEFAULT_CONTRAST_TYPE else f"_{contrast}"


def is_predictive_baseline(contrast_type: str, evidence_class_value: str) -> bool:
    """A candidate whose evidence is a baseline difference between groups that
    differ in a measured resilience outcome: the design for a marker that could
    be read before any stress. Its tier still says how well replicated it is."""
    return contrast_type == BASELINE_CONTRAST and evidence_class_value == "resilience_associated"


@lru_cache(maxsize=1)
def molecular_layers() -> dict:
    """feature_type -> molecular_layer, from the feature_types vocabulary."""
    return {
        term_id: term.get("molecular_layer", term_id)
        for term_id, term in load_vocab("feature_types").items()
    }


def layer_of(feature_type: str) -> str:
    return molecular_layers().get(feature_type, str(feature_type))


def evidence_class(meta_row: dict, relevance: str | None = None) -> str:
    """What kind of evidence this candidate rests on: the weakest class among
    its contributing records (`relevance`), else its phenotype's default."""
    return EVIDENCE_CLASS_BY_RELEVANCE.get(
        relevance or phenotype_relevance(meta_row["phenotype"]), DEFAULT_EVIDENCE_CLASS
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
        contrast_type_of(meta_row),
    )


class EvidenceIndex:
    """Positional indices into a normalized evidence frame, built once per run."""

    def __init__(self, evidence_df: pd.DataFrame):
        evidence = with_contrast_type(add_partition_columns(evidence_df))
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
        # Layers with at least one significant record, per feature + phenotype +
        # partition + contrast type. Expression differing at baseline and
        # methylation responding to heat are two questions, not convergence.
        significant = mapped[mapped["_significant"]]
        self._supported_layers = (
            significant.groupby(["feature_id_standardized", "phenotype", *PARTITION_COLUMNS, "contrast_type"])
            ["_layer"].agg(set).to_dict()
        )

    def _slice(self, index, key) -> pd.DataFrame:
        positions = index.get(key)
        if positions is None:
            return self.evidence.iloc[0:0]
        return self.evidence.iloc[positions]

    def for_candidate(self, meta_row: dict) -> pd.DataFrame:
        """All records for this feature / phenotype / feature type / partition / contrast type."""
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
        this candidate's own phenotype, partition and contrast type. The candidate's own layer
        also counts when its pooled estimate is significant, even if no single
        contributing record is."""
        key = candidate_key(meta_row)
        layers = set(self._supported_layers.get((key[0], key[1], key[3], key[4], key[5]), set()))
        if significance_p_value(meta_row) <= LAYER_SUPPORT_MAX_ADJUSTED_P:
            layers.add(layer_of(meta_row["feature_type"]))
        return layers

    def links_for(self, meta_row: dict) -> list:
        """Metabolite <-> gene links touching this candidate's feature, from the
        link table of its own partition (simulated demo or real)."""
        key = candidate_key(meta_row)
        table = links_for_partition(key[3], key[4])
        if not table.available:
            return []
        return table.links_for_feature(key[0], layer_of(meta_row["feature_type"]) == METABOLITE_LAYER)

    @staticmethod
    def linked_feature_id(link, from_metabolite: bool) -> str:
        return link.ncbi_gene_id if from_metabolite else link.metabolite

    def linked_records(self, meta_row: dict) -> pd.DataFrame:
        """Evidence for every feature linked to this candidate, under any
        phenotype, from layers other than the candidate's own, with the link
        that connects each record in `_link`. Callers use `phenotype`,
        `_significant` and `_link.credited` to tell support from context."""
        key = candidate_key(meta_row)
        own_layer = layer_of(meta_row["feature_type"])
        from_metabolite = own_layer == METABOLITE_LAYER
        frames = []
        for link in self.links_for(meta_row):
            linked_id = self.linked_feature_id(link, from_metabolite)
            records = self._slice(self._by_feature, (linked_id, key[3], key[4]))
            records = records[records["_layer"] != own_layer]
            if len(records):
                frames.append(records.assign(_link=[link] * len(records)))
        if not frames:
            return self.evidence.iloc[0:0].assign(_link=None)
        return pd.concat(frames)

    def linked_support(self, meta_row: dict, supporting: set) -> tuple:
        """(linked_layers, credited links behind them) for this candidate.

        A linked layer is one with a significant record for a CREDITED linked
        feature under this candidate's own phenotype, partition and contrast type, that is not
        already a shared-identifier supporting layer — shared identity is the
        stronger link and is never double-counted as a biochemical one.
        """
        key = candidate_key(meta_row)
        from_metabolite = layer_of(meta_row["feature_type"]) == METABOLITE_LAYER
        layers, links = set(), []
        for link in self.links_for(meta_row):
            if not link.credited:
                continue
            linked_id = self.linked_feature_id(link, from_metabolite)
            found = self._supported_layers.get((linked_id, key[1], key[3], key[4], key[5]), set()) - supporting
            found.discard(layer_of(meta_row["feature_type"]))
            if found:
                layers |= found
                links.append(link)
        return layers, links


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
        linked_layers, linked_links = index.linked_support(meta_dict, supporting)
        link_table = links_for_partition(_partition_value(meta_dict["simulated"]),
                                         _partition_value(meta_dict["species_taxid"]))

        classifications = (subset["resilience_classification"]
                           if "resilience_classification" in subset.columns else [])
        relevance = conservative_relevance(meta_dict["phenotype"], classifications)

        meta_dict["n_supporting_layers"] = n_supporting_layers
        orthogroups = (sorted(subset["orthogroup_id"].dropna().astype(str).unique())
                       if "orthogroup_id" in subset.columns else [])
        components = compute_components(meta_dict, mapping_confidences, quality_flags, relevance)
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
        # Same own-layer rule for the biochemical link: the candidate itself must
        # carry signal, and some other layer must carry it for a linked feature.
        is_metabolite_linked = own_layer in supporting and bool(linked_layers)

        if is_high_priority:
            tier = "high_priority_cross_study"
        elif is_multi_omics:
            tier = "multi_omics_convergence"
        elif is_metabolite_linked:
            tier = "metabolite_linked_convergence"
        else:
            tier = "emerging"

        if is_multi_omics:
            convergence_basis = "shared_identifier"
        elif is_metabolite_linked:
            convergence_basis = "metabolite_gene_link"
        else:
            convergence_basis = "none"

        klass = evidence_class(meta_dict, relevance)
        contrast = contrast_type_of(meta_dict)
        rows.append({
            **meta_dict,
            "contrast_type": contrast,
            "score": score,
            "tier": tier,
            "evidence_class": klass,
            "is_predictive_baseline": is_predictive_baseline(contrast, klass),
            "context_replication": context_replication(meta_dict),
            "pooled_stressors": "|".join(pooled_stressors(meta_dict)),
            "pooling_status": "pooled",
            "is_high_priority": is_high_priority,
            "is_multi_omics_convergence": is_multi_omics,
            "orthogroup_id": "|".join(orthogroups),
            "molecular_layer": own_layer,
            "n_supporting_layers": n_supporting_layers,
            "supporting_layers": "|".join(sorted(supporting)),
            "is_metabolite_linked_convergence": is_metabolite_linked,
            "convergence_basis": convergence_basis,
            "n_linked_layers": len(linked_layers),
            "linked_layers": "|".join(sorted(linked_layers)),
            "linked_features": "|".join(sorted({link.describe() for link in linked_links})),
            "link_confidences": "|".join(sorted({link.link_confidence for link in linked_links})),
            "link_table_status": link_table.status,
            "link_table_sha256": link_table.sha256 or "",
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
