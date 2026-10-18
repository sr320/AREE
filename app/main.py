"""AREE Streamlit interface.

Launch with:  streamlit run app/main.py

This is a read-only browser over artifacts produced by the `aree` CLI:
- registry/study_registry.csv         (aree build-registry)
- reports/evidence/evidence_table.tsv (aree harmonize)
- candidate synthesis                 (computed live from the evidence table)

It never writes to those files. If an artifact is missing, the relevant tab
tells you which command to run.
"""
from __future__ import annotations

import data_access as da
import streamlit as st

st.set_page_config(page_title="AREE — Aquaculture Resilience Evidence Engine", layout="wide")

st.title("AREE — Aquaculture Resilience Evidence Engine")
st.caption(
    "Harmonized, cross-study resilience biomarker evidence for the Pacific oyster (Magallana gigas), "
    "mussels and clams, each species against its own gene identifiers. Simulated demo data is hidden "
    "unless you include it from the sidebar, and is never pooled with real-study evidence."
)

st.warning(
    "Candidate scores reflect **associations across available evidence, not validated biomarkers**. "
    "A high score never overrides the study-count and direction-consistency gates required for the "
    "high-priority tier. Read the evidence card limitations before acting on any candidate.",
    icon="⚠️",
)

include_simulated = st.sidebar.checkbox(
    "Include simulated demo data", value=False, key="include_simulated",
    help="The simulated demo studies (registry/demo_studies/) are pipeline test fixtures, not findings. "
         "They are hidden unless this is checked. Simulated and real evidence are never pooled together.",
)

tab_studies, tab_evidence, tab_candidates, tab_search = st.tabs(
    ["Studies", "Harmonized evidence", "Candidate biomarkers", "Search a feature"]
)


# --------------------------------------------------------------------------- Studies
with tab_studies:
    st.header("Registered studies")
    st.subheader("Species")
    st.dataframe(da.load_species_summary(include_simulated), use_container_width=True, hide_index=True)
    st.caption("Gene identifiers come from each species' own crosswalk: `ncbi_gene` (NCBI Gene) or "
               "`gff_homology` (the authors' annotation plus a recorded homology search). Species in a "
               "hybridizing `complex` carry a species-identity basis on every study; see "
               "docs/species_identity.md.")
    registry = da.load_registry(include_simulated)
    if registry.empty:
        st.info("No studies registered yet. Run: `aree build-registry`")
    else:
        species_filter = st.multiselect("Species", sorted(registry["species_accepted"].dropna().unique()),
                                        key="studies_species")
        assay_filter = st.multiselect("Assay type", sorted(registry["assay_type"].unique()))
        mode_filter = st.multiselect("Analysis mode", sorted(registry["analysis_mode"].unique()))
        view = registry.copy()
        if species_filter:
            view = view[view["species_accepted"].isin(species_filter)]
        if assay_filter:
            view = view[view["assay_type"].isin(assay_filter)]
        if mode_filter:
            view = view[view["analysis_mode"].isin(mode_filter)]
        st.dataframe(view, use_container_width=True, hide_index=True)
        st.caption(f"{len(view)} of {len(registry)} studies shown. "
                   f"{int(registry['simulated'].sum()) if 'simulated' in registry else 0} are simulated.")
        st.download_button("Download filtered studies (CSV)", view.to_csv(index=False),
                           file_name="aree_studies_filtered.csv", mime="text/csv")


# --------------------------------------------------------------------------- Evidence
with tab_evidence:
    st.header("Harmonized evidence records")
    evidence = da.load_evidence(include_simulated)
    if evidence.empty:
        st.info("No evidence table yet. Run: `aree harmonize --study STUDY_ID` for each registered study.")
    else:
        evidence = evidence.assign(species_accepted=evidence["species_taxid"].map(da.species_label))
        col1, col2, col3 = st.columns(3)
        with col1:
            species_f = st.multiselect("Species", sorted(evidence["species_accepted"].dropna().unique()),
                                       key="evidence_species")
            stressor_f = st.multiselect("Stressor", sorted(evidence["stressor"].dropna().unique()))
            phenotype_f = st.multiselect("Phenotype", sorted(evidence["phenotype"].dropna().unique()))
        with col2:
            tissue_f = st.multiselect("Tissue", sorted(evidence["tissue"].dropna().unique()))
            life_stage_f = st.multiselect("Life stage", sorted(evidence["life_stage"].dropna().unique()))
        with col3:
            feature_type_f = st.multiselect("Feature type", sorted(evidence["feature_type"].dropna().unique()))
            mapping_f = st.multiselect("Mapping confidence", sorted(evidence["mapping_confidence"].dropna().unique()))
            data_origin = "Real studies only" if not include_simulated else st.radio(
                "Data origin", ["All", "Real studies only", "Simulated only"], index=0,
                help="Simulated demo evidence and evidence from real published studies are "
                     "never pooled together; this filter makes the split explicit.",
            )

        view = evidence.copy()
        if "simulated" in view.columns:
            flag = view["simulated"].astype(str).str.lower() == "true"
            if data_origin == "Real studies only":
                view = view[~flag]
            elif data_origin == "Simulated only":
                view = view[flag]
        for col, selected in [
            ("species_accepted", species_f),
            ("stressor", stressor_f), ("phenotype", phenotype_f), ("tissue", tissue_f),
            ("life_stage", life_stage_f), ("feature_type", feature_type_f), ("mapping_confidence", mapping_f),
        ]:
            if selected:
                view = view[view[col].isin(selected)]

        display_cols = [
            "study_id", "simulated", "species_accepted", "feature_id_standardized", "feature_type",
            "molecular_direction",
            "effect_size", "adjusted_p_value", "tissue", "life_stage", "stressor",
            "phenotype", "mapping_confidence",
        ]
        st.dataframe(view[display_cols], use_container_width=True, hide_index=True)
        st.caption(f"{len(view)} of {len(evidence)} evidence records shown.")
        st.download_button("Download filtered evidence (TSV)", view.to_csv(sep="\t", index=False),
                           file_name="aree_evidence_filtered.tsv", mime="text/tab-separated-values")


# --------------------------------------------------------------------------- Candidates
with tab_candidates:
    st.header("Candidate biomarkers")
    if not da.evidence_exists():
        st.info("No evidence table yet. Run `aree harmonize` first, then `aree build-evidence-cards`.")
    else:
        labels = da.phenotype_labels()
        evidence = da.load_evidence(include_simulated)
        phenotype_options = ["(all)"] + sorted(evidence["phenotype"].dropna().unique())
        chosen = st.selectbox("Phenotype", phenotype_options, format_func=lambda p: labels.get(p, p))
        phenotype = None if chosen == "(all)" else chosen

        try:
            candidates = da.load_candidates(phenotype=phenotype, include_simulated=include_simulated)
        except ValueError as exc:
            # The meta-analysis refuses evidence it cannot pool honestly (e.g. an
            # evidence table harmonized by older code); say why instead of crashing.
            st.error(f"Candidates could not be computed from the current evidence table: {exc}\n\n"
                     "Re-harmonize (`make demo`, or `aree harmonize` per study) and reload.")
            st.stop()
        if candidates.empty:
            st.info("No candidates for this filter.")
        else:
            candidates = candidates.assign(species=candidates["species_taxid"].map(da.species_label))
            species_c = st.multiselect("Species", sorted(candidates["species"].unique()), key="candidates_species",
                                       help="Candidates never pool across species; each is ranked within "
                                            "its own. Cross-species results are in "
                                            "reports/evidence_cards/cross_species/ (aree build-cross-species-cards).")
            tier_f = st.multiselect("Tier", sorted(candidates["tier"].unique()),
                                    default=sorted(candidates["tier"].unique()))
            view = candidates[candidates["tier"].isin(tier_f)] if tier_f else candidates
            if species_c:
                view = view[view["species"].isin(species_c)]
            if "evidence_class" in view.columns:
                class_f = st.multiselect(
                    "Evidence class", sorted(view["evidence_class"].dropna().unique()),
                    default=sorted(view["evidence_class"].dropna().unique()),
                )
                if class_f:
                    view = view[view["evidence_class"].isin(class_f)]
            show_cols = [
                "feature_id_standardized", "phenotype", "feature_type", "simulated",
                "species", "tier", "evidence_class", "context_replication", "score",
                "k_studies", "n_supporting_layers", "direction_consistency", "i_squared", "pooled_effect",
                "adjusted_p_value",
            ]
            # A candidates.tsv written before evidence_class existed still loads.
            st.dataframe(view[[c for c in show_cols if c in view.columns]],
                         use_container_width=True, hide_index=True)
            st.download_button("Download candidates (CSV)", view.to_csv(index=False),
                               file_name="aree_candidates.csv", mime="text/csv")

            st.subheader("Tier legend")
            st.markdown(
                "- **high_priority_cross_study** — ≥2 independent studies, interpretable phenotype, "
                "direction consistency ≥ 0.7, acceptable quality, and BH-adjusted pooled p ≤ 0.05.\n"
                "- **multi_omics_convergence** — ≥2 molecular layers (e.g. transcriptomics + DNA methylation) "
                "each with a significant record for the same standardized feature under the **same phenotype**, "
                "including the candidate's own layer.\n"
                "- **metabolite_linked_convergence** — the candidate's own layer is significant and a "
                "metabolite and a gene are linked only through KEGG reaction annotation "
                "(metabolite → reaction → KO → gene), under the same phenotype. Weaker than shared gene "
                "identity, ranked below it, and not counted in the score "
                "(see docs/metabolite_gene_links.md).\n"
                "- **emerging** — single-study or otherwise unconfirmed; **requires replication**."
            )

            st.subheader("Evidence class")
            st.markdown(
                "The tier says how well replicated a signal is; `evidence_class` says what it is "
                "evidence **of**, taken from the phenotype the contributing studies actually "
                "measured. Only `resilience_associated` candidates rest on a measured resilience "
                "outcome. A `disease_associated` or `stress_response` candidate in the top tier is "
                "well-replicated evidence of a *response*, not evidence that the feature marks a "
                "resilient animal — see the card's \"What this evidence does and does not show\"."
            )


# --------------------------------------------------------------------------- Search
with tab_search:
    st.header("Search a gene / protein / feature")
    evidence = da.load_evidence(include_simulated)
    if evidence.empty:
        st.info("No evidence table yet. Run `aree harmonize` first.")
    else:
        query = st.text_input("Feature identifier (standardized or original), e.g. LOC105333935 or hsp70")
        if query:
            q = query.strip().lower()
            mask = (
                evidence["feature_id_standardized"].astype(str).str.lower().str.contains(q, na=False)
                | evidence["feature_id_original"].astype(str).str.lower().str.contains(q, na=False)
                | evidence["orthogroup_id"].astype(str).str.lower().str.contains(q, na=False)
            )
            hits = evidence[mask]
            if hits.empty:
                st.warning("No matching evidence records.")
            else:
                st.success(f"{len(hits)} evidence records across {hits['study_id'].nunique()} studies.")
                st.dataframe(
                    hits[[
                        "study_id", "species", "feature_id_original", "feature_id_standardized", "feature_type",
                        "molecular_direction", "effect_size", "phenotype", "stressor", "mapping_confidence",
                    ]],
                    use_container_width=True, hide_index=True,
                )
