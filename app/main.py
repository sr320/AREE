"""AREE Streamlit interface.

Launch with:  streamlit run app/main.py

This is a read-only browser over artifacts produced by the `aree` CLI:
- registry/study_registry.csv         (aree register-study)
- reports/evidence/evidence_table.tsv (aree harmonize)
- candidate synthesis                 (computed live from the evidence table)

It never writes to those files. If an artifact is missing, the relevant tab
tells you which command to run.
"""
from __future__ import annotations

import streamlit as st

import data_access as da

st.set_page_config(page_title="AREE — Aquaculture Resilience Evidence Engine", layout="wide")

st.title("AREE — Aquaculture Resilience Evidence Engine")
st.caption(
    "Harmonized, cross-study resilience biomarker evidence for Pacific oyster (Crassostrea gigas) "
    "and other aquaculture species. All datasets shipped in this build are SIMULATED demo data."
)

st.warning(
    "Candidate scores reflect **associations across available evidence, not validated biomarkers**. "
    "A high score never overrides the study-count and direction-consistency gates required for the "
    "high-priority tier. Read the evidence card limitations before acting on any candidate.",
    icon="⚠️",
)

tab_studies, tab_evidence, tab_candidates, tab_search = st.tabs(
    ["Studies", "Harmonized evidence", "Candidate biomarkers", "Search a feature"]
)


# --------------------------------------------------------------------------- Studies
with tab_studies:
    st.header("Registered studies")
    registry = da.load_registry()
    if registry.empty:
        st.info("No studies registered yet. Run: `aree register-study registry/studies/STUDY_ID.yaml`")
    else:
        assay_filter = st.multiselect("Assay type", sorted(registry["assay_type"].unique()))
        mode_filter = st.multiselect("Analysis mode", sorted(registry["analysis_mode"].unique()))
        view = registry.copy()
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
    evidence = da.load_evidence()
    if evidence.empty:
        st.info("No evidence table yet. Run: `aree harmonize --study STUDY_ID` for each registered study.")
    else:
        col1, col2, col3 = st.columns(3)
        with col1:
            stressor_f = st.multiselect("Stressor", sorted(evidence["stressor"].dropna().unique()))
            phenotype_f = st.multiselect("Phenotype", sorted(evidence["phenotype"].dropna().unique()))
        with col2:
            tissue_f = st.multiselect("Tissue", sorted(evidence["tissue"].dropna().unique()))
            life_stage_f = st.multiselect("Life stage", sorted(evidence["life_stage"].dropna().unique()))
        with col3:
            feature_type_f = st.multiselect("Feature type", sorted(evidence["feature_type"].dropna().unique()))
            mapping_f = st.multiselect("Mapping confidence", sorted(evidence["mapping_confidence"].dropna().unique()))
            data_origin = st.radio(
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
            ("stressor", stressor_f), ("phenotype", phenotype_f), ("tissue", tissue_f),
            ("life_stage", life_stage_f), ("feature_type", feature_type_f), ("mapping_confidence", mapping_f),
        ]:
            if selected:
                view = view[view[col].isin(selected)]

        display_cols = [
            "study_id", "simulated", "feature_id_standardized", "feature_type", "molecular_direction",
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
        evidence = da.load_evidence()
        phenotype_options = ["(all)"] + sorted(evidence["phenotype"].dropna().unique())
        chosen = st.selectbox("Phenotype", phenotype_options, format_func=lambda p: labels.get(p, p))
        phenotype = None if chosen == "(all)" else chosen

        candidates = da.load_candidates(phenotype=phenotype)
        if candidates.empty:
            st.info("No candidates for this filter.")
        else:
            tier_f = st.multiselect("Tier", sorted(candidates["tier"].unique()),
                                    default=sorted(candidates["tier"].unique()))
            view = candidates[candidates["tier"].isin(tier_f)] if tier_f else candidates
            show_cols = [
                "feature_id_standardized", "phenotype", "feature_type", "tier", "score",
                "k_studies", "n_distinct_assays", "direction_consistency", "i_squared", "pooled_effect",
            ]
            st.dataframe(view[show_cols], use_container_width=True, hide_index=True)
            st.download_button("Download candidates (CSV)", view.to_csv(index=False),
                               file_name="aree_candidates.csv", mime="text/csv")

            st.subheader("Tier legend")
            st.markdown(
                "- **high_priority_cross_study** — ≥2 independent studies, interpretable phenotype, "
                "direction consistency ≥ 0.7, acceptable quality.\n"
                "- **multi_omics_convergence** — same standardized feature supported by ≥2 molecular layers.\n"
                "- **emerging** — single-study or otherwise unconfirmed; **requires replication**."
            )


# --------------------------------------------------------------------------- Search
with tab_search:
    st.header("Search a gene / protein / feature")
    evidence = da.load_evidence()
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
                        "study_id", "feature_id_original", "feature_id_standardized", "feature_type",
                        "molecular_direction", "effect_size", "phenotype", "stressor", "mapping_confidence",
                    ]],
                    use_container_width=True, hide_index=True,
                )
