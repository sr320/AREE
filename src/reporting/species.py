"""Per-species summary for the dashboard and the Streamlit app. WP8.

One row per registered species: its crosswalk (route, gene count, whether it
is built), its assemblies, how many real and simulated studies are registered
and harmonized, and how much evidence each origin holds. Species are taken
from the vocabulary, so a species registered without any study yet still
appears, with zeros, rather than being hidden.

A study's species is resolved through the vocabulary (accepted synonyms
included), so *Crassostrea gigas* and *Magallana gigas* count as one species.
"""
from __future__ import annotations

from collections import Counter, defaultdict

import pandas as pd

from common import load_assemblies, load_json, load_vocab, registered_crosswalks, resolve_species


def _crosswalk_facts(entry: dict) -> dict:
    path = entry["path"]
    facts = {"slug": path.name.replace("_gene_id_crosswalk.tsv", ""), "built": path.exists(),
             "route": None, "n_genes": None}
    provenance = path.with_name(f"{path.stem}.provenance.json")
    if provenance.exists():
        prov = load_json(provenance)
        facts["route"] = prov.get("route", "ncbi_gene")
        facts["n_genes"] = (prov.get("coverage") or {}).get("genes_total")
    return facts


def species_identity_basis(study: dict) -> str | None:
    """The basis recorded for a study's species identity, if any (WP1)."""
    identity = study.get("species_identity") or {}
    return identity.get("basis")


def species_summary(studies: list[dict], evidence: pd.DataFrame | None = None) -> list[dict]:
    """One row per registered species, in vocabulary order."""
    vocab = load_vocab("species")
    taxid_of_study = {}
    for study in studies:
        resolved = resolve_species(study.get("species"))
        if resolved is not None:
            taxid_of_study[study["study_id"]] = resolved.ncbi_taxid

    studies_by: dict[tuple, list[str]] = defaultdict(list)
    for study in studies:
        taxid = taxid_of_study.get(study["study_id"])
        studies_by[(taxid, bool(study.get("simulated", False)))].append(study["study_id"])

    records: Counter = Counter()
    harmonized: dict[tuple, set] = defaultdict(set)
    if evidence is not None and len(evidence):
        sim = evidence["simulated"].astype(str).str.lower() == "true"
        taxid = pd.to_numeric(evidence["species_taxid"], errors="coerce")
        for (t, s), n in pd.DataFrame({"t": taxid, "s": sim}).value_counts().items():
            if t == t:
                records[(int(t), bool(s))] = int(n)
        for (t, s, sid) in set(zip(taxid, sim, evidence["study_id"])):
            if t == t:
                harmonized[(int(t), bool(s))].add(sid)

    crosswalks = {int(e["ncbi_taxid"]): _crosswalk_facts(e) for e in registered_crosswalks()}
    assemblies: dict[int, list[str]] = defaultdict(list)
    for record in load_assemblies().values():
        assemblies[int(record["ncbi_taxid"])].append(record["assembly_id"])

    rows = []
    for term in vocab.values():
        t = int(term["ncbi_taxid"])
        cw = crosswalks.get(t, {})
        rows.append({
            "species_id": term["id"],
            "species": term["scientific_name"],
            "label": term["label"],
            "ncbi_taxid": t,
            "complex": term.get("complex"),
            "assemblies": sorted(assemblies.get(t, [])),
            "crosswalk": cw.get("slug"),
            "crosswalk_built": bool(cw.get("built")),
            "crosswalk_route": cw.get("route"),
            "crosswalk_genes": cw.get("n_genes"),
            "n_studies_real": len(studies_by.get((t, False), [])),
            "n_studies_simulated": len(studies_by.get((t, True), [])),
            "n_studies_real_harmonized": len(harmonized.get((t, False), set())),
            "n_evidence_records_real": records.get((t, False), 0),
            "n_evidence_records_simulated": records.get((t, True), 0),
            "studies_real": sorted(studies_by.get((t, False), [])),
        })
    return rows
