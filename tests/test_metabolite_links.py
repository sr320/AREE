"""Metabolite -> gene links: the weaker, separately labelled kind of convergence.

See docs/metabolite_gene_links.md.
"""
import json
import re
from pathlib import Path

import jsonschema
import pandas as pd
import pytest
from conftest import harmonize_all_demo_studies, oyster

import mappings.metabolite_links as ml
import prioritize.rank as rank
import reporting.evidence_cards as evidence_cards
from common import REPO_ROOT, SCHEMAS_DIR, load_json, load_vocab, sha256_file
from harmonize.core import load_evidence_table
from mappings.build_metabolite_links import build, equation_side, parse_entries
from prioritize.rank import TIER_ORDER, build_candidates

LDH, SDHA, GLUT1 = "LOC105342001", "LOC105342002", "LOC105342003"


@pytest.fixture(autouse=True)
def _fresh_link_caches():
    ml.clear_caches()
    yield
    ml.clear_caches()


# --------------------------------------------------------------------------- #
# Curated inputs
# --------------------------------------------------------------------------- #


def test_compound_map_is_well_formed():
    compound_map = ml.load_compound_map()
    assert compound_map, "compound map is empty"
    for name, compounds in compound_map.items():
        assert name == ml.normalize_metabolite_name(name)
        for c in compounds:
            assert re.fullmatch(r"C\d{5}", c["kegg_compound"]), c
            assert c["match_basis"] in {"exact_name", "synonym", "curator_judgement"}
            assert c["chebi_id"] is None or c["chebi_id"].startswith("CHEBI:")


def test_every_demo_metabolite_name_with_a_standardized_id_is_curated():
    """A level 1-2 demo metabolite missing from the map would silently never link."""
    table = pd.read_csv(
        REPO_ROOT / "data/demo/metabolomics/GIGAS_GROW06_hypoxia_vs_control_metabolite_features_demo.tsv",
        sep="\t",
    )
    names = table.loc[table["annotation_confidence_level"] <= 2, "putative_metabolite_name"]
    missing = {n for n in names if ml.normalize_metabolite_name(n) not in ml.load_compound_map()}
    assert not missing


def test_atp_is_a_currency_metabolite():
    assert "C00002" in ml.load_currency_compounds()


def test_demo_link_table_provenance_matches_the_file_and_is_simulated():
    prov = load_json(ml.provenance_path(ml.DEMO_LINKS_PATH))
    assert prov["simulated"] is True
    assert prov["artifact_sha256"] == sha256_file(ml.DEMO_LINKS_PATH)
    jsonschema.validate(prov, load_json(SCHEMAS_DIR / "metabolite_gene_links_provenance.schema.json"))


def test_real_kegg_link_table_is_not_tracked_by_git():
    """KEGG redistribution terms: the real table is built locally, never committed."""
    ignored = (REPO_ROOT / ".gitignore").read_text()
    assert "data/mappings/metabolite_gene_links_mgigas.*" in ignored


# --------------------------------------------------------------------------- #
# Mapping confidence
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize("basis,n_compounds,n_genes,expected", [
    ("exact_name", 1, 1, "inferred"),
    ("exact_name", 1, 4, "one_to_many_ortholog"),
    ("synonym", 1, 4, "one_to_many_ortholog"),
    ("curator_judgement", 1, 4, "inferred"),
    ("exact_name", 2, 4, "inferred"),
])
def test_link_confidence_is_never_exact(basis, n_compounds, n_genes, expected):
    confidence = ml.link_confidence(basis, n_compounds, n_genes)
    assert confidence == expected
    assert confidence not in {"exact", "one_to_one_ortholog"}
    assert confidence in load_vocab("mapping_confidence")


def _row(cpd, gene, reaction="R1", ko="K1", side="left"):
    return {"kegg_compound": cpd, "compound_name": cpd, "kegg_reaction": reaction,
            "equation_side": side, "ko": ko, "ncbi_gene_id": gene, "link_type": "enzyme"}


def test_hub_and_currency_compounds_are_shown_but_never_credited():
    compound_map = {
        "hubby": [{"kegg_compound": "C00099", "match_basis": "exact_name"}],
        "atp": [{"kegg_compound": "C00002", "match_basis": "exact_name"}],
        "rare": [{"kegg_compound": "C00050", "match_basis": "exact_name"}],
    }
    rows = [_row("C00099", f"G{i}") for i in range(4)] + [_row("C00002", "G1"), _row("C00050", "G9")]
    links, exclusions = ml.build_links(rows, compound_map, {"C00002": "currency"}, max_genes_per_compound=3)

    hub = [link for link in links if link.metabolite == "hubby"]
    assert len(hub) == 4 and not any(link.credited for link in hub)
    assert {(e.metabolite, e.reason) for e in exclusions} == {
        ("hubby", "hub_compound"), ("atp", "currency_metabolite")}
    assert not [link for link in links if link.metabolite == "atp"]
    (rare,) = [link for link in links if link.metabolite == "rare"]
    assert rare.credited and rare.link_confidence == "inferred"


def test_ambiguous_name_links_every_candidate_compound_as_inferred():
    """'lactate' from LC-MS cannot be told apart from its stereoisomer."""
    compound_map = {"lactate": [{"kegg_compound": "C00186", "match_basis": "curator_judgement"},
                                {"kegg_compound": "C00256", "match_basis": "curator_judgement"}]}
    links, exclusions = ml.build_links([_row("C00256", "G1"), _row("C00256", "G2", ko="K2")],
                                       compound_map, {})
    assert {link.link_confidence for link in links} == {"inferred"}
    assert [(e.kegg_compound, e.reason) for e in exclusions] == [("C00186", "no_reaction_in_table")]


# --------------------------------------------------------------------------- #
# Partition selection: demo and real tables are never mixed
# --------------------------------------------------------------------------- #


def test_simulated_oyster_partition_uses_the_demo_table_only():
    table = ml.links_for_partition("true", "29159")
    assert table.available and table.path == ml.DEMO_LINKS_PATH
    assert {link.metabolite for link in table.for_gene(LDH)} == {"lactate", "l_lactate", "pyruvate"}


def test_simulated_other_species_has_no_link_table():
    """The demo mussel reuses oyster demo IDs; it must not inherit oyster links."""
    assert ml.links_for_partition("true", "29158").status == ml.STATUS_NOT_REGISTERED


def test_real_partition_reports_an_unbuilt_table(monkeypatch, tmp_path):
    monkeypatch.setattr(ml, "registered_link_path", lambda taxid: tmp_path / "missing.tsv")
    table = ml.links_for_partition("false", "29159")
    assert table.status == ml.STATUS_NOT_BUILT
    assert table.label() == "not_built"


# --------------------------------------------------------------------------- #
# Ranking
# --------------------------------------------------------------------------- #


def test_linked_tier_sits_between_shared_identity_and_emerging():
    assert TIER_ORDER.index("multi_omics_convergence") < TIER_ORDER.index(
        "metabolite_linked_convergence") < TIER_ORDER.index("emerging")


def _candidate(feature, feature_type, evidence_id, study):
    return pd.DataFrame([{
        "feature_id_standardized": feature, "phenotype": "growth_under_stress",
        "feature_type": feature_type, "simulated": False, "species_taxid": 29159,
        "k_studies": 1, "studies": study, "n_evidence_records": 1, "total_sample_size": 10,
        "pooled_effect": 1.2, "pooled_se": 0.3, "ci_lower": 0.6, "ci_upper": 1.8, "z": 4.0,
        "p_value": 1e-4, "adjusted_p_value": 0.01, "n_tests_in_family": 10,
        "q_statistic": 0.0, "i_squared": 0.0, "tau_squared": 0.0, "direction_consistency": 1.0,
        "distinct_tissues": 1, "distinct_life_stages": 1, "distinct_stressors": "hypoxia",
        "contributing_evidence_ids": evidence_id,
    }])


def _evidence(rows):
    return pd.DataFrame([{
        "evidence_id": eid, "study_id": study, "feature_id_standardized": feature,
        "phenotype": phenotype, "feature_type": ftype, "simulated": False, "species_taxid": 29159,
        "mapping_confidence": "exact", "quality_flags": [], "adjusted_p_value": adj_p,
    } for eid, study, feature, ftype, phenotype, adj_p in rows])


@pytest.fixture
def synthetic_links(monkeypatch):
    """Real partition, taxid 29159, linked through an in-memory table."""
    links, exclusions = ml.build_links(
        [_row("C00186", "GENE_LDH", reaction="R00703", ko="K00016")],
        {"lactate": [{"kegg_compound": "C00186", "match_basis": "exact_name"}]}, {},
    )
    table = ml.MetaboliteGeneLinks(path=Path("synthetic.tsv"), status=ml.STATUS_AVAILABLE,
                                   sha256="0" * 64, source="test", links=links, exclusions=exclusions)
    monkeypatch.setattr(rank, "links_for_partition", lambda simulated, taxid: table)
    monkeypatch.setattr(evidence_cards, "links_for_partition", lambda simulated, taxid: table)
    return table


def test_gene_and_metabolite_converge_through_a_credited_link(synthetic_links):
    evidence = _evidence([
        ("rna", "RNA", "GENE_LDH", "gene", "growth_under_stress", 0.01),
        ("met", "MET", "lactate", "metabolite_feature", "growth_under_stress", 0.02),
    ])
    gene = build_candidates(_candidate("GENE_LDH", "gene", "rna", "RNA"), evidence).iloc[0]
    metabolite = build_candidates(_candidate("lactate", "metabolite_feature", "met", "MET"), evidence).iloc[0]

    for candidate, linked_layer in [(gene, "metabolomics"), (metabolite, "transcriptomics")]:
        assert candidate["tier"] == "metabolite_linked_convergence"
        assert candidate["convergence_basis"] == "metabolite_gene_link"
        assert candidate["linked_layers"] == linked_layer
        # A link is never counted as a shared-identifier layer.
        assert candidate["n_supporting_layers"] == 1
        assert candidate["linked_features"] == "lactate→GENE_LDH via C00186/R00703/K00016"
        assert candidate["link_table_sha256"] == "0" * 64


def test_link_does_not_change_the_score(synthetic_links, monkeypatch):
    evidence = _evidence([
        ("rna", "RNA", "GENE_LDH", "gene", "growth_under_stress", 0.01),
        ("met", "MET", "lactate", "metabolite_feature", "growth_under_stress", 0.02),
    ])
    linked = build_candidates(_candidate("GENE_LDH", "gene", "rna", "RNA"), evidence).iloc[0]
    monkeypatch.setattr(rank, "links_for_partition",
                        lambda simulated, taxid: ml.MetaboliteGeneLinks(path=None, status="not_registered"))
    unlinked = build_candidates(_candidate("GENE_LDH", "gene", "rna", "RNA"), evidence).iloc[0]
    assert unlinked["tier"] == "emerging"
    assert linked["score"] == unlinked["score"]


@pytest.mark.parametrize("phenotype,adj_p", [("thermal_tolerance", 0.01), ("growth_under_stress", 0.3)])
def test_linked_evidence_must_be_significant_and_same_phenotype(synthetic_links, phenotype, adj_p):
    evidence = _evidence([
        ("rna", "RNA", "GENE_LDH", "gene", "growth_under_stress", 0.01),
        ("met", "MET", "lactate", "metabolite_feature", phenotype, adj_p),
    ])
    candidate = build_candidates(_candidate("GENE_LDH", "gene", "rna", "RNA"), evidence).iloc[0]
    assert candidate["tier"] == "emerging"
    assert candidate["convergence_basis"] == "none"


def test_candidate_must_carry_its_own_signal(synthetic_links):
    evidence = _evidence([
        ("rna", "RNA", "GENE_LDH", "gene", "growth_under_stress", 0.4),
        ("met", "MET", "lactate", "metabolite_feature", "growth_under_stress", 0.01),
    ])
    meta = _candidate("GENE_LDH", "gene", "rna", "RNA").assign(p_value=0.3, adjusted_p_value=0.4)
    candidate = build_candidates(meta, evidence).iloc[0]
    assert candidate["tier"] == "emerging"


def test_shared_identity_outranks_a_link(synthetic_links):
    evidence = _evidence([
        ("rna", "RNA", "GENE_LDH", "gene", "growth_under_stress", 0.01),
        ("meth", "METH", "GENE_LDH", "methylation_region", "growth_under_stress", 0.01),
        ("met", "MET", "lactate", "metabolite_feature", "growth_under_stress", 0.01),
    ])
    candidate = build_candidates(_candidate("GENE_LDH", "gene", "rna", "RNA"), evidence).iloc[0]
    assert candidate["tier"] == "multi_omics_convergence"
    assert candidate["convergence_basis"] == "shared_identifier"
    # The metabolite layer is still reported as linked context, not as shared identity.
    assert candidate["supporting_layers"] == "dna_methylation|transcriptomics"
    assert candidate["linked_layers"] == "metabolomics"


# --------------------------------------------------------------------------- #
# Demo end to end: GIGAS_GROW06 (metabolomics) x GIGAS_HEAT01 chronic (RNA-seq)
# --------------------------------------------------------------------------- #


def _demo_growth_candidates():
    harmonize_all_demo_studies()
    from reporting.evidence_cards import rank_all_candidates
    candidates = rank_all_candidates(load_evidence_table(), "growth_under_stress", None)
    return oyster(candidates).set_index("feature_id_standardized")


def test_demo_lactate_and_ldh_converge_by_link_only(isolated_reports):
    candidates = _demo_growth_candidates()
    for feature in (LDH, "lactate"):
        assert candidates.loc[feature, "tier"] == "metabolite_linked_convergence"
        assert candidates.loc[feature, "link_confidences"] == "inferred"
    # sdha is linked to a significant succinate change, but is not significant itself.
    assert candidates.loc[SDHA, "tier"] == "emerging"
    assert candidates.loc["succinate", "tier"] == "emerging"
    # glucose <-> simulated glut1 converge through a transporter link, not a reaction.
    for feature in (GLUT1, "glucose"):
        assert candidates.loc[feature, "tier"] == "metabolite_linked_convergence"
        assert candidates.loc[feature, "linked_features"] == "glucose→LOC105342003 via C00031/transport/K07299"
        assert candidates.loc[feature, "link_confidences"] == "inferred"
    assert candidates.loc["ATP", "tier"] == "emerging"
    assert not (candidates["tier"] == "multi_omics_convergence").any()


def test_demo_link_examples_are_not_on_the_significance_boundary(isolated_reports):
    """A demo convergence case that sits at adj. p = 0.05 would flip on a change
    from <= to <, and silently stop demonstrating anything."""
    harmonize_all_demo_studies()
    evidence = oyster(load_evidence_table()).set_index("feature_id_standardized")
    for feature in ("lactate", "glucose", LDH, GLUT1):
        adjusted = float(evidence.loc[feature, "adjusted_p_value"])
        assert adjusted <= 0.04, (feature, adjusted)


def test_demo_cards_explain_the_link_and_its_exclusions(isolated_reports):
    harmonize_all_demo_studies()
    index = evidence_cards.build_evidence_cards(phenotype="growth_under_stress")
    cards = {(e["feature_id_standardized"], int(e["species_taxid"])): Path(e["card_file"]).read_text()
             for e in index}
    lactate = cards[("lactate", 29159)]
    assert "not by a shared identifier" in lactate
    assert "Metabolite–gene links (biochemical; weaker than shared gene identity)" in lactate
    assert "metabolite_gene_links_demo.tsv (sha256 " in lactate
    assert ("| GIGAS_HEAT01 | transcriptomics | LOC105342001 | C00186 | enzyme | R00703 | K00016 | left | "
            "kegg_reaction | inferred |") in lactate
    assert "Transporter links" not in lactate  # the demo table has no transporter rows
    assert "different stressor" in lactate
    assert "`C00256`: no enzyme or transporter link" in lactate

    glucose = cards[("glucose", 29159)]
    assert "**Transporter links** are weaker still" in glucose
    assert "SLC2: Facilitative GLUT transporter" in glucose
    assert ("| GIGAS_HEAT01 | transcriptomics | LOC105342003 | C00031 | transporter | transport | K07299 | — | "
            "ko_definition | inferred |") in glucose

    atp = cards[("ATP", 29159)]
    assert "`C00002`: currency metabolite" in atp
    assert "not by a shared identifier" not in atp

    candidates = pd.read_csv(isolated_reports["cards_dir"] / "candidates.tsv", sep="\t")
    for column in ("convergence_basis", "linked_layers", "linked_features", "link_confidences",
                   "link_table_status", "link_table_sha256"):
        assert column in candidates.columns


# --------------------------------------------------------------------------- #
# Builder (no network: a fake KEGG)
# --------------------------------------------------------------------------- #


def test_equation_side():
    assert equation_side("C00186 + C00003 <=> C00022 + C00004", "C00186") == "left"
    assert equation_side("C00186 + C00003 <=> C00022 + C00004", "C00022") == "right"
    assert equation_side("C00031 <=> C00031", "C00031") == "both"
    assert equation_side("C00031 <=> C00022", "C00186") is None


def test_parse_brite_records_every_family_path():
    from mappings.build_metabolite_links import parse_brite
    families = parse_brite(FAKE_KO02000 + "ABC transporters\nB  Other\nC    K14388  SMCT again\n")
    assert families["K14388"] == ("Solute carrier family (SLC) > SLC5: Sodium glucose cotransporter; "
                                  "BC transporters > Other")
    assert families["K08180"].endswith("SLC16: Monocarboxylate transporter")


def test_transporter_links_are_always_inferred_and_capped_separately():
    compound_map = {"glucose": [{"kegg_compound": "C00031", "match_basis": "exact_name"}]}
    rows = [_row("C00031", f"E{i}", ko=f"K{i}") for i in range(3)]
    rows += [{**_row("C00031", f"T{i}", reaction="", side=""), "link_type": "transporter",
              "link_basis": "ko_definition", "brite_family": "SLC2"} for i in range(5)]
    links, exclusions = ml.build_links(rows, compound_map, {}, max_genes_per_compound=4)

    enzyme = [link for link in links if link.link_type == "enzyme"]
    transport = [link for link in links if link.link_type == "transporter"]
    # Five transporters exceed the cap; the three enzyme links keep their credit.
    assert all(link.credited for link in enzyme)
    assert {link.link_confidence for link in enzyme} == {"one_to_many_ortholog"}
    assert not any(link.credited for link in transport)
    assert {link.link_confidence for link in transport} == {"inferred"}
    assert [(e.reason, e.link_type) for e in exclusions] == [("hub_compound", "transporter")]
    assert transport[0].describe() == "glucose→T0 via C00031/transport/K1"


def test_curated_transporter_rows_are_well_formed():
    rows = ml.load_transporter_substrates()
    curated_compounds = {c["kegg_compound"] for cs in ml.load_compound_map().values() for c in cs}
    assert rows
    for row in rows:
        assert re.fullmatch(r"K\d{5}", row["ko"]), row
        assert row["kegg_compound"] in curated_compounds, row
        assert row["kegg_compound"] not in ml.load_currency_compounds(), row
        assert row["substrate_basis"] in {"ko_definition", "substrate_class"}
        assert row["human_ortholog"] and row["note"]


def test_parse_entries_reads_first_value_per_entry():
    text = ("ENTRY       R00703                      Reaction\nNAME        (S)-lactate:NAD+ oxidoreductase\n"
            "EQUATION    C00186 + C00003 <=> C00022 + C00004 + C00080\n///\n"
            "ENTRY       C00186                      Compound\nNAME        (S)-Lactate;\n"
            "            L-Lactate;\n///\n")
    parsed = parse_entries(text, ("EQUATION", "NAME"))
    assert parsed["R00703"]["EQUATION"].startswith("C00186 + C00003")
    assert parsed["C00186"]["NAME"] == "(S)-Lactate;"


FAKE_KO02000 = """+D
!
ASolute carrier family (SLC)
B  SLC5: Sodium glucose cotransporter
C    K14388  SLC5A8_12, SMCT; solute carrier family 5 (sodium-coupled monocarboxylate transporter), member 8/12
B  SLC16: Monocarboxylate transporter
C    K08180  SLC16A3; MFS transporter, MCT family, solute carrier family 16 (monocarboxylic acid transporters), member 3
!
"""

FAKE_KEGG = {
    "info/kegg": "kegg\tKEGG (Kyoto Encyclopedia of Genes and Genomes)\n\tko          28,545  2026/10/16\n",
    "link/ko/crg": ("crg:1001\tko:K03778\ncrg:1002\tko:K03778\ncrg:9999\tko:K00234\n"
                    "crg:tRNA-Arg\tko:K00234\ncrg:1003\tko:K00873\ncrg:1005\tko:K14388\n"),
    "get/br:ko02000": FAKE_KO02000,
    "link/reaction/cpd:C00256": "cpd:C00256\trn:R00704\n",
    "link/reaction/cpd:C00042": "cpd:C00042\trn:R02164\n",
}


def _fake_fetch(calls):
    def fetch(operation):
        calls.append(operation)
        if operation in FAKE_KEGG:
            return FAKE_KEGG[operation]
        if operation.startswith("link/reaction/"):
            return ""
        if operation.startswith("link/ko/"):
            return "rn:R00704\tko:K03778\nrn:R02164\tko:K00234\n"
        if operation.startswith("get/rn:"):
            return ("ENTRY       R00704\nEQUATION    C00256 + C00003 <=> C00022 + C00004\n///\n"
                    "ENTRY       R02164\nEQUATION    C15602 + C00042 <=> C15603 + C00122\n///\n")
        if operation.startswith("get/cpd:"):
            return "ENTRY       C00256\nNAME        (R)-Lactate;\n///\nENTRY       C00042\nNAME        Succinate;\n///\n"
        raise AssertionError(f"unexpected KEGG call {operation}")
    return fetch


def test_builder_writes_a_resolved_table_with_provenance(tmp_path):
    calls = []
    out = tmp_path / "links.tsv"
    # 9999 is a retired GeneID replaced by 1004; "tRNA-Arg" resolves to nothing.
    resolver = ({"1001", "1002", "1003", "1004", "1005"}, {"9999": "1004"}, tmp_path / "crosswalk.tsv")
    build(taxid=29159, out_path=out, fetch=_fake_fetch(calls), cache_dir=tmp_path / "cache",
          gene_resolver=resolver)

    table = pd.read_csv(out, sep="\t", comment="#", dtype=str)
    assert set(table.columns) == set(ml.LINK_COLUMNS)
    enzyme = table[table["link_type"] == "enzyme"]
    assert set(zip(enzyme["kegg_compound"], enzyme["ncbi_gene_id"])) == {
        ("C00256", "1001"), ("C00256", "1002"), ("C00042", "1004")}
    assert set(enzyme.loc[enzyme["kegg_compound"] == "C00256", "equation_side"]) == {"left"}

    # Transporter links: the curated SMCT row, for every curated compound it carries.
    transport = table[table["link_type"] == "transporter"]
    assert set(zip(transport["kegg_compound"], transport["ko"], transport["ncbi_gene_id"])) == {
        ("C00186", "K14388", "1005"), ("C00256", "K14388", "1005"), ("C00022", "K14388", "1005")}
    assert set(transport["link_basis"]) == {"substrate_class"}
    assert set(transport["brite_family"]) == {
        "Solute carrier family (SLC) > SLC5: Sodium glucose cotransporter"}
    assert transport["kegg_reaction"].isna().all()
    # Currency compounds are never queried.
    assert "link/reaction/cpd:C00002" not in calls

    prov = load_json(ml.provenance_path(out))
    jsonschema.validate(prov, load_json(SCHEMAS_DIR / "metabolite_gene_links_provenance.schema.json"))
    assert prov["artifact_sha256"] == sha256_file(out)
    assert prov["coverage"]["kegg_genes_retired_remapped"] == 1
    assert prov["unresolved_kegg_gene_ids"] == ["tRNA-Arg"]
    assert "C00002" in prov["currency_metabolites"]["skipped"]
    assert all(len(s["sha256"]) == 64 for s in prov["sources"])
    transporters = prov["transporter_substrates"]
    # Curated KOs absent from (this fake) ko02000 are rejected and listed, never silently used.
    assert {"ko": "K07299", "kegg_compound": "C00031", "reason": "not in KEGG BRITE ko02000"} in \
        transporters["rejected_not_in_ko02000"]
    assert "K08180" in transporters["kos_without_crg_genes"]
    assert prov["coverage"]["rows_by_link_type"]["transporter"] == 3

    # Offline rebuild from the cache reproduces the same table.
    out2 = tmp_path / "links2.tsv"
    build(taxid=29159, out_path=out2, offline=True, cache_dir=tmp_path / "cache", gene_resolver=resolver,
          fetch=lambda op: (_ for _ in ()).throw(AssertionError("network used offline")))
    body = [line for line in out.read_text().splitlines() if not line.startswith("#")]
    assert body == [line for line in out2.read_text().splitlines() if not line.startswith("#")]
    assert json.loads(ml.provenance_path(out2).read_text())["artifact_sha256"] == sha256_file(out2)


def test_card_explains_a_transporter_link(monkeypatch):
    from prioritize.rank import EvidenceIndex, add_partition_columns

    transporter_row = {**_row("C00186", "GENE_MCT", reaction="", side="", ko="K08180"),
                       "link_type": "transporter", "link_basis": "substrate_class",
                       "brite_family": "Solute carrier family (SLC) > SLC16: Monocarboxylate transporter"}
    links, exclusions = ml.build_links(
        [transporter_row], {"lactate": [{"kegg_compound": "C00186", "match_basis": "exact_name"}]}, {})
    table = ml.MetaboliteGeneLinks(path=Path("synthetic.tsv"), status=ml.STATUS_AVAILABLE,
                                   sha256="0" * 64, source="test", links=links, exclusions=exclusions)
    monkeypatch.setattr(rank, "links_for_partition", lambda simulated, taxid: table)
    monkeypatch.setattr(evidence_cards, "links_for_partition", lambda simulated, taxid: table)

    evidence = _evidence([
        ("rna", "RNA", "GENE_MCT", "gene", "growth_under_stress", 0.01),
        ("met", "MET", "lactate", "metabolite_feature", "growth_under_stress", 0.02),
    ]).assign(stressor="hypoxia", molecular_direction="up", effect_size=1.0)
    candidate = build_candidates(_candidate("lactate", "metabolite_feature", "met", "MET"), evidence).iloc[0]
    candidate = candidate.to_dict()
    assert candidate["tier"] == "metabolite_linked_convergence"
    assert candidate["linked_features"] == "lactate→GENE_MCT via C00186/transport/K08180"
    assert candidate["link_confidences"] == "inferred"

    index = EvidenceIndex(add_partition_columns(evidence))
    content = evidence_cards._card_markdown(
        candidate, [], [], "Growth under stress", evidence_cards._linked_omics_context(index, candidate))
    assert "**Transporter links** are weaker still" in content
    assert "`Solute carrier family (SLC) > SLC16: Monocarboxylate transporter`" in content
    assert "| MET |" not in content  # the candidate's own record is not a linked row
    assert "| RNA | transcriptomics | GENE_MCT | C00186 | transporter | transport | K08180 | — | substrate_class | inferred |" in content
