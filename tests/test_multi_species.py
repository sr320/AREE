"""Tests for registering and harmonizing more than one species.

Covers WP0 and WP1 of docs/mussel_clam_integration_plan.md:

* WP0: the mussel and clam species terms, their reference assemblies, the new
  tissue/stressor/phenotype/quality-flag terms, and the validation rules for
  species in a hybridizing complex and for reads mapped to another species'
  genome.
* WP1: each real study resolves identifiers against the crosswalk registered
  for its own genome assembly, so studies of two species harmonized in one
  process never resolve against each other's identifiers.
"""
from __future__ import annotations

import json
import warnings

import pandas as pd
import pytest
import yaml

from common import (
    REPO_ROOT,
    STUDIES_DIR,
    assembly_crosswalk_path,
    load_assemblies,
    load_vocab,
    registered_crosswalks,
    resolve_assembly,
    resolve_species,
)
from harmonize.identifiers import (
    DEMO_CROSSWALK_PATH,
    crosswalk_for_study,
    crosswalk_taxid,
    set_crosswalk_path,
)

MGIGAS_CROSSWALK = REPO_ROOT / "data" / "reference" / "crosswalk" / "mgigas_gene_id_crosswalk.tsv"

NEW_SPECIES = {
    "Mytilus galloprovincialis": 29158,
    "Mytilus edulis": 6550,
    "Mytilus trossulus": 6551,
    "Mytilus californianus": 6549,
    "Mercenaria mercenaria": 6596,
    "Ruditapes philippinarum": 129788,
}


@pytest.fixture(autouse=True)
def _reset_crosswalk(monkeypatch):
    monkeypatch.delenv("AREE_CROSSWALK", raising=False)
    set_crosswalk_path(None)
    yield
    set_crosswalk_path(None)


# --------------------------------------------------------------------------- #
# WP0: species and reference data
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize("name,taxid", NEW_SPECIES.items())
def test_new_species_resolve_to_their_taxid(name, taxid):
    resolved = resolve_species(name)
    assert resolved is not None and resolved.ncbi_taxid == taxid
    assert not resolved.is_synonym


@pytest.mark.parametrize(
    "name,taxid",
    [
        ("Venerupis philippinarum", 129788),
        ("Tapes philippinarum", 129788),
        ("Ruditapes (Venerupis) philippinarum", 129788),
        ("Venus mercenaria", 6596),
        ("Cardium nuttallii", 410350),
    ],
)
def test_ncbi_listed_synonyms_resolve(name, taxid):
    resolved = resolve_species(name)
    assert resolved is not None and resolved.ncbi_taxid == taxid
    assert resolved.is_synonym


@pytest.mark.parametrize(
    "name",
    ["Mytilus edulis balthicus", "Mytilus edulis galloprovincialis", "Mytilus chilensis", "Panopea generosa"],
)
def test_names_not_deliberately_added_do_not_resolve(name):
    """A subspecies, an informal combination, and species with no RefSeq genome
    are not guessed onto an existing term."""
    assert resolve_species(name) is None


def test_mytilus_edulis_complex_members():
    members = {
        term["scientific_name"]
        for term in load_vocab("species").values()
        if term.get("complex") == "mytilus_edulis_complex"
    }
    assert members == {"Mytilus edulis", "Mytilus galloprovincialis", "Mytilus trossulus"}
    assert resolve_species("Mytilus californianus").complex is None


def test_every_complex_has_at_least_two_members():
    counts: dict[str, int] = {}
    for term in load_vocab("species").values():
        if term.get("complex"):
            counts[term["complex"]] = counts.get(term["complex"], 0) + 1
    assert all(n >= 2 for n in counts.values()), counts


def test_every_species_has_exactly_one_reference_assembly():
    """One NCBI reference assembly per species, except a species NCBI has not
    annotated (GFF-route crosswalk, WP7) or one with no annotated assembly at
    all yet (the basket cockle), which has none."""
    from harmonize.identifiers import crosswalk_route

    species_taxids = {t["ncbi_taxid"] for t in load_vocab("species").values()}
    gff_route = {e["ncbi_taxid"] for e in registered_crosswalks() if crosswalk_route(e["path"]) == "gff_homology"}
    annotated = {a["ncbi_taxid"] for a in load_assemblies().values() if a.get("annotation_release")}
    reference = {}
    for record in load_assemblies().values():
        if record.get("is_ncbi_reference"):
            reference.setdefault(record["ncbi_taxid"], []).append(record["assembly_id"])
    for taxid in species_taxids:
        expected = 0 if taxid in gff_route or taxid not in annotated else 1
        assert len(reference.get(taxid, [])) == expected, f"taxid {taxid}: {reference.get(taxid)}"


def test_every_assembly_belongs_to_a_registered_species():
    """An assembly's `species` must be the canonical name of the term with its taxid."""
    by_taxid = {t["ncbi_taxid"]: t for t in load_vocab("species").values()}
    for record in load_assemblies().values():
        term = by_taxid.get(record["ncbi_taxid"])
        assert term is not None, record["assembly_id"]
        assert record["species"] == term["scientific_name"], record["assembly_id"]


def test_every_assembly_records_a_verification_date_and_crosswalk_field():
    for record in load_assemblies().values():
        assert record.get("verified"), record["assembly_id"]
        assert "identifier_crosswalk" in record, record["assembly_id"]


@pytest.mark.parametrize(
    "value,expected",
    [
        ("GCF_965363235.1", "xbMytGall1.hap1.1"),
        ("GCA_963676685.2", "xbMytEdul2.2"),
        ("PNRI_Mtr1.1.1.hap1", "PNRI_Mtr1.1.1.hap1"),
        ("GCF_021730395.1", "MADL_Memer_1"),
        ("GCF_014805675.1", "ASM1480567v1.1"),
        ("GCF_026571515.1", "ASM2657151v2"),
        ("GCA_977929165.1", "xbCliNutt1.hap2.1"),
    ],
)
def test_new_assemblies_resolve(value, expected):
    assert resolve_assembly(value)["assembly_id"] == expected


def test_basket_cockle_is_registered_without_a_crosswalk():
    """Registered ahead of any annotation: a study can name it, but nothing can
    resolve identifiers against it until a crosswalk is built."""
    assert resolve_species("Clinocardium nuttallii").ncbi_taxid == 410350
    assembly = resolve_assembly("xbCliNutt1.hap2.1")
    assert assembly["is_ncbi_reference"] is False
    assert assembly["identifier_crosswalk"] is None
    assert assembly_crosswalk_path(assembly) is None


def test_registered_crosswalks_match_their_assemblys_taxon():
    """A crosswalk registered for an assembly must have been built for its taxid."""
    for entry in registered_crosswalks():
        if not entry["path"].exists():
            continue
        assert crosswalk_taxid(entry["path"]) == entry["ncbi_taxid"], entry


def test_both_oyster_assemblies_share_the_oyster_crosswalk():
    entries = [e for e in registered_crosswalks() if e["ncbi_taxid"] == 29159]
    assert len(entries) == 1
    assert entries[0]["path"] == MGIGAS_CROSSWALK
    assert set(entries[0]["assemblies"]) == {"xbMagGiga1.1", "cgigas_uk_roslin_v1"}


def test_every_registered_assembly_has_a_built_crosswalk():
    """WP2 built a crosswalk for every annotated species, so every annotated
    assembly can be harmonized. A registered path that does not exist would make
    validation and harmonization disagree about what is possible. An assembly
    with no annotation at all (the basket cockle) has nothing to build from."""
    for record in load_assemblies().values():
        if not record.get("annotation_release"):
            if record["ncbi_taxid"] not in {
                a["ncbi_taxid"] for a in load_assemblies().values() if a.get("annotation_release")
            }:
                assert record["identifier_crosswalk"] is None, record["assembly_id"]
                continue
        path = assembly_crosswalk_path(record)
        assert path is not None, record["assembly_id"]
        assert path.exists(), f"{record['assembly_id']}: {path} not built"


def test_assemblies_of_one_species_share_one_crosswalk():
    """The crosswalk is built per taxid, so two assemblies of one species must not
    point at two different files (both oyster and both hard clam assemblies)."""
    by_taxid: dict[int, set] = {}
    for record in load_assemblies().values():
        by_taxid.setdefault(record["ncbi_taxid"], set()).add(assembly_crosswalk_path(record))
    assert all(len(paths) == 1 for paths in by_taxid.values()), by_taxid


@pytest.fixture
def assembly_without_crosswalk(monkeypatch):
    """Pretend the Manila clam assembly has no crosswalk, as every new species
    starts out. Patched where each caller looks it up."""
    import harmonize.identifiers as identifiers
    import intake.schema_validate as schema_validate

    real = resolve_assembly

    def patched(value):
        record = real(value)
        if record and record["assembly_id"] == "ASM2657151v2":
            record = {**record, "identifier_crosswalk": None}
        return record

    monkeypatch.setattr(identifiers, "resolve_assembly", patched)
    monkeypatch.setattr(schema_validate, "resolve_assembly", patched)


@pytest.mark.parametrize(
    "vocab,term",
    [
        ("tissue_types", "hemocytes"),
        ("tissue_types", "foot"),
        ("tissue_types", "byssal_gland"),
        ("tissue_types", "siphon"),
        ("stressor_ontology", "emersion_desiccation"),
        ("phenotype_ontology", "byssal_attachment"),
        ("phenotype_ontology", "burrowing_performance"),
        ("quality_flags", "species_assignment_unverified"),
        ("quality_flags", "reference_species_mismatch"),
    ],
)
def test_new_vocabulary_terms_exist(vocab, term):
    assert term in load_vocab(vocab)


def test_new_phenotypes_carry_a_resilience_relevance():
    vocab = load_vocab("phenotype_ontology")
    for term in ("byssal_attachment", "burrowing_performance"):
        assert vocab[term]["resilience_relevance"] == "resilience"


# --------------------------------------------------------------------------- #
# WP0: validation rules
# --------------------------------------------------------------------------- #


def _study(tmp_path, study_id, **overrides):
    base = yaml.safe_load((STUDIES_DIR / "HESSER2024_VCOR.yaml").read_text())
    base.update(study_id=study_id, **overrides)
    path = tmp_path / f"{study_id}.yaml"
    path.write_text(yaml.safe_dump(base))
    return path


def _validate(path):
    from intake.schema_validate import validate_study_file

    return validate_study_file(path)


MUSSEL = {"species": "Mytilus galloprovincialis", "genome_assembly": "GCF_965363235.1"}


def test_complex_species_without_species_identity_fails(tmp_path):
    result = _validate(_study(tmp_path, "MUS_NOID", **MUSSEL))
    assert not result.valid
    assert any("declare species_identity" in e for e in result.errors)


def test_complex_species_assigned_by_site_needs_the_unverified_flag(tmp_path):
    result = _validate(_study(
        tmp_path, "MUS_SITE", **MUSSEL,
        species_identity={"basis": "collection_site", "notes": None},
        quality_flags=[],
    ))
    assert not result.valid
    assert any("species_assignment_unverified" in e for e in result.errors)


def test_complex_species_assigned_by_site_with_flag_validates(tmp_path):
    result = _validate(_study(
        tmp_path, "MUS_FLAG", **MUSSEL,
        species_identity={"basis": "collection_site", "notes": "Adriatic farm stock"},
        quality_flags=["species_assignment_unverified"],
    ))
    assert result.valid, result.errors


def test_genotyped_complex_species_needs_no_flag(tmp_path):
    result = _validate(_study(
        tmp_path, "MUS_GENO", **MUSSEL,
        species_identity={"basis": "genotyped", "notes": "Me15/16 marker"},
        quality_flags=[],
    ))
    assert result.valid, result.errors


def test_species_identity_basis_is_a_closed_vocabulary(tmp_path):
    result = _validate(_study(
        tmp_path, "MUS_BAD", **MUSSEL, species_identity={"basis": "vibes"},
    ))
    assert not result.valid


def test_species_outside_a_complex_does_not_need_species_identity(tmp_path):
    result = _validate(_study(
        tmp_path, "CLAM_OK", species="Ruditapes philippinarum", genome_assembly="ASM2657151v2",
    ))
    assert result.valid, result.errors


def test_real_study_on_an_assembly_without_a_crosswalk_warns(tmp_path, assembly_without_crosswalk):
    result = _validate(_study(
        tmp_path, "CLAM_WARN", species="Ruditapes philippinarum", genome_assembly="ASM2657151v2",
    ))
    assert result.valid
    assert any("no identifier_crosswalk" in w and "--taxid 129788" in w for w in result.warnings)


def test_reads_mapped_to_another_species_fail_without_the_flag(tmp_path):
    result = _validate(_study(
        tmp_path, "MUS_MISMATCH", species="Mytilus trossulus", genome_assembly="GCF_965363235.1",
        species_identity={"basis": "genotyped"},
    ))
    assert not result.valid
    assert any("reference_species_mismatch" in e for e in result.errors)


def test_reads_mapped_to_another_species_validate_with_the_flag(tmp_path):
    result = _validate(_study(
        tmp_path, "MUS_MISMATCH_OK", species="Mytilus trossulus", genome_assembly="GCF_965363235.1",
        species_identity={"basis": "genotyped"},
        quality_flags=["reference_species_mismatch"],
    ))
    assert result.valid, result.errors
    assert any("reference_species_mismatch" in w for w in result.warnings)


def test_mismatch_flag_without_a_mismatch_fails(tmp_path):
    result = _validate(_study(
        tmp_path, "OYS_FLAG", quality_flags=["reference_species_mismatch"],
    ))
    assert not result.valid


# --------------------------------------------------------------------------- #
# WP1: crosswalk selection follows the study's assembly
# --------------------------------------------------------------------------- #


def _real(study_id, assembly):
    return {"study_id": study_id, "simulated": False, "genome_assembly": assembly}


def test_real_study_on_an_assembly_with_no_crosswalk_fails_with_the_build_command(
    assembly_without_crosswalk,
):
    with pytest.raises(FileNotFoundError, match="build-crosswalk --taxid 129788"):
        crosswalk_for_study(_real("CLAM", "ASM2657151v2"))


@pytest.mark.parametrize(
    "assembly,slug",
    [
        ("xbMytGall1.hap1.1", "mgall"),
        ("xbMytEdul2.2", "medul"),
        ("PNRI_Mtr1.1.1.hap1", "mtros"),
        ("xbMytCali1.0.p", "mcali"),
        ("MADL_Memer_1", "mmerc"),
        ("ASM1480567v1.1", "mmerc"),
        ("ASM2657151v2", "rphil"),
    ],
)
def test_each_mussel_and_clam_assembly_selects_its_own_crosswalk(assembly, slug):
    path = crosswalk_for_study(_real("X", assembly))
    assert path.name == f"{slug}_gene_id_crosswalk.tsv"
    assert crosswalk_taxid(path) == resolve_assembly(assembly)["ncbi_taxid"]


def test_unknown_assembly_fails():
    with pytest.raises(ValueError, match="not in data/reference/genome_assemblies.yaml"):
        crosswalk_for_study(_real("X", "made_up_asm"))


def test_previous_study_does_not_decide_the_next_studys_crosswalk():
    """Regression: the crosswalk selected for one study used to stick for every
    later study in the process, so an oyster study harmonized first would have
    handed its crosswalk to a mussel study harmonized second."""
    if not MGIGAS_CROSSWALK.exists():
        pytest.skip("real crosswalk not built")
    set_crosswalk_path(crosswalk_for_study(_real("OYS", "cgigas_uk_roslin_v1")))
    mussel = crosswalk_for_study(_real("MUS", "xbMytGall1.hap1.1"))
    assert mussel != MGIGAS_CROSSWALK
    assert crosswalk_taxid(mussel) == 29158


def test_override_built_for_another_taxon_is_refused(monkeypatch):
    if not MGIGAS_CROSSWALK.exists():
        pytest.skip("real crosswalk not built")
    monkeypatch.setenv("AREE_CROSSWALK", str(MGIGAS_CROSSWALK))
    with pytest.raises(ValueError, match="built for taxid 29159"):
        crosswalk_for_study(_real("MUS", "xbMytGall1.hap1.1"))


def test_override_without_provenance_is_refused(monkeypatch):
    """The demo crosswalk has no provenance, so it cannot be vouched for."""
    monkeypatch.setenv("AREE_CROSSWALK", str(DEMO_CROSSWALK_PATH))
    with pytest.raises(ValueError, match="no provenance sidecar"):
        crosswalk_for_study(_real("OYS", "cgigas_uk_roslin_v1"))


# --------------------------------------------------------------------------- #
# WP1: two species in one evidence table
# --------------------------------------------------------------------------- #


def _write_crosswalk(directory, slug, taxid, genes):
    path = directory / f"{slug}_gene_id_crosswalk.tsv"
    rows = [
        {"ncbi_gene_id": gid, "ensembl_gene_id": "", "uniprot_accession": "",
         "uniprot_accessions_all": "", "locus_id": f"LOC{gid}", "gene_symbol": sym,
         "orthogroup_id": ""}
        for gid, sym in genes
    ]
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)
    path.with_name(f"{slug}_gene_id_crosswalk.provenance.json").write_text(json.dumps({
        "ncbi_taxid": taxid,
        "ncbi_reference_annotation_at_build": f"TEST_{slug.upper()}_ANNOTATION",
    }))
    return path


def _write_results(path, gene_ids):
    pd.DataFrame({
        "gene_id": gene_ids,
        "log2FoldChange": [1.0] * len(gene_ids),
        "lfcSE": [0.2] * len(gene_ids),
        "pvalue": [0.001] * len(gene_ids),
        "padj": [0.01] * len(gene_ids),
    }).to_csv(path, sep="\t", index=False)
    return path


def _comparison(results_file):
    return {
        "comparison_id": "heat_vs_control", "tissue": "gill", "life_stage": "adult",
        "stressor_original": "test", "stressor_standardized": "temperature",
        "control_condition": "ambient", "treatment_condition": "heat",
        "phenotype": "thermal_tolerance", "phenotype_direction": "decrease",
        "resilience_classification": "resilience", "contrast_type": "exposure_response", "sample_size": 6,
        "biological_replicates": 3, "results_file": str(results_file),
    }


def test_two_species_resolve_against_their_own_crosswalks(tmp_path, monkeypatch, isolated_reports):
    """Oyster and mussel crosswalks reuse the same LOC numbers for different genes.

    Each study must resolve only against its own assembly's crosswalk, the
    evidence must partition by species_taxid, and each harmonization manifest
    must name the crosswalk (and checksum) that study actually used.
    """
    import harmonize.core as core
    import harmonize.identifiers as identifiers

    oyster_cw = _write_crosswalk(tmp_path, "oystertest", 29159, [("1001", "oyster_gene_a")])
    mussel_cw = _write_crosswalk(tmp_path, "musseltest", 29158, [("2002", "mussel_gene_b")])

    assemblies = {
        "OYSTER_TEST_ASM": {"assembly_id": "OYSTER_TEST_ASM", "ncbi_taxid": 29159,
                            "species": "Magallana gigas", "identifier_crosswalk": str(oyster_cw)},
        "MUSSEL_TEST_ASM": {"assembly_id": "MUSSEL_TEST_ASM", "ncbi_taxid": 29158,
                            "species": "Mytilus galloprovincialis", "identifier_crosswalk": str(mussel_cw)},
    }
    monkeypatch.setattr(identifiers, "resolve_assembly", lambda value: assemblies.get(value))

    # Both studies report both IDs; each should resolve exactly one of them.
    studies = {
        "OYS_TEST": {
            "study_id": "OYS_TEST", "simulated": False, "species": "Crassostrea gigas",
            "genome_assembly": "OYSTER_TEST_ASM", "assay_type": ["rnaseq"],
            "analysis_mode": "raw_reanalysis", "analysis_status": "harmonized",
            "comparisons": [_comparison(_write_results(tmp_path / "oys.tsv", ["LOC1001", "LOC2002"]))],
        },
        "MUS_TEST": {
            "study_id": "MUS_TEST", "simulated": False, "species": "Mytilus galloprovincialis",
            "genome_assembly": "MUSSEL_TEST_ASM", "assay_type": ["rnaseq"],
            "analysis_mode": "raw_reanalysis", "analysis_status": "harmonized",
            "comparisons": [_comparison(_write_results(tmp_path / "mus.tsv", ["LOC1001", "LOC2002"]))],
        },
    }
    monkeypatch.setattr(core, "_load_study", lambda sid: studies[sid])

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        core.harmonize_study("OYS_TEST", date_generated="2026-10-11")
        core.harmonize_study("MUS_TEST", date_generated="2026-10-11")

    evidence = core.load_evidence_table()
    oys = evidence[evidence["study_id"] == "OYS_TEST"].set_index("feature_id_original")
    mus = evidence[evidence["study_id"] == "MUS_TEST"].set_index("feature_id_original")

    assert oys.loc["LOC1001", "mapping_confidence"] == "exact"
    assert oys.loc["LOC2002", "mapping_confidence"] == "unresolved"
    assert mus.loc["LOC2002", "mapping_confidence"] == "exact"
    assert mus.loc["LOC1001", "mapping_confidence"] == "unresolved"

    assert set(oys["species_taxid"].astype(int)) == {29159}
    assert set(mus["species_taxid"].astype(int)) == {29158}
    assert set(oys["identifier_annotation_release"]) == {"TEST_OYSTERTEST_ANNOTATION"}
    assert set(mus["identifier_annotation_release"]) == {"TEST_MUSSELTEST_ANNOTATION"}

    manifests = {
        json.loads(p.read_text())["study_id"]: json.loads(p.read_text())
        for p in isolated_reports["manifests_dir"].glob("*.json")
    }
    used = {sid: m["parameters"]["identifier_crosswalk"]["path"] for sid, m in manifests.items()}
    assert used["OYS_TEST"].endswith("oystertest_gene_id_crosswalk.tsv")
    assert used["MUS_TEST"].endswith("musseltest_gene_id_crosswalk.tsv")
    assert (
        manifests["OYS_TEST"]["parameters"]["identifier_crosswalk"]["sha256"]
        != manifests["MUS_TEST"]["parameters"]["identifier_crosswalk"]["sha256"]
    )



# --------------------------------------------------------------------------- #
# build-crosswalk refuses to mislabel or overwrite another species' crosswalk
# --------------------------------------------------------------------------- #


def _build_crosswalk(*args):
    from click.testing import CliRunner

    from aree.cli import main

    return CliRunner().invoke(main, ["build-crosswalk", *args])


def test_build_crosswalk_requires_a_slug_for_other_species():
    """Without this, --taxid 129788 would overwrite mgigas_gene_id_crosswalk.tsv."""
    result = _build_crosswalk("--taxid", "129788")
    assert result.exit_code == 1
    assert "--slug is required" in result.output


def test_build_crosswalk_refuses_an_unregistered_species():
    result = _build_crosswalk("--taxid", "6604", "--slug", "marenaria")
    assert result.exit_code == 1
    assert "not in registry/controlled_vocabularies/species.yaml" in result.output


def test_build_crosswalk_refuses_to_overwrite_another_species_file():
    if not MGIGAS_CROSSWALK.exists():
        pytest.skip("real crosswalk not built")
    result = _build_crosswalk("--taxid", "129788", "--slug", "mgigas")
    assert result.exit_code == 1
    assert "built for taxid 29159" in result.output


def test_no_gene_id_appears_in_two_species_crosswalks():
    """NCBI GeneIDs are unique across taxa. The dashboard searches every crosswalk
    together for gene symbols, which is only safe while this holds."""
    seen: dict[str, str] = {}
    for entry in registered_crosswalks():
        if not entry["path"].exists():
            continue
        ids = pd.read_csv(entry["path"], sep="\t", comment="#", dtype=str, usecols=["ncbi_gene_id"])
        for gid in ids["ncbi_gene_id"].dropna():  # a GFF-route crosswalk has no GeneIDs
            assert gid not in seen, f"GeneID {gid} in both {seen[gid]} and {entry['path'].name}"
            seen[gid] = entry["path"].name
