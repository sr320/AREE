"""Tests for real-crosswalk construction and crosswalk selection.

These tests are deliberately network-free: they build crosswalks from small
inline fixtures that reproduce the structural features of the real NCBI and
UniProt exports (multi-valued Ensembl xrefs, one gene with several UniProt
accessions, an accession shared between two genes, genes with no xrefs at all).
"""
import json

import pytest

from harmonize.identifiers import (
    DEMO_CROSSWALK_PATH,
    active_crosswalk_path,
    resolve_identifier,
    set_crosswalk_path,
)
from mappings.build_crosswalk import (
    build_rows,
    coverage_stats,
    parse_gene_history,
    parse_gene_info,
    parse_uniprot,
    write_crosswalk,
    write_provenance,
    write_retired_table,
)

GENE_INFO_FIXTURE = "\t".join(
    [
        "#tax_id", "GeneID", "Symbol", "LocusTag", "Synonyms", "dbXrefs", "chromosome",
        "map_location", "description", "type_of_gene", "Symbol_from_nomenclature_authority",
        "Full_name_from_nomenclature_authority", "Nomenclature_status", "Other_designations",
        "Modification_date", "Feature_type",
    ]
) + "\n" + "\n".join(
    [
        # named gene, single Ensembl xref
        "29159\t100\tCYTB\t-\tcob|cytb\tEnsembl:ENSMGIG00000000100\tMT\t-\tcytochrome b\tprotein-coding\t-\t-\t-\t-\t20240101\t-",
        # unnamed gene, two Ensembl xrefs, a legacy (oyster_v9-style) locus tag
        "29159\t200\tLOC200\tCGI_10000200\t-\tEnsembl:ENSMGIG00000000200|Ensembl:ENSMGIG00000000201\t1\t-\tuncharacterized LOC200\tprotein-coding\t-\t-\t-\t-\t20240101\t-",
        # no xrefs at all; a legacy locus tag
        "29159\t300\tLOC300\tCGI_10000300\t-\t-\t2\t-\tuncharacterized LOC300\tncRNA\t-\t-\t-\t-\t20240101\t-",
        # a different taxid must be ignored
        "6565\t999\tOTHER\t-\t-\t-\t1\t-\twrong species\tprotein-coding\t-\t-\t-\t-\t20240101\t-",
    ]
)

UNIPROT_FIXTURE = "\t".join(
    ["Entry", "Reviewed", "Protein names", "Gene Names", "Gene Names (primary)", "GeneID", "RefSeq"]
) + "\n" + "\n".join(
    [
        # reviewed entry for gene 100
        "P11111\treviewed\tCytochrome b\tCYTB\tCYTB\t100;\t",
        # two unreviewed accessions for gene 200
        "A0A111\tunreviewed\tUncharacterized\t\t\t200;\tXP_001",
        "A0A222\tunreviewed\tUncharacterized\t\t\t200;\tXP_002",
        # an accession shared between genes 100 and 200
        "A0A333\tunreviewed\tShared\t\t\t100;200;\t",
        # an entry with no GeneID xref at all -> contributes nothing
        "A0A444\tunreviewed\tOrphan\t\t\t\t",
    ]
)


GENE_HISTORY_FIXTURE = "\t".join(
    ["#tax_id", "GeneID", "Discontinued_GeneID", "Discontinued_Symbol", "Discontinue_Date"]
) + "\n" + "\n".join(
    [
        # retired id 150 was replaced by current gene 200
        "29159\t200\t150\tLOC150\t20200101",
        # retired id 160 has NO replacement -> must stay unresolvable
        "29159\t-\t160\tLOC160\t20200101",
        # retired id 170 points at a gene that is not in the crosswalk -> dropped
        "29159\t99999\t170\tLOC170\t20200101",
        # wrong taxid
        "6565\t1\t2\tLOC2\t20200101",
    ]
)


@pytest.fixture
def fixture_sources(tmp_path):
    gi = tmp_path / "gene_info.tsv"
    up = tmp_path / "uniprot.tsv"
    gh = tmp_path / "gene_history.tsv"
    gi.write_text(GENE_INFO_FIXTURE + "\n")
    up.write_text(UNIPROT_FIXTURE + "\n")
    gh.write_text(GENE_HISTORY_FIXTURE + "\n")
    return gi, up, gh


@pytest.fixture
def built_crosswalk(tmp_path, fixture_sources):
    gi, up, gh = fixture_sources
    genes = parse_gene_info(gi, 29159)
    rows = build_rows(genes, parse_uniprot(up))
    path = tmp_path / "test_gene_id_crosswalk.tsv"
    write_crosswalk(rows, path, 29159, "Test organism", "TEST_ASM")

    replacements, dead = parse_gene_history(gh, 29159)
    replacements = {k: v for k, v in replacements.items() if v["current_gene_id"] in genes}
    write_retired_table(replacements, dead, tmp_path / "test_retired_gene_ids.tsv", 29159)
    return path, rows


@pytest.fixture
def use_built_crosswalk(built_crosswalk):
    """Point the resolver at the fixture crosswalk, then restore the default."""
    path, rows = built_crosswalk
    set_crosswalk_path(path)
    yield path, rows
    set_crosswalk_path(None)


# --------------------------------------------------------------------------- #
# Source parsing
# --------------------------------------------------------------------------- #


def test_parse_gene_info_filters_to_requested_taxid(fixture_sources):
    genes = parse_gene_info(fixture_sources[0], 29159)
    assert set(genes) == {"100", "200", "300"}
    assert "999" not in genes


def test_parse_gene_info_normalises_ncbi_null_token(fixture_sources):
    genes = parse_gene_info(fixture_sources[0], 29159)
    # NCBI writes '-' for absent values; it must not leak into the crosswalk.
    assert genes["300"]["ensembl_gene_id"] == ""
    assert genes["300"]["gene_synonyms"] == ""


def test_parse_gene_info_collects_multiple_ensembl_xrefs(fixture_sources):
    genes = parse_gene_info(fixture_sources[0], 29159)
    assert genes["200"]["ensembl_gene_id"] == "ENSMGIG00000000200;ENSMGIG00000000201"


def test_locus_id_is_populated_for_every_gene(fixture_sources):
    """LOC<GeneID> is emitted even for genes that now carry an official symbol,
    because studies analysed against an older annotation still report the LOC form."""
    genes = parse_gene_info(fixture_sources[0], 29159)
    assert genes["100"]["locus_id"] == "LOC100"
    assert genes["100"]["gene_symbol"] == "CYTB"


def test_parse_gene_info_rejects_a_non_gene_info_file(tmp_path):
    bad = tmp_path / "bad.tsv"
    bad.write_text("col_a\tcol_b\n1\t2\n")
    with pytest.raises(ValueError, match="gene_info"):
        parse_gene_info(bad, 29159)


def test_parse_uniprot_groups_accessions_and_prefers_reviewed(fixture_sources):
    by_gene = parse_uniprot(fixture_sources[1])
    assert by_gene["100"][0] == "P11111"  # reviewed entry leads
    assert set(by_gene["200"]) == {"A0A111", "A0A222", "A0A333"}
    assert "A0A444" not in {a for accs in by_gene.values() for a in accs}


# --------------------------------------------------------------------------- #
# Assembly
# --------------------------------------------------------------------------- #


def test_crosswalk_has_exactly_one_row_per_gene(built_crosswalk):
    _, rows = built_crosswalk
    gene_ids = [r["ncbi_gene_id"] for r in rows]
    assert gene_ids == sorted(set(gene_ids), key=int)


def test_representative_and_full_uniprot_sets_are_both_recorded(built_crosswalk):
    _, rows = built_crosswalk
    gene200 = next(r for r in rows if r["ncbi_gene_id"] == "200")
    assert gene200["uniprot_accession"] in gene200["uniprot_accessions_all"].split(";")
    assert len(gene200["uniprot_accessions_all"].split(";")) == 3


def test_orthogroup_is_never_fabricated(built_crosswalk):
    _, rows = built_crosswalk
    assert all(r["orthogroup_id"] == "" for r in rows)


def test_coverage_stats_are_honest(built_crosswalk):
    _, rows = built_crosswalk
    stats = coverage_stats(rows)
    assert stats["genes_total"] == 3
    assert stats["with_ensembl_xref"]["n"] == 2
    assert stats["with_uniprot_accession"]["n"] == 2  # gene 300 has none
    assert stats["with_named_symbol"]["n"] == 1  # only CYTB; LOC-form is not a name
    assert stats["with_orthogroup"]["n"] == 0


def test_provenance_records_checksums_and_limitations(tmp_path, built_crosswalk):
    path, rows = built_crosswalk
    prov = tmp_path / "prov.json"
    write_provenance(prov, path, 29159, "Test organism", "TEST_ASM", [], coverage_stats(rows))
    doc = json.loads(prov.read_text())
    assert len(doc["artifact_sha256"]) == 64
    assert doc["ncbi_taxid"] == 29159
    assert doc["known_limitations"]


# --------------------------------------------------------------------------- #
# Resolution against a real-shaped crosswalk
# --------------------------------------------------------------------------- #


def test_numeric_gene_id_and_loc_form_resolve_to_the_same_gene(use_built_crosswalk):
    assert resolve_identifier("100", "ncbi_gene_id").feature_id_standardized == "100"
    loc = resolve_identifier("LOC100", "locus_id")
    assert loc.feature_id_standardized == "100"
    assert loc.mapping_confidence == "exact"


def test_every_uniprot_accession_resolves_not_just_the_representative(use_built_crosswalk):
    """Regression: with one representative accession per row, a study reporting a
    non-representative accession used to fall through to `unresolved`."""
    for acc in ("A0A111", "A0A222"):
        resolved = resolve_identifier(acc, "uniprot_accession")
        assert resolved.mapping_confidence == "exact"
        assert resolved.feature_id_standardized == "200"


def test_accession_shared_between_genes_is_downgraded(use_built_crosswalk):
    resolved = resolve_identifier("A0A333", "uniprot_accession")
    assert resolved.mapping_confidence == "many_to_one_ortholog"


def test_secondary_ensembl_xref_resolves(use_built_crosswalk):
    resolved = resolve_identifier("ENSMGIG00000000201", "ensembl_gene_id")
    assert resolved.mapping_confidence == "exact"
    assert resolved.feature_id_standardized == "200"


def test_empty_orthogroup_surfaces_as_none_not_empty_string(use_built_crosswalk):
    assert resolve_identifier("100", "ncbi_gene_id").orthogroup_id is None


def test_unknown_identifiers_stay_unresolved(use_built_crosswalk):
    assert resolve_identifier("LOC999999", "locus_id").mapping_confidence == "unresolved"
    assert resolve_identifier("Q99999", "uniprot_accession").mapping_confidence == "unresolved"


# --------------------------------------------------------------------------- #
# Retired GeneIDs (studies run against an older annotation)
# --------------------------------------------------------------------------- #


def test_legacy_locus_tags_are_carried_into_the_crosswalk(built_crosswalk):
    _, rows = built_crosswalk
    tags = {r["ncbi_gene_id"]: r["legacy_locus_tags"] for r in rows}
    assert tags == {"100": "", "200": "CGI_10000200", "300": "CGI_10000300"}


def test_legacy_locus_tag_resolves_as_inferred(use_built_crosswalk):
    resolved = resolve_identifier("CGI_10000200", "legacy_locus_tag")
    assert resolved.feature_id_standardized == "200"
    assert resolved.mapping_confidence == "inferred"


def test_legacy_protein_group_resolves_only_when_members_agree(use_built_crosswalk):
    assert resolve_identifier("CGI_10000200,CGI_10000200", "legacy_locus_tag").feature_id_standardized == "200"
    for group in ("CGI_10000200,CGI_10000300", "CGI_10000200,CGI_99999999"):
        resolved = resolve_identifier(group, "legacy_locus_tag")
        assert resolved.feature_id_standardized is None and resolved.mapping_confidence == "unresolved"


def test_unknown_legacy_locus_tag_stays_unresolved(use_built_crosswalk):
    assert resolve_identifier("CGI_10099999", "legacy_locus_tag").mapping_confidence == "unresolved"


def test_parse_gene_history_separates_replaced_from_dead(fixture_sources):
    replacements, dead = parse_gene_history(fixture_sources[2], 29159)
    assert replacements["150"]["current_gene_id"] == "200"
    assert "160" in dead  # discontinued with no replacement
    assert "160" not in replacements
    assert "2" not in replacements  # wrong taxid


def test_parse_gene_history_rejects_a_non_history_file(tmp_path):
    bad = tmp_path / "bad.tsv"
    bad.write_text("a\tb\n1\t2\n")
    with pytest.raises(ValueError, match="gene_history"):
        parse_gene_history(bad, 29159)


def test_retired_gene_id_resolves_to_its_replacement(use_built_crosswalk):
    """A study run on an older annotation reports GeneID 150, which NCBI retired
    in favour of 200. It must resolve, but as `inferred` rather than `exact`,
    because the remap crosses an annotation version."""
    for raw, id_type in (("150", "ncbi_gene_id"), ("LOC150", "locus_id")):
        resolved = resolve_identifier(raw, id_type)
        assert resolved.feature_id_standardized == "200", raw
        assert resolved.mapping_confidence == "inferred", raw


def test_retired_gene_id_without_replacement_stays_unresolved(use_built_crosswalk):
    assert resolve_identifier("160", "ncbi_gene_id").mapping_confidence == "unresolved"


def test_retired_gene_id_pointing_outside_the_crosswalk_stays_unresolved(use_built_crosswalk):
    assert resolve_identifier("170", "ncbi_gene_id").mapping_confidence == "unresolved"


def test_current_gene_ids_are_unaffected_by_the_retired_table(use_built_crosswalk):
    assert resolve_identifier("200", "ncbi_gene_id").mapping_confidence == "exact"


def test_absent_retired_sidecar_is_not_an_error(monkeypatch):
    """The demo crosswalk has no retired sidecar; resolution must still work."""
    monkeypatch.delenv("AREE_CROSSWALK", raising=False)
    set_crosswalk_path(None)
    assert resolve_identifier("LOC105333935", "ncbi_gene_id").mapping_confidence == "exact"
    assert resolve_identifier("LOC000000", "ncbi_gene_id").mapping_confidence == "unresolved"


# --------------------------------------------------------------------------- #
# Crosswalk selection
# --------------------------------------------------------------------------- #


def test_demo_crosswalk_is_the_default(monkeypatch):
    monkeypatch.delenv("AREE_CROSSWALK", raising=False)
    set_crosswalk_path(None)
    assert active_crosswalk_path() == DEMO_CROSSWALK_PATH


def test_environment_variable_selects_the_crosswalk(monkeypatch, built_crosswalk):
    path, _ = built_crosswalk
    set_crosswalk_path(None)
    monkeypatch.setenv("AREE_CROSSWALK", str(path))
    try:
        assert active_crosswalk_path() == path
        assert resolve_identifier("LOC200", "locus_id").feature_id_standardized == "200"
    finally:
        set_crosswalk_path(None)


def test_demo_and_real_crosswalks_are_kept_separate(use_built_crosswalk):
    """The demo LOC numbers collide with real GeneIDs but carry fabricated identities,
    so a demo symbol must not resolve against a real crosswalk."""
    assert resolve_identifier("hsp70", "gene_symbol").mapping_confidence == "unresolved"


def test_missing_crosswalk_raises_an_actionable_error(tmp_path):
    set_crosswalk_path(tmp_path / "does_not_exist.tsv")
    try:
        with pytest.raises(FileNotFoundError, match="build-crosswalk"):
            resolve_identifier("100", "ncbi_gene_id")
    finally:
        set_crosswalk_path(None)


# --------------------------------------------------------------------------- #
# Invariants of the committed real crosswalk (skipped if it has not been built)
# --------------------------------------------------------------------------- #

REAL_CROSSWALK = (
    DEMO_CROSSWALK_PATH.parents[1] / "reference" / "crosswalk" / "mgigas_gene_id_crosswalk.tsv"
)


@pytest.fixture
def real_crosswalk():
    if not REAL_CROSSWALK.exists():
        pytest.skip("real crosswalk not built; run `aree build-crosswalk`")
    set_crosswalk_path(REAL_CROSSWALK)
    yield REAL_CROSSWALK
    set_crosswalk_path(None)


def _committed_crosswalks():
    """Every real crosswalk registered in genome_assemblies.yaml, one test id per species."""
    from common import registered_crosswalks

    return [
        pytest.param(entry, id=entry["path"].name.replace("_gene_id_crosswalk.tsv", ""))
        for entry in registered_crosswalks()
    ]


@pytest.fixture(params=_committed_crosswalks())
def committed_crosswalk(request):
    entry = request.param
    if not entry["path"].exists():
        pytest.skip(f"{entry['path'].name} not built; run `aree build-crosswalk`")
    return entry


def _route(entry) -> str:
    from harmonize.identifiers import crosswalk_route

    return crosswalk_route(entry["path"])


def test_real_crosswalk_gene_ids_are_unique_and_numeric(committed_crosswalk):
    import pandas as pd

    if _route(committed_crosswalk) != "ncbi_gene":
        pytest.skip("GFF-route crosswalk: see test_gff_route_crosswalk_invariants")
    df = pd.read_csv(committed_crosswalk["path"], sep="\t", comment="#", dtype=str)
    assert len(df) > 1000
    assert df["ncbi_gene_id"].is_unique
    assert df["ncbi_gene_id"].str.fullmatch(r"\d+").all()
    assert (df["locus_id"] == "LOC" + df["ncbi_gene_id"]).all()


def test_gff_route_crosswalk_invariants(committed_crosswalk):
    """Anchored on the authors' gene IDs; UniProt only from the recorded homology search."""
    import pandas as pd

    if _route(committed_crosswalk) != "gff_homology":
        pytest.skip("NCBI-route crosswalk")
    path = committed_crosswalk["path"]
    df = pd.read_csv(path, sep="\t", comment="#", dtype=str).fillna("")
    assert len(df) > 1000
    assert df["gene_id"].is_unique and (df["gene_id"] != "").all()
    assert df["ncbi_gene_id"].eq("").all() or df["ncbi_gene_id"].str.fullmatch(r"\d*").all()
    assert set(df["uniprot_evidence"]) <= {"", "diamond_blastp"}
    assert ((df["uniprot_accession"] != "") == (df["uniprot_evidence"] == "diamond_blastp")).all()
    # A transcript or protein ID names exactly one gene.
    for column in ("transcript_ids", "protein_ids"):
        ids = [i for v in df[column] for i in v.split(";") if i]
        assert len(ids) == len(set(ids)), column
    prov = json.loads(path.with_name(f"{path.stem}.provenance.json").read_text())
    assert prov["anchor_column"] == "gene_id"
    thresholds = prov["homology_search"]["thresholds"]
    assert thresholds["max_evalue"] <= 1e-5 and thresholds["min_identity_pct"] >= 25
    assert {s["name"] for s in prov["sources"]} >= {"annotation_gff3", "annotation_proteins", "uniprot_database"}
    assert all(len(s["sha256"]) == 64 for s in prov["sources"])


def test_real_crosswalk_never_invents_orthogroups(committed_crosswalk):
    import pandas as pd

    df = pd.read_csv(committed_crosswalk["path"], sep="\t", comment="#", dtype=str)
    assert df["orthogroup_id"].isna().all()


def test_real_crosswalk_provenance_matches_the_artifact(committed_crosswalk):
    from common import sha256_file

    path = committed_crosswalk["path"]
    prov = json.loads(path.with_name(f"{path.stem}.provenance.json").read_text())
    assert prov["artifact_sha256"] == sha256_file(path), (
        f"{path.name} has changed since its provenance sidecar was written; rebuild it"
    )
    assert prov["ncbi_taxid"] == committed_crosswalk["ncbi_taxid"]


def test_real_crosswalk_provenance_names_the_species_reference_annotation(committed_crosswalk):
    """The annotation recorded at build time is what every evidence record reports as
    identifier_annotation_release, so it must be this species' reference assembly,
    not a default carried over from another species."""
    from common import load_assemblies

    path = committed_crosswalk["path"]
    prov = json.loads(path.with_name(f"{path.stem}.provenance.json").read_text())
    if _route(committed_crosswalk) == "gff_homology":
        # No NCBI reference annotation exists; the annotation named is the
        # registered assembly whose crosswalk this is.
        reference = [a for a in load_assemblies().values()
                     if a["assembly_id"] in committed_crosswalk["assemblies"]]
    else:
        reference = [
            a for a in load_assemblies().values()
            if a["ncbi_taxid"] == committed_crosswalk["ncbi_taxid"] and a.get("is_ncbi_reference")
        ]
    assert len(reference) == 1
    assert reference[0]["ncbi_assembly_accession"] in prov["ncbi_reference_annotation_at_build"]


def test_real_crosswalk_has_a_retired_gene_id_sidecar(committed_crosswalk):
    from harmonize.identifiers import retired_table_path

    if _route(committed_crosswalk) != "ncbi_gene":
        pytest.skip("retired GeneIDs exist only for NCBI-route crosswalks")
    sidecar = retired_table_path(committed_crosswalk["path"])
    assert sidecar is not None and sidecar.exists()


def test_real_crosswalk_resolves_a_known_gene(real_crosswalk):
    """LOC105331241 is a real M. gigas gene (coadhesin). It is also used by the demo
    crosswalk under a fabricated identity, which is why the two files stay separate."""
    from mappings.orthogroups import gene_orthogroup

    resolved = resolve_identifier("LOC105331241", "locus_id")
    assert resolved.mapping_confidence == "exact"
    assert resolved.feature_id_standardized == "105331241"
    # The crosswalk carries no orthogroup; resolution takes it from the active
    # orthogroup run (data/reference/orthogroups/active.yaml), if there is one.
    assert resolved.orthogroup_id == gene_orthogroup("105331241")
