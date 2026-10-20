"""GFF-route crosswalks for genomes with no NCBI Gene records, WP7 (src/mappings/gff_crosswalk.py).

DIAMOND is not run here: a hits table in the builder's output format stands in
for it, so the thresholds, best-hit choice and resolution rules are tested
without the search itself.
"""
from __future__ import annotations

import json

import pytest

from harmonize.identifiers import resolve_identifier, set_crosswalk_path
from mappings.gff_crosswalk import Thresholds, build, parse_gff, parse_hits

# An NCBI GenBank-style GFF (locus tags, rna-/cds- prefixes, protein_id on CDS)
GENBANK_GFF = """##gff-version 3
ctg1\tGenbank\tregion\t1\t9000\t.\t+\t.\tID=ctg1:1..9000
ctg1\tGenbank\tgene\t1\t1000\t.\t+\t.\tID=gene-MCOR_1;Name=MCOR_1;gbkey=Gene;gene_biotype=protein_coding;locus_tag=MCOR_1
ctg1\tGenbank\tmRNA\t1\t1000\t.\t+\t.\tID=rna-MCOR_1;Parent=gene-MCOR_1;locus_tag=MCOR_1
ctg1\tGenbank\tCDS\t1\t900\t.\t+\t0\tID=cds-CAC1.1;Parent=rna-MCOR_1;inference=protein motif:UniProtKB/Swiss-Prot:P11111.2;product=CAC1.1;protein_id=CAC1.1;standard_name=HSPA1A;locus_tag=MCOR_1
ctg1\tGenbank\tgene\t2000\t3000\t.\t+\t.\tID=gene-MCOR_2;Name=MCOR_2;gene_biotype=protein_coding;locus_tag=MCOR_2
ctg1\tGenbank\tmRNA\t2000\t3000\t.\t+\t.\tID=rna-MCOR_2;Parent=gene-MCOR_2;locus_tag=MCOR_2
ctg1\tGenbank\tCDS\t2000\t2600\t.\t+\t0\tID=cds-CAC2.1;Parent=rna-MCOR_2;protein_id=CAC2.1;standard_name=SOD1
ctg1\tGenbank\tmRNA\t2000\t3000\t.\t+\t.\tID=rna-MCOR_2-2;Parent=gene-MCOR_2;locus_tag=MCOR_2
ctg1\tGenbank\tCDS\t2000\t2900\t.\t+\t0\tID=cds-CAC3.1;Parent=rna-MCOR_2-2;protein_id=CAC3.1;standard_name=SOD1
ctg1\tGenbank\tgene\t4000\t5000\t.\t+\t.\tID=gene-MCOR_3;Name=MCOR_3;gene_biotype=protein_coding;locus_tag=MCOR_3
ctg1\tGenbank\tmRNA\t4000\t5000\t.\t+\t.\tID=rna-MCOR_3;Parent=gene-MCOR_3;locus_tag=MCOR_3
ctg1\tGenbank\tCDS\t4000\t4600\t.\t+\t0\tID=cds-CAC4.1;Parent=rna-MCOR_3;protein_id=CAC4.1;standard_name=HSPA1A
ctg1\tGenbank\tgene\t6000\t6200\t.\t+\t.\tID=gene-MCOR_4;Name=MCOR_4;gene_biotype=tRNA;locus_tag=MCOR_4
"""
GENBANK_FAA = (">CAC1.1 HSPA1A [Mytilus coruscus]\nMKKLLV\n>CAC2.1 SOD1\nMAAA\n"
               ">CAC3.1 SOD1\nMAAAAAAA\n>CAC4.1 HSPA1A\nMKKL\n")

# An author (EVidenceModeler-style) GFF: no locus tags; proteins named after mRNAs.
AUTHOR_GFF = """##gff-version 3
chr1\tEVM\tgene\t1\t500\t.\t+\t.\tID=evm.TU.chr1.1;Name=EVM%20prediction%20chr1.1
chr1\tEVM\tmRNA\t1\t500\t.\t+\t.\tID=evm.model.chr1.1;Parent=evm.TU.chr1.1
chr1\tEVM\tCDS\t1\t450\t.\t+\t0\tID=cds.evm.model.chr1.1;Parent=evm.model.chr1.1
"""
AUTHOR_FAA = ">evm.model.chr1.1\nMSTNPKPQRK\n"

HITS = [
    # gene, subject, pident, evalue, bitscore, qcov, scov, title
    ("MCOR_1", "sp|P11111|HSP7A_HUMAN", 72.0, 1e-120, 600.0, 95.0, 90.0,
     "sp|P11111|HSP7A_HUMAN Heat shock 70 kDa protein 1A OS=Homo sapiens OX=9606 GN=HSPA1A"),
    ("MCOR_1", "sp|P22222|HSP72_MOUSE", 70.0, 1e-110, 580.0, 95.0, 90.0,
     "sp|P22222|HSP72_MOUSE Heat shock-related 70 kDa protein 2 OS=Mus musculus OX=10090"),
    ("MCOR_3", "sp|P11111|HSP7A_HUMAN", 65.0, 1e-90, 450.0, 88.0, 85.0,
     "sp|P11111|HSP7A_HUMAN Heat shock 70 kDa protein 1A OS=Homo sapiens OX=9606 GN=HSPA1A"),
    ("MCOR_2", "sp|P33333|SODC_HUMAN", 25.0, 1e-20, 90.0, 80.0, 80.0,     # identity below 30
     "sp|P33333|SODC_HUMAN Superoxide dismutase [Cu-Zn] OS=Homo sapiens OX=9606"),
]


def _write_hits(path, rows):
    path.write_text("".join("\t".join(str(v) for v in r) + "\n" for r in rows))
    return path


@pytest.fixture
def built(tmp_path):
    (tmp_path / "a.gff").write_text(GENBANK_GFF)
    (tmp_path / "a.faa").write_text(GENBANK_FAA)
    (tmp_path / "sprot.fasta").write_text(">sp|P11111|X\nM\n")
    out = build(gff=tmp_path / "a.gff", proteins=tmp_path / "a.faa", uniprot_fasta=tmp_path / "sprot.fasta",
                uniprot_release="test release", uniprot_url=None, taxid=42192, organism="Mytilus coruscus",
                annotation="GCA_TEST (TEST1)", out_dir=tmp_path / "cw", slug="mtest", workdir=tmp_path / "work",
                hits_path=_write_hits(tmp_path / "hits.tsv", HITS), today="2026-10-17")
    set_crosswalk_path(out)
    yield out
    set_crosswalk_path(None)


def test_genbank_gff_genes_are_anchored_on_locus_tags():
    ann = parse_gff_text(GENBANK_GFF)
    assert set(ann.genes) == {"MCOR_1", "MCOR_2", "MCOR_3", "MCOR_4"}
    assert ann.transcript_gene["MCOR_2-2"] == "MCOR_2"
    assert ann.protein_gene["CAC3.1"] == "MCOR_2" and ann.protein_gene["CAC3"] == "MCOR_2"
    assert ann.genes["MCOR_1"]["gene_symbol"] == "HSPA1A"
    assert ann.genes["MCOR_1"]["author_uniprot_accession"] == "P11111"
    assert ann.genes["MCOR_4"]["gene_type"] == "tRNA"


def test_author_gff_without_locus_tags():
    ann = parse_gff_text(AUTHOR_GFF)
    assert set(ann.genes) == {"evm.TU.chr1.1"}
    assert ann.genes["evm.TU.chr1.1"]["gene_symbol"] == "EVM prediction chr1.1"
    assert ann.transcript_gene["evm.model.chr1.1"] == "evm.TU.chr1.1"


def parse_gff_text(text, tmp=None):
    import tempfile
    from pathlib import Path

    with tempfile.TemporaryDirectory() as d:
        p = Path(d) / "x.gff"
        p.write_text(text)
        return parse_gff(p)


def test_thresholds_are_reapplied_to_a_hits_table(tmp_path):
    hits = parse_hits(_write_hits(tmp_path / "h.tsv", HITS), Thresholds())
    assert [h["accession"] for h in hits["MCOR_1"]] == ["P11111", "P22222"]
    assert "MCOR_2" not in hits  # 25% identity fails the 30% threshold
    assert hits["MCOR_1"][0]["protein_name"] == "Heat shock 70 kDa protein 1A"


def test_crosswalk_rows_and_provenance(built):
    import pandas as pd

    df = pd.read_csv(built, sep="\t", comment="#", dtype=str).fillna("").set_index("gene_id")
    assert list(df.index) == ["MCOR_1", "MCOR_2", "MCOR_3", "MCOR_4"]
    assert (df["ncbi_gene_id"] == "").all()
    assert df.loc["MCOR_1", "uniprot_accession"] == "P11111"
    assert df.loc["MCOR_1", "uniprot_accessions_all"] == "P11111;P22222"
    assert df.loc["MCOR_1", "uniprot_evidence"] == "diamond_blastp"
    assert df.loc["MCOR_2", "uniprot_accession"] == ""
    assert df.loc["MCOR_2", "protein_ids"] == "CAC2;CAC2.1;CAC3;CAC3.1"
    assert df.loc["MCOR_4", "protein_ids"] == ""
    prov = json.loads(built.with_name(f"{built.stem}.provenance.json").read_text())
    assert prov["route"] == "gff_homology" and prov["anchor_column"] == "gene_id"
    assert prov["ncbi_taxid"] == 42192
    assert prov["homology_search"]["thresholds"]["min_identity_pct"] == 30.0
    assert prov["homology_search"]["run_by_builder"] is False
    cov = prov["coverage"]
    assert cov["genes_total"] == 4 and cov["genes_with_protein"]["n"] == 3
    assert cov["with_uniprot_accession"]["n"] == 2
    # MCOR_1's author accession P11111 is among our hits; MCOR_3 has no author accession.
    assert cov["author_uniprot_agreement"] == {"genes_with_both": 1, "author_accession_among_our_hits": 1,
                                               "pct": 100.0}


@pytest.mark.parametrize("raw,expected", [
    ("MCOR_2", "MCOR_2"),                 # the authors' gene ID / locus tag
    ("gene-MCOR_2", "MCOR_2"),            # with an NCBI GFF prefix
    ("MCOR_2-2", "MCOR_2"),               # a transcript ID
    ("rna-MCOR_2-2", "MCOR_2"),
    ("CAC3.1", "MCOR_2"),                 # a GenBank protein ID
    ("CAC3", "MCOR_2"),                   # unversioned
])
def test_authors_identifiers_resolve_exact(built, raw, expected):
    r = resolve_identifier(raw, "ncbi_gene_id")
    assert (r.feature_id_standardized, r.mapping_confidence) == (expected, "exact")


def test_a_uniprot_accession_reached_by_homology_is_inferred(built):
    r = resolve_identifier("P22222", "uniprot_accession")
    assert (r.feature_id_standardized, r.mapping_confidence) == ("MCOR_1", "inferred")


def test_a_uniprot_accession_shared_by_two_genes_is_many_to_one(built):
    # P11111 is the best hit of both MCOR_1 and MCOR_3 (two HSP70 copies).
    r = resolve_identifier("P11111", "uniprot_accession")
    assert r.mapping_confidence == "many_to_one_ortholog"


def test_a_shared_author_symbol_is_many_to_one(built):
    assert resolve_identifier("HSPA1A", "gene_symbol").mapping_confidence == "many_to_one_ortholog"
    r = resolve_identifier("SOD1", "gene_symbol")
    assert (r.feature_id_standardized, r.mapping_confidence) == ("MCOR_2", "inferred")


def test_unknown_and_numeric_ids_stay_unresolved(built):
    assert resolve_identifier("MCOR_999", "ncbi_gene_id").mapping_confidence == "unresolved"
    # An empty ncbi_gene_id column must never match anything.
    assert resolve_identifier("", "ncbi_gene_id").mapping_confidence == "unresolved"


def test_no_orthogroup_is_borrowed_from_the_geneid_keyed_run(built, monkeypatch):
    """The active orthogroup run is keyed by NCBI GeneID; an author ID must not consult it."""
    import mappings.orthogroups as og

    monkeypatch.setattr(og, "gene_orthogroup", lambda gene_id: "OG_SHOULD_NOT_APPEAR")
    assert resolve_identifier("MCOR_1", "ncbi_gene_id").orthogroup_id is None


def test_cli_requires_the_gff_inputs(tmp_path):
    from click.testing import CliRunner

    from aree.cli import main

    (tmp_path / "a.gff").write_text(GENBANK_GFF)
    out = CliRunner().invoke(main, ["build-crosswalk", "--taxid", "42192", "--slug", "mtest",
                                    "--from-gff", str(tmp_path / "a.gff")])
    assert out.exit_code == 1 and "--proteins" in out.output and "--uniprot-fasta" in out.output


def test_committed_mcoru_crosswalk_resolves_authors_ids():
    """The real M. coruscus MCOR1.1 crosswalk: a locus tag, a GenBank protein ID and a
    homology-assigned UniProt accession all reach MCOR_54798 (UHRF1)."""
    from common import REPO_ROOT

    path = REPO_ROOT / "data" / "reference" / "crosswalk" / "mcoru_gene_id_crosswalk.tsv"
    if not path.exists():
        pytest.skip("mcoru crosswalk not built")
    set_crosswalk_path(path)
    try:
        assert resolve_identifier("MCOR_54798", "ncbi_gene_id") == resolve_identifier("CAC5422767.1", "ncbi_gene_id")
        r = resolve_identifier("CAC5422767.1", "ncbi_gene_id")
        assert (r.feature_id_standardized, r.mapping_confidence, r.orthogroup_id) == ("MCOR_54798", "exact", None)
        u = resolve_identifier("B6CHA3", "uniprot_accession")
        assert u.feature_id_standardized == "MCOR_54798" and u.mapping_confidence in ("inferred",
                                                                                      "many_to_one_ortholog")
    finally:
        set_crosswalk_path(None)
