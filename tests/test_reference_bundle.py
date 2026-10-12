"""Tests for `aree fetch-reference` (src/intake/reference_bundle.py).

Network-free: NCBI's directory listing, md5checksums.txt and the compressed
files are served from a fake built in tmp_path.
"""
from __future__ import annotations

import gzip
import hashlib
import json
import shutil

import pytest

import intake.reference_bundle as rb

ACC = "GCF_000000001.1"
FOLDER = f"{ACC}_TestAsm1.0"
RELEASE = f"{ACC}-RS_2025_01"

GTF = f"""#gtf-version 2.2
#!genome-build TestAsm1.0
#!genome-build-accession NCBI_Assembly:{ACC}
#!annotation-source NCBI RefSeq {RELEASE}
chr1\tGnomon\tgene\t1\t900\t.\t+\t.\tgene_id "LOC100"; transcript_id ""; db_xref "GeneID:100";
chr1\tGnomon\ttranscript\t1\t900\t.\t+\t.\tgene_id "LOC100"; transcript_id "XM_2.1"; db_xref "GeneID:100";
chr1\tGnomon\texon\t1\t900\t.\t+\t.\tgene_id "LOC100"; transcript_id "XM_2.1"; db_xref "GeneID:100";
chr1\tGnomon\ttranscript\t1\t800\t.\t+\t.\tgene_id "LOC100"; transcript_id "XM_1.1"; db_xref "GeneID:100";
chr1\tGnomon\ttranscript\t1\t800\t.\t+\t.\tgene_id "LOC100"; transcript_id "XM_1.1"; db_xref "GeneID:100";
chr2\tGnomon\ttranscript\t1\t500\t.\t-\t.\tgene_id "LOC200"; transcript_id "XM_3.1"; db_xref "GeneID:200";
chr2\tRefSeq\ttranscript\t1\t500\t.\t-\t.\tgene_id "orphan"; transcript_id "XR_9.1";
"""
RNA = ">XM_1.1\nACGT\n>XM_2.1\nACGTACGT\n>XM_3.1\nTTTT\n"


def _md5(path):
    return hashlib.md5(path.read_bytes()).hexdigest()  # noqa: S324


@pytest.fixture
def fake_ncbi(tmp_path, monkeypatch):
    """A fake NCBI folder plus a registered RefSeq assembly with a crosswalk."""
    remote = tmp_path / "ncbi"
    remote.mkdir()
    for suffix, text in (("_rna.fna.gz", RNA), ("_genomic.gtf.gz", GTF),
                         ("_genomic.fna.gz", ">chr1\nACGT\n>chr2\nACGT\n")):
        with gzip.open(remote / f"{FOLDER}{suffix}", "wt") as fh:
            fh.write(text)
    md5_text = "".join(f"{_md5(p)}  ./{p.name}\n" for p in sorted(remote.iterdir()))

    crosswalk = tmp_path / "test_gene_id_crosswalk.tsv"
    crosswalk.write_text("# test\nncbi_gene_id\tgene_symbol\n100\tLOC100\n")

    record = {
        "assembly_id": "TestAsm1.0", "species": "Testus testus", "ncbi_taxid": 1,
        "ncbi_assembly_accession": ACC, "annotation_release": RELEASE,
        "identifier_crosswalk": str(crosswalk),
    }
    monkeypatch.setattr(rb, "resolve_assembly", lambda value: record if value in (ACC, "TestAsm1.0") else None)

    listing = f'<a href="{FOLDER}/">{FOLDER}/</a> <a href="../">../</a>'
    calls = {"downloads": 0}

    def get_text(url):
        if url.endswith("md5checksums.txt"):
            return md5_text
        assert url == rb.accession_dir_url(ACC), url
        return listing

    def download(url, dest):
        calls["downloads"] += 1
        shutil.copy(remote / url.rsplit("/", 1)[1], dest)

    return {"record": record, "remote": remote, "get_text": get_text, "download": download,
            "calls": calls, "dest": tmp_path / "bundle"}


def _fetch(fake, **kwargs):
    return rb.fetch_bundle(ACC, dest=fake["dest"], get_text=fake["get_text"],
                           download=fake["download"], today="2026-10-12", **kwargs)


# --------------------------------------------------------------------------- #
# Pure helpers
# --------------------------------------------------------------------------- #


def test_accession_dir_url():
    assert rb.accession_dir_url("GCF_965363235.1") == (
        "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/965/363/235/")
    with pytest.raises(rb.ReferenceBundleError):
        rb.accession_dir_url("xbMytGall1.hap1.1")


def test_folder_is_read_from_the_listing_not_reconstructed():
    html = '<a href="GCF_965363235.1_xbMytGall1.hap1.1/">x</a><a href="GCF_965363235.2_other/">y</a>'
    assert rb.find_assembly_folder(html, "GCF_965363235.1") == "GCF_965363235.1_xbMytGall1.hap1.1"
    with pytest.raises(rb.ReferenceBundleError, match="expected one folder"):
        rb.find_assembly_folder("<html></html>", "GCF_965363235.1")


def test_parse_md5_listing_strips_dot_slash():
    assert rb.parse_md5_listing("abc  ./x_rna.fna.gz\nbad line here\n") == {"x_rna.fna.gz": "abc"}


def test_tx2gene_is_sorted_unique_and_skips_transcripts_without_a_gene(tmp_path):
    gtf = tmp_path / "genomic.gtf"
    gtf.write_text(GTF)
    # XR_9.1 has no GeneID: the original awk would have given it gene 200
    # (its variables were never reset); it must be skipped instead.
    assert rb.build_tx2gene(gtf) == [("XM_1.1", "100"), ("XM_2.1", "100"), ("XM_3.1", "200")]


@pytest.mark.parametrize(
    "registered,source,expected",
    [
        (RELEASE, f"NCBI RefSeq {RELEASE}", True),
        ("NCBI Mytilus californianus Annotation Release 100",
         "NCBI Mytilus californianus Annotation Release 100", True),
        (RELEASE, f"NCBI RefSeq {ACC}-RS_2026_03", False),
        (RELEASE, None, False),
    ],
)
def test_annotation_release_matching(registered, source, expected):
    assert rb.annotation_release_matches(registered, source) is expected


# --------------------------------------------------------------------------- #
# Fetch, verify, record
# --------------------------------------------------------------------------- #


def test_fetch_writes_the_bundle_and_its_provenance(fake_ncbi):
    prov = _fetch(fake_ncbi)
    dest = fake_ncbi["dest"]

    assert (dest / "rna.fna").read_text() == RNA
    assert (dest / "tx2gene.tsv").read_text() == "XM_1.1\t100\nXM_2.1\t100\nXM_3.1\t200\n"
    assert not (dest / "genomic.fna").exists()

    on_disk = json.loads((dest / rb.PROVENANCE_NAME).read_text())
    assert on_disk == prov
    assert prov["annotation_release_matches_registry"] is True
    assert prov["tx2gene"] == {"transcripts": 3, "genes": 2}
    assert {f["file"] for f in prov["files"]} == {"rna.fna", "genomic.gtf", "tx2gene.tsv"}
    assert all(f["source_md5_ncbi"] for f in prov["files"] if "source_file" in f)
    # GeneID 200 is not in the crosswalk: half the bundle's genes would be unresolved.
    cov = prov["crosswalk_coverage"]
    assert (cov["found_in_crosswalk"], cov["pct"], cov["missing_examples"]) == (1, 50.0, ["200"])
    assert rb.check_bundle(dest) == []


def test_genome_is_fetched_only_when_asked(fake_ncbi):
    prov = _fetch(fake_ncbi, with_genome=True)
    assert (fake_ncbi["dest"] / "genomic.fna").exists()
    assert "genomic.fna" in {f["file"] for f in prov["files"]}


def test_verified_files_are_reused_not_downloaded_again(fake_ncbi):
    _fetch(fake_ncbi)
    first = fake_ncbi["calls"]["downloads"]
    _fetch(fake_ncbi, with_genome=True)
    assert fake_ncbi["calls"]["downloads"] == first + 1  # only the genome


def test_md5_mismatch_is_refused_and_the_file_deleted(fake_ncbi):
    def corrupt(url, dest):
        dest.write_bytes(b"not what NCBI published")

    with pytest.raises(rb.ReferenceBundleError, match="does not match NCBI"):
        rb.fetch_bundle(ACC, dest=fake_ncbi["dest"], get_text=fake_ncbi["get_text"],
                        download=corrupt)
    assert not list(fake_ncbi["dest"].glob("*.gz"))


def test_check_detects_a_modified_file(fake_ncbi):
    _fetch(fake_ncbi)
    (fake_ncbi["dest"] / "tx2gene.tsv").write_text("XM_1.1\t999\n")
    assert any("tx2gene.tsv" in p for p in rb.check_bundle(fake_ncbi["dest"]))


def test_check_without_provenance_fails(tmp_path):
    assert any("missing" in p for p in rb.check_bundle(tmp_path))


def test_a_reannotated_assembly_is_recorded_and_flagged(fake_ncbi):
    """NCBI can replace an assembly's annotation in place. The bundle must say so
    rather than let evidence claim the registered release."""
    fake_ncbi["record"]["annotation_release"] = f"{ACC}-RS_2023_01"
    prov = _fetch(fake_ncbi)
    assert prov["annotation_release_matches_registry"] is False
    assert any("does not name the registered release" in p for p in rb.check_bundle(fake_ncbi["dest"]))


def test_refuses_an_assembly_without_refseq_annotation():
    """The real registry: the Roslin oyster assembly has no RefSeq annotation."""
    with pytest.raises(rb.ReferenceBundleError, match="no NCBI RefSeq annotation"):
        rb.plan_bundle("cgigas_uk_roslin_v1", get_text=lambda url: "")
    with pytest.raises(rb.ReferenceBundleError, match="not in data/reference"):
        rb.plan_bundle("made_up_assembly", get_text=lambda url: "")


def test_every_species_reference_assembly_can_have_a_bundle():
    """Raw reanalysis needs a RefSeq transcriptome and GTF, which NCBI publishes only
    for RefSeq-annotated (GCF_) assemblies with an annotation release."""
    from common import load_assemblies

    for record in load_assemblies().values():
        if record.get("is_ncbi_reference"):
            assert record["ncbi_assembly_accession"].startswith("GCF_"), record["assembly_id"]
            assert record.get("annotation_release"), record["assembly_id"]
            rb.accession_dir_url(record["ncbi_assembly_accession"])
