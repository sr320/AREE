"""Offline tests for per-run design labels and the methylKit sample sheet.

Some BioProjects deposit every library under one BioSample, so the design is
only recoverable from a run-level string such as experiment_title. These pin
that parse failing loudly instead of defaulting, and the methylKit coding
refusing levels that are not in the sheet. No network access.
"""
from __future__ import annotations

import pytest


def _runs(*titles):
    return [{"run_accession": f"SRR{i}", "experiment_title": t} for i, t in enumerate(titles)]


PATTERN = r"(?P<phenotype_class>RR|SS)(?P<family>\d+)-R(?P<replicate>\d+)$"


def test_parse_run_labels_extracts_named_groups():
    from intake.ena_samplesheet import parse_run_labels

    out = parse_run_labels(_runs("HiSeq: RR52-R1", "HiSeq: SS05-R3"), "experiment_title", PATTERN)
    assert out["SRR0"] == {"phenotype_class": "RR", "family": "52", "replicate": "1", "_label": "RR52-R1"}
    assert out["SRR1"]["family"] == "05"


def test_parse_run_labels_rejects_unmatched_duplicate_and_unnamed():
    from intake.ena_samplesheet import ENAError, parse_run_labels

    with pytest.raises(ENAError, match="do not match"):
        parse_run_labels(_runs("HiSeq: RR52-R1", "HiSeq: pool"), "experiment_title", PATTERN)
    with pytest.raises(ENAError, match="more than one run"):
        parse_run_labels(_runs("a RR52-R1", "b RR52-R1"), "experiment_title", PATTERN)
    with pytest.raises(ENAError, match="named group"):
        parse_run_labels(_runs("RR52-R1"), "experiment_title", r"RR\d+")


def test_methylkit_sheet_codes_levels_by_run_accession(tmp_path):
    from intake.ena_samplesheet import write_methylkit_samplesheet

    rows = [{"run_accession": "SRR2", "condition": "RR"},
            {"run_accession": "SRR1", "condition": "SS"},
            {"run_accession": "SRR3", "condition": "other"}]
    out = write_methylkit_samplesheet(rows, tmp_path / "mk.csv", "SS", "RR")
    assert (tmp_path / "mk.csv").read_text() == "sample_id,treatment\nSRR1,0\nSRR2,1\n"
    assert out == {"control_level": "SS", "treatment_level": "RR", "n_control": 1, "n_treatment": 1}


def test_methylkit_sheet_rejects_unknown_level(tmp_path):
    from intake.ena_samplesheet import ENAError, write_methylkit_samplesheet

    rows = [{"run_accession": "SRR1", "condition": "RR"}, {"run_accession": "SRR2", "condition": "SS"}]
    with pytest.raises(ENAError, match="not found"):
        write_methylkit_samplesheet(rows, tmp_path / "mk.csv", "control", "RR")


def _design(tmp_path, rows):
    path = tmp_path / "design.tsv"
    path.write_text("isolate\tploidy\ttreatment\n" + "".join("\t".join(r) + "\n" for r in rows))
    return path


ATTRS = {"SAMN1": {"isolate": "D54", "tissue": "ctenidium"},
         "SAMN2": {"isolate": "T59", "tissue": "ctenidium"}}


def test_design_table_joins_on_a_deposited_attribute(tmp_path):
    from intake.ena_samplesheet import read_design_table

    path = _design(tmp_path, [("D54", "diploid", "control"), ("T59", "triploid", "heat")])
    out = read_design_table(path, "isolate", ATTRS)
    assert out == {"SAMN1": {"ploidy": "diploid", "treatment": "control"},
                   "SAMN2": {"ploidy": "triploid", "treatment": "heat"}}


def test_design_table_rejects_partial_or_ambiguous_joins(tmp_path):
    from intake.ena_samplesheet import ENAError, read_design_table

    with pytest.raises(ENAError, match="no row"):
        read_design_table(_design(tmp_path, [("D54", "diploid", "control")]), "isolate", ATTRS)
    with pytest.raises(ENAError, match="match no sample"):
        read_design_table(_design(tmp_path, [("D54", "diploid", "control"), ("T59", "triploid", "heat"),
                                             ("X99", "triploid", "heat")]), "isolate", ATTRS)
    with pytest.raises(ENAError, match="repeats"):
        read_design_table(_design(tmp_path, [("D54", "diploid", "control"), ("D54", "triploid", "heat")]),
                          "isolate", ATTRS)
    with pytest.raises(ENAError, match="no 'sample_name' column"):
        read_design_table(_design(tmp_path, [("D54", "diploid", "control")]), "sample_name", ATTRS)


def test_design_table_overrides_deposited_attributes_and_is_recorded(tmp_path, monkeypatch):
    import json

    import intake.ena_samplesheet as ena

    runs = [{"run_accession": "SRR1", "sample_accession": "SAMN1", "sample_alias": "D54",
             "library_strategy": "RNA-Seq", "library_layout": "SINGLE", "instrument_model": "NovaSeq",
             "read_count": "10", "fastq_ftp": "ftp.x/SRR1.fastq.gz", "fastq_md5": "a" * 32, "fastq_bytes": "5"},
            {"run_accession": "SRR2", "sample_accession": "SAMN2", "sample_alias": "T59",
             "library_strategy": "RNA-Seq", "library_layout": "SINGLE", "instrument_model": "NovaSeq",
             "read_count": "10", "fastq_ftp": "ftp.x/SRR2.fastq.gz", "fastq_md5": "b" * 32, "fastq_bytes": "5"}]
    attrs = {"SAMN1": {"isolate": "D54", "treatment": "deposited"}, "SAMN2": {"isolate": "T59"}}
    monkeypatch.setattr(ena, "fetch_run_report", lambda bp, extra=None: runs)
    monkeypatch.setattr(ena, "fetch_sample_attributes", lambda accs: attrs)

    table = _design(tmp_path, [("D54", "diploid", "control"), ("T59", "triploid", "heat")])
    prov = ena.build_samplesheet("PRJNA1", study_id="X", out_dir=tmp_path / "out",
                                 condition_attributes=["ploidy", "treatment"], extra_attributes=["isolate"],
                                 design_table=table, design_join="isolate", design_source="https://example.org/s")
    assert prov["group_sizes"] == {"diploid_control": 1, "triploid_heat": 1}
    recorded = json.loads((tmp_path / "out" / "ena_provenance.json").read_text())["design_table"]
    assert recorded["join_attribute"] == "isolate"
    assert recorded["source"] == "https://example.org/s"
    assert len(recorded["sha256"]) == 64
    with pytest.raises(ena.ENAError, match="together"):
        ena.build_samplesheet("PRJNA1", study_id="X", out_dir=tmp_path / "o2",
                              condition_attributes=["ploidy"], design_table=table)
