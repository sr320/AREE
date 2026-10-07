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
