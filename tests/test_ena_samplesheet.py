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


# --------------------------------------------------------------------------- #
# Mixed library strategies (RNA-seq + 16S amplicon from the same animals)
# --------------------------------------------------------------------------- #


def _mixed_project(monkeypatch):
    """Two groups of 3 animals, each with one RNA-seq and one amplicon run, as in
    PRJNA1132140 (Manila clam priming), where counting both doubles every group."""
    import intake.ena_samplesheet as es

    runs, attrs = [], {}
    for i, (group, strategy, layout) in enumerate(
        [(g, s, lay) for g in ("primed", "naive") for s, lay in (("RNA-Seq", "SINGLE"), ("AMPLICON", "PAIRED"))
         for _ in range(3)]
    ):
        sample = f"SAMN{i}"
        runs.append({"run_accession": f"SRR{i}", "sample_accession": sample, "sample_alias": f"{group}{i}",
                     "library_strategy": strategy, "library_layout": layout,
                     "instrument_model": "NovaSeq", "read_count": "1000",
                     "fastq_ftp": f"ftp.example/SRR{i}.fastq.gz", "fastq_md5": "abc", "fastq_bytes": "10"})
        attrs[sample] = {"Treatment": group, "isolate": f"I{i}"}
    monkeypatch.setattr(es, "fetch_run_report", lambda *a, **k: runs)
    monkeypatch.setattr(es, "fetch_sample_attributes", lambda accs: {a: attrs[a] for a in accs})
    return es


def test_mixed_strategies_are_reported_per_strategy_without_filtering(monkeypatch, tmp_path):
    es = _mixed_project(monkeypatch)
    report = es.build_samplesheet("PRJNA1", study_id="X", out_dir=tmp_path,
                                  condition_attributes=["Treatment"])
    assert report["group_sizes"] == {"naive": 6, "primed": 6}  # inflated: both assays counted
    assert report["group_sizes_by_strategy"] == {
        "AMPLICON": {"naive": 3, "primed": 3}, "RNA-Seq": {"naive": 3, "primed": 3}}
    assert "library_strategy_filter" not in report


def test_library_strategy_filter_keeps_only_the_named_assay_and_records_it(monkeypatch, tmp_path):
    import json

    es = _mixed_project(monkeypatch)
    report = es.build_samplesheet("PRJNA1", study_id="X", out_dir=tmp_path,
                                  condition_attributes=["Treatment"], library_strategies=["RNA-Seq"])
    assert report["group_sizes"] == {"naive": 3, "primed": 3}
    assert report["library_layouts"] == ["SINGLE"]
    sheet = (tmp_path / "samplesheet.tsv").read_text().splitlines()
    assert len(sheet) == 7
    manifest = (tmp_path / "fastq_manifest.tsv").read_text().splitlines()
    assert len(manifest) == 7
    prov = json.loads((tmp_path / "ena_provenance.json").read_text())
    assert prov["library_strategy_filter"] == {
        "kept": ["RNA-Seq"], "runs_by_strategy_before_filter": {"AMPLICON": 6, "RNA-Seq": 6}}


def test_unknown_library_strategy_is_an_error(monkeypatch, tmp_path):
    es = _mixed_project(monkeypatch)
    with pytest.raises(es.ENAError, match="not in this project"):
        es.build_samplesheet("PRJNA1", study_id="X", out_dir=tmp_path,
                             condition_attributes=["Treatment"], library_strategies=["WGS"])


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


def test_library_strategy_filter_applies_before_the_design_table_join(monkeypatch, tmp_path):
    """A design table needs rows only for the runs kept: the filter runs first."""
    es = _mixed_project(monkeypatch)
    rnaseq = [i for i in range(12) if (i // 3) % 2 == 0]  # SAMN0-2 (primed), SAMN6-8 (naive)
    table = tmp_path / "design.tsv"
    table.write_text("isolate\tarm\n" + "".join(f"I{i}\t{'P' if i < 6 else 'N'}\n" for i in rnaseq))

    report = es.build_samplesheet("PRJNA1", study_id="X", out_dir=tmp_path / "ok",
                                  condition_attributes=["arm"], library_strategies=["RNA-Seq"],
                                  design_table=table, design_join="isolate")
    assert report["group_sizes"] == {"N": 3, "P": 3}
    with pytest.raises(es.ENAError, match="no row"):
        es.build_samplesheet("PRJNA1", study_id="X", out_dir=tmp_path / "bad",
                             condition_attributes=["arm"], design_table=table, design_join="isolate")
