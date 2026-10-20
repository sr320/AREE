"""Tests for the replicate-QC scripts: the duplicate-sample_id fix in
check_replicate_dispersion.py and the read streaming in prescreen_replicates.py.
"""
from __future__ import annotations

import gzip
import importlib.util
import shutil

import numpy as np
import pandas as pd
import pytest

from common import REPO_ROOT


def _load(name: str):
    spec = importlib.util.spec_from_file_location(name, REPO_ROOT / "scripts" / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


crd = _load("check_replicate_dispersion")


def _sheet(sample_ids: list[str]) -> pd.DataFrame:
    runs = [f"SRR{i}" for i in range(len(sample_ids))]
    return pd.DataFrame({"run_accession": runs, "sample_id": sample_ids,
                         "condition": ["a"] * 3 + ["b"] * 3}).set_index("run_accession")


def _write_quant(root, run: str, reads: dict[str, float]):
    d = root / run
    d.mkdir(parents=True)
    pd.DataFrame({"Name": list(reads), "NumReads": list(reads.values())}).to_csv(d / "quant.sf", sep="\t", index=False)


def test_unique_sample_ids_are_kept():
    sheet = _sheet(["A1", "A2", "A3", "B1", "B2", "B3"])
    assert list(crd.library_ids(sheet)["sample_id"]) == ["A1", "A2", "A3", "B1", "B2", "B3"]


def test_shared_sample_ids_fall_back_to_run_accessions(tmp_path):
    """PRJNA341808 gives all three runs of a group one sample_alias. Keyed by
    it, the three replicates collapsed into a single column."""
    sheet = _sheet(["Ctrl"] * 3 + ["Heat"] * 3)
    fixed = crd.library_ids(sheet)
    assert list(fixed["sample_id"]) == list(sheet.index)
    assert list(sheet["sample_id"]) == ["Ctrl"] * 3 + ["Heat"] * 3   # input not mutated

    tx2gene = pd.DataFrame({"tx": ["t1", "t2"], "gene": ["g1", "g2"]})
    for i, run in enumerate(sheet.index):
        _write_quant(tmp_path, run, {"t1": 1000 + i, "t2": 2000 + i})
    counts = crd.gene_counts([tmp_path], sheet, tx2gene)
    assert counts.shape == (2, 6)


def test_identical_libraries_have_near_zero_dispersion():
    rng = np.random.default_rng(0)
    base = rng.uniform(600, 5000, size=200)
    counts = pd.DataFrame({f"s{i}": rng.poisson(base) for i in range(6)}, index=[f"g{i}" for i in range(200)])
    groups = pd.Series(["a"] * 3 + ["b"] * 3, index=counts.columns)
    table = crd.group_dispersion(counts, groups)
    assert (table["median_dispersion"].abs() < 0.005).all()


@pytest.mark.skipif(shutil.which("curl") is None, reason="curl not available")
def test_stream_head_takes_only_the_first_reads(tmp_path):
    pre = _load("prescreen_replicates")
    src = tmp_path / "full.fastq.gz"
    with gzip.open(src, "wt") as fh:
        for i in range(1000):
            fh.write(f"@r{i}\nACGT\n+\nIIII\n")
    dest = tmp_path / "head.fastq.gz"
    assert pre.stream_head(f"file://{src}", 250, dest) == 250
    with gzip.open(dest, "rt") as fh:
        lines = fh.read().splitlines()
    assert len(lines) == 1000 and lines[0] == "@r0" and lines[-4] == "@r249"
    # A file shorter than the request is taken whole.
    assert pre.stream_head(f"file://{src}", 5000, dest) == 1000


def test_ftp_urls_are_streamed_over_https():
    pre = _load("prescreen_replicates")
    assert pre.https("ftp://ftp.sra.ebi.ac.uk/vol1/x.fastq.gz") == "https://ftp.sra.ebi.ac.uk/vol1/x.fastq.gz"
    assert pre.https("file:///tmp/x.fastq.gz") == "file:///tmp/x.fastq.gz"
