"""Guards on the RNA-seq workflow's library protocols (full_length, three_prime_tag).

3' Tag-Seq (Lexogen QuantSeq and similar) needs no length correction anywhere:
not in Salmon, not in DESeq2. These tests read the workflow source, so they
run in the Python CI job without Nextflow, Salmon or R; they guard the wiring
against being undone. The behaviour itself was checked by running the workflow
on simulated Tag-seq reads in both modes (workflows/rnaseq/README.md), and CI
stub-runs the three_prime_tag DAG and rejects an unknown protocol.
"""
from __future__ import annotations

import re

from common import REPO_ROOT

WORKFLOW = REPO_ROOT / "workflows" / "rnaseq"
MODULES = REPO_ROOT / "modules" / "rnaseq"


def _read(path):
    return path.read_text(encoding="utf-8")


def test_full_length_is_the_default_so_existing_runs_are_unchanged():
    assert re.search(r"^\s*library_protocol\s*=\s*'full_length'", _read(WORKFLOW / "nextflow.config"), re.M)


def test_workflow_rejects_an_unknown_protocol():
    main = _read(WORKFLOW / "main.nf")
    assert "params.library_protocol in ['full_length', 'three_prime_tag']" in main
    assert "error \"--library_protocol must be" in main


def test_salmon_drops_length_correction_only_for_three_prime_tag():
    salmon = _read(MODULES / "align_quant_salmon.nf")
    assert "params.library_protocol == 'three_prime_tag' ? '--noLengthCorrection' : ''" in salmon
    assert "${length_model}" in salmon


def test_fastp_trims_poly_a_only_for_three_prime_tag():
    fastp = _read(MODULES / "trim_fastp.nf")
    assert "params.library_protocol == 'three_prime_tag' ? '--trim_poly_x' : ''" in fastp
    assert "${poly_x}" in fastp


def test_deseq2_uses_raw_counts_without_length_offsets_for_three_prime_tag():
    deseq2 = _read(MODULES / "differential_expression_deseq2.nf")
    assert 'library_protocol = "${params.library_protocol}"' in deseq2
    # Per-library fit: raw counts into a plain count matrix, not FromTximport.
    assert re.search(
        r"else if \(three_prime_tag\) \{\s*DESeqDataSetFromMatrix\(countData = round\(txi\\\$counts\)", deseq2)
    # Per-unit fit: raw counts summed, not length-scaled TPM.
    assert 'countsFromAbundance = if (three_prime_tag) "no" else "lengthScaledTPM"' in deseq2
    # The length model actually used is written for the manifest.
    assert '"length_model": ' in deseq2 and '"library_protocol": ' in deseq2


def test_manifest_records_the_protocol():
    assert '"library_protocol": "${params.library_protocol}"' in _read(MODULES / "emit_manifest.nf")
