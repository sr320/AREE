"""Guards on the RNA-seq workflow's 3' Tag-seq path (--tag_seq).

3' Tag-Seq (Lexogen QuantSeq and similar) needs no length correction anywhere:
not in Salmon, not in DESeq2. Applying it can manufacture fold changes from
isoform switching (workflows/rnaseq/README.md, "3' Tag-seq libraries"). These
tests read the workflow source, so they run in the Python CI job without
Nextflow, Salmon or R, and guard the wiring against being undone. CI's
Nextflow job stub-runs the --tag_seq DAG.
"""
from __future__ import annotations

import re

from common import REPO_ROOT

WORKFLOW = REPO_ROOT / "workflows" / "rnaseq"
MODULES = REPO_ROOT / "modules" / "rnaseq"


def _read(path):
    return path.read_text(encoding="utf-8")


def test_tag_seq_is_off_by_default():
    assert re.search(r"^\s*tag_seq\s*=\s*false\b", _read(WORKFLOW / "nextflow.config"), re.M)


def test_salmon_drops_length_correction_and_gc_bias_for_tag_seq():
    salmon = _read(MODULES / "align_quant_salmon.nf")
    assert "params.tag_seq ? '--noLengthCorrection' : ''" in salmon
    assert "(params.tag_seq ? '' : ' --gcBias')" in salmon


def test_fastp_trims_poly_a_for_tag_seq():
    assert "params.tag_seq ? '--trim_poly_x' : ''" in _read(MODULES / "trim_fastp.nf")


def test_deseq2_uses_counts_without_length_offsets_for_tag_seq():
    deseq2 = _read(MODULES / "differential_expression_deseq2.nf")
    # Per-library fit: counts into a plain matrix, not DESeqDataSetFromTximport.
    assert re.search(r"else if \(args\\\$tag_seq\) \{.*?DESeqDataSetFromMatrix\(countData = round\(txi\\\$counts\)",
                     deseq2, re.S)
    # Per-unit fit: raw counts summed, not length-scaled TPM.
    assert 'countsFromAbundance = if (args\\$tag_seq) "no" else "lengthScaledTPM"' in deseq2


def test_the_length_model_used_is_written_for_the_manifest():
    deseq2 = _read(MODULES / "differential_expression_deseq2.nf")
    assert '"none: counts without length offsets (tag_seq)"' in deseq2
    assert "paste0('  \"length_model\": ', json_str(length_model)" in deseq2
    assert "paste0('  \"tag_seq\": '" in deseq2


def test_manifest_records_tag_seq():
    assert '"tag_seq": ${params.tag_seq ?' in _read(MODULES / "emit_manifest.nf")
