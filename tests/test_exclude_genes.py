"""Tests for --exclude_genes support: the gene-list builder and its wiring.

The R side (dropping genes before DESeq2) is exercised end to end on real data
in the PR that added it; here the list format, the committed rRNA list, and
the workflow plumbing are pinned.
"""
from __future__ import annotations

import importlib.util
import re

from common import REPO_ROOT

RRNA_LIST = REPO_ROOT / "data/reference/gene_lists/GCF_963853765.1-RS_2024_06_rRNA.txt"


def _builder():
    spec = importlib.util.spec_from_file_location("lgb", REPO_ROOT / "scripts/list_genes_by_biotype.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _ids(path):
    """Parse a list exactly as the DESeq2 step does: strip '#' comments, blanks."""
    ids = [re.sub(r"#.*$", "", line).strip() for line in path.read_text().splitlines()]
    return [i for i in ids if i]


GTF = "\n".join([
    "#!genome-build test",
    'chr1\tRefSeq\tgene\t1\t10\t.\t+\t.\tgene_id "LOC1"; db_xref "GeneID:111"; gene "LOC1"; gene_biotype "rRNA";',
    'chr1\tRefSeq\ttranscript\t1\t10\t.\t+\t.\tgene_id "LOC1"; transcript_id "XR_1.1"; db_xref "GeneID:111"; gene_biotype "rRNA";',
    'chr1\tRefSeq\tgene\t20\t90\t.\t+\t.\tgene_id "LOC2"; db_xref "GeneID:222"; gene "RPL3"; gene_biotype "protein_coding";',
    'chr1\tRefSeq\tgene\t95\t99\t.\t-\t.\tgene_id "LOC3"; db_xref "GeneID:33"; gene "LOC3"; gene_biotype "rRNA";',
    'chr1\tRefSeq\tgene\t100\t110\t.\t-\t.\tgene_id "LOC4"; db_xref "GeneID:444"; gene "trnL"; gene_biotype "tRNA";',
]) + "\n"


def test_builder_lists_only_the_requested_biotypes(tmp_path):
    gtf = tmp_path / "genomic.gtf"
    gtf.write_text(GTF)
    tx2gene = tmp_path / "tx2gene.tsv"
    tx2gene.write_text("XR_1.1\t111\nXM_2.1\t222\n")
    out = tmp_path / "rrna.txt"
    lgb = _builder()
    assert lgb.genes_by_biotype(gtf, {"rRNA"}) == {"111": "LOC1", "33": "LOC3"}
    import sys
    argv = sys.argv
    sys.argv = ["x", "--gtf", str(gtf), "--biotype", "rRNA", "--tx2gene", str(tx2gene), "--out", str(out)]
    try:
        assert lgb.main() == 0
    finally:
        sys.argv = argv
    assert _ids(out) == ["33", "111"]          # numeric order, comments stripped
    text = out.read_text()
    assert "1 of 2 are in tx2gene.tsv" in text   # LOC3 is not in this map
    assert "sha256" in text


def test_committed_rrna_list_matches_the_reference_annotation():
    ids = _ids(RRNA_LIST)
    assert len(ids) == len(set(ids)) == 68
    assert all(i.isdigit() for i in ids)
    head = RRNA_LIST.read_text()
    assert 'gene_biotype rRNA' in head
    assert "68 of 68 are in tx2gene.tsv" in head


def test_workflow_wires_the_option_and_its_placeholder():
    wf = REPO_ROOT / "workflows/rnaseq"
    assert (wf / "assets/NO_EXCLUDE_GENES").exists()
    assert re.search(r"^\s*exclude_genes\s*=\s*null", (wf / "nextflow.config").read_text(), re.M)
    main = (wf / "main.nf").read_text()
    assert "params.exclude_genes ? Channel.fromPath(params.exclude_genes" in main
    module = (REPO_ROOT / "modules/rnaseq/differential_expression_deseq2.nf").read_text()
    assert "path exclude_genes" in module
    # Applied after tximport in both the per-library and the pseudo-bulk branch.
    assert module.count("txi <- drop_genes(txi)") == 2
    assert '"exclude_genes"' in (REPO_ROOT / "modules/rnaseq/emit_manifest.nf").read_text()
