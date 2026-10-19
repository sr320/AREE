"""Blocked (paired) limma fit in modules/methylation/dmr_pooled_tiles.nf.

Runs the module's embedded R script directly on synthetic cytosine reports:
5 matched blocks (populations) x 2 arms, a large block effect shared by both
arms of each block, and a treatment effect on the first 20 of 400 tiles. The
unpaired fit cannot see the effect through the block variation; the fit with
block_by = 'population' should, while staying calibrated on null tiles.

Skipped when Rscript, limma or data.table is not installed.
"""

from __future__ import annotations

import gzip
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

MODULE = Path(__file__).resolve().parents[1] / "modules" / "methylation" / "dmr_pooled_tiles.nf"
N_BLOCKS, N_TILES, N_CPG, COV, N_TRUE = 5, 400, 12, 15, 20


def _r_ready() -> bool:
    if shutil.which("Rscript") is None:
        return False
    probe = subprocess.run(
        ["Rscript", "-e", "library(limma); library(data.table)"], capture_output=True
    )
    return probe.returncode == 0


pytestmark = pytest.mark.skipif(not _r_ready(), reason="Rscript with limma and data.table not available")


def _render_r(values: dict[str, str]) -> str:
    """Extract the module's R script and fill its Groovy placeholders."""
    text = MODULE.read_text()
    script = text.split("cat <<-'EOF_R' > run_pooled_tiles.R\n", 1)[1].split("    EOF_R", 1)[0]
    for key, val in values.items():
        script = script.replace("${" + key + "}", val)
    assert "${" not in script
    # Undo Groovy GString escaping.
    return script.replace("\\$", "$").replace("\\\\", "\\")


def _simulate(workdir: Path, seed: int = 1) -> dict[str, str]:
    rng = np.random.default_rng(seed)
    blocks = [f"GT{i}" for i in range(1, N_BLOCKS + 1)]
    sample_ids = [f"S_{b}" for b in blocks] + [f"R_{b}" for b in blocks]
    pops = blocks + blocks
    trts = [0] * N_BLOCKS + [1] * N_BLOCKS
    mu = rng.normal(0, 1.5, N_TILES)
    block_eff = rng.normal(0, 0.8, (N_TILES, N_BLOCKS))
    eff = np.r_[np.tile([1.0, -1.0], N_TRUE // 2), np.zeros(N_TILES - N_TRUE)]
    offsets = np.arange(N_CPG) * 50 + 10
    pos = (np.repeat(np.arange(N_TILES) * 1000, N_CPG) + np.tile(offsets, N_TILES)).astype(int)
    for sid, pop, trt in zip(sample_ids, pops, trts):
        logit = mu + block_eff[:, blocks.index(pop)] + trt * eff + rng.normal(0, 0.3, N_TILES)
        p = 1 / (1 + np.exp(-np.repeat(logit, N_CPG)))
        m = rng.binomial(COV, p)
        with gzip.open(workdir / f"{sid}.CX_report.txt.gz", "wt") as fh:
            for ps, mm in zip(pos, m):
                fh.write(f"chr1\t{ps}\t+\t{mm}\t{COV - mm}\tCG\tCGA\n")
    pd.DataFrame({"sample_id": sample_ids, "population": pops}).to_csv(
        workdir / "design.tsv", sep="\t", index=False
    )
    return {
        "sample_ids.join(',')": ",".join(sample_ids),
        "treatment_labels.join(',')": ",".join(map(str, trts)),
        "cx_reports.join(',')": ",".join(f"{s}.CX_report.txt.gz" for s in sample_ids),
    }


def _run(workdir: Path, inputs: dict[str, str], block_by: str) -> subprocess.CompletedProcess:
    values = {
        **inputs,
        "dmr_test": "limma",
        "unit_min_coverage": "5",
        "qvalue_cutoff": "0.05",
        "meth_diff_cutoff": "0",
        "tile_size": "1000",
        "tile_min_cpgs": "10",
        "feature_units": "tiles",
        "annotation_gtf": "NA",
        "promoter_upstream_bp": "1000",
        "replicate_unit": "",
        "block_by": block_by,
        "design_sheet": "design.tsv",
        "study_id": "SIM",
        "comparison_id": block_by or "unblocked",
    }
    script = workdir / f"run_{block_by or 'unblocked'}.R"
    script.write_text(_render_r(values))
    return subprocess.run(["Rscript", script.name], cwd=workdir, capture_output=True, text=True)


def test_blocked_fit_recovers_paired_effect_and_stays_calibrated(tmp_path):
    inputs = _simulate(tmp_path)
    for block in ("", "population"):
        res = _run(tmp_path, inputs, block)
        assert res.returncode == 0, res.stderr

    truth = np.arange(N_TILES) < N_TRUE
    unblocked = pd.read_csv(tmp_path / "SIM_unblocked_dmr_raw.tsv", sep="\t")
    blocked = pd.read_csv(tmp_path / "SIM_population_dmr_raw.tsv", sep="\t")

    assert (blocked["analysis_method"] == "limma_blocked_pooled_tiles").all()
    assert (unblocked["analysis_method"] == "limma_pooled_tiles").all()

    # The block effect swamps the treatment effect in the unpaired fit.
    assert (blocked["qvalue"][truth] < 0.05).sum() >= 10
    assert (blocked["qvalue"][truth] < 0.05).sum() > (unblocked["qvalue"][truth] < 0.05).sum()
    # Calibrated on null tiles.
    assert (blocked["pvalue"][~truth] < 0.05).mean() < 0.10
    assert (blocked["qvalue"][~truth] < 0.05).sum() <= 5

    # Balanced complete blocks: mean within-block difference equals the
    # difference of arm means, which the unblocked fit reports.
    np.testing.assert_allclose(blocked["meth_diff_percent"], unblocked["meth_diff_percent"], atol=0.002)

    qc = pd.read_csv(tmp_path / "SIM_population_methylkit_qc.tsv", sep="\t")
    assert qc.loc[0, "block_by"] == "population"
    assert qc.loc[0, "n_blocks"] == N_BLOCKS
    assert "~ population + treatment" in qc.loc[0, "test"]


def test_block_missing_an_arm_is_an_error(tmp_path):
    inputs = _simulate(tmp_path)
    design = pd.read_csv(tmp_path / "design.tsv", sep="\t")
    design.loc[design["sample_id"] == "R_GT5", "population"] = "GT1"  # GT5 loses its treated unit
    design.to_csv(tmp_path / "design.tsv", sep="\t", index=False)
    res = _run(tmp_path, inputs, "population")
    assert res.returncode != 0
    assert "missing an arm" in res.stderr
