"""A second comparison of the same study reuses the per-library RNA-seq tasks.

Two `-stub-run`s of the raw DAG share one launch directory (so one cache
session) and one work directory, and differ only in comparison-level params:
--comparison_id, --control_level/--treatment_level (swapped),
--exclude_samples and --outdir. With `-resume`, the second run must take
every per-library and index task from the cache and rerun only the
comparison-specific tail.

On 2026-10-12 three DUPERRET2025_STARVE_RNASEQ comparisons each recomputed
every library. `-dump-hashes` showed identical task hashes across
comparisons. The cause was the launch: `java -jar nextflow-<ver>-one.jar`
skips the launcher's --add-opens flags, so Nextflow could not write any cache
entry. main.nf now warns about that, and errors under -resume. This test
guards both halves: a param leaking into a per-library task hash, and a
Nextflow launch that cannot cache.

Needs Nextflow, so it skips where `nextflow` is not on PATH (the Python CI
job) and runs in CI's Nextflow job. Set AREE_NEXTFLOW to use another launch
command, e.g. a wrapper script.
"""
from __future__ import annotations

import csv
import os
import re
import shlex
import shutil
import subprocess
from collections import Counter

import pytest

from common import REPO_ROOT

NEXTFLOW = shlex.split(os.environ.get("AREE_NEXTFLOW", "nextflow"))
FIXTURES = REPO_ROOT / "tests" / "fixtures"
VIRAL_REFERENCE = REPO_ROOT / "data" / "reference" / "pathogens" / "OsHV-1" / "NC_005881.2.fasta"

def _nextflow_version():
    """(major, minor) of a runnable Nextflow, or None."""
    if shutil.which(NEXTFLOW[0]) is None:
        return None
    try:
        proc = subprocess.run([*NEXTFLOW, "-version"], capture_output=True, text=True, timeout=300)
    except (OSError, subprocess.TimeoutExpired):
        return None
    match = re.search(r"version (\d+)\.(\d+)", proc.stdout)
    return (int(match[1]), int(match[2])) if proc.returncode == 0 and match else None


# workflows/rnaseq/nextflow.config: manifest.nextflowVersion = '>=23.10.0'
pytestmark = pytest.mark.skipif(
    (_nextflow_version() or (0, 0)) < (23, 10),
    reason=f"no runnable Nextflow >=23.10 from {' '.join(NEXTFLOW)!r} (set AREE_NEXTFLOW)",
)

# Comparison-independent: one task per library, or one per study reference.
PER_LIBRARY = ["FASTQC", "TRIM_FASTP", "SALMON_QUANT", "VIRAL_QUANT"]
PER_STUDY = ["SALMON_INDEX", "VIRAL_INDEX", "SAMPLE_QC"]
# Comparison-specific: must rerun, or the second run did not change anything.
PER_COMPARISON = [
    "DIFFERENTIAL_EXPRESSION_DESEQ2",
    "VIRAL_SUMMARY",
    "STANDARDIZE_OUTPUT",
    "EMIT_MANIFEST",
    "RENDER_REPORT",
]

LAYOUTS = {
    "paired_end": {
        "reads": str(FIXTURES / "rnaseq_raw" / "*_R{1,2}.fastq"),
        "sample_sheet": FIXTURES / "rnaseq_raw" / "samples.tsv",
        "extra": [],
    },
    "single_end": {
        "reads": str(FIXTURES / "rnaseq_raw_se" / "ERR*.fastq"),
        "sample_sheet": FIXTURES / "rnaseq_raw_se" / "samples.tsv",
        "extra": ["--single_end"],
    },
}


def _run(launch_dir, work_dir, layout, comparison_id, control, treatment, extra=()):
    outdir = launch_dir / f"out_{comparison_id}"
    cmd = [
        *NEXTFLOW, "run", str(REPO_ROOT / "workflows" / "rnaseq"),
        "-profile", "local,ci_stub", "-stub-run", "-resume", "-ansi-log", "false",
        "-work-dir", str(work_dir),
        "--mode", "raw_reanalysis",
        "--study_id", "CI_RESUME",
        "--comparison_id", comparison_id,
        "--reads", layout["reads"],
        "--transcript_fasta", str(FIXTURES / "rnaseq_raw" / "reference.fa"),
        "--gtf", str(FIXTURES / "rnaseq_raw" / "annotation.gtf"),
        "--tx2gene", str(FIXTURES / "rnaseq_raw" / "tx2gene.tsv"),
        "--sample_sheet", str(layout["sample_sheet"]),
        "--viral_reference", str(VIRAL_REFERENCE),
        "--control_level", control,
        "--treatment_level", treatment,
        "--outdir", str(outdir),
        "--tracedir", str(outdir / "pipeline_info"),
        *layout["extra"], *extra,
    ]
    proc = subprocess.run(cmd, cwd=launch_dir, capture_output=True, text=True)
    assert proc.returncode == 0, f"nextflow failed:\n{proc.stdout}\n{proc.stderr}"
    with open(outdir / "pipeline_info" / "execution_trace.txt", newline="") as fh:
        return [(row["name"].split(" (")[0], row["status"])
                for row in csv.DictReader(fh, delimiter="\t")]


@pytest.mark.parametrize("layout_name", sorted(LAYOUTS))
def test_second_comparison_reuses_per_library_tasks(tmp_path, layout_name):
    layout = LAYOUTS[layout_name]
    work_dir = tmp_path / "work"

    first = _run(tmp_path, work_dir, layout, "cmp_a", "control", "treatment")
    assert all(status == "COMPLETED" for _, status in first), first

    second = _run(tmp_path, work_dir, layout, "cmp_b", "treatment", "control",
                  extra=["--exclude_samples", "treatment2"])
    cached = Counter(name for name, status in second if status == "CACHED")
    rerun = Counter(name for name, status in second if status != "CACHED")

    assert cached, (
        "The second run reused no task at all. If every hash matches (-dump-hashes), "
        "Nextflow is not writing its cache: launch it with the `nextflow` script, "
        "not `java -jar`."
    )
    for name in PER_LIBRARY:
        assert cached[name] == 4 and not rerun[name], (name, cached[name], rerun[name])
    for name in PER_STUDY:
        assert cached[name] == 1 and not rerun[name], (name, cached[name], rerun[name])
    for name in PER_COMPARISON:
        assert rerun[name] == 1 and not cached[name], (name, cached[name], rerun[name])
