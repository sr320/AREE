"""workflows/methylation/bin/sum_bismark_pe_reports.sh: per-chunk Bismark
PE reports summed into one library report in Bismark's layout."""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "workflows" / "methylation" / "bin" / "sum_bismark_pe_reports.sh"

pytestmark = pytest.mark.skipif(shutil.which("bash") is None or shutil.which("awk") is None, reason="needs bash and awk")


def _report(pairs, unique, none, multi, meth_cpg, unmeth_cpg, meth_chh, unmeth_chh):
    t = "\t"
    return "\n".join([
        "Bismark report for: X_val_1.part_001.fq.gz and X_val_2.part_001.fq.gz (version: v0.24.2)",
        "Bismark was run with Bowtie 2 against the bisulfite genome of /g/ with the specified options: -q --score-min L,0,-0.6",
        "",
        "Final Alignment report",
        "======================",
        f"Sequence pairs analysed in total:{t}{pairs}",
        f"Number of paired-end alignments with a unique best hit:{t}{unique}",
        f"Mapping efficiency:{t}{100 * unique / pairs:.1f}% ",
        f"Sequence pairs with no alignments under any condition:{t}{none}",
        f"Sequence pairs did not map uniquely:{t}{multi}",
        "",
        f"CT/GA/CT:{t}{unique // 2}{t}((converted) top strand)",
        f"CT/GA/GA:{t}{unique - unique // 2}{t}((converted) bottom strand)",
        "",
        "Final Cytosine Methylation Report",
        "=================================",
        f"Total methylated C's in CpG context:{t}{meth_cpg}",
        f"Total methylated C's in CHH context:{t}{meth_chh}",
        "",
        f"Total unmethylated C's in CpG context:{t}{unmeth_cpg}",
        f"Total unmethylated C's in CHH context:{t}{unmeth_chh}",
        "",
        f"C methylated in CpG context:{t}{100 * meth_cpg / (meth_cpg + unmeth_cpg):.1f}%",
        f"C methylated in CHH context:{t}{100 * meth_chh / (meth_chh + unmeth_chh):.1f}%",
        "",
        "",
        "Bismark completed in 0d 1h 2m 3s",
        "",
    ])


def test_counts_summed_and_percentages_recomputed(tmp_path):
    a = tmp_path / "S1.chunk0001_bismark_bt2_PE_report.txt"
    b = tmp_path / "S1.chunk0002_bismark_bt2_PE_report.txt"
    a.write_text(_report(1000, 600, 300, 100, 30, 170, 5, 995))
    b.write_text(_report(500, 200, 250, 50, 70, 30, 15, 485))

    out = subprocess.run(["bash", str(SCRIPT), "S1", str(a), str(b)], capture_output=True, text=True, check=True).stdout
    lines = dict(line.split("\t", 1) for line in out.splitlines() if "\t" in line)

    assert out.splitlines()[0] == "Bismark report for: S1 (sum of 2 chunk reports; version: v0.24.2)"
    assert lines["Sequence pairs analysed in total:"] == "1500"
    assert lines["Number of paired-end alignments with a unique best hit:"] == "800"
    # Recomputed from sums (800 / 1500), not the mean of 60.0% and 40.0%.
    assert lines["Mapping efficiency:"].strip() == "53.3%"
    assert lines["CT/GA/CT:"] == "400\t((converted) top strand)"
    assert lines["Total methylated C's in CpG context:"] == "100"
    assert lines["C methylated in CpG context:"] == "33.3%"  # 100 / 300
    assert lines["C methylated in CHH context:"] == "1.3%"   # 20 / (20 + 1480)
    assert "Bismark completed" not in out


def test_no_reports_is_an_error():
    res = subprocess.run(["bash", str(SCRIPT), "S1"], capture_output=True, text=True)
    assert res.returncode != 0
