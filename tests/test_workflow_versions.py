"""A workflow manifest's software_versions must come from what actually ran.

Until 2026-10-11 every EMIT_MANIFEST wrote a hard-coded copy of the container
tags, so native runs (-profile local) reported e.g. Salmon 1.10.3 while Salmon
2.3.4 ran. These tests pin the replacement: each process's versions.yml is
wired into EMIT_MANIFEST, merged by modules/shared/collect_versions.py, and
scripts/check_manifest_versions.py (run in CI after the Nextflow smoke runs)
fails a manifest that disagrees with the work dir.
"""
import importlib.util
import json
import re
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
ASSAYS = ["rnaseq", "methylation", "proteomics", "metabolomics"]
# Run after (RENDER_REPORT) or are (EMIT_MANIFEST) the manifest step.
NOT_UPSTREAM = {"EMIT_MANIFEST", "RENDER_REPORT"}


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


collect_versions = _load("collect_versions", REPO / "modules/shared/collect_versions.py")
checker = _load("check_manifest_versions", REPO / "scripts/check_manifest_versions.py")

# data/studies/IOCAS2024_OA_TRANSGEN/process_versions.yml (sr320/AREE#17):
# what each process reported on the 2026-10-11 native run.
IOCAS_NATIVE = {
    "TRIM_FASTP": "TRIM_FASTP:\n    fastp: 1.3.6\n",
    "FASTQC": "FASTQC:\n    fastqc: 0.12.1\n",
    "SALMON_QUANT": "SALMON_QUANT:\n    salmon: 2.3.4\n",
    "SAMPLE_QC": "SAMPLE_QC:\n    multiqc: 1.9\n",
    "DIFFERENTIAL_EXPRESSION_DESEQ2": (
        "DIFFERENTIAL_EXPRESSION_DESEQ2:\n    r-base: 4.3.2\n"
        "    bioconductor-deseq2: 1.42.1\n    bioconductor-tximport: 1.30.0\n"),
    "STANDARDIZE_OUTPUT": "STANDARDIZE_OUTPUT:\n    python: 3.9.6\n    pandas: 2.3.3\n",
}


def _work_dir(tmp_path, files):
    """Fake Nextflow work dir: one task dir per (process, versions.yml text)."""
    for i, text in enumerate(files):
        task = tmp_path / "work" / f"{i:02x}" / f"task{i:04d}"
        task.mkdir(parents=True)
        (task / "versions.yml").write_text(text)
    return tmp_path / "work"


# ---- parsing and merging ---------------------------------------------------

def test_parse_nf_core_versions_yml_with_quotes_and_comments():
    text = 'DIFFERENTIAL_ABUNDANCE:\n    r_base: "4.3.2"\n    limma: \'3.58.1\'  # pinned\n'
    versions, warnings = collect_versions.parse_versions_yml(text)
    assert versions == {"DIFFERENTIAL_ABUNDANCE": {"r_base": "4.3.2", "limma": "3.58.1"}}
    assert warnings == []


def test_stub_versions_files_are_recorded_not_dropped():
    # `echo 'PROC:' > versions.yml` stubs name the process with no tools;
    # `touch versions.yml` stubs name nothing and get a warning.
    assert collect_versions.parse_versions_yml("FASTQC:\n") == ({"FASTQC": {}}, [])
    versions, warnings = collect_versions.parse_versions_yml("")
    assert versions == {} and "no process block" in warnings[0]


def test_identical_per_sample_files_collapse_and_disagreements_are_kept():
    merged, warnings = collect_versions.merge_versions([
        {"SALMON_QUANT": {"salmon": "2.3.4"}},
        {"SALMON_QUANT": {"salmon": "2.3.4"}},
        {"FASTQC": {"fastqc": "0.12.1"}},
        {"FASTQC": {"fastqc": "0.11.9"}},
    ])
    assert merged["SALMON_QUANT"] == {"salmon": "2.3.4"}
    assert merged["FASTQC"] == {"fastqc": ["0.11.9", "0.12.1"]}
    assert any("software_version_inconsistent: FASTQC" in w for w in warnings)


def test_empty_version_string_is_flagged():
    merged, warnings = collect_versions.merge_versions([{"TRIM_FASTP": {"fastp": ""}}])
    assert merged == {"TRIM_FASTP": {"fastp": ""}}
    assert any("software_version_unresolved: TRIM_FASTP" in w for w in warnings)


def test_iocas_native_run_reports_the_tools_that_ran(tmp_path):
    files = []
    for process, text in IOCAS_NATIVE.items():
        path = tmp_path / f"{process}.yml"
        path.write_text(text)
        files.append(str(path))
    merged, warnings = collect_versions.collect(files)
    assert warnings == []
    assert merged["SALMON_QUANT"] == {"salmon": "2.3.4"}
    assert merged["TRIM_FASTP"] == {"fastp": "1.3.6"}
    assert merged["SAMPLE_QC"] == {"multiqc": "1.9"}


def test_cli_writes_software_versions_json(tmp_path):
    src = tmp_path / "v.yml"
    src.write_text(IOCAS_NATIVE["SALMON_QUANT"])
    out = tmp_path / "out.json"
    assert collect_versions.main(["collect_versions.py", str(out), str(src)]) == 0
    assert json.loads(out.read_text()) == {
        "software_versions": {"SALMON_QUANT": {"salmon": "2.3.4"}}, "warnings": []}


# ---- manifest vs. work dir --------------------------------------------------

def test_checker_accepts_a_manifest_built_from_the_run(tmp_path):
    work = _work_dir(tmp_path, list(IOCAS_NATIVE.values()) + [
        IOCAS_NATIVE["SALMON_QUANT"],             # a second sample's task
        "RENDER_REPORT:\n    quarto: 1.6.40\n",   # runs after the manifest
    ])
    ran, _ = collect_versions.collect(checker.task_versions_files(work))
    manifest = {
        "software_versions": {p: v for p, v in ran.items() if p != "RENDER_REPORT"},
        "software_versions_not_captured": {"RENDER_REPORT": "runs after the manifest"},
    }
    assert checker.compare(manifest, ran) == []


def test_checker_rejects_container_pins_on_a_native_run(tmp_path):
    work = _work_dir(tmp_path, list(IOCAS_NATIVE.values()))
    ran, _ = collect_versions.collect(checker.task_versions_files(work))
    # The pre-fix manifest block: flat container tags, not what ran.
    pinned = {"software_versions": {"fastp": "0.23.4", "salmon": "1.10.3", "multiqc": "1.21"}}
    problems = checker.compare(pinned, ran)
    assert any("SALMON_QUANT ran but is missing" in p for p in problems)
    assert any("salmon is in software_versions but left no versions.yml" in p for p in problems)

    wrong = {"software_versions": dict(ran, SALMON_QUANT={"salmon": "1.10.3"})}
    assert checker.compare(wrong, ran) == [
        'SALMON_QUANT: manifest says {"salmon": "1.10.3"}, the process reported {"salmon": "2.3.4"}']


def test_checker_ignores_versions_staged_into_emit_manifest(tmp_path):
    work = _work_dir(tmp_path, [IOCAS_NATIVE["FASTQC"]])
    staged = work / "ff" / "emit" / "process_versions"
    staged.mkdir(parents=True)
    (staged / "versions_1.yml").write_text(IOCAS_NATIVE["FASTQC"])
    assert len(checker.task_versions_files(work)) == 1


# ---- workflow wiring ---------------------------------------------------------

def _processes_with_versions(module_path):
    """Process names in a module file that declare `emit: versions`."""
    text = module_path.read_text()
    blocks = re.split(r"^process\s+", text, flags=re.M)[1:]
    return {b.split()[0] for b in blocks if re.search(r"emit:\s*versions\b", b)}


def _included_processes(main_nf):
    included = {}
    for names, module in re.findall(r"include\s*\{([^}]*)\}\s*from\s*'([^']+)'", main_nf.read_text()):
        for name in re.split(r"[;\s]+", names.strip()):
            if name:
                included[name] = (main_nf.parent / module).resolve()
    return included


@pytest.mark.parametrize("assay", ASSAYS)
def test_every_process_versions_channel_reaches_emit_manifest(assay):
    main_nf = REPO / "workflows" / assay / "main.nf"
    text = main_nf.read_text()
    missing = []
    for name, module in _included_processes(main_nf).items():
        if name in NOT_UPSTREAM:
            continue
        if name in _processes_with_versions(module) and f"{name}.out.versions" not in text:
            missing.append(name)
    assert missing == [], (
        f"{main_nf}: these processes write versions.yml but it never reaches "
        f"EMIT_MANIFEST, so the manifest would omit them: {missing}")
    assert "collect_versions_py" in text


@pytest.mark.parametrize("assay", ASSAYS)
def test_emit_manifest_builds_software_versions_from_collected_files(assay):
    text = (REPO / "modules" / assay / "emit_manifest.nf").read_text()
    assert "stageAs: 'process_versions/" in text
    assert "collect_versions_py" in text
    assert '"software_versions": software_versions' in text
    # No literal {"tool": "1.2.3"} block standing in for what ran.
    assert not re.search(r'software_versions\s*=\s*\{\s*"', text)
    assert not re.search(r'"software_versions":\s*\{\s*(#[^\n]*\n\s*)*"', text)
