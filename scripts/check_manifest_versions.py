#!/usr/bin/env python3
"""Check a workflow manifest's software_versions against what actually ran.

Every task directory in a Nextflow work dir holds the versions.yml that task
wrote. This script merges them with the same code EMIT_MANIFEST uses
(modules/shared/collect_versions.py) and fails when the manifest disagrees:

  * a process that ran is missing from software_versions (its versions
    channel was not wired into EMIT_MANIFEST), unless the manifest lists it
    under software_versions_not_captured;
  * a version in the manifest differs from the one the process reported
    (e.g. hard-coded container pins on a native run);
  * the manifest names a process that left no versions.yml.

Use a work dir that holds only the run being checked: tasks cached from an
earlier run with other inputs would be counted too.

Usage:
    python scripts/check_manifest_versions.py MANIFEST.json WORK_DIR
"""
import importlib.util
import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]


def _load_collector():
    spec = importlib.util.spec_from_file_location(
        "collect_versions", REPO / "modules" / "shared" / "collect_versions.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


collect_versions = _load_collector()


def task_versions_files(work_dir):
    """versions.yml at the top of each task dir (work/ab/cdef.../versions.yml).

    Files EMIT_MANIFEST staged under process_versions/ are one level deeper,
    so they are not double counted.
    """
    return sorted(str(p) for p in Path(work_dir).glob("*/*/versions.yml"))


def compare(manifest, ran):
    """Return a list of problems; empty means the manifest matches the run."""
    reported = manifest.get("software_versions")
    if not isinstance(reported, dict):
        return ["manifest has no software_versions mapping"]
    not_captured = manifest.get("software_versions_not_captured") or {}
    problems = []
    for process in sorted(ran):
        if process in reported:
            if reported[process] != ran[process]:
                said = json.dumps(reported[process], sort_keys=True)
                saw = json.dumps(ran[process], sort_keys=True)
                problems.append(f"{process}: manifest says {said}, the process reported {saw}")
        elif process not in not_captured:
            problems.append(
                f"{process} ran but is missing from software_versions (wire its "
                ".out.versions into EMIT_MANIFEST)")
    for process in sorted(set(reported) - set(ran)):
        problems.append(
            f"{process} is in software_versions but left no versions.yml in the work dir")
    return problems


def main(argv):
    if len(argv) != 3:
        sys.stderr.write(__doc__)
        return 2
    manifest = json.loads(Path(argv[1]).read_text())
    files = task_versions_files(argv[2])
    if not files:
        sys.stderr.write(f"no */*/versions.yml under {argv[2]}\n")
        return 1
    ran, _ = collect_versions.collect(files)
    problems = compare(manifest, ran)
    for problem in problems:
        print("FAIL " + problem)
    if problems:
        return 1
    print(f"OK {argv[1]}: software_versions matches {len(ran)} processes in {argv[2]}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
