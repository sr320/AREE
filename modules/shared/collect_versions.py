#!/usr/bin/env python3
"""Merge the versions.yml files written by a workflow's processes.

Each AREE Nextflow process writes a versions.yml of the nf-core shape:

    PROCESS_NAME:
        tool: 1.2.3
        other-tool: "4.5"

EMIT_MANIFEST stages every one it receives and calls this script, so the
manifest's ``software_versions`` reports the tools that actually ran (native
Homebrew/conda tools or a container, whichever the profile used) instead of
the declared container pins. A 2026-10-11 native run reported Salmon 1.10.3,
fastp 0.23.4 and MultiQC 1.21 from those pins while Salmon 2.3.4, fastp 1.3.6
and MultiQC 1.9 ran.

Standard library only and Python 3.6+, because it runs inside
python:3.11-slim and under whatever python3 a native run finds. The parser
covers the flat two-level shape above only; it is not a general YAML reader.

Usage:
    collect_versions.py OUT.json [versions.yml ...]

OUT.json holds {"software_versions": {...}, "warnings": [...]}.
"""
import json
import sys


def parse_versions_yml(text, source="versions.yml"):
    """Return ({process: {tool: version}}, warnings) for one versions.yml."""
    versions = {}
    warnings = []
    current = None
    for lineno, raw in enumerate(text.splitlines(), 1):
        line = raw.split(" #", 1)[0].rstrip()
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        if not line[0].isspace():
            if not line.endswith(":"):
                warnings.append(f"{source}:{lineno}: expected 'PROCESS:' line, got {raw!r}")
                current = None
                continue
            current = line[:-1].strip().strip("'\"")
            versions.setdefault(current, {})
            continue
        if current is None:
            warnings.append(f"{source}:{lineno}: tool line outside a process block: {raw!r}")
            continue
        if ":" not in line:
            warnings.append(f"{source}:{lineno}: expected 'tool: version', got {raw!r}")
            continue
        tool, value = line.strip().split(":", 1)
        versions[current][tool.strip().strip("'\"")] = value.strip().strip("'\"")
    if not versions:
        warnings.append(f"{source}: no process block (empty versions.yml, e.g. a stub task)")
    return versions, warnings


def merge_versions(parsed):
    """Merge per-file results into one {process: {tool: version}} mapping.

    A per-sample process (SALMON_QUANT, FASTQC, ...) writes one identical
    file per task, so identical values collapse. Tasks of one process that
    disagree keep every distinct value, sorted, and get a warning: that is a
    mixed environment, not something to paper over by picking one.
    """
    seen = {}
    for versions in parsed:
        for process, tools in versions.items():
            slot = seen.setdefault(process, {})
            for tool, value in tools.items():
                slot.setdefault(tool, set()).add(value)
    merged = {}
    warnings = []
    for process in sorted(seen):
        merged[process] = {}
        for tool in sorted(seen[process]):
            values = sorted(seen[process][tool])
            if len(values) == 1:
                merged[process][tool] = values[0]
            else:
                merged[process][tool] = values
                warnings.append(
                    f"software_version_inconsistent: {process} reported {tool} as "
                    f"{', '.join(values)} across its tasks")
            if values == [""]:
                warnings.append(
                    f"software_version_unresolved: {process} wrote an empty version for {tool} "
                    "(version command failed or tool missing)")
    return merged, warnings


def collect(paths):
    parsed = []
    warnings = []
    for path in sorted(paths):
        with open(path) as fh:
            versions, file_warnings = parse_versions_yml(fh.read(), source=path)
        parsed.append(versions)
        warnings.extend(file_warnings)
    merged, merge_warnings = merge_versions(parsed)
    return merged, warnings + merge_warnings


def main(argv):
    if len(argv) < 2:
        sys.stderr.write(__doc__)
        return 2
    merged, warnings = collect(argv[2:])
    with open(argv[1], "w") as fh:
        json.dump({"software_versions": merged, "warnings": warnings}, fh, indent=2, sort_keys=True)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
