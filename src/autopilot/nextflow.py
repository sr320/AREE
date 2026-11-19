"""Read a Nextflow run's progress (execution_trace.txt) and why it stopped (.nextflow.log).

The trace file gets one row per finished task, so its modification time is the
time of the run's last finished task: the stall signal.
"""
from __future__ import annotations

import csv
import re
from collections import Counter
from pathlib import Path

TRACE_NAME = "execution_trace.txt"
LOG_NAME = ".nextflow.log"

# Order matters: the first match wins. Each entry: class, pattern, recoverable,
# and what Phase 2 would do (or what a person needs to do).
LOG_CLASSES = [
    ("task_walltime", re.compile(r"exceeded running time limit|exit status \(140\)", re.I), False,
     "lower methylation.align_chunk_pairs or raise the process time in the run-only config"),
    ("task_oom", re.compile(r"exit status \(137\)|OUT_OF_MEMORY|OutOfMemory", re.I), True,
     "raise that process's memory one step in a run-only override, then -resume"),
    ("container", re.compile(r"Failed to pull singularity image|apptainer.*(pull|FATAL)|"
                             r"singularity.*FATAL|Cannot find (an )?image", re.I), False,
     "fix the container or image cache, then -resume"),
    ("input_missing", re.compile(r"Missing output file|No such file or directory|"
                                 r"Cannot find any reads matching|does not exist|checksum|md5", re.I), False,
     "restore or re-stage the missing input, then -resume"),
    ("slurm_controller", re.compile(r"Socket timed out|Unable to contact slurm controller", re.I), True,
     "-resume"),
]
SUCCESS_RE = re.compile(r"Execution complete -- Goodbye")
_ERROR_MARK = "ERROR ~"


def trace_summary(outdir: Path | str) -> dict:
    """Task counts by status, the failed task names and the trace file's mtime (or None)."""
    path = Path(outdir) / "pipeline_info" / TRACE_NAME
    summary = {"trace_path": str(path), "trace_mtime": None, "status_counts": {}, "failed_tasks": []}
    if not path.exists():
        return summary
    summary["trace_mtime"] = path.stat().st_mtime
    counts = Counter()
    with path.open(newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            status = (row.get("status") or "").strip()
            counts[status] += 1
            if status in ("FAILED", "ABORTED"):
                summary["failed_tasks"].append(row.get("name", ""))
    summary["status_counts"] = dict(counts)
    return summary


def read_log_tail(launch_dir: Path | str, max_bytes: int = 200_000) -> str:
    path = Path(launch_dir) / LOG_NAME
    if not path.exists():
        return ""
    with path.open("rb") as fh:
        fh.seek(0, 2)
        size = fh.tell()
        fh.seek(max(0, size - max_bytes))
        return fh.read().decode("utf-8", errors="replace")


def error_section(log_text: str, lines: int = 80) -> str:
    """The text from the last `ERROR ~` on, or the last `lines` lines when there is none.

    Slurm socket timeouts are logged and retried throughout a healthy run, so
    only the final error says why a run stopped.
    """
    idx = log_text.rfind(_ERROR_MARK)
    if idx >= 0:
        return log_text[idx:]
    return "\n".join(log_text.splitlines()[-lines:])


def log_succeeded(log_text: str) -> bool:
    return bool(SUCCESS_RE.search(log_text))


def classify_log(log_text: str) -> tuple[str, bool, str]:
    """(class, recoverable, next step) for a run that stopped."""
    section = error_section(log_text)
    for name, pattern, recoverable, action in LOG_CLASSES:
        if pattern.search(section):
            return name, recoverable, action
    return "unknown", False, "read the end of .nextflow.log"
