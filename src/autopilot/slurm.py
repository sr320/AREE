"""Read-only Slurm queries: squeue for live jobs, sacct for how a job ended.

Plain `--parsable2`/`--format` output rather than `--json`, which depends on the
Slurm build. Every query takes a `run` callable (args -> stdout) so tests can
feed recorded output.
"""
from __future__ import annotations

import os
import subprocess
from typing import Callable

Runner = Callable[[list], str]

SQUEUE_FIELDS = ["job_id", "name", "state", "elapsed", "time_limit", "partition", "work_dir", "reason"]
SACCT_FIELDS = ["job_id", "name", "state", "exit_code", "start", "end", "elapsed"]


def run_command(args: list) -> str:
    return subprocess.run(args, capture_output=True, text=True, check=True, timeout=120).stdout


def parse_squeue(text: str) -> list[dict]:
    jobs = []
    for line in text.splitlines():
        if not line.strip():
            continue
        parts = line.split("|")
        if len(parts) != len(SQUEUE_FIELDS):
            continue
        jobs.append(dict(zip(SQUEUE_FIELDS, (p.strip() for p in parts))))
    return jobs


def squeue(user: str | None = None, run: Runner = run_command) -> list[dict]:
    user = user or os.environ.get("USER", "")
    fmt = "%i|%j|%T|%M|%l|%P|%Z|%r"
    return parse_squeue(run(["squeue", "-u", user, "--noheader", f"--format={fmt}"]))


def parse_sacct(text: str) -> dict | None:
    """The allocation line of a job (`-X`), with `state` reduced to its first word
    ("CANCELLED by 1234" -> "CANCELLED")."""
    for line in text.splitlines():
        parts = line.strip().split("|")
        if len(parts) < len(SACCT_FIELDS):
            continue
        row = dict(zip(SACCT_FIELDS, (p.strip() for p in parts)))
        row["state_detail"] = row["state"]
        row["state"] = row["state"].split()[0] if row["state"] else ""
        return row
    return None


def sacct(job_id: str, run: Runner = run_command) -> dict | None:
    return parse_sacct(run(["sacct", "-j", str(job_id), "-X", "--noheader", "--parsable2",
                            "--format=JobID,JobName,State,ExitCode,Start,End,Elapsed"]))


def elapsed_hours(text: str) -> float:
    """Slurm elapsed time ("D-HH:MM:SS", "HH:MM:SS" or "MM:SS") in hours."""
    if not text or text in ("INVALID", "UNLIMITED"):
        return 0.0
    days = 0
    if "-" in text:
        day_part, text = text.split("-", 1)
        days = int(day_part)
    fields = [int(x) for x in text.split(":")]
    while len(fields) < 3:
        fields.insert(0, 0)
    hours, minutes, seconds = fields[-3:]
    return days * 24 + hours + minutes / 60 + seconds / 3600


def sbatch(argv: list, run: Runner = run_command) -> str:
    """Submit with `sbatch --parsable ...` and return the job id."""
    if "--parsable" not in argv:
        argv = [argv[0], "--parsable", *argv[1:]]
    out = run(argv).strip()
    job_id = out.splitlines()[-1].split(";")[0].strip() if out else ""
    if not job_id.isdigit():
        raise RuntimeError(f"sbatch returned no job id: {out!r}")
    return job_id


def scancel(job_id: str, run: Runner = run_command) -> None:
    run(["scancel", str(job_id)])
