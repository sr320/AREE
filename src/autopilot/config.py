"""Autopilot settings: config/autopilot.yaml over built-in defaults, with $USER expanded."""
from __future__ import annotations

import os
from pathlib import Path

import yaml

from common import REPO_ROOT

DEFAULT_CONFIG_PATH = REPO_ROOT / "config" / "autopilot.yaml"

DEFAULTS = {
    # Where the ledger (state.json) and the local status copy live. Not on
    # scrubbed: the ledger must outlive a 21-day purge.
    "state_dir": "~/.local/share/aree-autopilot",
    # Head jobs are found in squeue by this job-name prefix: aree-<STUDY_ID>.
    "head_job_prefix": "aree-",
    # Nextflow names its Slurm task jobs nf-<PROCESS>.
    "task_job_prefix": "nf-",
    # Default --outdir of a run; {study_id} is filled in. pipeline_info/ is under it.
    "results_template": "/gscratch/scrubbed/$USER/aree/results/{study_id}",
    "scrubbed_root": "/gscratch/scrubbed",
    "min_scrubbed_free_tb": 8.0,
    # A running head job whose trace file has not changed for this long is stalled.
    "stall_hours": 6.0,
    # Scrubbed deletes files untouched for 21 days; warn well before that.
    "scratch_alert_days": 14,
    "github_repo": "sr320/AREE",
    "status_issue_title": "Klone autopilot status",
    "status_label": "autopilot-status",
    "alert_label": "autopilot-alert",
    "token_file": "~/.config/aree/github_token",
    # "report": observe and report only (Phase 1). "resume": also resume runs
    # that stopped for a recoverable reason (Phase 2). A PAUSE file in
    # state_dir, or the paused_label on the status issue, forces "report".
    "mode": "report",
    "paused_label": "autopilot-paused",
    # Resume budget per run; past it the run needs a person.
    "max_resumes_per_day": 3,
    "max_resumes_total": 10,
    # New memory for a process that ran out = its peak RSS x this, rounded up
    # to whole GB, capped at oom_max_memory_gb so Slurm never queues an
    # unsatisfiable request. One escalation per process.
    "oom_memory_factor": 2.0,
    "oom_max_memory_gb": 192,
    # Head jobs submitted by `aree autopilot launch` (scripts/klone/head.sbatch).
    "head_script": str(REPO_ROOT / "scripts" / "klone" / "head.sbatch"),
    "head_account": "srlab",
    "head_partition": "cpu-g2",
    "head_time": "4-00:00:00",
    "head_cpus": 2,
    "head_mem": "8G",
    # Exported to the head job (see docs/running_on_klone.md#running).
    "miniforge": "~/miniforge3",
    "nextflow_env": "nextflow",
    "nxf_root": "/gscratch/scrubbed/$USER/aree",
    "task_account": "srlab",
}


def _expand(value):
    if isinstance(value, str):
        return os.path.expanduser(os.path.expandvars(value))
    return value


def load_config(path: Path | str | None = None) -> dict:
    """Return the settings, reading `path` (default config/autopilot.yaml) if it exists."""
    cfg = dict(DEFAULTS)
    path = Path(path) if path else DEFAULT_CONFIG_PATH
    if path.exists():
        cfg.update(yaml.safe_load(path.read_text()) or {})
    return {key: _expand(value) for key, value in cfg.items()}
