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
