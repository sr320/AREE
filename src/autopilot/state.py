"""The run ledger (state.json) and the rules that turn observations into a run state.

States: pending, running, stalled, stopped, finished, unknown. `assess` is a pure
function of what Slurm, the trace file and the log say, so every rule is tested
on recorded output.
"""
from __future__ import annotations

import json
import os
import tempfile
from pathlib import Path

from autopilot.nextflow import classify_log, log_succeeded
from autopilot.slurm import elapsed_hours

# Head-job ends that say nothing about the pipeline: the job was interrupted.
INTERRUPTED_STATES = {"PREEMPTED", "NODE_FAIL", "TIMEOUT", "BOOT_FAIL", "DEADLINE"}
ACTIVE_SQUEUE = {"RUNNING", "COMPLETING", "CONFIGURING", "SUSPENDED", "STAGE_OUT", "RESIZING"}


def load_state(path: Path) -> dict:
    if not path.exists():
        return {"runs": {}}
    data = json.loads(path.read_text())
    data.setdefault("runs", {})
    return data


def save_state(path: Path, data: dict) -> None:
    """Write atomically, so a tick killed mid-write leaves the old ledger intact."""
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(dir=path.parent, prefix=".state.", suffix=".json")
    with os.fdopen(fd, "w") as fh:
        json.dump(data, fh, indent=2, sort_keys=True)
        fh.write("\n")
    os.replace(tmp, path)


def assess(head: dict | None, acct: dict | None, trace: dict, log_text: str,
           now: float, stall_hours: float) -> dict:
    """Classify one run.

    head  -- its squeue row, or None when the head job is not in the queue
    acct  -- its sacct row, or None
    trace -- autopilot.nextflow.trace_summary()
    Returns {"state", "reason", "recoverable", "next_step"}.
    """
    if head is not None:
        if head["state"] == "PENDING":
            return _result("pending", f"head job pending ({head.get('reason') or 'no reason'})")
        if head["state"] in ACTIVE_SQUEUE:
            mtime = trace.get("trace_mtime")
            elapsed = elapsed_hours(head.get("elapsed", ""))
            # A fresh head job (e.g. just resumed) is not stalled by an old trace.
            idle = min((now - mtime) / 3600, elapsed) if mtime else elapsed
            if idle >= stall_hours:
                since = "no task has finished" if mtime else "no trace file yet"
                return _result("stalled", f"head job running but {since} for {idle:.1f} h",
                               recoverable=True, next_step="scancel the head job, then -resume")
            return _result("running", f"head job {head['state'].lower()}")
        return _result("running", f"head job {head['state'].lower()}")

    if acct is None:
        if log_succeeded(log_text):
            return _result("finished", "log reports Execution complete")
        return _result("unknown", "head job not in squeue and no sacct record",
                       next_step="check the run by hand, or untrack it")

    state = acct["state"]
    if state == "PENDING" or state in ACTIVE_SQUEUE:
        # Just submitted and not yet visible in squeue.
        return _result("pending" if state == "PENDING" else "running", f"head job {state.lower()} (sacct)")
    if state == "COMPLETED" and acct.get("exit_code", "0:0") == "0:0":
        return _result("finished", "head job completed, exit 0")
    if state in INTERRUPTED_STATES:
        return _result("stopped", f"head job {state} ({acct.get('end') or 'end unknown'})",
                       recoverable=True, next_step="-resume", cls="head_interrupted")
    if state == "CANCELLED" and not log_text.strip():
        return _result("stopped", f"head job {acct['state_detail']}",
                       next_step="resume or untrack by hand", cls="cancelled")
    cls, recoverable, step = classify_log(log_text)
    return _result("stopped", f"head job {acct['state_detail']}, exit {acct.get('exit_code')}",
                   recoverable=recoverable, next_step=step, cls=cls)


def _result(state, reason, recoverable=False, next_step="", cls=""):
    return {"state": state, "reason": reason, "recoverable": recoverable,
            "next_step": next_step, "failure_class": cls}
