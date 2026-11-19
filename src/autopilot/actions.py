"""Phase 2: resume runs that stopped or stalled for a recoverable reason.

Only runs started by `aree autopilot launch` can be resumed, because only those
have a recorded head-job submission (`rec["launch"]`). Their head job always
runs `nextflow run ... -resume`, so resubmitting the same command continues the
run from its cache. A run tracked any other way stays report-only.

`plan_action` is pure (what would be done, and why); `apply_action` does it.
"""
from __future__ import annotations

import math
from pathlib import Path

from autopilot.nextflow import memory_gb, process_name
from autopilot.slurm import sbatch, scancel

DAY = 86400


def override_path(cfg: dict, study_id: str) -> Path:
    """Run-only Nextflow config the autopilot writes (memory escalations); passed last with -c."""
    return Path(cfg["state_dir"]) / "overrides" / f"{study_id}.config"


def launch_argv(cfg: dict, study_id: str, workflow: str, config: str, outdir: str,
                launch_dir: str, run_configs: list) -> list:
    """The sbatch command for a head job. Only the listed variables are exported
    (not ALL): a TMPDIR inherited from the caller would break the run
    (docs/running_on_klone.md#chunked-alignment)."""
    export = ",".join([
        f"AREE_MINIFORGE={cfg['miniforge']}",
        f"AREE_NF_ENV={cfg['nextflow_env']}",
        f"AREE_NXF_ROOT={cfg['nxf_root']}",
        f"AREE_SLURM_ACCOUNT={cfg['task_account']}",
    ])
    return ["sbatch", "--parsable",
            f"--job-name={cfg['head_job_prefix']}{study_id}",
            f"--account={cfg['head_account']}", f"--partition={cfg['head_partition']}",
            f"--time={cfg['head_time']}", f"--cpus-per-task={cfg['head_cpus']}", f"--mem={cfg['head_mem']}",
            f"--chdir={launch_dir}", f"--output={launch_dir}/aree-head-%j.log",
            f"--export={export}",
            cfg["head_script"], workflow, config, outdir, str(override_path(cfg, study_id)), *run_configs]


def resumes_left(rec: dict, cfg: dict, now: float) -> tuple[bool, str]:
    history = rec.get("resumes", [])
    if len(history) >= int(cfg["max_resumes_total"]):
        return False, f"resume budget spent ({len(history)} in total)"
    recent = [r for r in history if now - r["at"] < DAY]
    if len(recent) >= int(cfg["max_resumes_per_day"]):
        return False, f"resume budget spent ({len(recent)} in the last 24 h)"
    return True, ""


def _oom_override(rec: dict, cfg: dict) -> tuple[dict | None, str]:
    """New memory settings for the processes that ran out, or why there are none."""
    escalated = rec.get("memory_escalations", {})
    new = {}
    for task in rec.get("failed_details", []):
        if task.get("exit") not in ("137", "-1", ""):
            continue
        proc = process_name(task["name"])
        if proc in escalated or proc in new:
            continue
        peak = memory_gb(task.get("peak_rss", ""))
        if not peak:
            continue
        new[proc] = int(math.ceil(peak * float(cfg["oom_memory_factor"])))
    if new:
        return new, ""
    if escalated:
        return None, "memory already raised once for " + ", ".join(sorted(escalated))
    return None, "no failed task with a recorded peak memory to scale from"


def plan_action(rec: dict, cfg: dict, now: float) -> dict | None:
    """{"kind": "resume"|"cancel_and_resume"|"blocked", "why": ..., ...} or None (nothing to do)."""
    state = rec.get("state")
    if state not in ("stopped", "stalled") or not rec.get("recoverable"):
        return None
    if not rec.get("launch"):
        return {"kind": "blocked", "why": "not started by `aree autopilot launch`; resume by hand"}
    ok, why = resumes_left(rec, cfg, now)
    if not ok:
        return {"kind": "blocked", "why": why}
    if state == "stalled":
        return {"kind": "cancel_and_resume", "why": rec.get("reason", "stalled")}
    if rec.get("failure_class") == "task_oom":
        memory, why = _oom_override(rec, cfg)
        if memory is None:
            return {"kind": "blocked", "why": why}
        return {"kind": "resume", "why": "out of memory: " + ", ".join(f"{p} -> {gb} GB"
                                                                      for p, gb in sorted(memory.items())),
                "memory": memory}
    return {"kind": "resume", "why": rec.get("reason", "")}


def _write_override(cfg: dict, rec: dict, memory: dict) -> None:
    escalated = rec.setdefault("memory_escalations", {})
    escalated.update(memory)
    lines = ["// Written by aree autopilot after out-of-memory failures. Passed last with -c.",
             "process {"]
    lines += [f"    withName: '{proc}' {{ memory = '{gb} GB' }}" for proc, gb in sorted(escalated.items())]
    lines.append("}")
    path = override_path(cfg, rec["study_id"])
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")


def apply_action(action: dict, rec: dict, cfg: dict, now: float, run) -> str:
    """Carry out a planned action and update the record. Returns a one-line log entry."""
    sid = rec["study_id"]
    if action["kind"] == "blocked":
        rec["blocked"] = action["why"]
        return f"{sid}: needs a person ({action['why']})"
    if action["kind"] == "cancel_and_resume" and rec.get("job_id"):
        scancel(rec["job_id"], run=run)
    if action.get("memory"):
        _write_override(cfg, rec, action["memory"])
    old = rec.get("job_id")
    new = sbatch(rec["launch"]["argv"], run=run)
    rec.setdefault("resumes", []).append({"at": now, "from_job": old, "to_job": new, "why": action["why"]})
    rec.update(job_id=new, state="pending", state_since=now, reason=f"resumed as job {new}",
               recoverable=False, next_step="", failure_class="")
    rec.pop("blocked", None)
    return f"{sid}: resumed job {old} as {new} ({action['why']})"
