"""One autopilot pass: observe -> classify -> act -> record -> report.

With `mode: report` (the default) nothing is submitted or cancelled. With
`mode: resume` (Phase 2), runs started by `aree autopilot launch` that stopped
or stalled for a recoverable reason are resubmitted, within a resume budget
(autopilot.actions). A PAUSE file in state_dir, or the paused label on the
status issue, forces report mode.

Runs are tracked in the ledger either by `aree autopilot track` or by finding a
head job named <head_job_prefix><STUDY_ID> in squeue. A tracked run stays in the
ledger, so a head job that has left the queue is still seen (and classified)
until a person untracks it.
"""
from __future__ import annotations

import datetime as dt
import socket
import time
from pathlib import Path

from autopilot import github as gh_mod
from autopilot.actions import apply_action, launch_argv, plan_action
from autopilot.nextflow import read_log_tail, trace_summary
from autopilot.slurm import run_command, sacct, sbatch, squeue
from autopilot.state import assess, load_state, save_state
from autopilot.storage import free_tb

REPORT_NOTE = "Report-only: the autopilot submits and cancels nothing."
RESUME_NOTE = ("Resume mode: runs started with `aree autopilot launch` that stop or stall for a "
               "recoverable reason are resubmitted with -resume, within the resume budget.")
ALERT_STATES = {"stalled", "stopped", "unknown"}
# The watchdog's own job (scripts/klone/) shares the head-job prefix.
WATCHDOG_JOB_NAME = "aree-autopilot"


def state_path(cfg: dict) -> Path:
    return Path(cfg["state_dir"]) / "state.json"


def track(cfg: dict, study_id: str, launch_dir: str, outdir: str | None = None,
          job_id: str | None = None, now: float | None = None) -> dict:
    path = state_path(cfg)
    data = load_state(path)
    run = data["runs"].setdefault(study_id, {"study_id": study_id, "first_seen": now or time.time()})
    run["launch_dir"] = str(launch_dir)
    run["outdir"] = str(outdir or cfg["results_template"].format(study_id=study_id))
    if job_id:
        run["job_id"] = str(job_id)
    save_state(path, data)
    return run


def launch(cfg: dict, study_id: str, workflow: str, config: str, launch_dir: str,
           outdir: str | None = None, run_configs: list | None = None, run=run_command,
           now: float | None = None) -> dict:
    """Submit a study's head job (scripts/klone/head.sbatch) and record how, so it can be resumed."""
    now = now or time.time()
    path = state_path(cfg)
    data = load_state(path)
    rec = data["runs"].get(study_id)
    if rec and rec.get("state") in ("pending", "running", "stalled"):
        raise RuntimeError(f"{study_id} already has head job {rec.get('job_id')} ({rec['state']})")
    outdir = str(outdir or cfg["results_template"].format(study_id=study_id))
    argv = launch_argv(cfg, study_id, str(workflow), str(config), outdir, str(launch_dir),
                       [str(c) for c in run_configs or []])
    job_id = sbatch(argv, run=run)
    rec = data["runs"].setdefault(study_id, {"study_id": study_id, "first_seen": now})
    rec.update(launch_dir=str(launch_dir), outdir=outdir, job_id=job_id, state="pending",
               state_since=now, reason=f"launched as job {job_id}", launch={"argv": argv, "at": now})
    save_state(path, data)
    return rec


def paused(cfg: dict, gh=None) -> str:
    """Why actions are off this tick ("" when they are on)."""
    if cfg.get("mode") != "resume":
        return f"mode is {cfg.get('mode')!r}"
    if (Path(cfg["state_dir"]) / "PAUSE").exists():
        return "PAUSE file present"
    if gh is not None and gh_mod.status_paused(gh, cfg):
        return f"status issue labelled {cfg['paused_label']}"
    return ""


def untrack(cfg: dict, study_id: str) -> bool:
    path = state_path(cfg)
    data = load_state(path)
    found = data["runs"].pop(study_id, None) is not None
    save_state(path, data)
    return found


def tick(cfg: dict, now: float | None = None, run=run_command, gh=None, free=free_tb,
         host: str | None = None) -> dict:
    """Run one pass. `gh` is an autopilot.github.GitHub (None: report locally only)."""
    now = now or time.time()
    host = host or socket.gethostname()
    path = state_path(cfg)
    data = load_state(path)
    runs = data["runs"]

    hold = paused(cfg, gh)
    actions_log = []

    jobs = squeue(run=run)
    prefix, task_prefix = cfg["head_job_prefix"], cfg["task_job_prefix"]
    heads = {j["name"][len(prefix):]: j for j in jobs
             if j["name"].startswith(prefix) and j["name"] != WATCHDOG_JOB_NAME}
    tasks = [j for j in jobs if j["name"].startswith(task_prefix)]

    for study_id, job in heads.items():
        rec = runs.setdefault(study_id, {"study_id": study_id, "first_seen": now})
        rec["job_id"] = job["job_id"]
        rec.setdefault("launch_dir", job["work_dir"])
        rec.setdefault("outdir", cfg["results_template"].format(study_id=study_id))

    alerts = []
    for study_id, rec in sorted(runs.items()):
        head = heads.get(study_id)
        acct = None
        if head is None and rec.get("job_id"):
            acct = sacct(rec["job_id"], run=run)
        trace = trace_summary(rec["outdir"])
        log_text = read_log_tail(rec["launch_dir"]) if rec.get("launch_dir") else ""
        result = assess(head, acct, trace, log_text, now, float(cfg["stall_hours"]))
        if result["state"] != rec.get("state"):
            rec["state_since"] = now
        rec.update(result)
        rec["trace_mtime"] = trace["trace_mtime"]
        rec["status_counts"] = trace["status_counts"]
        rec["failed_tasks"] = trace["failed_tasks"][-5:]
        rec["failed_details"] = trace["failed_details"][-20:]
        if result["state"] not in ("stopped", "stalled"):
            rec.pop("blocked", None)

        rec.pop("would", None)
        action = plan_action(rec, cfg, now)
        resumed = False
        if action is None:
            pass
        elif action["kind"] == "blocked":
            if rec.get("blocked") != action["why"]:
                actions_log.append(apply_action(action, rec, cfg, now, run))
        elif hold:
            rec["would"] = f"{action['kind'].replace('_', ' ')} ({action['why']})"
        else:
            try:
                actions_log.append(apply_action(action, rec, cfg, now, run))
                resumed = True
            except Exception as exc:  # one failed sbatch must not stop the tick for other runs
                rec["blocked"] = f"{action['kind'].replace('_', ' ')} failed: {exc}"
                actions_log.append(f"{study_id}: needs a person ({rec['blocked']})")

        if resumed:
            pass
        elif rec.get("blocked"):
            alerts.append((f"{study_id}: needs a person (job {rec.get('job_id', '?')})",
                           _run_alert_body(rec)))
        elif result["state"] in ALERT_STATES:
            alerts.append((f"{study_id}: {result['state']} (job {rec.get('job_id', '?')})",
                           _run_alert_body(rec)))
        age_days = (now - rec["first_seen"]) / 86400
        if age_days >= float(cfg["scratch_alert_days"]):
            alerts.append((f"{study_id}: on scrubbed for {int(cfg['scratch_alert_days'])}+ days",
                           f"`{study_id}` was first seen {age_days:.0f} days ago. Scrubbed deletes files "
                           f"untouched for 21 days: copy what you keep out of `{rec['outdir']}`, "
                           f"then `aree autopilot untrack {study_id}`.\n"))

    free = free(cfg["scrubbed_root"])
    floor = float(cfg["min_scrubbed_free_tb"])
    if free is not None and free < floor:
        alerts.append(("scrubbed free space below floor",
                       f"`{cfg['scrubbed_root']}` has {free:.1f} TB free (floor {floor:g} TB). "
                       f"Phase 3 will not start new studies below the floor.\n"))

    now_dt = dt.datetime.fromtimestamp(now, dt.timezone.utc)
    body = render_status(cfg, runs, tasks, free, now_dt, host, alerts, hold, actions_log)
    data["last_tick"] = now
    save_state(path, data)
    (path.parent / "status.md").write_text(body)

    issue = opened = None
    if gh is not None:
        issue = gh_mod.publish_status(gh, cfg, body)
        opened = [key for key, text in alerts if gh_mod.open_alert(gh, cfg, key, text + _footer(host))]
    return {"runs": runs, "alerts": [k for k, _ in alerts], "opened": opened, "actions": actions_log,
            "status_issue": issue, "body": body}


def _fmt_age(now: float, then: float | None) -> str:
    if not then:
        return "—"
    hours = (now - then) / 3600
    return f"{hours:.1f} h ago" if hours < 48 else f"{hours / 24:.0f} d ago"


def _tasks_done(counts: dict) -> str:
    if not counts:
        return "—"
    done = counts.get("COMPLETED", 0) + counts.get("CACHED", 0)
    failed = counts.get("FAILED", 0) + counts.get("ABORTED", 0)
    return f"{done} / {failed}"


def render_status(cfg, runs, tasks, free, now_dt, host, alerts, hold="", actions_log=()) -> str:
    now = now_dt.timestamp()
    running = sum(1 for t in tasks if t["state"] == "RUNNING")
    pending = sum(1 for t in tasks if t["state"] == "PENDING")
    free_txt = f"{free:.1f} TB" if free is not None else "n/a"
    lines = [
        gh_mod.heartbeat_marker(now_dt),
        "## Klone autopilot status",
        "",
        f"_{RESUME_NOTE if not hold else REPORT_NOTE}_" + (f" ({hold})" if hold and cfg.get("mode") == "resume" else ""),
        "",
        f"- Last tick: **{now_dt.strftime('%Y-%m-%d %H:%M UTC')}** on `{host}`",
        f"- Scrubbed free: {free_txt} (floor {float(cfg['min_scrubbed_free_tb']):g} TB)",
        f"- Nextflow task jobs in queue: {len(tasks)} ({running} running, {pending} pending)",
        f"- Open alert conditions: {len(alerts)}",
        "",
    ]
    if not runs:
        lines.append("No runs tracked. Name head jobs `"
                     f"{cfg['head_job_prefix']}<STUDY_ID>` or use `aree autopilot track`.")
    else:
        lines += ["| Study | State | Head job | Tasks done / failed | Last task finished | "
                  "Tracked for | Why | Next step |",
                  "|---|---|---|---|---|---|---|---|"]
        for sid, rec in sorted(runs.items()):
            step = rec.get("next_step") or ""
            if rec.get("blocked"):
                step = f"needs a person: {rec['blocked']}"
            elif rec.get("would"):
                step = f"would {rec['would']}"
            elif step and rec.get("recoverable"):
                step = f"recoverable: {step}"
            resumes = len(rec.get("resumes", []))
            if resumes:
                step = f"{step} (resumed {resumes}x)".strip()
            lines.append(f"| `{sid}` | **{rec.get('state', '?')}** | {rec.get('job_id', '—')} | "
                         f"{_tasks_done(rec.get('status_counts') or {})} | "
                         f"{_fmt_age(now, rec.get('trace_mtime'))} | "
                         f"{_fmt_age(now, rec.get('first_seen')).replace(' ago', '')} | "
                         f"{rec.get('reason', '')} | {step} |")
    if actions_log:
        lines += ["", "**This tick:**", ""] + [f"- {a}" for a in actions_log]
    lines.append("")
    return "\n".join(lines)


def _run_alert_body(rec: dict) -> str:
    failed = ", ".join(f"`{t}`" for t in rec.get("failed_tasks") or []) or "none in trace"
    return (f"**{rec['study_id']}** is **{rec['state']}**: {rec['reason']}.\n\n"
            f"- Failure class: `{rec.get('failure_class') or '—'}`"
            f" ({'recoverable' if rec.get('recoverable') else 'needs a person'})\n"
            f"- Next step: {rec.get('next_step') or '—'}\n"
            f"- Head job: {rec.get('job_id', '—')}\n"
            f"- Launch dir: `{rec.get('launch_dir', '—')}`\n"
            f"- Results: `{rec.get('outdir', '—')}`\n"
            f"- Recently failed tasks: {failed}\n\n"
            f"Close this issue once handled.\n")


def _footer(host: str) -> str:
    return f"\n_Opened by `aree autopilot tick` on `{host}`._\n"
