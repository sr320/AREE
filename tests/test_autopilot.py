"""Klone autopilot, Phase 1 (report-only): parsing, run classification, the tick, GitHub reporting.

Slurm output is recorded text in the formats `autopilot.slurm` asks for; GitHub
is an in-memory fake. Nothing here calls Slurm or the network.
"""
import datetime as dt
import os

import pytest

from autopilot import github as gh_mod
from autopilot.config import load_config
from autopilot.nextflow import classify_log, trace_summary
from autopilot.slurm import elapsed_hours, parse_sacct, parse_squeue
from autopilot.state import assess, load_state
from autopilot.tick import tick, track, untrack

NOW = dt.datetime(2026, 11, 19, 12, 0, tzinfo=dt.timezone.utc).timestamp()
HOUR = 3600

SQUEUE = """\
2001|aree-VALDIVIESO2025_POMS_RS_EMSEQ|RUNNING|1-02:03:04|4-00:00:00|cpu-g2|/gscratch/srlab/u/github/AREE|None
2002|aree-STALLED_STUDY|RUNNING|10:00:00|4-00:00:00|cpu-g2|/gscratch/srlab/u/stalled|None
2003|aree-autopilot|RUNNING|00:01|00:15:00|cpu-g2|/gscratch/srlab/u/github/AREE|None
3001|nf-TRIM_ALIGN_CHUNK_(1)|RUNNING|01:00:00|04:45:00|ckpt-all|/gscratch/scrubbed/u/aree/work|None
3002|nf-TRIM_ALIGN_CHUNK_(2)|PENDING|00:00|04:45:00|ckpt-all|/gscratch/scrubbed/u/aree/work|Priority
"""

TRACE_HEADER = "task_id\thash\tnative_id\tname\tstatus\texit\tsubmit\tduration\trealtime\t%cpu\tpeak_rss\tpeak_vmem\n"

OOM_LOG = """\
Nov-18 10:00:00.000 [Task monitor] DEBUG nextflow.executor.GridTaskHandler - Socket timed out on send/recv, retrying
Nov-18 11:00:00.000 [Task monitor] ERROR nextflow.processor.TaskProcessor - Error executing process > 'MERGE_DEDUP_CHUNKS (ERR1)'
ERROR ~ Error executing process > 'MERGE_DEDUP_CHUNKS (ERR1)'
Caused by:
  Process `MERGE_DEDUP_CHUNKS (ERR1)` terminated with an error exit status (137)
"""


def _runner(squeue_text=SQUEUE, sacct_by_job=None):
    sacct_by_job = sacct_by_job or {}
    calls = []

    def run(args):
        calls.append(args)
        if args[0] == "squeue":
            return squeue_text
        if args[0] == "sacct":
            return sacct_by_job.get(args[2], "")
        raise AssertionError(f"unexpected command {args}")
    run.calls = calls
    return run


class FakeGitHub:
    def __init__(self, issues=None):
        self.issues = issues or []
        self.labels = set()
        self.next = 100

    def open_issues(self, label):
        return [i for i in self.issues if label in i["labels"] and i["state"] == "open"]

    def ensure_label(self, name, color):
        self.labels.add(name)

    def create_issue(self, title, body, labels):
        self.next += 1
        issue = {"number": self.next, "title": title, "body": body, "labels": labels, "state": "open"}
        self.issues.append(issue)
        return issue

    def update_issue_body(self, number, body):
        next(i for i in self.issues if i["number"] == number)["body"] = body


@pytest.fixture
def cfg(tmp_path):
    c = load_config(tmp_path / "absent.yaml")
    c.update(state_dir=str(tmp_path / "state"),
             results_template=str(tmp_path / "results" / "{study_id}"),
             scrubbed_root=str(tmp_path))
    return c


def _write_trace(outdir, rows, mtime):
    info = outdir / "pipeline_info"
    info.mkdir(parents=True, exist_ok=True)
    path = info / "execution_trace.txt"
    path.write_text(TRACE_HEADER + "".join(
        f"{i}\th\t{9000 + i}\t{name}\t{status}\t0\t-\t-\t-\t-\t-\t-\n" for i, (name, status) in enumerate(rows)))
    os.utime(path, (mtime, mtime))


# --- parsing ---------------------------------------------------------------

def test_parse_squeue_and_elapsed():
    jobs = parse_squeue(SQUEUE + "garbage line\n")
    assert [j["job_id"] for j in jobs] == ["2001", "2002", "2003", "3001", "3002"]
    assert jobs[0]["work_dir"] == "/gscratch/srlab/u/github/AREE"
    assert elapsed_hours("1-02:03:04") == pytest.approx(26 + 3 / 60 + 4 / 3600)
    assert elapsed_hours("10:00") == pytest.approx(10 / 60)
    assert elapsed_hours("") == 0.0


def test_parse_sacct_reduces_state():
    row = parse_sacct("2001|aree-X|CANCELLED by 12345|0:15|2026-11-18T00:00:00|2026-11-18T05:00:00|05:00:00\n")
    assert row["state"] == "CANCELLED"
    assert row["state_detail"] == "CANCELLED by 12345"
    assert parse_sacct("") is None


def test_trace_summary_counts(tmp_path):
    _write_trace(tmp_path, [("A", "COMPLETED"), ("B", "CACHED"), ("C", "FAILED")], NOW - HOUR)
    s = trace_summary(tmp_path)
    assert s["status_counts"] == {"COMPLETED": 1, "CACHED": 1, "FAILED": 1}
    assert s["failed_tasks"] == ["C"]
    assert s["trace_mtime"] == pytest.approx(NOW - HOUR)
    assert trace_summary(tmp_path / "none")["trace_mtime"] is None


@pytest.mark.parametrize("log,cls,recoverable", [
    (OOM_LOG, "task_oom", True),
    ("ERROR ~ Process `X` exceeded running time limit (4h 45m)", "task_walltime", False),
    ("ERROR ~ Failed to pull singularity image\n FATAL: ...", "container", False),
    ("ERROR ~ Cannot find any reads matching: /x/*_{1,2}.fastq.gz", "input_missing", False),
    ("ERROR ~ Unable to contact slurm controller (connect failure)", "slurm_controller", True),
    ("ERROR ~ something nobody has seen before", "unknown", False),
])
def test_classify_log(log, cls, recoverable):
    got_cls, got_rec, step = classify_log(log)
    assert (got_cls, got_rec) == (cls, recoverable)
    assert step


def test_retried_socket_timeouts_before_the_error_do_not_decide_the_class():
    # OOM_LOG has a benign, retried socket timeout before the final error.
    assert classify_log(OOM_LOG)[0] == "task_oom"


# --- assess -----------------------------------------------------------------

def _head(state="RUNNING", elapsed="10:00:00"):
    return {"state": state, "elapsed": elapsed, "reason": "None"}


def test_assess_running_and_stalled():
    fresh = {"trace_mtime": NOW - HOUR}
    stale = {"trace_mtime": NOW - 7 * HOUR}
    assert assess(_head(), None, fresh, "", NOW, 6)["state"] == "running"
    stalled = assess(_head(), None, stale, "", NOW, 6)
    assert stalled["state"] == "stalled" and stalled["recoverable"]
    # No trace yet: fall back to the head job's elapsed time.
    assert assess(_head(elapsed="02:00:00"), None, {"trace_mtime": None}, "", NOW, 6)["state"] == "running"
    assert assess(_head(elapsed="07:00:00"), None, {"trace_mtime": None}, "", NOW, 6)["state"] == "stalled"
    assert assess(_head("PENDING"), None, {}, "", NOW, 6)["state"] == "pending"


@pytest.mark.parametrize("sacct_state,exit_code,log,state,cls,recoverable", [
    ("COMPLETED", "0:0", "", "finished", "", False),
    ("PREEMPTED", "0:0", "", "stopped", "head_interrupted", True),
    ("TIMEOUT", "0:0", "", "stopped", "head_interrupted", True),
    ("NODE_FAIL", "0:0", "", "stopped", "head_interrupted", True),
    ("FAILED", "1:0", OOM_LOG, "stopped", "task_oom", True),
    ("CANCELLED", "0:15", "", "stopped", "cancelled", False),
])
def test_assess_ended_head_job(sacct_state, exit_code, log, state, cls, recoverable):
    acct = {"state": sacct_state, "state_detail": sacct_state, "exit_code": exit_code, "end": "x"}
    r = assess(None, acct, {}, log, NOW, 6)
    assert (r["state"], r["failure_class"], r["recoverable"]) == (state, cls, recoverable)


def test_assess_without_any_slurm_record():
    assert assess(None, None, {}, "Execution complete -- Goodbye", NOW, 6)["state"] == "finished"
    assert assess(None, None, {}, "", NOW, 6)["state"] == "unknown"


# --- tick -------------------------------------------------------------------

def test_tick_report_only(cfg, tmp_path):
    _write_trace(tmp_path / "results" / "VALDIVIESO2025_POMS_RS_EMSEQ", [("A", "COMPLETED")] * 3, NOW - HOUR)
    _write_trace(tmp_path / "results" / "STALLED_STUDY", [("A", "COMPLETED")], NOW - 8 * HOUR)
    # A run tracked earlier whose head job has since been preempted.
    launch = tmp_path / "launch_gone"
    launch.mkdir()
    (launch / ".nextflow.log").write_text("Nov-18 [main] DEBUG nextflow.Session - Session aborted\n")
    track(cfg, "GONE_STUDY", str(launch), job_id="1999", now=NOW - 15 * 24 * HOUR)

    run = _runner(sacct_by_job={"1999": "1999|aree-GONE_STUDY|PREEMPTED|0:0|a|b|3-00:00:00\n"})
    gh = FakeGitHub()
    res = tick(cfg, now=NOW, run=run, gh=gh, free=lambda p: 50.0, host="klone-test")

    runs = res["runs"]
    assert set(runs) == {"VALDIVIESO2025_POMS_RS_EMSEQ", "STALLED_STUDY", "GONE_STUDY"}  # not the watchdog
    assert runs["VALDIVIESO2025_POMS_RS_EMSEQ"]["state"] == "running"
    assert runs["VALDIVIESO2025_POMS_RS_EMSEQ"]["launch_dir"] == "/gscratch/srlab/u/github/AREE"
    assert runs["STALLED_STUDY"]["state"] == "stalled"
    assert runs["GONE_STUDY"]["state"] == "stopped"
    assert runs["GONE_STUDY"]["failure_class"] == "head_interrupted"

    # Report-only: only squeue and sacct, never sbatch or scancel.
    assert {c[0] for c in run.calls} == {"squeue", "sacct"}

    # Status issue with a heartbeat, and one alert per condition.
    status = gh.open_issues("autopilot-status")
    assert len(status) == 1
    assert gh_mod.parse_heartbeat(status[0]["body"]).timestamp() == NOW
    assert "Phase 2 would: -resume" in status[0]["body"]
    titles = sorted(i["title"] for i in gh.open_issues("autopilot-alert"))
    assert titles == [
        "[autopilot-alert] GONE_STUDY: on scrubbed for 14+ days",
        "[autopilot-alert] GONE_STUDY: stopped (job 1999)",
        "[autopilot-alert] STALLED_STUDY: stalled (job 2002)",
    ]

    # The ledger and the local status page are written.
    assert load_state(tmp_path / "state" / "state.json")["last_tick"] == NOW
    assert (tmp_path / "state" / "status.md").read_text() == res["body"]

    # A second tick edits the same status issue and opens no duplicate alerts.
    res2 = tick(cfg, now=NOW + 1800, run=run, gh=gh, free=lambda p: 50.0, host="klone-test")
    assert res2["opened"] == []
    assert len(gh.open_issues("autopilot-status")) == 1
    assert len(gh.open_issues("autopilot-alert")) == 3


def test_tick_low_storage_and_untrack(cfg):
    gh = FakeGitHub()
    res = tick(cfg, now=NOW, run=_runner(squeue_text=""), gh=gh, free=lambda p: 2.0, host="h")
    assert res["alerts"] == ["scrubbed free space below floor"]
    assert "No runs tracked" in res["body"]

    track(cfg, "X", "/nowhere", now=NOW)
    assert untrack(cfg, "X") is True
    assert untrack(cfg, "X") is False


def test_tick_local_only_needs_no_github(cfg):
    res = tick(cfg, now=NOW, run=_runner(), gh=None, free=lambda p: None, host="h")
    assert res["status_issue"] is None
    assert "Scrubbed free: n/a" in res["body"]


# --- heartbeat ---------------------------------------------------------------

def _status_issue(last_tick):
    body = gh_mod.heartbeat_marker(dt.datetime.fromtimestamp(last_tick, dt.timezone.utc)) + "\nstatus"
    return {"number": 7, "title": "Klone autopilot status", "body": body,
            "labels": ["autopilot-status"], "state": "open"}


def test_heartbeat(cfg):
    now = dt.datetime.fromtimestamp(NOW, dt.timezone.utc)
    assert gh_mod.check_heartbeat(FakeGitHub(), cfg, 2, now)[0] == "absent"

    gh = FakeGitHub([_status_issue(NOW - HOUR)])
    assert gh_mod.check_heartbeat(gh, cfg, 2, now)[0] == "ok"
    assert gh.open_issues("autopilot-alert") == []

    gh = FakeGitHub([_status_issue(NOW - 3 * HOUR)])
    assert gh_mod.check_heartbeat(gh, cfg, 2, now)[0] == "stale"
    assert gh_mod.check_heartbeat(gh, cfg, 2, now)[0] == "stale"
    assert [i["title"] for i in gh.open_issues("autopilot-alert")] == ["[autopilot-alert] heartbeat missing"]
