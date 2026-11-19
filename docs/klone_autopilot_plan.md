# Plan: Klone autopilot (continuous monitoring, restart and backfill)

Status: **Phase 1 (report-only) implemented** (2026-11-19); Phases 2-4 planned.
Setup and use: [running_on_klone.md#autopilot](running_on_klone.md#autopilot).

## Goal

Keep Klone busy with AREE reanalyses without someone watching it:

1. Know, at all times, which study runs are active, stalled, failed or done.
2. Resume a run that stopped for a recoverable reason (preempted head job,
   Slurm controller timeout, walltime, node failure).
3. When fewer runs are active than the target, start the next study that is
   ready to run.
4. Report state somewhere visible off-cluster, and alert a person when a run
   needs a decision.

What it must not do: invent study metadata, mark a study PR ready for review,
merge anything, delete results, or work around Hyak's scrubbed-storage purge.
Curation stays a reviewed human (or human-approved Claude) step, as
[CLAUDE.md](../CLAUDE.md) and [parallel_study_runs.md](parallel_study_runs.md)
require.

## Constraints that shape the design

| Constraint | Consequence |
|---|---|
| Klone login needs Duo; laptop `ControlPersist` lasts 12 h | The monitor cannot live on a laptop or in GitHub Actions. It runs **on Klone** and pushes state out. |
| ckpt tasks are preempted; Nextflow already retries them (`-profile klone`, six attempts) | The autopilot watches the **head job**, not individual tasks. Task-level retry stays Nextflow's job. |
| The head job runs on a non-preemptible partition with a time limit (`cpu-g2`, ~4 days) | A run longer than that limit needs a scheduled `-resume`, not a failure alert. |
| `/gscratch/scrubbed` deletes files untouched for 21 days and is often >97% full | Gate new starts on free space; alert at day 14 of a run; copy kept outputs to group storage. Never `touch` files to dodge the purge. |
| One study per `study/<ACCESSION>/<STUDY_ID>` branch; the push is the claim | The queue and the run state map one-to-one onto study branches and their draft PRs. |
| Curation (phenotype, stressor, comparisons) needs judgment | Only studies whose YAML + config are already reviewed are started automatically. |

## Architecture

```
          Klone                                     GitHub (sr320/AREE)
┌────────────────────────────────────┐       ┌───────────────────────────────┐
│ watchdog job (every 30 min)        │ API   │ issues labelled klone-queue   │
│  aree autopilot tick               │──────▶│  (the run queue)              │
│   ├ squeue / sacct     (Slurm)     │◀──────│ study/ branches + draft PRs   │
│   ├ trace.txt, .nextflow.log (NF)  │       │ "Klone autopilot" status issue│
│   ├ hyakstorage        (storage)   │       │  (heartbeat comment, edited)  │
│   ├ decide (rules, no LLM)         │       └───────────────┬───────────────┘
│   ├ act: sbatch / scancel / resume │                       │
│   └ state.json ledger              │       GitHub Actions cron (hourly):
│ head jobs (cpu-g2) ─▶ tasks (ckpt) │       heartbeat older than 2 h → open
└────────────────────────────────────┘       an alert issue (dead-man switch)
                │ optional
                ▼
   Anthropic API (claude -p, headless): failure-log triage notes and
   draft study YAMLs for the "needs curation" tier — never on the
   restart path, never merged without review.
```

The **APIs** used, and for what:

- **Slurm CLI** (`squeue --json`, `sacct --json`, `sbatch --parsable`,
  `scancel`): job state and actions. If Hyak exposes `slurmrestd` to users,
  the same module can switch to it; the CLI is the baseline.
- **Nextflow**: `pipeline_info/trace.txt` (task progress, last update time),
  `.nextflow.log` (error class), and `nextflow log` for run names to resume.
  Optional: launch with `-with-tower` to Seqera Platform for a hosted run UI
  and API, if the lab wants one; not required.
- **GitHub REST API** (`gh api`, fine-grained token limited to `sr320/AREE`:
  issues read/write, contents read, pull requests read/write): the queue, the
  status issue, comments on each study's draft PR.
- **Anthropic API** (optional, Phase 4): judgment steps only, see below.

## Run state machine

One record per study run in `state.json`
(`/gscratch/<group>/$USER/aree-autopilot/`, not scrubbed), mirrored to the
status issue every tick.

```
queued ─▶ staging ─▶ running ─▶ finishing ─▶ done (PR comment, awaits review)
             │          │  ▲
             │          ▼  │ resume (bounded)
             │       stalled / stopped
             ▼          │
          blocked ◀─────┘ (needs a person: alert)
```

| State | Detected by | Autopilot action |
|---|---|---|
| `queued` | issue labelled `klone-queue` + `ready-to-run`, study branch exists, config validates | start when capacity and storage allow |
| `staging` | download job in `squeue` | wait; on failure resubmit (resumable, MD5-checked) up to 3× |
| `running` | head job `RUNNING`, trace updated within `stall_hours` | nothing |
| `stalled` | head job running, no trace update for `stall_hours` (default 6) **and** no AREE tasks in `squeue` | `scancel` head job, then resume |
| `stopped` | head job gone, no `pipeline_info/*manifest*` / success marker | classify (below), resume if recoverable |
| `finishing` | Nextflow exited 0, manifest present | run `aree harmonize` on the branch, commit standardized outputs, push, comment on the PR |
| `done` | outputs pushed | stop tracking; a person reviews and marks the PR ready |
| `blocked` | unrecoverable class, or resume budget spent | alert; no further action until the label `autopilot-retry` is set |

**Stop classification** (from `sacct` state/exit code and the last lines of
`.nextflow.log`):

| Class | Signal | Action |
|---|---|---|
| head job preempted / node fail / `TIMEOUT` | `sacct` `PREEMPTED`, `NODE_FAIL`, `TIMEOUT` | resume with `-resume <run-name>` |
| Slurm controller timeout | `Socket timed out` in log | resume |
| task exhausted retries, OOM | exit 137 / `OUT_OF_MEMORY` on a task | write a run-only override raising that process's memory one step (max one escalation), resume |
| task exhausted retries, walltime | task `TIMEOUT` at 4 h 45 m | `blocked`: needs smaller `align_chunk_pairs` (a reviewed config change) |
| missing input / checksum mismatch / container pull failure | log patterns | `blocked` |
| anything else | — | `blocked`, with the last 50 log lines in the alert |

**Resume budget**: at most 3 resumes per run per 24 h and 10 in total; then
`blocked`. This prevents a restart loop burning allocation.

## Keeping Klone full (backfill)

Each tick, after handling existing runs:

1. `active = runs in staging|running|finishing`.
2. Read limits from `config/autopilot.yaml`, e.g.
   `max_active: 3`, `max_active_methylation: 1`, `min_scrubbed_free_tb: 8`,
   `max_queue_age_days`.
3. If `active < max_active`, storage is above the floor, and the group's
   `hyakstorage` quota has room for the next study's FASTQ manifest size × 2.8
   (raw + ~1.8× intermediates), pop the highest-priority `queued` item and
   submit its download job, then its head job (dependency `afterok`).
4. If the queue is empty, update the status issue ("queue empty") so a person
   adds studies. Optionally (Phase 4) trigger the curation drafter.

### The queue

[candidate_studies.md](candidate_studies.md) is prose and stays the screening
record. The machine-readable queue is GitHub issues, one per study:

- title `[klone-queue] <STUDY_ID> (<ACCESSION>)`
- labels: `klone-queue`, `assay:<rnaseq|methylation|…>`, `priority:<1-3>`,
  and `ready-to-run` once the study branch has a reviewed YAML and config
- body: a small YAML block (`study_id`, `accession`, `branch`, `config`,
  `est_raw_gb`)

Two tiers:

- **ready-to-run**: claimed with `aree start-study`, YAML and config reviewed.
  The autopilot downloads and runs these with no further input.
- **needs-curation**: accession screened, no reviewed YAML yet. The autopilot
  never runs these. Phase 4 can draft them.

`aree autopilot enqueue --from-branch study/<ACC>/<ID>` creates the issue from
an existing claim, so the queue never drifts from the branches.

## Visibility and alerting

- **Status issue** "Klone autopilot status": one comment edited every tick
  with a table (study, state, head job id, % tasks done, last trace update,
  resumes used, scratch age in days) and a `last_tick` timestamp.
- **Study PR comments** on state changes only (started, resumed and why,
  blocked and why, finished).
- **Alerts**: a new issue labelled `autopilot-alert` (GitHub notifies by
  email/mobile) for `blocked`, scratch age ≥14 days, scrubbed free space
  below the floor, or a tick error.
- **Dead-man switch**: `.github/workflows/autopilot-heartbeat.yml`, hourly
  cron, reads `last_tick` from the status issue; older than 2 h opens an
  alert. This catches the watchdog itself dying, which the watchdog cannot
  report.
- Optional: a daily Claude Code scheduled routine (`/schedule`) that reads the
  status issue and posts a short summary to Slack. It reads GitHub only; it
  never touches Klone.

## Where the watchdog runs

Preferred: **`scrontab`** on Klone (`*/30 * * * *`, account `srlab`,
partition `cpu-g2`, 1 CPU, 2 GB, 10 min). Verify first that Hyak allows
scrontab for the group.

Fallback: a self-resubmitting job, `scripts/klone/autopilot.sbatch`, that
runs one tick and then `sbatch --begin=now+30minutes` itself. A lock file
(`flock`) prevents two ticks overlapping. Either way, the tick is a short job,
not a process sitting on a login node.

## Use of the Anthropic API (optional, Phase 4)

Kept off the restart path, which is deterministic rules only. Two uses:

1. **Blocked-run triage**: `claude -p` reads the tail of `.nextflow.log` and
   the failing task's `.command.err` and writes a short diagnosis into the
   alert issue. Advice only; it takes no action.
2. **Curation drafter**: when the ready queue is empty, take the top
   `needs-curation` item, run `aree start-study` and `aree fetch-samplesheet`,
   and have Claude draft the study YAML (comparisons, phenotype and stressor
   terms from the controlled vocabularies) into the draft PR with every
   inferred field flagged. A person reviews and adds `ready-to-run`. Never
   auto-labelled.

Key on Klone in `~/.config/aree/anthropic_key` (mode 600) with a monthly spend
cap set in the Anthropic console.

## Code to add

| Path | Contents |
|---|---|
| `src/autopilot/state.py` | ledger load/save, state transitions, resume budget |
| `src/autopilot/slurm.py` | `squeue`/`sacct` JSON parsing, `sbatch`/`scancel` wrappers |
| `src/autopilot/nextflow.py` | trace parsing (progress, last update), log classification, run-name lookup |
| `src/autopilot/storage.py` | free space (statvfs now; `hyakstorage` parsing if Phase 0 shows it is needed) |
| `src/autopilot/github.py` | queue issues, status comment, PR comments, alerts (via `gh api`) |
| `src/autopilot/tick.py` | one pass: observe → classify → act → report |
| `src/aree/cli.py` | `aree autopilot tick [--dry-run]`, `status`, `enqueue`, `pause`/`unpause` |
| `config/autopilot.yaml` | limits, thresholds, account/partition, paths |
| `scripts/klone/head.sbatch.tmpl` (Phase 2) | head job (the command from [running_on_klone.md](running_on_klone.md)), with `-resume` and run name |
| `scripts/klone/download.sbatch.tmpl` (Phase 3) | resumable, MD5-checked ENA download job |
| `scripts/klone/autopilot.sbatch`, `autopilot.scrontab`, `autopilot_tick.sh` | watchdog (done, Phase 1) |
| `.github/workflows/autopilot-heartbeat.yml` | dead-man switch |
| `tests/test_autopilot_*.py` | fixtures of real `squeue`/`sacct`/trace/log output; each classification; resume budget; backfill gating; dry-run makes no calls |
| `docs/running_on_klone.md` | new "Autopilot" section; link this plan |

Safety switches: a `PAUSE` file in the autopilot dir or a `autopilot-paused`
label on the status issue makes every tick report-only. `--dry-run` logs the
actions it would take.

## Rollout

| Phase | What runs | Exit criterion |
|---|---|---|
| 0. Verify | scrontab allowed? head-job partition and time limit; GitHub token on Klone; `hyakstorage` output format | each answered and recorded in this doc |
| 1. Observe | tick in report-only mode: status issue, heartbeat, alerts; no `sbatch`/`scancel` | 1 week of correct status for the runs already on Klone (e.g. `VALDIVIESO2025_POMS_RS_EMSEQ`) |
| 2. Resume | resume of `stopped`/`stalled` runs, OOM escalation, budget | one real preempted or timed-out head job resumed without a person |
| 3. Backfill | download + start from `ready-to-run` queue; finishing step (harmonize, push, PR comment) | two studies run end to end from queue to "awaiting review" |
| 4. Assist | Anthropic API triage notes and curation drafts | drafts judged useful in review of ~3 studies |

Phases 1–3 are a few days of work each; Phase 1 alone already answers "is
anything running?".

## Open questions

1. Is `scrontab` enabled for `srlab` on Klone, and is a 30-minute 1-CPU tick
   acceptable to Hyak support?
2. Which non-preemptible partition and maximum time for head jobs (`cpu-g2`
   and 4 days assumed)?
3. Target concurrency: how many studies at once, and how many deep
   methylation runs?
4. Where kept outputs go after a run (group `/gscratch/srlab`, the
   workstation, or only the committed standardized results)?
5. Who receives alerts, and is GitHub notification enough or should alerts
   also go to Slack?
