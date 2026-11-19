# Running a study on Klone (UW Hyak)

A workstation is enough for most studies. A deep methylation study is not:
40 EM-seq libraries of ~224 M read pairs each (`VALDIVIESO2025_POMS_RS_EMSEQ`)
would take about 33 days of Bismark alignment on a 10-core machine. Klone's
checkpoint (ckpt) partitions offer hundreds of idle cores. The catch is that
ckpt jobs can be preempted and requeued at any time, and Hyak documents a cap
of a few hours per CPU job. The `klone` profile and chunked alignment are
built around that limit.

## What runs where

| Piece | Where | Why |
|---|---|---|
| Nextflow head job | a non-preemptible partition you can use (e.g. `cpu-g2` under your group's account), 2 CPUs, 8 GB, a few days | if the head job is preempted, the whole run stops |
| Every task | `ckpt-all` (`-profile klone`) | idle capacity; each task is short and retried |
| Read-sized temporary files | the compute node's `/scr` (~350 GB per node) | chunk FASTQs, trimmed reads, Bismark temp files, merged BAMs and extractor context files never touch shared storage |
| Raw FASTQ, work dir, results | `/gscratch/scrubbed/$USER/aree/` | with chunked alignment, about 1.8x the raw reads beyond the raw reads themselves. **Scrubbed is a shared fileset (153.6 TB) that is often >97% full, and deletes files untouched for 21 days.** Check `hyakstorage /gscratch/scrubbed` before a run, and copy what you keep before 21 days |
| Repository checkout | `/gscratch/<group>/$USER/github/AREE` | small, and not scrubbed |

Do not run Nextflow, downloads or tools on a login node. Submit them as jobs.

## One-time setup

1. **SSH without a Duo prompt every time.** In `~/.ssh/config` on your own
   machine:

   ```
   Host klone klone.hyak.uw.edu
     HostName klone.hyak.uw.edu
     User <netid>
     ControlMaster auto
     ControlPath ~/.ssh/cm/%r@%h:%p
     ControlPersist 12h
   ```

   Run `mkdir -p ~/.ssh/cm`, then `ssh klone` once with your password and
   Duo. Later `ssh`/`scp` calls reuse that connection for 12 hours.

2. **Nextflow and Java** in a conda environment, e.g.
   `mamba create -n nextflow -c bioconda -c conda-forge nextflow openjdk=17`.
   Apptainer is on the compute nodes at `/usr/bin/apptainer`.

3. **Reference bundle.** Copy `genomic.fna` and `genomic.gtf` for the study's
   assembly to `/gscratch/scrubbed/$USER/aree/reference/<assembly>/`, and
   check the MD5s against the source.

4. **Checkout.** Clone the repository and check out the study branch (see
   [parallel_study_runs.md](parallel_study_runs.md)).

Compute nodes reach ENA over HTTPS through the cluster proxy. In October 2026
aggregate throughput was about 12 MB/s whether 12 or 24 files were fetched at
once: plan on ~15 h per 600 GB. Use the same resumable, MD5-checking download
pattern as on the workstation, submitted as a job.

## Running

```bash
#SBATCH --account=<group> --partition=cpu-g2 --cpus-per-task=2 --mem=8G --time=4-00:00:00
source <miniforge>/etc/profile.d/conda.sh && conda activate nextflow
export NXF_HOME=/gscratch/scrubbed/$USER/aree/nxf-home
export NXF_TEMP=/gscratch/scrubbed/$USER/aree/tmp   # not TMPDIR: see below
export AREE_SLURM_ACCOUNT=<group>          # account for the ckpt jobs (default srlab)
nextflow run workflows/methylation \
  -config config/<STUDY_ID>.config \
  -c run_resources_klone.config \
  -profile apptainer,klone \
  --outdir /gscratch/scrubbed/$USER/aree/results/<STUDY_ID> \
  --tracedir /gscratch/scrubbed/$USER/aree/results/<STUDY_ID>/pipeline_info
```

`-profile klone` (in `config/base.config`) sets the Slurm executor, the
`ckpt-all` queue, a 4 h 45 min limit on every task, retry on any failure
(six attempts), and per-process CPUs and memory. Override the account,
partition, work dir or Apptainer cache with `AREE_SLURM_ACCOUNT`,
`AREE_SLURM_QUEUE`, `AREE_WORKDIR` and `AREE_APPTAINER_CACHEDIR`. (Not `APPTAINER_CACHEDIR`: Hyak sets it to node-local `/tmp` inside jobs, which compute nodes cannot see.)

A study config's own `process` block is read after the profile, so its
longer time limits (written for a workstation) would win. Put Klone run
settings in a run-only config passed last with `-c`, as the workstation runs
do, and set at least:

```groovy
params.methylation.align_chunk_pairs = 20000000   // ~1-1.5 h Bismark tasks at 12 CPUs
process {
    withName: 'FASTQC|COUNT_READ_PAIRS|TRIM_ALIGN_CHUNK|MERGE_DEDUP_CHUNKS|BISMARK_.*|COVERAGE_FILTER|DMR_.*' { time = '4h 45m' }
}
```

## Chunked alignment

`methylation.align_chunk_pairs` trims, aligns and deduplicates each library
as chunks of read pairs, with every read-sized intermediate on the node's
`/scr` (see
[workflows/methylation/README.md](../workflows/methylation/README.md#chunked-alignment-methylation-align_chunk_pairs)).
A preempted chunk costs one chunk's work, not a library's, and the chunks of
all libraries align in parallel. The chunk task uses a pinned conda
environment, which `-profile klone` enables (cache in
`/gscratch/scrubbed/$USER/aree/conda-cache`, or `AREE_CONDA_CACHEDIR`).

**Do not export `TMPDIR` in the head job, and pin scratch to `/scr`.** Slurm
passes the head job's environment to every task. A `TMPDIR` on gscratch then
(1) puts `scratch true` back on shared storage, and (2) is forwarded into
containers, where `/gscratch/...` is not mounted (only `/mmfs1/gscratch`).
Quarto in `RENDER_REPORT` then fails with `NotFound ... tmpdir`. Use
`NXF_TEMP` for Nextflow's own temporary files. The klone profile also sets
`scratch = '/scr'` for the chunk, merge and extractor tasks.

**Select a subset of libraries with a directory of links, not a brace
pattern.** `fromFilePairs` names each sample from the part of the file name
before the first wildcard. A pattern such as
`ERR13913{46[2-9],47[0-9],48[01]}_{1,2}.fastq.gz` gave all 20 libraries the
sample name `ERR13913`. Link the wanted files into their own directory and
use `<dir>/ERR*_{1,2}.fastq.gz`.

## Provenance

Record the machine, partitions and storage in the study PR body (*Machine and
storage used*), as for any run. The workflow manifest lists the tool
versions, which come from the same pinned containers as the workstation runs.
Report a profile change that alters results in the study YAML.

## Autopilot

`aree autopilot` watches study runs on Klone and reports their state to a
GitHub status issue, opening `autopilot-alert` issues when a run needs a
person. **Phase 1 is report-only: it never submits or cancels a job.** The
plan for the later phases (automatic resume, starting queued studies) is in
[klone_autopilot_plan.md](klone_autopilot_plan.md).

### What it reports

Each tick (every 30 minutes) reads `squeue`, `sacct`, each run's
`pipeline_info/execution_trace.txt` and `.nextflow.log`, and free space on
scrubbed, and classifies every tracked run:

| State | Meaning | Alert |
|---|---|---|
| `pending` / `running` | head job queued / running, tasks finishing | no |
| `stalled` | head job running, but no task has finished for `stall_hours` (6) | yes |
| `stopped` | head job ended without success; the alert gives the failure class (preempted or timed-out head job, out-of-memory task, task walltime, missing input, container, Slurm controller, unknown) and the next step | yes |
| `finished` | head job completed with exit 0 | no |
| `unknown` | head job not in `squeue` and not in `sacct` | yes |

It also alerts when a run has been tracked for 14 days (scrubbed deletes files
untouched for 21) and when scrubbed free space falls below
`min_scrubbed_free_tb`. An open alert is not reopened while it is open; close
it once handled.

A GitHub Action (`autopilot-heartbeat.yml`) checks every hour that the status
issue was updated in the last 2 hours, and opens a `heartbeat missing` alert
if not, since the watchdog cannot report its own death. Closing the status
issue turns that check off.

### Setup

1. **GitHub token.** Create a fine-grained token for `sr320/AREE` only, with
   *Issues: read and write*. On Klone:
   `mkdir -p ~/.config/aree && umask 077 && cat > ~/.config/aree/github_token`
   (paste, then Ctrl-D).
2. **Environment.** A conda environment named `aree` with the repository
   installed (`pip install -e .`) in the checkout under `/gscratch/<group>`.
3. **Settings.** Check `config/autopilot.yaml`: `results_template` must match
   the `--outdir` your runs use. The ledger lives in `state_dir`
   (`~/.local/share/aree-autopilot`, not scrubbed).
4. **Name head jobs `aree-<STUDY_ID>`** (`sbatch --job-name=aree-<STUDY_ID>`)
   so the tick finds them. A run started otherwise, or one that already ended,
   is added with
   `aree autopilot track --study <STUDY_ID> --launch-dir <dir nextflow ran in> [--outdir ...] [--job-id ...]`.
5. **Try one tick** from an interactive job (`salloc`), not the login node:
   `aree autopilot tick --local-only`, then without `--local-only` to create
   the status issue.
6. **Schedule it.** Preferred, if scrontab is enabled for the account:
   edit the placeholders in `scripts/klone/autopilot.scrontab`, `mkdir -p
   /gscratch/srlab/$USER/aree-autopilot`, then `scrontab scripts/klone/autopilot.scrontab`.
   Otherwise, from the repository root:
   `sbatch scripts/klone/autopilot.sbatch <miniforge prefix>`. That job queues
   its successor 30 minutes out; stop the chain with
   `scancel -n aree-autopilot`.

`aree autopilot status` prints the last tick's report. When a run's outputs
are copied off scrubbed, or the run is abandoned,
`aree autopilot untrack <STUDY_ID>` removes it.
