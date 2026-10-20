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
| Raw FASTQ, work dir, results | `/gscratch/scrubbed/$USER/aree/` | room for ~2–3 TB of intermediates. **Scrubbed deletes files untouched for 21 days:** copy what you keep before then |
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
export TMPDIR=/gscratch/scrubbed/$USER/aree/tmp NXF_TEMP=$TMPDIR
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
    withName: 'TRIM_GALORE|BISMARK_.*|SPLIT_FASTQ_PAIRS|FASTQC|DMR_.*' { time = '4h 45m' }
}
```

## Chunked alignment

`methylation.align_chunk_pairs` splits each trimmed library into chunks of
read pairs. It aligns each chunk as its own task, then concatenates the chunk
BAMs in order before deduplication (see
[workflows/methylation/README.md](../workflows/methylation/README.md#chunked-alignment-methylation-align_chunk_pairs)).
A preempted chunk costs one chunk's work, not a library's, and the chunks of
all libraries align in parallel.

## Provenance

Record the machine, partitions and storage in the study PR body (*Machine and
storage used*), as for any run. The workflow manifest lists the tool
versions, which come from the same pinned containers as the workstation runs.
Report a profile change that alters results in the study YAML.
