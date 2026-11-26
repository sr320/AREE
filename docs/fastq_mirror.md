# Mirroring study FASTQ to gannet

Downloading is the slowest step of a raw reanalysis and needs no compute. So
it runs ahead of the runs: `aree mirror-fastq`, from raven, a little every
night, until the FASTQ of every study still to be analysed is on gannet. A
study then stages from a lab copy with `aree fetch-fastq` instead of a
multi-day download from ENA.

This is the same design as project-white-whale's mirror, ported. AREE keeps
its own area on gannet and its own environment variables; the two mirrors do
not share files, ledgers or staging directories.

## What gets mirrored

With no arguments, every real registered study that has a committed
`data/studies/<STUDY_ID>/fastq_manifest.tsv`, a BioProject accession, and an
`analysis_status` of `in_progress` or `not_started`. Studies in progress go
first. The status is read from the study YAMLs in the checkout the command
runs in, so keep that checkout on an up-to-date `main`.

| To mirror | Use |
|---|---|
| the default set | `aree mirror-fastq` |
| one study, whatever its status | `--study STUDY_ID` (repeatable) |
| harmonized studies too (an archive copy) | `--all-studies` |
| a candidate that is not registered yet | `--accession PRJ…`, with `--library-strategy RNA-Seq` to leave out other assays |

The file list always comes from ENA's read-run report at the time of the run,
restricted to the runs in the study's manifest. Several studies on one
BioProject share one folder.

## Where the reads come from

Measured from raven on 2026-10-23 for White Whale, ENA served FASTQ at
0.04–1 MB/s per connection and NCBI's AWS open-data bucket
(`sra-pub-run-odp`) served the same run at about 25 MB/s. So by default
(`--source ncbi`) each run is:

1. located with NCBI's locator (`locate.ncbi.nlm.nih.gov/sdl/2/retrieve`),
   which gives the S3 link, size and MD5 of the full-quality `.sra`;
2. downloaded into staging (resuming a partial file) and checked against that
   MD5;
3. converted on raven with `fasterq-dump --split-3` and gzipped (pigz if
   available, else gzip with one process per mate) under the file names ENA
   lists (`<RUN>_1.fastq.gz`, `<RUN>_2.fastq.gz`), so the study config's
   `reads` glob needs no change;
4. checked: both mates have the same number of reads, and reads and bases
   equal ENA's `read_count` and `base_count` for the run. An unpaired leftover
   file counts towards those totals and is then discarded, as the study
   manifests leave ENA's third file out.

A run falls back to ENA's FASTQ (downloaded over HTTPS and checked against
ENA's MD5) when NCBI has no full-quality copy, a check fails, or only some of
its files are missing. `--source ena` uses ENA for everything.

One case to know: a run with no ENA-generated FASTQ (`source` =
`submitted_ftp` in the manifest, e.g. PRJEB105019's ERR15996548) is mirrored
under the submitter's file names, which a `<RUN>_{1,2}.fastq.gz` glob will
not match. Rename or symlink those after staging.

## On gannet

```
raw/<ACCESSION>/<file>                  every FASTQ, named as ENA lists it
raw/<ACCESSION>/mirror_manifest.tsv     the ledger, one row per verified file
```

under `/volume1/v1_web/aree`, served at
`https://gannet.fish.washington.edu/v1_web/aree/raw/`
(`config/data_host.yaml`). These are public reads; nothing private is
mirrored. Ledger columns:

| Column | |
|---|---|
| `run_accession`, `mate`, `file` | as in ENA's report |
| `url` | where it came from: the S3 `.sra` or ENA's FASTQ |
| `md5`, `bytes` | **of the file on gannet**, recomputed there after upload |
| `source` | `ncbi_sra`, `fastq_ftp` or `submitted_ftp` |
| `ena_md5` | ENA's MD5 for the file at the time it was mirrored |
| `sra_md5`, `read_count` | for converted runs: the `.sra` MD5 and reads per mate |
| `mirrored_on` | date |

A converted FASTQ holds the same reads and qualities as ENA's, but its read
names and gzip bytes differ, so **its MD5 is not ENA's**. Check copies taken
from gannet against the ledger's `md5`, not the study's `fastq_manifest.tsv`.
`aree fetch-fastq` does this for you and records which MD5 it checked.

A file is skipped when the ledger has it against ENA's current MD5 and it is
on gannet at its recorded size. If ENA reprocesses a run (new MD5), the next
run fetches it again. A second copy started while one is running with the
same staging directory exits at once, so overlapping cron runs are harmless.
The exit code is 1 when any run failed.

Each ledger update takes a lock on gannet (the directory
`raw/<ACCESSION>/.ledger.lock`), re-reads the ledger, adds its row and writes
it back, so mirrors on two machines never drop each other's rows. If a mirror
is killed in that moment the lock stays behind; after 5 minutes the next
writer stops and names the directory to remove.

## How far along it is

```bash
aree mirror-fastq --plan
```

lists, per accession, the studies, the files wanted, how many are mirrored
and verified, and the GB left. It downloads nothing. On 2026-11-25 the
default set was seven studies and about 1.9 TB, none of it mirrored yet.

## Running it on raven

Once, in a checkout on `main` with the package installed (`pip install -e .`
into a venv):

```bash
export AREE_SRATOOLKIT=/home/shared/sratoolkit.3.0.2-ubuntu64
export AREE_PIGZ=/home/shared/pigz-2.8/bin/pigz
export AREE_MIRROR_STAGING=/home/shared/8TB_HDD_03/$USER/aree/mirror-staging
aree mirror-fastq --plan
```

Disk: staging holds the `.sra` and the gzipped FASTQ of the runs in flight.
`--tmp` holds fasterq-dump's scratch and the uncompressed FASTQ, roughly 3–4×
the `.sra`, so 100–200 GB for one of the deepest EM-seq runs. Put `--tmp` on
fast local disk and check `df -h` there before raising `--jobs`: White
Whale's mirror and the workflows use the same disks.

Nightly from cron: `scripts/raven/mirror_fastq.cron` has an entry to adapt
(`crontab -e`). It caps one night at 500 GB and 3 runs at a time. Cron needs
raven's SSH key for gannet to work without a passphrase (the commands use
`BatchMode`). Or run the same command by hand in `tmux`, without `--max-gb`.

The first run creates `/volume1/v1_web/aree` on gannet.

| Option | Meaning |
|---|---|
| `--study ID`, `--accession PRJ…` | Only these (repeatable), e.g. the next study to be run |
| `--all-studies` | The default set plus harmonized studies |
| `--source ncbi\|ena` | `ncbi` (default) converts NCBI's `.sra`, with ENA as fallback per run; `ena` uses ENA's FASTQ only |
| `--jobs N` | Runs at once (default 4) |
| `--threads N` | fasterq-dump / pigz threads per run (default 4) |
| `--tmp DIR` | fasterq-dump scratch (default `<staging>/tmp`) |
| `--limit-rate MB/s` | Total download cap across all jobs |
| `--max-gb G` | Stop starting new runs once G GB (ENA's sizes) are queued in this run |
| `--no-verify-remote` | Skip the MD5 on gannet (rsync already checks what it transfers) |
| `--dest DIR` | Mirror into a local directory instead (a mounted share, or testing) |

## Staging a study from the mirror

```bash
aree fetch-fastq --study STUDY_ID --dry-run
aree fetch-fastq --study STUDY_ID
```

stages the study's files into `$AREE_RAW/<STUDY_ID>/` (default
`data/raw/<STUDY_ID>/`, the folder the study config reads). It works on any
machine that can reach gannet over HTTPS, Klone included. For each file in
the study's `fastq_manifest.tsv`:

- **mirrored against the manifest's current ENA MD5:** copied from gannet over
  HTTPS (resumable) and checked against the ledger's `md5`;
- **not mirrored, or mirrored against an older ENA MD5:** downloaded from ENA
  and checked against ENA's MD5;
- **mirror copy fails** (download or MD5): deleted and fetched from ENA.

A partial file is only resumed from the source it was started from. Files
already verified from either source are kept, and each gets a `<file>.md5`
sidecar so a re-run does not re-hash it. Free space is checked first.

`data/studies/<STUDY_ID>/fastq_staged.tsv` records, per file, the MD5
actually checked, the `source` (`mirror` or `ena`), ENA's MD5, the machine
and the date. Commit it on the study branch with the run's results: it is
the integrity record for the reads that run used.

`--source mirror` refuses to start unless every file is usable from the
mirror; `--source ena` ignores the mirror; `--mirror-url` points at another
copy; `--run` limits it to some runs.

On Klone, `/gscratch/scrubbed` deletes files untouched for 21 days. A study
staged there and then left can lose its FASTQ; re-running `aree fetch-fastq`
restores only what is missing, from gannet.
