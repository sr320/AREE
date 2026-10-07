# Running studies on several machines

AREE's backlog ([candidate_studies.md](candidate_studies.md)) is worked one
study at a time, and a raw reanalysis can tie up a machine for a day or more.
To run several at once, each job is **one study on its own branch**, claimed
from whichever machine is going to run it, and merged back as a pull request.

```
study/<ACCESSION>/<STUDY_ID>          e.g. study/PRJNA690951/CIBNOR2021_HEAT_RRSS
```

The accession names the job, because that is what you pick from the backlog.
The `study_id` is still the key everywhere inside AREE.

## Why study PRs don't conflict

A study branch changes only files that belong to its study:

| Path | What |
|---|---|
| `registry/studies/<STUDY_ID>.yaml` | the registration |
| `config/<STUDY_ID>.config` | the Nextflow run configuration |
| `data/studies/<STUDY_ID>/` | sample sheet, FASTQ manifest, provenance, standardized results, workflow manifests, reports |
| `tests/test_<study>_study.py` | curation tests |

Nothing else, including `registry/study_registry.csv`. That file is a generated
index (`aree build-registry`) and is not committed. Shared outputs (evidence
table, meta-analysis, evidence cards, dashboard) are rebuilt on `main` after
merges and are never committed from a study branch either. So two study PRs
never touch the same file, and they merge in any order.

CI enforces this. On any PR from a `study/` branch, the `study-pr-scope` job
runs `aree check-study-scope`. It fails if:

- the PR changes a file outside that study's own paths;
- the YAML is missing, or its `study_id` differs from the branch;
- the branch's accession is not in the YAML's `accessions` block.

It also warns when another registered study already lists the same accession.

**When a study needs a code, vocabulary or docs change**, for example a new
quality flag or a parser option, put that change on its own branch and open its
own PR. Merge that PR first, then rebase the study branch onto `main`. List the
waiting study branches in the code PR so whoever merges it knows to rebase them.

## Claiming a study

From an up-to-date checkout on the machine that will do the run:

```bash
aree start-study --accession PRJNA690951 --study-id CIBNOR2021_HEAT_RRSS
```

This:

1. fetches `origin/main` and refuses the claim if the `study_id` or the
   accession is already registered there, or if any pushed `study/` branch
   already uses either;
2. creates the branch in its own git worktree, `../AREE-<STUDY_ID>` by default
   (`--dir` to change it), so one machine can hold several jobs, and your main
   checkout is untouched;
3. writes a stub study YAML (with `study_id`, the accession, `registered_by`
   and the date filled in) and a stub `config/<STUDY_ID>.config`;
4. links this checkout's `data/raw` and `data/reference/GCF_*` drive symlinks
   into the worktree, so the run reads from and writes to the same drives;
5. commits, **pushes the branch (this is the claim)**, and opens a draft PR
   titled `[study] <STUDY_ID> (<ACCESSION>)`. The PR body names the machine and
   includes the study checklist.

If two machines claim the same study at the same moment, both may pass step 1,
but only one push can create the branch. The other is rejected, and
`start-study` stops before anything is downloaded and prints the two commands
that remove the local worktree and branch.

`--no-push` builds the worktree and stubs without claiming anything (useful for
a dry look). `--no-pr` claims without opening a PR. Without the `gh` CLI, no
PR is opened; open one from the pushed branch.

## Seeing what is claimed

```bash
git ls-remote --heads origin 'refs/heads/study/*'
gh pr list --search "head:study/" --state open
```

Every open draft PR is a study in progress. The PR body says which machine has it.

## Doing the work

In the worktree, follow [adding_a_study.md](adding_a_study.md) and, for raw
data, [first_raw_reanalysis.md](first_raw_reanalysis.md):

```bash
cd ../AREE-CIBNOR2021_HEAT_RRSS
aree fetch-samplesheet --bioproject PRJNA690951 --study CIBNOR2021_HEAT_RRSS ...
# stage and MD5-check FASTQs onto this machine's data drive, then
nextflow run workflows/<assay> -config config/CIBNOR2021_HEAT_RRSS.config -profile local,workstation
```

Commit as you go and push. The draft PR's CI shows whether the branch is still
valid. When the run and curation are done, work through the checklist in the
PR body and mark the PR ready for review.

## Setting up an additional machine

Each machine needs:

- a clone of the repository and `pip install -e ".[dev,intake]"`;
- Nextflow, and the container engine the profile uses;
- the reference bundle (`data/reference/GCF_963853765.1/`), staged once per
  machine, as a directory or a symlink to a data drive. It is not in git
  (see [handling_genome_versions.md](handling_genome_versions.md));
- a `data/raw` directory or symlink on a drive with room for the study's FASTQs
  (the FASTQ manifest gives the total size);
- `export AREE_WORKDIR=/path/on/fast/local/disk` so `-profile workstation`
  puts the Nextflow work directory there instead of the default path,
  which only exists on the original workstation;
- `gh auth login`, if `start-study` should open the draft PR.

Avoid exFAT for any of these (see the storage notes in
[first_raw_reanalysis.md](first_raw_reanalysis.md#storage-layout)).

Results are comparable across machines because each run's workflow manifest
records tool versions, parameters and input checksums. Check them in review: a
study run with a different container or reference than the rest of the pool
should say why in its YAML.

## After merging

On `main`, rebuild what the merged study feeds:

```bash
aree build-registry
aree harmonize --study <STUDY_ID> --comparison <id> --input data/studies/<STUDY_ID>/...
aree meta-analyze --feature-type gene
aree build-evidence-cards
```

Update shared pages such as [candidate_studies.md](candidate_studies.md)
and [implementation_status.md](implementation_status.md) in a small PR on
`main`, not on the study branch. Then remove the worktree:

```bash
git worktree remove ../AREE-<STUDY_ID>
```
