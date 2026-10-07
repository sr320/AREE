"""Claim one study for curation or reanalysis: one branch per accession.

Several machines can work through the study backlog at once. Each job is one
study on its own branch, named

    study/<ACCESSION>/<STUDY_ID>       e.g. study/PRJNA690951/CIBNOR2021_HEAT_RRSS

and opened as a draft pull request against main. Pushing the branch is the
claim. Two machines that race for the same study push different commits to
the same branch name, so git rejects the second push; the loser stops before
any data is downloaded. Before pushing, `start_study` also refuses an
accession or study_id that is already registered on main or already claimed
by another `study/` branch.

A study branch only touches files the study owns (see
`validation.study_scope`), so study PRs merge in any order without conflicts.
Shared code, vocabularies and docs change in their own PRs.
"""
from __future__ import annotations

import re
import shutil
import socket
import subprocess
from datetime import date
from pathlib import Path

import yaml

BRANCH_PREFIX = "study/"
ACCESSION_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")
STUDY_ID_RE = re.compile(r"^[A-Za-z0-9_-]+$")
_BRANCH_RE = re.compile(r"^study/(?P<accession>[^/]+)/(?P<study_id>[^/]+)$")

# Prefix -> field under `accessions:` in the study YAML. First match wins.
ACCESSION_FIELDS = [
    ("PRJ", "bioproject"),
    ("GSE", "geo"),
    ("PXD", "proteomexchange"),
    ("SRP", "sra"),
    ("ERP", "ena"),
    ("DRP", "sra"),
]

STUDY_PR_TEMPLATE = Path(".github/PULL_REQUEST_TEMPLATE/study.md")

# Gitignored symlinks that point a checkout at this machine's data drives.
# A new worktree does not get them from git, so they are copied across.
DRIVE_LINK_GLOBS = ["data/raw", "data/reference/GCF_*"]


class ClaimError(RuntimeError):
    pass


def infer_accession_field(accession: str) -> str:
    for prefix, field in ACCESSION_FIELDS:
        if accession.upper().startswith(prefix):
            return field
    return "other"


def branch_name(accession: str, study_id: str) -> str:
    if not ACCESSION_RE.match(accession):
        raise ClaimError(
            f"accession {accession!r} must be a repository accession (letters, digits, '.', '_', '-'). "
            "A DOI cannot name a branch; claim the study under its data accession instead."
        )
    if not STUDY_ID_RE.match(study_id):
        raise ClaimError(f"study_id {study_id!r} must match {STUDY_ID_RE.pattern}.")
    return f"{BRANCH_PREFIX}{accession}/{study_id}"


def parse_branch(name: str) -> tuple[str, str] | None:
    """(accession, study_id) for a study branch, or None for any other branch."""
    m = _BRANCH_RE.match(name)
    return (m["accession"], m["study_id"]) if m else None


def accessions_in(study: dict) -> set[str]:
    """Every accession-like token in a study's `accessions` block.

    Values are free text in places (e.g. "BioSample SAMN17269691 (...)"), so
    each value is split into tokens as well as kept whole.
    """
    out = set()
    for value in (study.get("accessions") or {}).values():
        if not value:
            continue
        out.add(str(value))
        out.update(t for t in re.split(r"[^A-Za-z0-9._]+", str(value)) if t)
    return out


def find_conflicts(accession: str, study_id: str, registered: dict[str, set[str]],
                   claimed_branches: list[str]) -> list[str]:
    """Reasons this (accession, study_id) cannot be claimed. Empty means free.

    `registered` maps each study_id on main to its accessions; `claimed_branches`
    are the remote branch names currently pushed.
    """
    problems = []
    if study_id in registered:
        problems.append(f"study_id {study_id} is already registered on main.")
    for sid, accs in sorted(registered.items()):
        if sid != study_id and accession in accs:
            problems.append(f"accession {accession} is already registered on main as {sid}.")
    for branch in sorted(claimed_branches):
        parsed = parse_branch(branch)
        if parsed and (parsed[0] == accession or parsed[1] == study_id):
            problems.append(f"already claimed by branch {branch}.")
    return problems


def _sub_once(text: str, pattern: str, repl: str, what: str) -> str:
    new, n = re.subn(pattern, repl, text, count=1, flags=re.MULTILINE)
    if n != 1:
        raise ClaimError(f"registry/studies/_TEMPLATE.yaml has no {what} line to fill; "
                         "update intake/study_claims.py to match the template.")
    return new


def render_study_yaml(template: str, *, study_id: str, accession: str, registered_by: str,
                      today: str) -> str:
    """The template with identity and provenance filled; everything else left to the curator."""
    field = infer_accession_field(accession)
    text = _sub_once(template, r"^study_id: .*$", f"study_id: {study_id}", "study_id")
    text = _sub_once(text, rf"^(  {field}:) null\s*$", rf'\1 "{accession}"', f"accessions.{field}")
    text = _sub_once(text, r"^(  registered_by:) .*$", rf'\1 "{registered_by}"', "registered_by")
    text = _sub_once(text, r"^(  date_registered:) .*$", rf'\1 "{today}"', "date_registered")
    return text


def render_config(study_id: str, accession: str) -> str:
    return f"""/*
 * AREE :: run configuration for {study_id} ({accession})
 *
 * Claimed with `aree start-study`. Fill in params for the workflow this
 * study runs (see config/IOCAS2022_OA_ENERGY.config for RNA-seq and
 * config/CIBNOR2021_HEAT_RRSS.config for methylation), then
 *
 *   nextflow run workflows/<assay> -config config/{study_id}.config -profile local,workstation
 *
 * Generate the sample sheet and FASTQ manifest rather than writing them:
 *
 *   aree fetch-samplesheet --bioproject {accession} --study {study_id} --condition-attribute <attr>
 */

includeConfig 'base.config'

params {{
    mode = 'raw_reanalysis'

    sample_sheet = 'data/studies/{study_id}/samplesheet.tsv'
    reads        = 'data/raw/{study_id}/*_{{1,2}}.fastq.gz'
    outdir       = 'results/{study_id}'
}}
"""


# --- git / GitHub ------------------------------------------------------------

def _run(args: list[str], cwd) -> str:
    proc = subprocess.run(args, cwd=cwd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise ClaimError(f"`{' '.join(args)}` failed:\n{proc.stderr.strip() or proc.stdout.strip()}")
    return proc.stdout


def registered_on(ref: str, repo) -> dict[str, set[str]]:
    """study_id -> accessions for every study YAML at `ref` (e.g. origin/main)."""
    names = _run(["git", "ls-tree", "--name-only", ref, "registry/studies/"], repo).split()
    out = {}
    for name in names:
        if not name.endswith(".yaml") or Path(name).name.startswith("_"):
            continue
        study = yaml.safe_load(_run(["git", "show", f"{ref}:{name}"], repo)) or {}
        out[study.get("study_id") or Path(name).stem] = accessions_in(study)
    return out


def remote_study_branches(remote: str, repo) -> list[str]:
    out = _run(["git", "ls-remote", "--heads", remote, f"refs/heads/{BRANCH_PREFIX}*"], repo)
    return [line.split("refs/heads/", 1)[1] for line in out.splitlines() if "refs/heads/" in line]


def link_drives(repo: Path, worktree: Path) -> list[str]:
    """Recreate the repo's data-drive symlinks in a new worktree. Returns the paths linked."""
    linked = []
    for pattern in DRIVE_LINK_GLOBS:
        for src in sorted(repo.glob(pattern)):
            if not src.is_symlink():
                continue
            rel = src.relative_to(repo)
            dest = worktree / rel
            if dest.exists() or dest.is_symlink():
                continue
            dest.parent.mkdir(parents=True, exist_ok=True)
            dest.symlink_to(src.readlink())
            linked.append(str(rel))
    return linked


def local_study_branches(repo) -> list[str]:
    """Study branches on this machine, including unpushed (`--no-push`) ones."""
    out = _run(["git", "branch", "--list", f"{BRANCH_PREFIX}*", "--format=%(refname:short)"], repo)
    return [b.strip() for b in out.splitlines() if b.strip()]


def start_study(accession: str, study_id: str, *, repo, worktree_dir=None, registered_by: str,
                machine: str | None = None, base: str = "main", remote: str = "origin",
                push: bool = True, open_pr: bool = True) -> dict:
    """Create the study branch in its own worktree, commit the stubs, and claim it.

    Returns {branch, worktree, files, linked, pushed, pr_url}.
    """
    repo = Path(repo)
    branch = branch_name(accession, study_id)
    worktree = Path(worktree_dir) if worktree_dir else repo.parent / f"AREE-{study_id}"
    if worktree.exists():
        raise ClaimError(f"{worktree} already exists; pass --dir to put the worktree elsewhere.")

    _run(["git", "fetch", "--quiet", remote, base], repo)
    base_ref = f"{remote}/{base}"
    claimed = sorted(set(remote_study_branches(remote, repo)) | set(local_study_branches(repo)))
    problems = find_conflicts(accession, study_id, registered_on(base_ref, repo), claimed)
    if problems:
        raise ClaimError("Cannot claim this study:\n  - " + "\n  - ".join(problems))

    _run(["git", "worktree", "add", "--quiet", "-b", branch, str(worktree), base_ref], repo)

    today = date.today().isoformat()
    template = (worktree / "registry/studies/_TEMPLATE.yaml").read_text()
    files = {
        f"registry/studies/{study_id}.yaml": render_study_yaml(
            template, study_id=study_id, accession=accession,
            registered_by=registered_by, today=today),
        f"config/{study_id}.config": render_config(study_id, accession),
    }
    for rel, text in files.items():
        (worktree / rel).parent.mkdir(parents=True, exist_ok=True)
        (worktree / rel).write_text(text)
    _run(["git", "add", *files], worktree)
    _run(["git", "commit", "--quiet", "-m", f"Claim {study_id} ({accession})"], worktree)

    result = {"branch": branch, "worktree": str(worktree), "files": sorted(files),
              "linked": link_drives(repo, worktree), "pushed": False, "pr_url": None}
    if not push:
        return result
    try:
        _run(["git", "push", "--quiet", "-u", remote, branch], worktree)
    except ClaimError as exc:
        raise ClaimError(
            f"{exc}\n\nThe push was rejected, most likely because another machine claimed "
            f"{branch} first. Nothing was downloaded. Remove the local claim with:\n"
            f"  git worktree remove {worktree} && git branch -D {branch}"
        ) from None
    result["pushed"] = True

    if open_pr and shutil.which("gh"):
        machine = machine or socket.gethostname().split(".")[0]
        body = (f"Claimed by `aree start-study` on **{machine}**, {today}. "
                f"Draft until the run is finished.\n\n"
                + (worktree / STUDY_PR_TEMPLATE).read_text())
        result["pr_url"] = _run(
            ["gh", "pr", "create", "--draft", "--base", base, "--head", branch,
             "--title", f"[study] {study_id} ({accession})", "--body", body], worktree).strip()
    return result
