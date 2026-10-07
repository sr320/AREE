"""Keep a study PR inside the files its study owns.

Study branches (`study/<ACCESSION>/<STUDY_ID>`, see `intake.study_claims`) are
developed on different machines at the same time. They can only merge in any
order without conflicts if none of them edits a shared file, so a study PR
may touch only:

    registry/studies/<STUDY_ID>.yaml
    config/<STUDY_ID>.config
    data/studies/<STUDY_ID>/**
    tests/test_<study_id>.py  or  tests/test_<first token of study_id>_study.py

A code fix, vocabulary term, or docs change the study needs goes in its own
PR, which merges first; the study branch then rebases onto main.
"""
from __future__ import annotations

from intake.study_claims import accessions_in, parse_branch


def owned_paths(study_id: str) -> tuple[set[str], str]:
    """(exact files, directory prefix) a study branch may change."""
    first = study_id.split("_")[0].lower()
    files = {
        f"registry/studies/{study_id}.yaml",
        f"config/{study_id}.config",
        f"tests/test_{study_id.lower()}.py",
        f"tests/test_{first}_study.py",
    }
    return files, f"data/studies/{study_id}/"


def check_study_scope(branch: str, changed: list[str], study: dict | None,
                      registered: dict[str, set[str]]) -> tuple[list[str], list[str]]:
    """(errors, warnings) for a study PR.

    `changed` are paths changed relative to the base branch, `study` is the
    branch's parsed study YAML (None if it is missing), and `registered` maps
    study_id -> accessions for every study on the base branch.
    """
    parsed = parse_branch(branch)
    if parsed is None:
        return [f"{branch!r} is not a study branch (expected study/<ACCESSION>/<STUDY_ID>)."], []
    accession, study_id = parsed
    errors, warnings = [], []

    files, prefix = owned_paths(study_id)
    outside = sorted(p for p in changed if p not in files and not p.startswith(prefix))
    if outside:
        errors.append(
            f"changes files outside {study_id}'s own paths; move these to a separate PR:\n    "
            + "\n    ".join(outside))

    if study is None:
        errors.append(f"registry/studies/{study_id}.yaml is missing.")
        return errors, warnings
    if study.get("study_id") != study_id:
        errors.append(f"branch names study {study_id} but the YAML has study_id {study.get('study_id')!r}.")
    if accession not in accessions_in(study):
        errors.append(f"branch accession {accession} does not appear in the YAML's `accessions` block.")

    for sid, accs in sorted(registered.items()):
        if sid != study_id and accession in accs:
            warnings.append(f"{accession} is also an accession of {sid} on main. That is fine only "
                            "if one deposit genuinely backs two studies; say so in the YAML.")
    return errors, warnings
