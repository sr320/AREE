"""Branch-per-study claims, the study-PR scope check, and the generated registry.

The end-to-end claim tests run real git against a bare repository in a temp
directory standing in for GitHub, so they need no network.
"""
from __future__ import annotations

import shutil
import subprocess

import pytest
import yaml

from common import REPO_ROOT, STUDIES_DIR

TEMPLATE = (STUDIES_DIR / "_TEMPLATE.yaml").read_text()


# --- naming and conflicts ------------------------------------------------------

def test_branch_name_round_trips():
    from intake.study_claims import branch_name, parse_branch

    b = branch_name("PRJNA690951", "CIBNOR2021_HEAT_RRSS")
    assert b == "study/PRJNA690951/CIBNOR2021_HEAT_RRSS"
    assert parse_branch(b) == ("PRJNA690951", "CIBNOR2021_HEAT_RRSS")
    assert parse_branch("intake/run-label-parsing") is None
    assert parse_branch("study/PRJNA1") is None


def test_branch_name_rejects_doi_and_bad_study_id():
    from intake.study_claims import ClaimError, branch_name

    with pytest.raises(ClaimError, match="DOI"):
        branch_name("10.1038/s41437-020-0351-7", "WANG2021")
    with pytest.raises(ClaimError, match="study_id"):
        branch_name("PRJNA1", "bad id")


@pytest.mark.parametrize("acc,field", [
    ("PRJNA690951", "bioproject"), ("PRJEB81880", "bioproject"), ("GSE12345", "geo"),
    ("PXD002316", "proteomexchange"), ("SRP000001", "sra"), ("ERP000001", "ena"),
    ("iProX123", "other"),
])
def test_accession_field_inferred_from_prefix(acc, field):
    from intake.study_claims import infer_accession_field

    assert infer_accession_field(acc) == field


def test_accessions_in_tokenizes_free_text():
    from intake.study_claims import accessions_in

    study = yaml.safe_load((STUDIES_DIR / "CALLA2026_OSHV.yaml").read_text())
    accs = accessions_in(study)
    assert any(a.startswith("PRJ") for a in accs)
    assert accessions_in({"accessions": {"other": "BioSample SAMN17269691 (single)", "geo": None}}) >= {
        "SAMN17269691", "BioSample"}


def test_find_conflicts_reports_registered_and_claimed():
    from intake.study_claims import find_conflicts

    registered = {"OLD_STUDY": {"PRJNA1"}}
    branches = ["study/PRJNA2/OTHER", "intake/x", "study/PRJNA9/NEW_STUDY"]
    assert find_conflicts("PRJNA3", "FRESH", registered, branches) == []
    assert "already registered on main as OLD_STUDY" in find_conflicts("PRJNA1", "FRESH", registered, [])[0]
    assert "already registered on main" in find_conflicts("PRJNA3", "OLD_STUDY", registered, [])[0]
    assert find_conflicts("PRJNA2", "FRESH", registered, branches) == [
        "already claimed by branch study/PRJNA2/OTHER."]
    assert find_conflicts("PRJNA3", "NEW_STUDY", registered, branches) == [
        "already claimed by branch study/PRJNA9/NEW_STUDY."]


def test_render_study_yaml_fills_identity_and_provenance_only():
    from intake.study_claims import render_study_yaml

    text = render_study_yaml(TEMPLATE, study_id="NEW_STUDY", accession="PXD002316",
                             registered_by="sr320", today="2026-10-07")
    study = yaml.safe_load(text)
    assert study["study_id"] == "NEW_STUDY"
    assert study["accessions"]["proteomexchange"] == "PXD002316"
    assert study["accessions"]["bioproject"] is None
    assert study["provenance"] == {**yaml.safe_load(TEMPLATE)["provenance"],
                                   "registered_by": "sr320", "date_registered": "2026-10-07"}
    # The curator still has to fill the science; the stub must not validate.
    assert study["comparisons"][0]["sample_size"] == 0


def test_render_study_yaml_fails_loudly_on_template_drift():
    from intake.study_claims import ClaimError, render_study_yaml

    with pytest.raises(ClaimError, match="study_id"):
        render_study_yaml("id: x\n", study_id="A", accession="PRJNA1", registered_by="a", today="b")


# --- scope check -----------------------------------------------------------------

STUDY = {"study_id": "NEW_STUDY", "accessions": {"bioproject": "PRJNA1"}}
OWNED = ["registry/studies/NEW_STUDY.yaml", "config/NEW_STUDY.config",
         "data/studies/NEW_STUDY/samplesheet.tsv", "tests/test_new_study.py"]


def test_scope_accepts_owned_paths_only():
    from validation.study_scope import check_study_scope

    errors, warnings = check_study_scope("study/PRJNA1/NEW_STUDY", OWNED, STUDY, {})
    assert errors == warnings == []
    errors, _ = check_study_scope("study/PRJNA1/NEW_STUDY", ["tests/test_new_study.py"], STUDY, {})
    assert errors == []


def test_scope_rejects_shared_and_other_study_files():
    from validation.study_scope import check_study_scope

    changed = [*OWNED, "src/aree/cli.py", "docs/candidate_studies.md",
               "data/studies/NEW_STUDY_2/x.tsv", "registry/study_registry.csv"]
    errors, _ = check_study_scope("study/PRJNA1/NEW_STUDY", changed, STUDY, {})
    assert len(errors) == 1
    for p in changed[len(OWNED):]:
        assert p in errors[0]


def test_scope_checks_yaml_identity_and_accession():
    from validation.study_scope import check_study_scope

    errs, _ = check_study_scope("study/PRJNA1/NEW_STUDY", OWNED, {**STUDY, "study_id": "OTHER"}, {})
    assert any("study_id 'OTHER'" in e for e in errs)
    errs, _ = check_study_scope("study/PRJNA2/NEW_STUDY", OWNED, STUDY, {})
    assert any("PRJNA2 does not appear" in e for e in errs)
    errs, _ = check_study_scope("study/PRJNA1/NEW_STUDY", OWNED, None, {})
    assert any("missing" in e for e in errs)
    errs, _ = check_study_scope("main", [], STUDY, {})
    assert "not a study branch" in errs[0]


def test_scope_warns_on_accession_shared_with_registered_study():
    from validation.study_scope import check_study_scope

    errors, warnings = check_study_scope("study/PRJNA1/NEW_STUDY", OWNED, STUDY, {"OLD": {"PRJNA1"}})
    assert errors == []
    assert "also an accession of OLD" in warnings[0]


def test_existing_real_study_files_fit_their_own_scope():
    """Every real study's committed files must already sit in the paths a study branch owns."""
    from validation.study_scope import owned_paths

    for sid in ("CALLA2026_OSHV", "DELISLE2020_OSHV_TEMP", "IOCAS2022_OA_ENERGY"):
        files, prefix = owned_paths(sid)
        assert (REPO_ROOT / f"registry/studies/{sid}.yaml").exists()
        assert (REPO_ROOT / prefix).is_dir()
        assert any((REPO_ROOT / f).exists() for f in files if f.startswith("tests/"))


# --- generated registry ------------------------------------------------------------

def test_rebuild_registry_indexes_valid_studies_and_reports_invalid(tmp_path, isolated_registry):
    from intake.registry import list_studies, rebuild_registry

    studies = tmp_path / "studies"
    studies.mkdir()
    shutil.copy(STUDIES_DIR / "GIGAS_HEAT01.yaml", studies)
    shutil.copy(STUDIES_DIR / "_TEMPLATE.yaml", studies)
    (studies / "BROKEN.yaml").write_text("study_id: BROKEN\n")

    rows, invalid = rebuild_registry(studies)
    assert [r["study_id"] for r in rows] == ["GIGAS_HEAT01"]
    assert list(invalid) == [str(studies / "BROKEN.yaml")]
    assert [r["study_id"] for r in list_studies()] == ["GIGAS_HEAT01"]

    rebuild_registry(studies)  # replaces rather than appends
    assert len(list_studies()) == 1


def test_registry_csv_is_not_tracked():
    out = subprocess.run(["git", "ls-files", "registry/study_registry.csv"], cwd=REPO_ROOT,
                         capture_output=True, text=True)
    if out.returncode != 0:
        pytest.skip("not a git checkout")
    assert out.stdout.strip() == ""


# --- end to end against a local "GitHub" --------------------------------------------

def _git(*args, cwd):
    subprocess.run(["git", *args], cwd=cwd, check=True, capture_output=True, text=True)


@pytest.fixture
def remote_and_clones(tmp_path, monkeypatch):
    for k, v in {"GIT_AUTHOR_NAME": "t", "GIT_AUTHOR_EMAIL": "t@example.org",
                 "GIT_COMMITTER_NAME": "t", "GIT_COMMITTER_EMAIL": "t@example.org"}.items():
        monkeypatch.setenv(k, v)
    seed = tmp_path / "seed"
    (seed / "registry/studies").mkdir(parents=True)
    (seed / "config").mkdir()
    (seed / ".github/PULL_REQUEST_TEMPLATE").mkdir(parents=True)
    (seed / "registry/studies/_TEMPLATE.yaml").write_text(TEMPLATE)
    (seed / "registry/studies/OLD_STUDY.yaml").write_text(
        "study_id: OLD_STUDY\naccessions:\n  bioproject: PRJNA1\n")
    (seed / ".github/PULL_REQUEST_TEMPLATE/study.md").write_text("checklist\n")
    _git("init", "-q", "-b", "main", cwd=seed)
    _git("add", ".", cwd=seed)
    _git("commit", "-q", "-m", "seed", cwd=seed)

    remote = tmp_path / "remote.git"
    _git("clone", "-q", "--bare", str(seed), str(remote), cwd=tmp_path)
    clones = []
    for name in ("machine_a", "machine_b"):
        _git("clone", "-q", str(remote), name, cwd=tmp_path)
        clones.append(tmp_path / name)
    return remote, clones


def _remote_branches(remote):
    out = subprocess.run(["git", "branch", "--list"], cwd=remote, capture_output=True, text=True)
    return {b.strip("* ").strip() for b in out.stdout.splitlines()}


def test_start_study_claims_branch_and_blocks_second_machine(remote_and_clones, tmp_path):
    from intake.study_claims import ClaimError, start_study

    remote, (a, b) = remote_and_clones
    drive = tmp_path / "drive_raw"
    drive.mkdir()
    (a / "data").mkdir()
    (a / "data/raw").symlink_to(drive)
    res = start_study("PRJNA2", "NEW_STUDY", repo=a, worktree_dir=tmp_path / "wt_a",
                      registered_by="curator", open_pr=False)
    assert res["linked"] == ["data/raw"]
    assert (tmp_path / "wt_a/data/raw").resolve() == drive.resolve()
    assert res["pushed"] and res["branch"] == "study/PRJNA2/NEW_STUDY"
    assert "study/PRJNA2/NEW_STUDY" in _remote_branches(remote)

    study = yaml.safe_load((tmp_path / "wt_a/registry/studies/NEW_STUDY.yaml").read_text())
    assert study["accessions"]["bioproject"] == "PRJNA2"
    assert study["provenance"]["registered_by"] == "curator"
    assert (tmp_path / "wt_a/config/NEW_STUDY.config").exists()

    # Same accession under another name, from another machine.
    with pytest.raises(ClaimError, match="already claimed by branch study/PRJNA2/NEW_STUDY"):
        start_study("PRJNA2", "RENAMED", repo=b, worktree_dir=tmp_path / "wt_b",
                    registered_by="other", open_pr=False)
    # Accession already registered on main.
    with pytest.raises(ClaimError, match="registered on main as OLD_STUDY"):
        start_study("PRJNA1", "AGAIN", repo=b, worktree_dir=tmp_path / "wt_c",
                    registered_by="other", open_pr=False)
    assert not (tmp_path / "wt_b").exists()


def test_simultaneous_claim_loses_at_push(remote_and_clones, tmp_path, monkeypatch):
    """Machine B checks before A's push lands, so the pre-check passes; git rejects B's push."""
    import intake.study_claims as claims

    remote, (a, b) = remote_and_clones
    claims.start_study("PRJNA3", "RACE", repo=a, worktree_dir=tmp_path / "wt_a",
                       registered_by="a", open_pr=False)
    monkeypatch.setattr(claims, "remote_study_branches", lambda *args: [])
    monkeypatch.setattr(claims, "local_study_branches", lambda *args: [])
    with pytest.raises(claims.ClaimError, match="another machine claimed study/PRJNA3/RACE first"):
        claims.start_study("PRJNA3", "RACE", repo=b, worktree_dir=tmp_path / "wt_b",
                           registered_by="b", open_pr=False)
    # A's claim is untouched on the remote.
    log = subprocess.run(["git", "log", "-1", "--format=%s", "study/PRJNA3/RACE"], cwd=remote,
                         capture_output=True, text=True).stdout.strip()
    assert log == "Claim RACE (PRJNA3)"


def test_unpushed_local_claim_blocks_the_same_machine(remote_and_clones, tmp_path):
    from intake.study_claims import ClaimError, start_study

    _, (a, _b) = remote_and_clones
    start_study("PRJNA4", "LOCAL", repo=a, worktree_dir=tmp_path / "wt_1",
                registered_by="a", push=False, open_pr=False)
    with pytest.raises(ClaimError, match="already claimed by branch study/PRJNA4/LOCAL"):
        start_study("PRJNA4", "LOCAL", repo=a, worktree_dir=tmp_path / "wt_2",
                    registered_by="a", push=False, open_pr=False)
