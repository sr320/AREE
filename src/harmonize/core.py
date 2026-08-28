"""Dispatch harmonization by assay type and maintain the master evidence table."""
from __future__ import annotations

import ast
import warnings
from pathlib import Path

import pandas as pd

from common import EVIDENCE_TABLE_PATH, REPO_ROOT, REPORTS_DIR, STUDIES_DIR, load_yaml, sha256_file

from .identifiers import (
    active_crosswalk_path,
    crosswalk_for_study,
    retired_table_path,
    set_crosswalk_path,
)
from .metabolomics import harmonize_metabolomics
from .methylation import harmonize_methylation
from .proteomics import harmonize_proteomics
from .rnaseq import harmonize_rnaseq
from .schema import EVIDENCE_COLUMNS

MANIFESTS_DIR = REPORTS_DIR / "manifests"

ASSAY_HARMONIZERS = {
    "rnaseq": harmonize_rnaseq,
    "methylation": harmonize_methylation,
    "proteomics": harmonize_proteomics,
    "metabolomics": harmonize_metabolomics,
}

WORKFLOW_VERSION = "aree-harmonize-0.1.0"


def _generated_by() -> str:
    return f"automated:{WORKFLOW_VERSION}"


def _load_study(study_id: str) -> dict:
    path = STUDIES_DIR / f"{study_id}.yaml"
    if not path.exists():
        raise FileNotFoundError(f"No registered study file at {path}")
    return load_yaml(path)


def _harmonize_comparison(
    study: dict, comparison: dict, date_generated: str, results_override=None
) -> pd.DataFrame:
    """Harmonize one comparison.

    `results_override` is the file the caller actually wants harmonized. It
    takes precedence over the comparison's registered `results_file`, which is
    what makes `aree harmonize --input` mean what it says and is the only way a
    raw_reanalysis study — whose `results_file` is null until a workflow run
    produces one — can be harmonized at all.
    """
    assay_type = study["assay_type"][0]
    harmonizer = ASSAY_HARMONIZERS.get(assay_type)
    if harmonizer is None:
        raise ValueError(f"No harmonizer registered for assay_type {assay_type!r}")

    if results_override is not None:
        results_path = Path(results_override)
    else:
        results_file = comparison.get("results_file")
        if not results_file:
            raise ValueError(
                f"comparison {comparison['comparison_id']!r} in study {study['study_id']!r} "
                "has no results_file to harmonize. For a raw_reanalysis study, pass the "
                "workflow output with --input and name the comparison with --comparison."
            )
        results_path = REPO_ROOT / results_file
    if not results_path.exists():
        raise FileNotFoundError(f"results file not found: {results_path}")

    # The study decides which crosswalk it resolves against, so that simulated and
    # real studies can coexist in one evidence table without cross-contamination.
    set_crosswalk_path(crosswalk_for_study(bool(study.get("simulated", False))))

    evidence_df = harmonizer(
        study, comparison, results_path,
        workflow_version=WORKFLOW_VERSION,
        date_generated=date_generated,
        generated_by=_generated_by(),
    )
    _write_harmonize_manifest(study, comparison, results_path, evidence_df, date_generated)
    return evidence_df


def _repo_relative(path: Path) -> str:
    """Repo-relative where possible, absolute otherwise."""
    return str(path.relative_to(REPO_ROOT)) if path.is_relative_to(REPO_ROOT) else str(path)


def _crosswalk_ref() -> dict:
    """Path + checksum of the identifier crosswalk currently in force."""
    path = active_crosswalk_path()
    ref = {
        "path": _repo_relative(path),
        "sha256": sha256_file(path) if path.exists() else None,
    }
    retired = retired_table_path(path)
    if retired is not None and retired.exists():
        ref["retired_gene_id_table"] = (
            _repo_relative(retired)
        )
        ref["retired_gene_id_table_sha256"] = sha256_file(retired)
    return ref


def _write_harmonize_manifest(study: dict, comparison: dict, results_path, evidence_df: pd.DataFrame, date_generated: str) -> Path:
    MANIFESTS_DIR.mkdir(parents=True, exist_ok=True)
    manifest = {
        "study_id": study["study_id"],
        "comparison_id": comparison["comparison_id"],
        "assay_type": study["assay_type"][0],
        "analysis_mode": study["analysis_mode"],
        # Workflow output legitimately lives outside the repo (results/, or a
        # staging volume), so fall back to the absolute path rather than raising.
        "input_file": _repo_relative(Path(results_path)),
        "input_checksum": f"sha256:{sha256_file(results_path)}",
        "parameters": {
            "genome_assembly": study["genome_assembly"],
            "annotation_version": study.get("annotation_version"),
            # Record the crosswalk that was actually used, with its checksum. A
            # hardcoded path here silently misreported provenance for real studies.
            "identifier_crosswalk": _crosswalk_ref(),
        },
        "software_versions": {"aree-harmonize": WORKFLOW_VERSION},
        "workflow_version": WORKFLOW_VERSION,
        "date_generated": date_generated,
        "generated_by": _generated_by(),
        "n_evidence_records": int(len(evidence_df)),
        "mapping_confidence_counts": evidence_df["mapping_confidence"].value_counts().to_dict() if len(evidence_df) else {},
        "quality_flags_declared_on_study": study.get("quality_flags") or [],
        "warnings": (
            ["raw_data_not_available: study analysis_mode is processed_results_harmonization; "
             "raw QC/normalization were not independently re-verified by AREE."]
            if study["analysis_mode"] == "processed_results_harmonization" else []
        ),
    }
    out_path = MANIFESTS_DIR / f"{study['study_id']}_{comparison['comparison_id']}_manifest.json"
    with open(out_path, "w") as fh:
        import json
        json.dump(manifest, fh, indent=2, sort_keys=True, default=str)
    return out_path


def harmonize_study(study_id: str, date_generated: str) -> pd.DataFrame:
    """Harmonize every comparison in a registered study and upsert into the
    master evidence table (reports/evidence/evidence_table.tsv)."""
    study = _load_study(study_id)
    frames = [
        _harmonize_comparison(study, comparison, date_generated)
        for comparison in study["comparisons"]
    ]
    new_rows = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(columns=EVIDENCE_COLUMNS)
    _upsert_evidence_table(new_rows)
    _warn_if_analysis_status_stale(study)
    return new_rows



def _warn_if_analysis_status_stale(study: dict) -> None:
    """Flag a registered `analysis_status` that harmonization has just outdated.

    `analysis_status` is curator-maintained, so it drifts: HESSER2024_VCOR sat at
    `not_started` through several harmonization runs. This does not rewrite the
    curator's YAML — a registration file is theirs to edit — but it makes the
    staleness visible at the moment it is created.
    """
    status = study.get("analysis_status", "not_started")
    if status in ("not_started", "in_progress"):
        warnings.warn(
            f"Study {study['study_id']} is registered with analysis_status="
            f"{status!r} but has just been harmonized. Update "
            f"registry/studies/{study['study_id']}.yaml to 'harmonized' and re-run "
            "`aree register-study --update`.",
            stacklevel=2,
        )

def harmonize_processed_table(
    study_id: str, input_path, date_generated: str, comparison_id: str | None = None
) -> pd.DataFrame:
    """CLI entry point for `aree harmonize --study STUDY_ID --input path`.

    The comparison is identified either explicitly by `comparison_id`, or —
    when it is omitted — by matching the input filename against a comparison's
    declared `results_file`.

    Filename matching alone is not enough for a raw_reanalysis study: the
    Nextflow workflow names its output for the pipeline stage that produced it
    (`..._dge_standardized.tsv`), which is deliberately not the same as any
    path the curator wrote into the registry, and a raw-mode study has no
    `results_file` to match against at registration time. Passing
    `--comparison` is the supported way to connect a freshly produced workflow
    output to the comparison it belongs to.
    """
    study = _load_study(study_id)
    input_path = Path(input_path)
    known = [c["comparison_id"] for c in study["comparisons"]]

    if comparison_id:
        matches = [c for c in study["comparisons"] if c["comparison_id"] == comparison_id]
        if not matches:
            raise ValueError(
                f"Study {study_id!r} has no comparison {comparison_id!r}. Known: {known}"
            )
    else:
        matches = [
            c for c in study["comparisons"]
            if c.get("results_file") and Path(c["results_file"]).name == input_path.name
        ]
        if not matches:
            raise ValueError(
                f"No comparison in study {study_id!r} declares results_file matching "
                f"{input_path.name!r}, so the comparison cannot be inferred from the "
                f"filename. Pass --comparison explicitly. Known comparisons: {known}"
            )

    frames = [
        _harmonize_comparison(study, comparison, date_generated, results_override=input_path)
        for comparison in matches
    ]
    new_rows = pd.concat(frames, ignore_index=True)
    _upsert_evidence_table(new_rows)
    return new_rows


def _read_evidence_table() -> pd.DataFrame:
    if not EVIDENCE_TABLE_PATH.exists():
        return pd.DataFrame(columns=EVIDENCE_COLUMNS)
    df = pd.read_csv(EVIDENCE_TABLE_PATH, sep="\t", dtype=str)
    if "quality_flags" in df.columns:
        df["quality_flags"] = df["quality_flags"].apply(
            lambda v: ast.literal_eval(v) if isinstance(v, str) and v.startswith("[") else []
        )
    for numeric_col in ("effect_size", "standard_error", "ci_lower", "ci_upper", "p_value", "adjusted_p_value", "sample_size"):
        if numeric_col in df.columns:
            df[numeric_col] = pd.to_numeric(df[numeric_col], errors="coerce")
    return df


def _upsert_evidence_table(new_rows: pd.DataFrame) -> None:
    EVIDENCE_TABLE_PATH.parent.mkdir(parents=True, exist_ok=True)
    existing = _read_evidence_table()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", FutureWarning)
        if len(existing) and len(new_rows):
            existing = existing[~existing["evidence_id"].isin(new_rows["evidence_id"])]
            combined = pd.concat([existing, new_rows], ignore_index=True) if len(existing) else new_rows
        elif len(existing):
            combined = existing
        else:
            combined = new_rows
    combined = combined[EVIDENCE_COLUMNS]
    combined.to_csv(EVIDENCE_TABLE_PATH, sep="\t", index=False)


def load_evidence_table() -> pd.DataFrame:
    return _read_evidence_table()
