"""Run a declarative intake config: published supplementary table -> AREE result files.

`supplementary_table.py` holds the mechanical reshaping. This module makes that
reshaping *reproducible*: the source file, its checksum, the per-sheet column
maps, and the curator's transformation notes all live in a committed YAML
config, so the step from a published artifact to the TSVs AREE harmonizes is a
command anyone can re-run rather than an undocumented one-off.

Two modes:

* **write** (default) — regenerate the result files and the intake provenance.
* **check** (``--check``) — regenerate into a temporary directory and compare
  checksums against the committed provenance, without touching the tree. This
  is what CI runs: it proves the committed derived files still follow from the
  committed source, and fails loudly if either has drifted.
"""
from __future__ import annotations

import tempfile
from pathlib import Path

import pandas as pd

from common import REPO_ROOT, load_json, load_yaml, sha256_file

from .supplementary_table import RNASEQ_COLUMNS, convert_de_table, write_intake_provenance

SPREADSHEET_SUFFIXES = {".xls", ".xlsx", ".xlsm"}


class IntakeError(RuntimeError):
    """A config, source file, or checksum problem that must stop the intake."""


def _resolve(path_str: str) -> Path:
    path = Path(path_str)
    return path if path.is_absolute() else REPO_ROOT / path


def _display(path: Path) -> str:
    """Repo-relative where possible, absolute otherwise (e.g. a temp dir in tests)."""
    return str(path.relative_to(REPO_ROOT)) if path.is_relative_to(REPO_ROOT) else str(path)


def _read_source_table(source_path: Path, sheet: str | None) -> pd.DataFrame:
    """Read one table out of the published artifact.

    Spreadsheet engines (xlrd for .xls, openpyxl for .xlsx) are optional
    dependencies, so a missing engine reports how to install it rather than
    surfacing pandas' generic ImportError.
    """
    suffix = source_path.suffix.lower()
    if suffix in SPREADSHEET_SUFFIXES:
        if sheet is None:
            raise IntakeError(
                f"{source_path.name} is a spreadsheet, so each conversion must name a "
                f"`source_sheet`."
            )
        try:
            return pd.read_excel(source_path, sheet_name=sheet)
        except ImportError as exc:
            raise IntakeError(
                f"Reading {suffix} files needs an optional spreadsheet engine. "
                f'Install it with: pip install -e ".[intake]"  ({exc})'
            ) from exc
        except ValueError as exc:
            raise IntakeError(f"Sheet {sheet!r} not found in {source_path.name}: {exc}") from exc
    if suffix in {".tsv", ".txt"}:
        return pd.read_csv(source_path, sep="\t")
    if suffix == ".csv":
        return pd.read_csv(source_path)
    raise IntakeError(f"Unsupported source format {suffix!r} for {source_path.name}")


def _verify_source_checksum(source_path: Path, declared: str | None) -> str:
    actual = sha256_file(source_path)
    if declared and declared != actual:
        raise IntakeError(
            f"Source checksum mismatch for {source_path.name}.\n"
            f"  config declares: {declared}\n"
            f"  file on disk:    {actual}\n"
            "The published artifact has changed, or the wrong file is staged. "
            "Resolve this deliberately — do not update the config to match without "
            "confirming which version the study was curated from."
        )
    return actual


def _validate_config(cfg: dict, config_path: Path) -> None:
    """Fail loudly on a malformed config rather than writing degraded provenance."""
    for required in ("study_id", "source", "conversions"):
        if required not in cfg:
            raise IntakeError(f"{config_path.name}: missing required key {required!r}")
    if "local_copy" not in cfg["source"]:
        raise IntakeError(f"{config_path.name}: source is missing required key 'local_copy'")
    if not cfg["conversions"]:
        raise IntakeError(f"{config_path.name}: `conversions` is empty")

    for i, conv in enumerate(cfg["conversions"]):
        for required in ("output_file", "column_map"):
            if required not in conv:
                raise IntakeError(f"{config_path.name}: conversion {i} is missing {required!r}")
        unknown = set(conv["column_map"]) - set(RNASEQ_COLUMNS)
        if unknown:
            raise IntakeError(
                f"{config_path.name}: conversion {i} maps unknown AREE column(s) "
                f"{sorted(unknown)}; expected a subset of {RNASEQ_COLUMNS}"
            )
        if "gene_id" not in conv["column_map"]:
            raise IntakeError(
                f"{config_path.name}: conversion {i} must map 'gene_id' — without an "
                "identifier the rows cannot be harmonized"
            )

    # An unquoted note containing ': ' parses as a YAML mapping, which would land
    # in the provenance as a dict and quietly corrupt the curator's caveats.
    for note in cfg.get("transformation_notes", []):
        if not isinstance(note, str):
            raise IntakeError(
                f"{config_path.name}: transformation_notes must be strings, got "
                f"{type(note).__name__}: {note!r}. A note containing ': ' needs quoting."
            )


def run_intake(config_path: Path, *, check: bool = False) -> dict:
    """Execute one intake config.

    Returns a report dict. In `check` mode nothing on disk is modified and the
    report carries a `mismatches` list; in write mode the result files and the
    intake provenance JSON are (re)written.
    """
    config_path = Path(config_path)
    if not config_path.exists():
        raise IntakeError(f"No intake config at {config_path}")
    cfg = load_yaml(config_path)

    _validate_config(cfg, config_path)

    study_id = cfg["study_id"]
    source = cfg["source"]
    source_path = _resolve(source["local_copy"])
    if not source_path.exists():
        raise IntakeError(
            f"Source artifact not found: {source_path}\n"
            f"Download it from {source.get('url', '(no url in config)')} and place it there."
        )
    _verify_source_checksum(source_path, source.get("sha256"))

    output_dir = _resolve(cfg.get("output_dir", source_path.parent.parent))
    provenance_path = _resolve(
        cfg.get("provenance_file", str(Path(output_dir) / "intake_provenance.json"))
    )

    with tempfile.TemporaryDirectory() as tmp:
        write_dir = Path(tmp) if check else output_dir
        conversions = []
        for conv in cfg["conversions"]:
            frame = _read_source_table(source_path, conv.get("source_sheet"))
            out_name = conv["output_file"]
            report = convert_de_table(
                frame,
                conv["column_map"],
                out_path=write_dir / Path(out_name).name,
            )
            # In check mode the temp path is meaningless; always record the
            # config-declared destination so provenance is stable across modes.
            report["output_file"] = _display(output_dir / Path(out_name).name)
            if "source_sheet" in conv:
                report["source_sheet"] = conv["source_sheet"]
            conversions.append(report)

        if check:
            return {
                "study_id": study_id,
                "mode": "check",
                "conversions": conversions,
                "mismatches": _compare_to_committed(provenance_path, conversions, output_dir),
            }

    # A re-run that reproduces byte-identical outputs keeps the original derivation
    # date, so verifying reproducibility never shows up as a spurious diff.
    unchanged = not _compare_to_committed(provenance_path, conversions, output_dir)
    prior_date = load_json(provenance_path).get("date_generated") if provenance_path.exists() else None

    write_intake_provenance(
        provenance_path,
        study_id=study_id,
        source_description=source.get("description", ""),
        source_url=source.get("url", ""),
        source_license=source.get("license", ""),
        source_file=source_path,
        citation=source.get("citation", ""),
        conversions=conversions,
        notes=cfg.get("transformation_notes", []),
        date_generated=prior_date if unchanged else None,
    )
    return {
        "study_id": study_id,
        "mode": "write",
        "conversions": conversions,
        "provenance_file": _display(provenance_path),
        "mismatches": [],
    }


def _compare_to_committed(
    provenance_path: Path, conversions: list, output_dir: Path
) -> list[str]:
    """Compare freshly derived outputs against what is committed to the repo.

    Checks the regenerated checksum against both the committed provenance record
    and the committed file itself, so that a hand-edited TSV is caught even if
    the provenance JSON was updated to match it.
    """
    mismatches = []
    committed = {}
    if provenance_path.exists():
        for entry in load_json(provenance_path).get("conversions", []):
            committed[Path(entry["output_file"]).name] = entry
    else:
        mismatches.append(f"no committed provenance at {provenance_path}")

    for report in conversions:
        name = Path(report["output_file"]).name
        record = committed.get(name)
        if record is None:
            if committed:
                mismatches.append(f"{name}: not recorded in committed provenance")
            continue
        if record.get("output_sha256") != report["output_sha256"]:
            mismatches.append(
                f"{name}: regenerated checksum {report['output_sha256'][:12]}… "
                f"!= provenance {str(record.get('output_sha256'))[:12]}…"
            )
        if record.get("rows_written") != report["rows_written"]:
            mismatches.append(
                f"{name}: regenerated {report['rows_written']} rows, "
                f"provenance records {record.get('rows_written')}"
            )

        on_disk = output_dir / name
        if not on_disk.exists():
            mismatches.append(f"{name}: committed result file is missing from {output_dir}")
        elif sha256_file(on_disk) != report["output_sha256"]:
            mismatches.append(
                f"{name}: the committed file on disk does not match what the config regenerates"
            )

    return mismatches
