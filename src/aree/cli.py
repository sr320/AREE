"""AREE command-line interface.

See README.md for the full command reference.
"""
from __future__ import annotations

import datetime
import sys

import click
from tabulate import tabulate

from aree import __version__
from harmonize.core import harmonize_processed_table, harmonize_study
from intake.registry import DuplicateStudyError, list_studies, register_study
from intake.schema_validate import validate_study_file
from meta_analysis.run import write_meta_analysis
from reporting.evidence_cards import build_evidence_cards


def _today() -> str:
    return datetime.date.today().isoformat()


@click.group()
@click.version_option(version=__version__, prog_name="aree")
def main():
    """AREE — Aquaculture Resilience Evidence Engine."""


@main.command("validate-study")
@click.argument("path", type=click.Path(exists=True))
def validate_study_cmd(path):
    """Validate a study registration YAML file against the schema and controlled vocabularies."""
    result = validate_study_file(path)
    if result.warnings:
        click.echo(click.style("Warnings:", fg="yellow"))
        for w in result.warnings:
            click.echo(f"  - {w}")
    if result.valid:
        click.echo(click.style(f"VALID: {path}", fg="green"))
    else:
        click.echo(click.style(f"INVALID: {path}", fg="red"))
        for e in result.errors:
            click.echo(f"  - {e}")
        sys.exit(1)


@main.command("register-study")
@click.argument("path", type=click.Path(exists=True))
@click.option("--update", "allow_update", is_flag=True, help="Overwrite an existing registry entry for this study_id.")
def register_study_cmd(path, allow_update):
    """Validate and add a study to registry/study_registry.csv."""
    try:
        row = register_study(path, allow_update=allow_update)
    except DuplicateStudyError as exc:
        click.echo(click.style(str(exc), fg="red"))
        sys.exit(1)
    except ValueError as exc:
        click.echo(click.style(str(exc), fg="red"))
        sys.exit(1)
    click.echo(click.style(f"Registered {row['study_id']}", fg="green"))


@main.command("list-studies")
def list_studies_cmd():
    """List all studies currently in the registry."""
    rows = list_studies()
    if not rows:
        click.echo("No studies registered yet. Run `aree register-study` first.")
        return
    click.echo(tabulate(
        [[r["study_id"], r["assay_type"], r["analysis_mode"], r["qc_status"], r["analysis_status"]] for r in rows],
        headers=["study_id", "assay_type", "analysis_mode", "qc_status", "analysis_status"],
    ))


@main.command("harmonize")
@click.option("--study", "study_id", required=True, help="Registered study_id.")
@click.option("--input", "input_path", required=False, type=click.Path(exists=True),
              help="Processed results file to harmonize. If omitted, harmonizes every comparison in the study.")
@click.option("--date", "date_generated", default=None, help="Override date_generated (ISO 8601). Defaults to today.")
def harmonize_cmd(study_id, input_path, date_generated):
    """Harmonize a processed study result table (or an entire study) into the shared evidence table."""
    date_generated = date_generated or _today()
    try:
        if input_path:
            df = harmonize_processed_table(study_id, input_path, date_generated=date_generated)
        else:
            df = harmonize_study(study_id, date_generated=date_generated)
    except (FileNotFoundError, ValueError) as exc:
        click.echo(click.style(str(exc), fg="red"))
        sys.exit(1)
    click.echo(click.style(f"Harmonized {len(df)} evidence records for {study_id}", fg="green"))
    click.echo("Evidence table: reports/evidence/evidence_table.tsv")


@main.command("meta-analyze")
@click.option("--phenotype", default=None, help="Phenotype ontology term id to filter on (default: all).")
@click.option("--feature-type", default=None, help="Feature type to filter on, e.g. gene, protein (default: all).")
def meta_analyze_cmd(phenotype, feature_type):
    """Run a random-effects meta-analysis over the harmonized evidence table."""
    try:
        result, out_path = write_meta_analysis(phenotype, feature_type)
    except FileNotFoundError as exc:
        click.echo(click.style(str(exc), fg="red"))
        sys.exit(1)
    if len(result) == 0:
        click.echo("No poolable evidence found for the given filters.")
        return
    click.echo(tabulate(
        result[["feature_id_standardized", "phenotype", "feature_type", "k_studies", "pooled_effect", "p_value", "i_squared"]].head(20),
        headers="keys", showindex=False, floatfmt=".3g",
    ))
    click.echo(f"\nFull results ({len(result)} rows) written to {out_path}")


@main.command("build-evidence-cards")
@click.option("--phenotype", default=None, help="Phenotype ontology term id to filter on (default: all).")
@click.option("--feature-type", default=None, help="Feature type to filter on (default: all).")
def build_evidence_cards_cmd(phenotype, feature_type):
    """Generate one evidence card per candidate biomarker under reports/evidence_cards/."""
    index = build_evidence_cards(phenotype=phenotype, feature_type=feature_type)
    if not index:
        click.echo("No candidates found for the given filters. Run `aree meta-analyze` first if needed.")
        return
    click.echo(click.style(f"Wrote {len(index)} evidence cards to reports/evidence_cards/", fg="green"))
    tiers = {}
    for row in index:
        tiers[row["tier"]] = tiers.get(row["tier"], 0) + 1
    for tier, count in sorted(tiers.items()):
        click.echo(f"  {tier}: {count}")


if __name__ == "__main__":
    main()
