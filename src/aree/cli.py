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
from intake.run_intake import IntakeError, run_intake
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


@main.command("intake-supplementary")
@click.argument("config", type=click.Path(exists=True))
@click.option("--check", is_flag=True,
              help="Regenerate into a temporary directory and verify against the committed "
                   "files and provenance, without modifying the repository.")
def intake_supplementary_cmd(config, check):
    """Convert a published supplementary table into AREE result files.

    CONFIG is an intake YAML (see data/studies/HESSER2024_VCOR/intake.yaml).
    Reshaping is mechanical only: no statistic is ever imputed, and identifiers
    are preserved verbatim for `harmonize` to resolve.
    """
    try:
        report = run_intake(config, check=check)
    except IntakeError as exc:
        click.echo(click.style(str(exc), fg="red"))
        sys.exit(1)

    for conv in report["conversions"]:
        absent = ", ".join(conv["columns_absent_from_source"]) or "none"
        click.echo(
            f"  {conv['output_file']}: {conv['rows_written']} rows "
            f"(dropped {conv['rows_dropped_missing_id_or_effect']} blank, "
            f"{conv['rows_dropped_non_numeric']} non-numeric); "
            f"not reported by source: {absent}"
        )

    if check:
        if report["mismatches"]:
            click.echo(click.style(
                f"\n{report['study_id']}: intake is NOT reproducible from the committed source.",
                fg="red",
            ))
            for m in report["mismatches"]:
                click.echo(click.style(f"  - {m}", fg="red"))
            sys.exit(1)
        click.echo(click.style(
            f"\n{report['study_id']}: committed result files reproduce exactly from the "
            "committed source.", fg="green",
        ))
        return

    click.echo(click.style(f"\nWrote {report['provenance_file']}", fg="green"))


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


@main.command("build-crosswalk")
@click.option("--taxid", default=29159, show_default=True, type=int,
              help="NCBI taxonomy id to build the crosswalk for.")
@click.option("--slug", default="mgigas", show_default=True,
              help="Output filename prefix, e.g. 'mgigas' -> mgigas_gene_id_crosswalk.tsv.")
@click.option("--organism", default="Magallana gigas (Crassostrea gigas), Pacific oyster",
              show_default=False, help="Human-readable organism label recorded in provenance.")
@click.option("--gene-info", "gene_info", default=None, type=click.Path(exists=True),
              help="Pre-downloaded taxid-filtered NCBI gene_info TSV.")
@click.option("--uniprot", default=None, type=click.Path(exists=True),
              help="Pre-downloaded UniProt TSV export.")
@click.option("--download", is_flag=True, help="Force re-download of reference sources.")
def build_crosswalk_cmd(taxid, slug, organism, gene_info, uniprot, download):
    """Build a real identifier crosswalk from NCBI Gene and UniProtKB.

    Downloads reference data from public sources (~230 MB streamed from NCBI) and
    writes a crosswalk plus a JSON provenance sidecar to data/reference/crosswalk/.
    """
    from mappings.build_crosswalk import build

    try:
        path = build(taxid=taxid, organism=organism, gene_info_path=gene_info,
                     uniprot_path=uniprot, download=download, slug=slug)
    except (RuntimeError, ValueError, FileNotFoundError) as exc:
        click.echo(click.style(str(exc), fg="red"))
        sys.exit(1)

    prov = path.with_suffix("").with_suffix(".provenance.json")
    click.echo(click.style(f"Wrote {path}", fg="green"))
    click.echo(f"Provenance: {prov}")
    click.echo("")
    click.echo("To harmonize real studies against it, set:")
    click.echo(click.style(f"  export AREE_CROSSWALK={path}", fg="cyan"))


if __name__ == "__main__":
    main()
