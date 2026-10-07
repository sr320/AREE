"""AREE command-line interface.

See README.md for the full command reference.
"""
from __future__ import annotations

import datetime
import subprocess
import sys
from pathlib import Path

import click
from tabulate import tabulate

from aree import __version__
from common import REPO_ROOT
from harmonize.core import harmonize_processed_table, harmonize_study
from intake.ena_samplesheet import ENAError, build_samplesheet
from intake.registry import DuplicateStudyError, list_studies, rebuild_registry, register_study
from intake.run_intake import IntakeError, run_intake
from intake.schema_validate import validate_study_file
from intake.study_claims import ClaimError, registered_on, start_study
from meta_analysis.run import write_meta_analysis
from prioritize.rank import TIER_ORDER
from reporting.dashboard import DEFAULT_OUT_PATH as DEFAULT_DASHBOARD_PATH
from reporting.dashboard import build_dashboard_data
from reporting.evidence_cards import DEFAULT_MAX_ADJUSTED_P, build_evidence_cards_report
from reporting.top_candidates import (
    DEFAULT_CANDIDATES_PATH,
    DEFAULT_OUT_PATH,
    build_top_candidates_summary,
)
from validation.study_scope import check_study_scope


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
    """Validate and add a study to registry/study_registry.csv (a generated, gitignored index)."""
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
        click.echo("No studies in the index yet. Run `aree build-registry` to build it from registry/studies/.")
        return
    click.echo(tabulate(
        [[r["study_id"], r["assay_type"], r["analysis_mode"], r["qc_status"], r["analysis_status"]] for r in rows],
        headers=["study_id", "assay_type", "analysis_mode", "qc_status", "analysis_status"],
    ))


@main.command("build-registry")
@click.option("--strict", is_flag=True, help="Exit non-zero if any study YAML fails validation.")
def build_registry_cmd(strict):
    """Rebuild registry/study_registry.csv from every study YAML.

    The CSV is a generated index and is not committed, so study branches
    developed on different machines never conflict on it.
    """
    rows, invalid = rebuild_registry()
    click.echo(click.style(f"Indexed {len(rows)} studies.", fg="green"))
    for path, errors in invalid.items():
        click.echo(click.style(f"INVALID, left out: {path}", fg="red"))
        for e in errors:
            click.echo(f"  - {e}")
    if invalid and strict:
        sys.exit(1)


@main.command("start-study")
@click.option("--accession", required=True,
              help="Data accession this job covers (e.g. PRJNA690951, GSE12345, PXD002316). Names the branch.")
@click.option("--study-id", "study_id", required=True, help="AREE study_id to register it under.")
@click.option("--dir", "worktree_dir", default=None, type=click.Path(),
              help="Where to create the study's git worktree (default: ../AREE-<STUDY_ID>).")
@click.option("--registered-by", default=None, help="Curator name for the YAML (default: git user.name).")
@click.option("--machine", default=None, help="Machine name recorded in the PR (default: this host's short name).")
@click.option("--no-push", is_flag=True, help="Create the branch and stubs locally only; nothing is claimed.")
@click.option("--no-pr", is_flag=True, help="Push the claim but do not open a draft PR.")
def start_study_cmd(accession, study_id, worktree_dir, registered_by, machine, no_push, no_pr):
    """Claim one study for this machine: branch study/<ACCESSION>/<STUDY_ID>, stub YAML + config, draft PR.

    Refuses an accession or study_id already registered on main or already
    claimed by another study/ branch. Pushing the branch is the claim; if
    another machine pushed it first, the push is rejected and nothing is
    downloaded. See docs/parallel_study_runs.md.
    """
    if registered_by is None:
        registered_by = subprocess.run(["git", "config", "user.name"], capture_output=True,
                                       text=True, cwd=REPO_ROOT).stdout.strip() or "unknown"
    try:
        res = start_study(accession, study_id, repo=REPO_ROOT, worktree_dir=worktree_dir,
                          registered_by=registered_by, machine=machine,
                          push=not no_push, open_pr=not no_pr)
    except ClaimError as exc:
        click.echo(click.style(str(exc), fg="red"))
        sys.exit(1)
    click.echo(click.style(f"Branch   : {res['branch']}", fg="green"))
    click.echo(f"Worktree : {res['worktree']}")
    for f in res["files"]:
        click.echo(f"  stub   : {f}")
    for f in res["linked"]:
        click.echo(f"  linked : {f} (same drive as this checkout)")
    if not res["linked"]:
        click.echo(click.style("  No data/raw or data/reference/GCF_* symlinks found to copy; link this "
                               "machine's drives into the worktree before running.", fg="yellow"))
    if res["pushed"]:
        click.echo(click.style("Claimed  : branch pushed", fg="green"))
    else:
        click.echo(click.style("NOT claimed: --no-push given; other machines cannot see this study yet.",
                               fg="yellow"))
    if res["pr_url"]:
        click.echo(f"Draft PR : {res['pr_url']}")
    click.echo(f"\nNext: cd {res['worktree']} and follow docs/adding_a_study.md.")


@main.command("check-study-scope")
@click.option("--branch", required=True, help="Head branch name, study/<ACCESSION>/<STUDY_ID>.")
@click.option("--base", default="origin/main", show_default=True, help="Ref the PR merges into.")
def check_study_scope_cmd(branch, base):
    """Fail if a study PR changes files its study does not own (run by CI on study/ branches)."""
    from common import STUDIES_DIR, load_yaml
    from intake.study_claims import parse_branch

    diff = subprocess.run(["git", "diff", "--name-only", f"{base}...HEAD"], capture_output=True,
                          text=True, cwd=REPO_ROOT)
    if diff.returncode != 0:
        click.echo(click.style(diff.stderr.strip(), fg="red"))
        sys.exit(1)
    changed = [p for p in diff.stdout.splitlines() if p]
    parsed = parse_branch(branch)
    yaml_path = STUDIES_DIR / f"{parsed[1]}.yaml" if parsed else None
    study = load_yaml(yaml_path) if yaml_path and yaml_path.exists() else None
    errors, warnings = check_study_scope(branch, changed, study, registered_on(base, REPO_ROOT))
    for w in warnings:
        click.echo(click.style(f"WARNING: {w}", fg="yellow"))
    for e in errors:
        click.echo(click.style(f"ERROR: {e}", fg="red"))
    if errors:
        sys.exit(1)
    click.echo(click.style(f"{branch}: {len(changed)} changed files, all owned by the study.", fg="green"))


@main.command("fetch-samplesheet")
@click.option("--bioproject", required=True, help="ENA/SRA BioProject accession, e.g. PRJNA1329250.")
@click.option("--study", "study_id", required=True, help="AREE study_id this project will be registered under.")
@click.option("--condition-attribute", "condition_attributes", multiple=True, required=True,
              help="Sample attribute defining an experimental group. Repeatable; "
                   "values are combined in order to form the `condition` column.")
@click.option("--attribute", "extra_attributes", multiple=True,
              help="Additional sample attribute to carry into the sheet (repeatable).")
@click.option("--out", "out_dir", default=None, type=click.Path(),
              help="Output directory (default: data/studies/<study_id>/).")
def fetch_samplesheet_cmd(bioproject, study_id, condition_attributes, extra_attributes, out_dir):
    """Build a sample sheet + FASTQ manifest from a BioProject's deposited metadata.

    Downloads no sequence data. Reads ENA's run report and sample attributes and
    writes the design, the FASTQ URLs with ENA's MD5s, and a provenance record.
    """
    out_dir = Path(out_dir) if out_dir else REPO_ROOT / "data" / "studies" / study_id
    try:
        report = build_samplesheet(
            bioproject,
            study_id=study_id,
            out_dir=out_dir,
            condition_attributes=list(condition_attributes),
            extra_attributes=list(extra_attributes),
        )
    except ENAError as exc:
        click.echo(click.style(str(exc), fg="red"))
        sys.exit(1)

    click.echo(f"{report['bioproject']}: {report['n_runs']} runs / {report['n_samples']} samples")
    click.echo(f"  strategies : {', '.join(report['library_strategies'])}")
    click.echo(f"  layouts    : {', '.join(report['library_layouts'])}")
    click.echo(f"  download   : {report['total_fastq_bytes'] / 1e9:.1f} GB of FASTQ")
    click.echo("  groups:")
    for group, n in report["group_sizes"].items():
        click.echo(f"    {group:28} n={n}")

    smallest = report["min_group_size"]
    if smallest < 3:
        click.echo(click.style(
            f"\nWARNING: smallest group has n={smallest}. Differential expression needs "
            "at least 3 biological replicates per group; this design may not support it.",
            fg="yellow",
        ))
    click.echo(click.style(f"\nWrote samplesheet, FASTQ manifest, and provenance to {out_dir}", fg="green"))


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
@click.option("--comparison", "comparison_id", default=None,
              help="comparison_id the --input file belongs to. Required when the filename "
                   "does not match a declared results_file, e.g. for workflow output.")
@click.option("--date", "date_generated", default=None, help="Override date_generated (ISO 8601). Defaults to today.")
def harmonize_cmd(study_id, input_path, comparison_id, date_generated):
    """Harmonize a processed study result table (or an entire study) into the shared evidence table."""
    date_generated = date_generated or _today()
    try:
        if input_path:
            df = harmonize_processed_table(
                study_id, input_path, date_generated=date_generated, comparison_id=comparison_id
            )
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
        result[["feature_id_standardized", "phenotype", "feature_type", "k_studies", "pooled_effect",
                "p_value", "adjusted_p_value", "n_tests_in_family", "i_squared"]].head(20),
        headers="keys", showindex=False, floatfmt=".3g",
    ))
    click.echo(f"\nFull results ({len(result)} rows) written to {out_path}")
    click.echo("adjusted_p_value is Benjamini-Hochberg within each phenotype / feature-type / "
               "origin / species family (n_tests_in_family).")


@main.command("build-evidence-cards")
@click.option("--phenotype", default=None, help="Phenotype ontology term id to filter on (default: all).")
@click.option("--feature-type", default=None, help="Feature type to filter on (default: all).")
@click.option("--max-adjusted-p", default=DEFAULT_MAX_ADJUSTED_P, show_default=True, type=float,
              help="Render a card only for candidates whose BH-adjusted pooled p, or at least one "
                   "contributing study's adjusted p, is at or below this value. Every candidate is "
                   "still listed in reports/evidence_cards/candidates.tsv.")
@click.option("--all-cards", is_flag=True,
              help="Render a card for every candidate regardless of significance. On a genome-wide "
                   "pool this means tens of thousands of files.")
def build_evidence_cards_cmd(phenotype, feature_type, max_adjusted_p, all_cards):
    """Rank candidate biomarkers and write evidence cards under reports/evidence_cards/."""
    report = build_evidence_cards_report(
        phenotype=phenotype, feature_type=feature_type,
        max_adjusted_p=max_adjusted_p, all_cards=all_cards,
    )
    if report.n_candidates == 0:
        click.echo("No candidates found for the given filters. Run `aree meta-analyze` first if needed.")
        return
    index = report.index
    click.echo(click.style(
        f"Ranked {report.n_candidates} candidates ({report.candidates_path}); "
        f"Wrote {len(index)} evidence cards to {report.candidates_path.parent}/", fg="green",
    ))
    if not all_cards:
        click.echo(f"  cards rendered for candidates with adjusted p <= {max_adjusted_p} "
                   "(pooled, or in any contributing study); pass --all-cards to render every candidate")
    tiers = {}
    for row in index:
        tiers[row["tier"]] = tiers.get(row["tier"], 0) + 1
    for tier in TIER_ORDER:
        if tier in tiers:
            click.echo(f"  {tier}: {tiers[tier]}")


@main.command("top-candidates")
@click.option("--n", "n", default=10, show_default=True, type=int,
              help="Number of top-ranked candidates to show per phenotype.")
@click.option("--candidates", "candidates_path", default=DEFAULT_CANDIDATES_PATH, show_default=True,
              type=click.Path(exists=True, path_type=Path),
              help="Path to a ranked candidates.tsv (from `aree build-evidence-cards`).")
@click.option("--out", "out_path", default=DEFAULT_OUT_PATH, show_default=True, type=click.Path(path_type=Path),
              help="Where to write the Markdown summary.")
def top_candidates_cmd(n, candidates_path, out_path):
    """Summarize the top-N ranked candidates per phenotype from an existing candidates.tsv."""
    try:
        out = build_top_candidates_summary(n=n, candidates_path=candidates_path, out_path=out_path)
    except FileNotFoundError as exc:
        click.echo(click.style(str(exc), fg="red"))
        sys.exit(1)
    click.echo(click.style(f"Wrote top-{n}-per-phenotype summary to {out}", fg="green"))


@main.command("build-dashboard")
@click.option("--out", "out_path", default=DEFAULT_DASHBOARD_PATH, show_default=True,
              type=click.Path(path_type=Path), help="Where to write the dashboard JSON.")
def build_dashboard_cmd(out_path):
    """Summarize registry, pipeline, and findings state into docs/dashboard/data.json.

    Reads whatever the pipeline has produced so far; missing outputs are reported
    as missing on the page rather than hidden. `quarto render docs` then renders
    docs/index.qmd (the site home page) from it.
    """
    data = build_dashboard_data(out_path)
    present = data["artifacts_present"]
    missing = [k for k, v in present.items() if not v]
    k = data["kpis"]
    click.echo(
        f"{k['n_studies_real']} real + {k['n_studies_simulated']} simulated studies, "
        f"{k['n_evidence_records_real']} real evidence records, "
        f"{k['n_candidates_real_high_priority']} real high-priority candidates"
    )
    if missing:
        click.echo(click.style(f"Pipeline outputs not found: {', '.join(missing)} "
                               "(run `make demo` to populate them)", fg="yellow"))
    click.echo(click.style(f"Wrote dashboard data to {out_path}", fg="green"))


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
