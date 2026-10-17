"""AREE command-line interface.

See README.md for the full command reference.
"""
from __future__ import annotations

import datetime
import json
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


@main.command("study-ids")
@click.option("--simulated", "which", flag_value="simulated", help="Only simulated (demo) studies.")
@click.option("--real", "which", flag_value="real", help="Only real studies.")
def study_ids_cmd(which):
    """Print registered study_ids, one per line, read from registry/studies/.

    `aree study-ids --simulated` is the demo study list the Makefile, CI and the
    Pages build iterate over.
    """
    from common import study_ids

    for sid in study_ids({"simulated": True, "real": False}.get(which)):
        click.echo(sid)


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
                   "values are combined in order to form the `condition` column. "
                   "May also name a group from --label-pattern.")
@click.option("--attribute", "extra_attributes", multiple=True,
              help="Additional sample attribute to carry into the sheet (repeatable).")
@click.option("--label-field", default=None,
              help="ENA run-report field holding per-library design labels (e.g. experiment_title), "
                   "for projects that deposit all libraries under one BioSample.")
@click.option("--label-pattern", default=None,
              help="Regex with named groups applied to --label-field; each group becomes an attribute "
                   "and the matched text becomes sample_id. Every run must match.")
@click.option("--methylkit-levels", nargs=2, default=None, metavar="CONTROL TREATMENT",
              help="Also write methylkit_samplesheet.csv coding these two condition levels as 0 and 1.")
@click.option("--library-strategy", "library_strategies", multiple=True,
              help="Keep only runs with this ENA library_strategy (repeatable), e.g. RNA-Seq when a "
                   "project also deposits 16S amplicon runs from the same animals. Check that the "
                   "labels are right first: they are submitter-supplied and sometimes wrong.")
@click.option("--design-table", default=None, type=click.Path(exists=True, dir_okay=False),
              help="TSV of design attributes, for projects whose BioSamples do not carry the design "
                   "(e.g. the authors' submission sheet, committed with the study). Its columns become "
                   "attributes; every sample must match exactly one row.")
@click.option("--design-join", default=None,
              help="Deposited sample attribute (e.g. isolate) that matches the design table's column of "
                   "the same name.")
@click.option("--design-source", default=None,
              help="Where the design table came from (URL or citation); recorded in the provenance file.")
@click.option("--out", "out_dir", default=None, type=click.Path(),
              help="Output directory (default: data/studies/<study_id>/).")
def fetch_samplesheet_cmd(bioproject, study_id, condition_attributes, extra_attributes,
                          label_field, label_pattern, methylkit_levels, library_strategies, design_table,
                          design_join, design_source, out_dir):
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
            label_field=label_field,
            label_pattern=label_pattern,
            methylkit_levels=tuple(methylkit_levels) if methylkit_levels else None,
            library_strategies=list(library_strategies) or None,
            design_table=Path(design_table) if design_table else None,
            design_join=design_join,
            design_source=design_source,
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

    by_strategy = report.get("group_sizes_by_strategy") or {}
    if len(by_strategy) > 1:
        click.echo(click.style(
            "\nWARNING: this sheet mixes library strategies, so every group size above counts "
            "runs of each. Per strategy:", fg="yellow"))
        for strategy, groups in by_strategy.items():
            click.echo(click.style(
                f"    {strategy:10} " + ", ".join(f"{g} n={n}" for g, n in groups.items()), fg="yellow"))
        click.echo(click.style(
            "  If some runs are another assay (e.g. 16S amplicon from the same animals), re-run "
            "with --library-strategy. Check layouts and read counts first: the label itself can "
            "be wrong (LUTIER2022_OA_TIPPING's AMPLICON runs are RNA-seq).", fg="yellow"))

    smallest = report["min_group_size"]
    if smallest < 3:
        click.echo(click.style(
            f"\nWARNING: smallest group has n={smallest}. Differential expression needs "
            "at least 3 biological replicates per group; this design may not support it.",
            fg="yellow",
        ))
    if report["n_samples"] < report["n_runs"]:
        click.echo(click.style(
            f"\nWARNING: {report['n_runs']} runs share {report['n_samples']} BioSample(s). The "
            "deposit does not establish that libraries are distinct animals; confirm "
            "from the paper or the submitter before treating them as biological replicates.",
            fg="yellow",
        ))
    if report.get("submitted_ftp_runs"):
        runs = report["submitted_ftp_runs"]
        click.echo(click.style(
            f"\nNOTE: ENA has no generated FASTQ for {len(runs)} run(s); the manifest uses the "
            f"submitter's own FASTQ (source=submitted_ftp): {', '.join(runs)}. Their file names "
            "and read headers are the submitter's, not ENA's.", fg="yellow"))
    if report.get("unpaired_files_dropped"):
        click.echo(click.style(
            f"\nNOTE: {len(report['unpaired_files_dropped'])} paired run(s) also carry an unpaired "
            "FASTQ, left out of the manifest and listed in ena_provenance.json.", fg="yellow"))
    if report.get("methylkit_samplesheet"):
        mk = report["methylkit_samplesheet"]
        click.echo(f"  methylKit  : {mk['control_level']}=0 (n={mk['n_control']}), "
                   f"{mk['treatment_level']}=1 (n={mk['n_treatment']})")
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
@click.option("--level", type=click.Choice(["gene", "orthogroup"]), default="gene", show_default=True,
              help="gene: pool within each species. orthogroup: pool species' gene-level estimates "
                   "across species through orthogroups (single-copy orthogroups only by default).")
@click.option("--include-multi-copy", is_flag=True,
              help="With --level orthogroup, also pool multi-copy orthogroups where each species has "
                   "evidence for exactly one gene. Off by default: paralogs are not like for like.")
def meta_analyze_cmd(phenotype, feature_type, level, include_multi_copy):
    """Run a random-effects meta-analysis over the harmonized evidence table."""
    if level == "orthogroup":
        from meta_analysis.run import write_orthogroup_meta_analysis

        try:
            result, _species, summary, out_path = write_orthogroup_meta_analysis(
                phenotype, feature_type, include_multi_copy=include_multi_copy)
        except (FileNotFoundError, ValueError) as exc:
            click.echo(click.style(str(exc), fg="red"))
            sys.exit(1)
        click.echo(f"Orthogroup groups seen: {summary['orthogroup_groups']}; pooled across species: "
                   f"{summary['pooled']}; one species only: {summary['single_species']}; set aside: "
                   f"{summary['multi_copy_excluded']} multi-copy, {summary['paralogs_in_species_excluded']} "
                   f"with paralogs in one species, {summary['not_in_active_run_excluded']} not in the "
                   "active orthogroup run.")
        if len(result) == 0:
            click.echo("No orthogroup has gene-level estimates from two or more species.")
            return
        click.echo(tabulate(
            result[["orthogroup_id", "copy_class", "phenotype", "n_species", "k_studies_by_species",
                    "pooled_effect", "p_value", "adjusted_p_value", "i_squared_between_species",
                    "meets_replication_gate"]].head(20),
            headers="keys", showindex=False, floatfmt=".3g",
        ))
        click.echo(f"\nFull results ({len(result)} rows) written to {out_path}; per-species "
                   "estimates beside it. meets_replication_gate: at least two species each with k >= 2.")
        return
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


@main.command("build-cross-species-cards")
@click.option("--phenotype", default=None, help="Phenotype ontology term id to filter on (default: all).")
@click.option("--feature-type", default=None, help="Feature type to filter on (default: all).")
@click.option("--max-adjusted-p", default=DEFAULT_MAX_ADJUSTED_P, show_default=True, type=float,
              help="Render a card only when the pooled cross-species adjusted p, or any one species' "
                   "gene-level adjusted p, is at or below this value. Every row is still listed in "
                   "reports/evidence_cards/cross_species/candidates.tsv.")
@click.option("--all-cards", is_flag=True, help="Render a card for every orthogroup row.")
@click.option("--include-multi-copy", is_flag=True,
              help="Also pool multi-copy orthogroups where each species has evidence for exactly one gene.")
def build_cross_species_cards_cmd(phenotype, feature_type, max_adjusted_p, all_cards, include_multi_copy):
    """Tier orthogroup-level (cross-species) results and write cross-species evidence cards."""
    from reporting.cross_species_cards import TIER_ORDER as CROSS_TIERS
    from reporting.cross_species_cards import build_cross_species_cards_report

    try:
        report = build_cross_species_cards_report(
            phenotype=phenotype, feature_type=feature_type, max_adjusted_p=max_adjusted_p,
            all_cards=all_cards, include_multi_copy=include_multi_copy)
    except (FileNotFoundError, ValueError) as exc:
        click.echo(click.style(str(exc), fg="red"))
        sys.exit(1)
    if report.n_candidates == 0:
        click.echo("No orthogroup has gene-level estimates from two or more species, so there are no "
                   "cross-species candidates. Harmonize a study of a second species first.")
        return
    click.echo(click.style(
        f"Tiered {report.n_candidates} cross-species rows ({report.candidates_path}); "
        f"wrote {len(report.index)} cards to {report.candidates_path.parent}/", fg="green"))
    counts = {}
    for row in report.index:
        counts[row["tier"]] = counts.get(row["tier"], 0) + 1
    for tier in CROSS_TIERS:
        if tier in counts:
            click.echo(f"  {tier}: {counts[tier]}")


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


@main.command("fetch-reference")
@click.option("--assembly", "assembly", required=True,
              help="Registered assembly: assembly_id or accession, e.g. GCF_965363235.1.")
@click.option("--dest", default=None, type=click.Path(),
              help="Bundle directory. Defaults to data/reference/<ACCESSION> (often a symlink "
                   "to external storage; see docs/reference_bundles.md).")
@click.option("--genome", "with_genome", is_flag=True,
              help="Also fetch the genome FASTA (genomic.fna), needed for methylation.")
@click.option("--proteins", "with_proteins", is_flag=True,
              help="Also fetch the RefSeq proteins and feature table and write primary_proteins.faa "
                   "(one protein per gene), the input to `aree build-orthogroups`.")
@click.option("--plan", "plan_only", is_flag=True,
              help="Only resolve the NCBI folder and checksums; download no sequence data.")
@click.option("--check", "check_only", is_flag=True,
              help="Verify an existing bundle against its reference_provenance.json.")
def fetch_reference_cmd(assembly, dest, with_genome, with_proteins, plan_only, check_only):
    """Fetch a RefSeq reference bundle (transcriptome, GTF, tx2gene) for raw reanalysis.

    Every file is checked against NCBI's MD5s and recorded with SHA-256 in
    reference_provenance.json, along with whether the GTF's annotation release
    matches the registry and how many of its genes the species' crosswalk knows.
    """
    from intake.reference_bundle import (
        ReferenceBundleError,
        check_bundle,
        default_bundle_dir,
        fetch_bundle,
        plan_bundle,
    )

    try:
        if check_only:
            from common import resolve_assembly

            record = resolve_assembly(assembly)
            if record is None:
                raise ReferenceBundleError(f"{assembly!r} is not a registered assembly")
            target = Path(dest) if dest else default_bundle_dir(record["ncbi_assembly_accession"])
            problems = check_bundle(target)
            for problem in problems:
                click.echo(click.style(f"  - {problem}", fg="red"))
            click.echo(click.style(f"FAIL: {target}", fg="red") if problems
                       else click.style(f"OK: {target}", fg="green"))
            sys.exit(1 if problems else 0)

        if plan_only:
            plan = plan_bundle(assembly, with_genome=with_genome, with_proteins=with_proteins)
            click.echo(f"{plan['assembly']['assembly_id']}: {plan['base_url']}")
            for f in plan["files"]:
                click.echo(f"  {f['remote']}  md5 {f['ncbi_md5']}  -> {f['local']}")
            return

        prov = fetch_bundle(assembly, dest=Path(dest) if dest else None, with_genome=with_genome,
                            with_proteins=with_proteins)
    except (ReferenceBundleError, OSError) as exc:
        click.echo(click.style(str(exc), fg="red"))
        sys.exit(1)

    click.echo(click.style(f"Bundle {prov['bundle']} for {prov['assembly_id']} ({prov['species']})",
                           fg="green"))
    click.echo(f"  tx2gene: {prov['tx2gene']['transcripts']} transcripts, "
               f"{prov['tx2gene']['genes']} genes")
    cov = prov["crosswalk_coverage"]
    if cov.get("checked"):
        colour = "green" if cov["pct"] >= 99 else "yellow"
        click.echo(click.style(f"  genes in {cov['crosswalk']}: {cov['found_in_crosswalk']} "
                               f"({cov['pct']}%)", fg=colour))
    else:
        click.echo(click.style("  no crosswalk registered for this assembly; gene IDs unchecked",
                               fg="yellow"))
    if prov.get("primary_proteome"):
        pp = prov["primary_proteome"]
        click.echo(f"  primary proteome: {pp['genes']} genes (one protein each)"
                   + (f", {pp['crosswalk_coverage']['pct']}% in the crosswalk"
                      if pp["crosswalk_coverage"].get("checked") else ""))
    if not prov["annotation_release_matches_registry"]:
        click.echo(click.style(
            f"  WARNING: GTF annotation-source {prov['gtf_header'].get('annotation-source')!r} "
            f"does not name the registered release {prov['registered_annotation_release']!r}. "
            "Update genome_assemblies.yaml (and rebuild the crosswalk) before using it.",
            fg="red"))
        sys.exit(1)


@main.command("build-orthogroups")
@click.option("--orthofinder", "orthofinder_bin", default="orthofinder", show_default=True,
              help="OrthoFinder executable; its directory is put first on PATH so its DIAMOND and "
                   "MCL are used (e.g. ~/miniforge3/envs/orthofinder/bin/orthofinder).")
@click.option("--threads", default=4, show_default=True, type=int, help="OrthoFinder -t and -a.")
@click.option("--run-id", default=None,
              help="Name of the run directory (default: orthofinder<version>_<date>).")
@click.option("--work-dir", default=None, type=click.Path(),
              help="Where proteomes are staged and OrthoFinder writes its intermediates "
                   "(default: data/reference/orthogroups/_work/<run_id>, gitignored).")
@click.option("--from-results", "from_results", default=None, type=click.Path(exists=True, file_okay=False),
              help="Read an existing OrthoFinder Results_* directory instead of running OrthoFinder.")
@click.option("--activate", is_flag=True,
              help="Make this run the one harmonization attaches to evidence (active.yaml).")
def build_orthogroups_cmd(orthofinder_bin, threads, run_id, work_dir, from_results, activate):
    """Infer orthogroups across every species with a crosswalk (OrthoFinder on RefSeq proteomes).

    Needs each species' primary proteome: `aree fetch-reference --assembly <acc> --proteins`.
    Writes data/reference/orthogroups/<run_id>/ (gene_orthogroups.tsv, orthogroups.tsv,
    provenance.json). See src/mappings/orthogroups.py.
    """
    from datetime import date as _date

    from mappings import orthogroups as og

    try:
        version = og.orthofinder_version(orthofinder_bin)
        run_id = run_id or f"orthofinder{version.split()[-1]}_{_date.today().isoformat()}"
        work = Path(work_dir) if work_dir else og.ORTHOGROUPS_DIR / "_work" / run_id
        species = og.run_species()
        staged = og.stage_proteomes(species, work / "proteomes")
        click.echo(f"Staged {len(staged)} proteomes: "
                   + ", ".join(f"{s['slug']} ({s['proteins']})" for s in staged))
        if from_results:
            results, command = Path(from_results), None
        else:
            click.echo(f"Running {version} with {threads} threads (hours for this many proteomes)...")
            results, command = og.run_orthofinder(work / "proteomes", orthofinder_bin, threads)
        slug_taxid = {s["slug"]: s["ncbi_taxid"] for s in staged}
        rows = og.parse_orthogroups(results, slug_taxid)
        summary = og.summarize(rows, [s["ncbi_taxid"] for s in staged])
        prov = og.write_run(og.ORTHOGROUPS_DIR / run_id, rows, summary, staged, version=version,
                            command=command, results_dir=results)
    except og.OrthogroupError as exc:
        click.echo(click.style(str(exc), fg="red"))
        sys.exit(1)

    click.echo(click.style(f"Wrote {og.ORTHOGROUPS_DIR / run_id}", fg="green"))
    click.echo(f"  orthogroups: {prov['orthogroups']['total']} "
               + ", ".join(f"{k} {v}" for k, v in prov["orthogroups"]["by_copy_class"].items()))
    for slug, a in prov["per_species"].items():
        click.echo(f"  {slug:7} {a['in_multi_species_orthogroups']}/{a['genes']} genes in "
                   "orthogroups shared with another species")
    if activate:
        og.activate(run_id)
        click.echo(f"Active orthogroup run: {run_id}")


@main.command("build-crosswalk")
@click.option("--taxid", default=29159, show_default=True, type=int,
              help="NCBI taxonomy id to build the crosswalk for. Must be a registered species.")
@click.option("--slug", default=None,
              help="Output filename prefix, e.g. 'mgigas' -> mgigas_gene_id_crosswalk.tsv. "
                   "Defaults to 'mgigas' for the Pacific oyster; required for any other species.")
@click.option("--organism", default=None,
              help="Human-readable organism label recorded in provenance. "
                   "Defaults to the species vocabulary's name, synonyms and label.")
@click.option("--assembly", default=None,
              help="Reference annotation recorded in provenance (and on every evidence record "
                   "as identifier_annotation_release). Defaults to the species' NCBI reference "
                   "assembly in data/reference/genome_assemblies.yaml.")
@click.option("--gene-info", "gene_info", default=None, type=click.Path(exists=True),
              help="Pre-downloaded taxid-filtered NCBI gene_info TSV.")
@click.option("--uniprot", default=None, type=click.Path(exists=True),
              help="Pre-downloaded UniProt TSV export.")
@click.option("--download", is_flag=True, help="Force re-download of reference sources.")
def build_crosswalk_cmd(taxid, slug, organism, assembly, gene_info, uniprot, download):
    """Build a real identifier crosswalk from NCBI Gene and UniProtKB.

    Downloads reference data from public sources (~230 MB streamed from NCBI) and
    writes a crosswalk plus a JSON provenance sidecar to data/reference/crosswalk/.
    The species and its reference assembly must already be registered (see
    docs/adding_a_species.md), so that provenance names the right annotation.
    """
    from common import load_assemblies, load_vocab
    from harmonize.identifiers import crosswalk_taxid
    from mappings.build_crosswalk import CROSSWALK_DIR

    def fail(message):
        click.echo(click.style(message, fg="red"))
        sys.exit(1)

    term = next((t for t in load_vocab("species").values() if t["ncbi_taxid"] == taxid), None)
    if term is None:
        fail(f"taxid {taxid} is not in registry/controlled_vocabularies/species.yaml; "
             "register the species first (docs/adding_a_species.md).")
    reference = [a for a in load_assemblies().values()
                 if a.get("ncbi_taxid") == taxid and a.get("is_ncbi_reference")]
    if assembly is None:
        if not reference:
            fail(f"no NCBI reference assembly for taxid {taxid} in "
                 "data/reference/genome_assemblies.yaml; register it, or pass --assembly.")
        ref = reference[0]
        assembly = (f"{ref['ncbi_assembly_accession']} ({ref['assembly_id']}), "
                    f"annotation release {ref.get('annotation_release')}")
    if organism is None:
        synonyms = term.get("accepted_synonyms") or []
        organism = term["scientific_name"]
        if synonyms:
            organism += f" ({', '.join(synonyms)})"
        organism += f", {term['label']}"
    if slug is None:
        if taxid != 29159:
            fail("--slug is required for any species other than the Pacific oyster, "
                 "so its crosswalk cannot overwrite another species' file.")
        slug = "mgigas"

    # Never overwrite a crosswalk built for a different species.
    existing = CROSSWALK_DIR / f"{slug}_gene_id_crosswalk.tsv"
    if existing.exists():
        existing_taxid = crosswalk_taxid(existing)
        if existing_taxid is not None and existing_taxid != taxid:
            fail(f"{existing} was built for taxid {existing_taxid}; choose a different --slug "
                 f"for taxid {taxid}.")

    from mappings.build_crosswalk import build

    try:
        path = build(taxid=taxid, organism=organism, assembly=assembly, gene_info_path=gene_info,
                     uniprot_path=uniprot, download=download, slug=slug)
    except (RuntimeError, ValueError, FileNotFoundError) as exc:
        click.echo(click.style(str(exc), fg="red"))
        sys.exit(1)

    prov = path.with_suffix("").with_suffix(".provenance.json")
    click.echo(click.style(f"Wrote {path}", fg="green"))
    click.echo(f"Provenance: {prov}")
    click.echo("")
    try:
        shown = path.relative_to(REPO_ROOT)
    except ValueError:
        shown = path
    click.echo("Real studies resolve against it once it is registered for this species'")
    click.echo("assemblies in data/reference/genome_assemblies.yaml:")
    click.echo(click.style(f"  identifier_crosswalk: \"{shown}\"", fg="cyan"))



@main.command("build-metabolite-links")
@click.option("--taxid", default=29159, show_default=True, type=int,
              help="NCBI taxid; must have a KEGG organism code in mappings.build_metabolite_links.")
@click.option("--offline", is_flag=True,
              help="Rebuild only from cached KEGG responses in data/mappings/_sources/kegg/.")
def build_metabolite_links_cmd(taxid, offline):
    """Build the real metabolite -> enzyme / transporter gene link table from KEGG.

    Queries KEGG only for the compounds curated in
    data/mappings/metabolite_compound_map.tsv; transporter links come from the
    curated data/mappings/transporter_substrates.tsv, each KO verified against
    KEGG BRITE ko02000. Writes
    data/mappings/metabolite_gene_links_<slug>.tsv plus a .provenance.json
    sidecar. The output is git-ignored and must not be committed (KEGG terms).
    See docs/metabolite_gene_links.md.
    """
    from mappings.build_metabolite_links import build
    from mappings.metabolite_links import clear_caches

    try:
        path = build(taxid=taxid, offline=offline)
    except (OSError, ValueError) as exc:
        click.echo(click.style(str(exc), fg="red"))
        sys.exit(1)
    clear_caches()
    prov = json.loads(path.with_suffix(".provenance.json").read_text())
    cov = prov["coverage"]
    by_type = cov["rows_by_link_type"]
    click.echo(click.style(f"Wrote {path} ({cov['rows']} links: {by_type['enzyme']} enzyme, "
                           f"{by_type['transporter']} transporter)", fg="green"))
    transporters = prov["transporter_substrates"]
    click.echo(f"  curated transporter KOs: {len(transporters['kos_used'])} in KEGG BRITE ko02000, "
               f"{len(transporters['rejected_not_in_ko02000'])} rejected (not in ko02000), "
               f"{len(transporters['kos_without_crg_genes'])} with no crg gene")
    click.echo(f"  compounds with links: {cov['compounds_with_links']} of {cov['compounds_queried']} queried"
               f" (currency skipped: {', '.join(prov['currency_metabolites']['skipped']) or 'none'})")
    click.echo(f"  hub compounds (not credited): {', '.join(cov['hub_compounds']) or 'none'}")
    click.echo(f"  KEGG genes: {cov['kegg_genes_current']} current, {cov['kegg_genes_retired_remapped']} "
               f"remapped from retired IDs, {cov['kegg_genes_unresolved']} unresolved")
    click.echo("  Not for redistribution: do not commit this file (docs/metabolite_gene_links.md).")

if __name__ == "__main__":
    main()
