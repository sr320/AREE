# Contributing to AREE

Thank you for helping build an open, reproducible resilience-biomarker resource.
There are two main contribution paths: **code** and **curated datasets**.

## Ground rules

AREE has non-negotiable scientific principles (see
[docs/design.md](docs/design.md) and [docs/governance_and_provenance.md](docs/governance_and_provenance.md)):

- Preserve provenance for every result (source, checksum, parameters, versions).
- Never overstate causality — AREE reports associations, not confirmed mechanism.
- Never silently drop missing metadata or contradictory findings; surface them.
- Preserve original identifiers and treatment descriptions; map to controlled
  vocabularies additively, never destructively.

Any contribution that weakens these guarantees will be asked to change.

## Development setup

```bash
python -m venv .venv && source .venv/bin/activate
pip install -e ".[dev,app]"
pytest -q
ruff check src tests
```

## Contributing a curated dataset

1. Copy `registry/studies/_TEMPLATE.yaml` to
   `registry/studies/<STUDY_ID>.yaml` (or use `registry/studies/_batch_template.csv`
   for several studies at once).
2. Fill in every field. Use controlled-vocabulary term IDs from
   `registry/controlled_vocabularies/` for phenotype, stressor, tissue,
   life stage, and assay type. Preserve the exact source wording in
   `stressor_original`.
3. Place the processed results table (if any) under `data/` and point
   `comparisons[].results_file` at it. Match the column schema of the relevant
   `data/demo/<assay>/` example.
4. Validate: `aree validate-study registry/studies/<STUDY_ID>.yaml`.
5. Register and harmonize locally to confirm it flows end-to-end:
   `aree register-study ...` then `aree harmonize --study <STUDY_ID>`.
6. Open a pull request. Do **not** commit generated `reports/` artifacts
   (they are gitignored).

See [docs/adding_a_study.md](docs/adding_a_study.md) for a full walkthrough,
and [docs/resilience_vs_exposure.md](docs/resilience_vs_exposure.md) for how to
classify a study's phenotype relevance correctly.

## Contributing code

1. Open an issue describing the change first for anything non-trivial.
2. Add or update tests under `tests/` — new logic without a test will not be
   merged. The suite must stay green (`pytest -q`) and lint-clean
   (`ruff check src tests`).
3. Keep the scoring formula transparent: any change to candidate scoring must
   go through the named, documented components in
   `src/prioritize/scoring.py` — no opaque or per-candidate special-casing.
4. Update the relevant `docs/` page in the same PR.

## Adding a new species or genome version

See [docs/adding_a_species.md](docs/adding_a_species.md) and
[docs/handling_genome_versions.md](docs/handling_genome_versions.md).

## Code of conduct

By participating you agree to abide by [CODE_OF_CONDUCT.md](CODE_OF_CONDUCT.md).
