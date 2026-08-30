# Installation

AREE is a standard Python `src`-layout package plus Nextflow workflows and a
Streamlit app. This page covers the Python package, which is
what the `aree` CLI, harmonization, meta-analysis, and prioritization code
depend on.

## Requirements

- Python >= 3.9
- `pip`
- (optional) [Quarto](https://quarto.org) if you plan to render `docs/` as a
  static site
- (optional) Nextflow + Docker/Apptainer if you plan to execute the raw-data
  workflows under `workflows/` — see
  [roadmap.md](roadmap.md) for the current status of that path

## Install

Create and activate a virtual environment, then install AREE in editable mode
with the `dev` and `app` extras:

```bash
python3 -m venv .venv
source .venv/bin/activate

# Python 3.9's bundled pip is too old for editable installs — upgrade first.
pip install --upgrade pip setuptools wheel

pip install -e ".[dev,app]"
```

This installs the packages declared in `pyproject.toml`'s `[tool.setuptools]`
table — `aree`, `common`, `intake`, `harmonize`, `meta_analysis`, `prioritize`,
`reporting`, `validation` — all under the `src/` layout, plus:

- core runtime dependencies: `click`, `pyyaml`, `pandas`, `numpy`, `scipy`,
  `jsonschema`, `tabulate`
- `dev` extra: `pytest`, `ruff`
- `app` extra: `streamlit`

The `aree` console script is registered via `[project.scripts]` and becomes
available on your `PATH` once the package is installed.

## Verify the install

```bash
aree --help
aree list-studies
```

`list-studies` reads `registry/study_registry.csv`. On a freshly cloned repo
this is header-only, so `list-studies` shows no studies until you register the
demo studies once:

```bash
for f in registry/studies/GIGAS_*.yaml; do aree register-study "$f"; done
```

(Re-running plain `register-study` on an already-registered study fails by
design — pass `--update` to overwrite, or use `make demo`, which is
idempotent. See [adding_a_study.md](adding_a_study.md).)

## Running the demo end to end

Once installed, the full demo pipeline runs against the committed demo data
with no external downloads:

```bash
aree validate-study registry/studies/GIGAS_HEAT01.yaml
aree harmonize --study GIGAS_HEAT01
aree meta-analyze --phenotype thermal_tolerance --feature-type gene
aree build-evidence-cards --phenotype thermal_tolerance
```

See [adding_a_study.md](adding_a_study.md), [interpreting_meta_analysis.md](interpreting_meta_analysis.md),
and [interpreting_candidate_scores.md](interpreting_candidate_scores.md) for
what each step produces.

## Launching the interface

```bash
streamlit run app/main.py
```

or, for the static documentation/report site:

```bash
quarto render docs/
```

Neither requires login or a cloud deployment — both run against local files
(`registry/`, `reports/`). See [architecture.md](architecture.md) for how the
interface layer relates to the rest of the system.

## Related documentation

- [design.md](design.md) — data model and architecture this installation is built on
- [architecture.md](architecture.md) — how the installed packages fit together
- [roadmap.md](roadmap.md) — what is fully functional today vs. scaffolded
