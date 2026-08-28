"""Shared paths, IO helpers, and controlled-vocabulary loaders used across all AREE packages."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import yaml

REPO_ROOT = Path(__file__).resolve().parents[2]
SCHEMAS_DIR = REPO_ROOT / "schemas"
REGISTRY_DIR = REPO_ROOT / "registry"
STUDIES_DIR = REGISTRY_DIR / "studies"
VOCAB_DIR = REGISTRY_DIR / "controlled_vocabularies"
STUDY_REGISTRY_CSV = REGISTRY_DIR / "study_registry.csv"
DATA_DIR = REPO_ROOT / "data"
DEMO_DIR = DATA_DIR / "demo"
MAPPINGS_DIR = DATA_DIR / "mappings"
REPORTS_DIR = REPO_ROOT / "reports"
EVIDENCE_TABLE_PATH = REPORTS_DIR / "evidence" / "evidence_table.tsv"


class VocabError(ValueError):
    """Raised when a value does not belong to a controlled vocabulary."""


def load_yaml(path) -> dict:
    with open(path) as fh:
        return yaml.safe_load(fh)


def load_json(path) -> dict:
    with open(path) as fh:
        return json.load(fh)


def sha256_file(path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def load_vocab(name: str) -> dict:
    """Load a controlled vocabulary YAML file by stem name (e.g. 'phenotype_ontology').

    Returns a dict keyed by term id -> full term record.
    """
    path = VOCAB_DIR / f"{name}.yaml"
    doc = load_yaml(path)
    return {term["id"]: term for term in doc["terms"]}


def vocab_ids(name: str) -> set:
    return set(load_vocab(name).keys())


def check_in_vocab(value: str, name: str, field_label: str) -> None:
    ids = vocab_ids(name)
    if value not in ids:
        raise VocabError(
            f"{field_label!r} value {value!r} is not a valid term in "
            f"controlled vocabulary '{name}'. Valid ids: {sorted(ids)}"
        )


# --------------------------------------------------------------------------- #
# Species identity
# --------------------------------------------------------------------------- #

GENOME_ASSEMBLIES_PATH = DATA_DIR / "reference" / "genome_assemblies.yaml"


class SpeciesRecord:
    """A resolved species: what the study called it, and what it canonically is."""

    __slots__ = ("term_id", "scientific_name", "ncbi_taxid", "as_reported", "is_synonym")

    def __init__(self, term: dict, as_reported: str, is_synonym: bool):
        self.term_id = term["id"]
        self.scientific_name = term["scientific_name"]
        self.ncbi_taxid = term["ncbi_taxid"]
        self.as_reported = as_reported
        self.is_synonym = is_synonym

    def __repr__(self) -> str:  # pragma: no cover - debugging aid
        return f"SpeciesRecord({self.scientific_name!r}, taxid={self.ncbi_taxid})"


def resolve_species(name: str) -> SpeciesRecord | None:
    """Resolve a scientific name (or accepted synonym) to a canonical species.

    Returns None if the name is not in the vocabulary. Matching ignores case and
    surrounding whitespace but nothing else — a genuinely different name must be
    added to the vocabulary deliberately, not guessed at.
    """
    if not name:
        return None
    wanted = str(name).strip().casefold()
    for term in load_vocab("species").values():
        if wanted == term["scientific_name"].casefold():
            return SpeciesRecord(term, name, is_synonym=False)
        for synonym in term.get("accepted_synonyms") or []:
            if wanted == synonym.casefold():
                return SpeciesRecord(term, name, is_synonym=True)
    return None


# --------------------------------------------------------------------------- #
# Reference assemblies
# --------------------------------------------------------------------------- #


def load_assemblies() -> dict:
    """Assembly records keyed by assembly_id."""
    doc = load_yaml(GENOME_ASSEMBLIES_PATH)
    return {a["assembly_id"]: a for a in doc["assemblies"]}


def resolve_assembly(value: str) -> dict | None:
    """Resolve a study's `genome_assembly` string to a known assembly record.

    Accepts the assembly_id, either accession, or any declared alias. Studies
    sometimes write the accession and the name together
    ("GCA_902806645.1 (cgigas_uk_roslin_v1)"), so each whitespace/parenthesis
    separated token is tried before giving up.
    """
    if not value:
        return None
    assemblies = load_assemblies()
    tokens = [str(value).strip()]
    tokens += [t.strip(" ()") for t in str(value).replace("(", " ").replace(")", " ").split()]

    for token in tokens:
        for record in assemblies.values():
            candidates = {
                record["assembly_id"],
                record.get("ncbi_assembly_accession"),
                record.get("genbank_accession"),
                *(record.get("aliases") or []),
            }
            if token in {c for c in candidates if c}:
                return record
    return None
