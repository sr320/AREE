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
