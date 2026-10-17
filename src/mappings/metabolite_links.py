"""Metabolite -> enzyme gene links, the weaker second kind of multi-omics convergence.

Gene-level layers (transcriptomics, proteomics, DNA methylation) converge by
sharing a standardized gene identifier. A metabolite is not a gene product, so
it can only be *linked* to genes biochemically:

    reported metabolite name
      -> KEGG compound      (curated: data/mappings/metabolite_compound_map.tsv)
      -> KEGG reaction      (the compound appears on either side of the equation)
      -> KEGG orthology     (the KO annotated to catalyse that reaction)
      -> crg gene           (KEGG's KO assignment for C. gigas = an NCBI GeneID)

Every step is annotation, not measurement, so a link is never `exact` and never
`one_to_one_ortholog` (see docs/metabolite_gene_links.md section 7). Currency
metabolites (ATP, NAD+, H2O, ...) and hub compounds reaching more than
MAX_GENES_PER_COMPOUND genes are kept visible but never credited: without that
filter a change in ATP would "link" to most of the metabolic genome.

Link tables are chosen per evidence partition, like the identifier crosswalk:
simulated evidence resolves against the simulated demo table and real evidence
against the table registered for its species' assembly in
data/reference/genome_assemblies.yaml (`metabolite_gene_links`). The two are
never merged. The real KEGG-derived table is built locally by
`aree build-metabolite-links` and is not committed (KEGG redistribution terms).
"""
from __future__ import annotations

import csv
import re
from collections import defaultdict
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path

from common import (
    MAPPINGS_DIR,
    REPO_ROOT,
    load_assemblies,
    load_json,
    load_yaml,
    sha256_file,
)

COMPOUND_MAP_PATH = MAPPINGS_DIR / "metabolite_compound_map.tsv"
CURRENCY_PATH = MAPPINGS_DIR / "currency_metabolites.yaml"
DEMO_LINKS_PATH = MAPPINGS_DIR / "metabolite_gene_links_demo.tsv"

# A compound annotated to more crg genes than this is a hub: crediting it would
# turn "this metabolite changed" into "much of the metabolic genome is linked".
MAX_GENES_PER_COMPOUND = 25

LINK_COLUMNS = [
    "kegg_compound", "compound_name", "kegg_reaction", "equation_side", "ko", "ncbi_gene_id", "link_type",
]

# Link-table availability for a partition, as reported on every candidate.
STATUS_AVAILABLE = "available"
STATUS_NOT_BUILT = "not_built"            # registered for the assembly, not built on this machine
STATUS_NOT_REGISTERED = "not_registered"  # no table exists for this species / data origin


def normalize_metabolite_name(name) -> str:
    """Lowercase, trim, and fold spaces/hyphens to underscores: 'Glycine betaine' -> 'glycine_betaine'."""
    return re.sub(r"[\s\-]+", "_", str(name).strip().casefold())


def _read_tsv(path: Path) -> list:
    with open(path) as fh:
        return list(csv.DictReader((line for line in fh if not line.startswith("#")), delimiter="\t"))


@lru_cache(maxsize=4)
def load_compound_map(path: Path = COMPOUND_MAP_PATH) -> dict:
    """normalized metabolite name -> list of {kegg_compound, chebi_id, match_basis, note}."""
    out = defaultdict(list)
    for row in _read_tsv(path):
        out[normalize_metabolite_name(row["metabolite_name"])].append({
            "kegg_compound": row["kegg_compound"],
            "chebi_id": row.get("chebi_id") or None,
            "match_basis": row["match_basis"],
            "note": row.get("note") or "",
        })
    return dict(out)


@lru_cache(maxsize=4)
def load_currency_compounds(path: Path = CURRENCY_PATH) -> dict:
    """KEGG compound -> reason it is never used as a link."""
    doc = load_yaml(path)
    return {c["kegg_compound"]: c["reason"] for c in doc["compounds"]}


@dataclass(frozen=True)
class MetaboliteGeneLink:
    """One metabolite -> compound -> reaction -> KO -> gene path."""

    metabolite: str
    kegg_compound: str
    compound_match_basis: str
    kegg_reaction: str
    equation_side: str
    ko: str
    ncbi_gene_id: str
    link_type: str
    n_genes_for_compound: int
    link_confidence: str
    credited: bool
    exclusion_reason: str | None = None

    def describe(self) -> str:
        return (f"{self.metabolite}→{self.ncbi_gene_id} via {self.kegg_compound}/"
                f"{self.kegg_reaction}/{self.ko}")


@dataclass(frozen=True)
class CompoundExclusion:
    """A compound a metabolite resolves to that contributes no credited link."""

    metabolite: str
    kegg_compound: str
    reason: str  # currency_metabolite | hub_compound | no_reaction_in_table


def link_confidence(match_basis: str, n_compounds_for_name: int, n_genes_for_compound: int) -> str:
    """Mapping confidence of a metabolite -> gene link (never 'exact' or one-to-one).

    A curator's judgement call on the name, or a name that resolves to more
    than one compound (e.g. 'lactate' -> L- and D-lactate), is the weakest
    rung, `inferred`. Otherwise a compound reaching several genes is
    `one_to_many_ortholog`, and one reaching a single gene is still `inferred`:
    the enzyme assignment is KEGG orthology transfer, not oyster biochemistry.
    """
    if match_basis == "curator_judgement" or n_compounds_for_name > 1:
        return "inferred"
    if n_genes_for_compound > 1:
        return "one_to_many_ortholog"
    return "inferred"


@dataclass
class MetaboliteGeneLinks:
    """A loaded link table, indexed in both directions."""

    path: Path | None
    status: str
    sha256: str | None = None
    source: str | None = None
    links: list = field(default_factory=list)
    exclusions: list = field(default_factory=list)

    def __post_init__(self):
        self._by_metabolite = defaultdict(list)
        self._by_gene = defaultdict(list)
        self._exclusions_by_metabolite = defaultdict(list)
        for link in self.links:
            self._by_metabolite[link.metabolite].append(link)
            self._by_gene[link.ncbi_gene_id].append(link)
        for exclusion in self.exclusions:
            self._exclusions_by_metabolite[exclusion.metabolite].append(exclusion)

    @property
    def available(self) -> bool:
        return self.status == STATUS_AVAILABLE

    def for_metabolite(self, name) -> list:
        return self._by_metabolite.get(normalize_metabolite_name(name), [])

    def for_gene(self, gene_id) -> list:
        return self._by_gene.get(str(gene_id), [])

    def excluded_for_metabolite(self, name) -> list:
        return self._exclusions_by_metabolite.get(normalize_metabolite_name(name), [])

    def links_for_feature(self, feature_id, is_metabolite: bool) -> list:
        """Links touching this standardized feature, from whichever side it is on."""
        return self.for_metabolite(feature_id) if is_metabolite else self.for_gene(feature_id)

    def label(self) -> str:
        if not self.available:
            return self.status
        rel = self.path.relative_to(REPO_ROOT) if self.path.is_relative_to(REPO_ROOT) else self.path
        return f"{rel} (sha256 {self.sha256[:12]}…, {self.source})"


def build_links(
    link_rows: list,
    compound_map: dict,
    currency: dict,
    max_genes_per_compound: int = MAX_GENES_PER_COMPOUND,
) -> tuple:
    """Expand link-table rows into per-metabolite links plus visible exclusions.

    Returns (links, exclusions). Pure: tests build tables from plain rows.
    """
    rows_by_compound = defaultdict(list)
    for row in link_rows:
        rows_by_compound[row["kegg_compound"]].append(row)
    genes_by_compound = {
        cpd: {str(r["ncbi_gene_id"]) for r in rows} for cpd, rows in rows_by_compound.items()
    }

    links, exclusions = [], []
    for metabolite, compounds in compound_map.items():
        n_compounds = len(compounds)
        for compound in compounds:
            cpd = compound["kegg_compound"]
            if cpd in currency:
                exclusions.append(CompoundExclusion(metabolite, cpd, "currency_metabolite"))
                continue
            rows = rows_by_compound.get(cpd)
            if not rows:
                exclusions.append(CompoundExclusion(metabolite, cpd, "no_reaction_in_table"))
                continue
            n_genes = len(genes_by_compound[cpd])
            hub = n_genes > max_genes_per_compound
            if hub:
                exclusions.append(CompoundExclusion(metabolite, cpd, "hub_compound"))
            confidence = link_confidence(compound["match_basis"], n_compounds, n_genes)
            seen = set()
            for row in rows:
                key = (row["kegg_reaction"], row["ko"], str(row["ncbi_gene_id"]))
                if key in seen:
                    continue
                seen.add(key)
                links.append(MetaboliteGeneLink(
                    metabolite=metabolite,
                    kegg_compound=cpd,
                    compound_match_basis=compound["match_basis"],
                    kegg_reaction=row["kegg_reaction"],
                    equation_side=row["equation_side"],
                    ko=row["ko"],
                    ncbi_gene_id=str(row["ncbi_gene_id"]),
                    link_type=row.get("link_type") or "enzyme",
                    n_genes_for_compound=n_genes,
                    link_confidence=confidence,
                    credited=not hub,
                    exclusion_reason="hub_compound" if hub else None,
                ))
    return links, exclusions


def provenance_path(table_path: Path) -> Path:
    return table_path.with_suffix(".provenance.json")


def load_link_table(path: Path) -> MetaboliteGeneLinks:
    """Load and index a built link table, or report why there is none."""
    if not path.exists():
        return MetaboliteGeneLinks(path=path, status=STATUS_NOT_BUILT)
    provenance = load_json(provenance_path(path)) if provenance_path(path).exists() else {}
    links, exclusions = build_links(_read_tsv(path), load_compound_map(), load_currency_compounds())
    return MetaboliteGeneLinks(
        path=path,
        status=STATUS_AVAILABLE,
        sha256=sha256_file(path),
        source=provenance.get("source_summary", "no provenance sidecar"),
        links=links,
        exclusions=exclusions,
    )


def _demo_taxid() -> str | None:
    sidecar = provenance_path(DEMO_LINKS_PATH)
    return str(load_json(sidecar)["ncbi_taxid"]) if sidecar.exists() else None


def registered_link_path(taxid) -> Path | None:
    """The real link table registered for any assembly of this taxid, or None."""
    for record in load_assemblies().values():
        if str(record.get("ncbi_taxid")) != str(taxid):
            continue
        value = record.get("metabolite_gene_links")
        if value:
            path = Path(value)
            return path if path.is_absolute() else REPO_ROOT / path
    return None


@lru_cache(maxsize=16)
def links_for_partition(simulated: str, taxid: str) -> MetaboliteGeneLinks:
    """The link table for one evidence partition (`_simulated`, `_species_taxid`).

    Simulated evidence only ever sees the simulated demo table, and only for the
    taxid it was fabricated for. Real evidence only ever sees a real table.
    """
    if simulated == "true":
        if taxid == _demo_taxid():
            return load_link_table(DEMO_LINKS_PATH)
        return MetaboliteGeneLinks(path=None, status=STATUS_NOT_REGISTERED)
    path = registered_link_path(taxid)
    if path is None:
        return MetaboliteGeneLinks(path=None, status=STATUS_NOT_REGISTERED)
    return load_link_table(path)


def clear_caches() -> None:
    """Forget loaded tables (after a rebuild, or in tests that swap files)."""
    load_compound_map.cache_clear()
    load_currency_compounds.cache_clear()
    links_for_partition.cache_clear()
