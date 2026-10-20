"""Build the real C. gigas metabolite -> gene link table from the KEGG REST API.

    aree build-metabolite-links            # fetch (or reuse cached) KEGG responses
    aree build-metabolite-links --offline  # rebuild from the cached responses only

Only the compounds named in data/mappings/metabolite_compound_map.tsv are
queried — never a mirror of KEGG — and currency compounds
(data/mappings/currency_metabolites.yaml) are not queried at all. The chain is:

    link/ko/crg              every C. gigas gene -> KO assignment (one call)
    link/reaction/cpd:C...   reactions each curated compound takes part in
    link/ko/R...+R...        KOs annotated to catalyse those reactions
    get/R...+R...            reaction equations, for the compound's side
    get/cpd:C...+C...        compound names
    get/br:ko02000           KEGG BRITE Transporters: which KOs are transporters,
                             and their family

Transporter links come from data/mappings/transporter_substrates.tsv, a curated
KO -> compound table: ko02000 classifies transporters but records no substrate
compounds. Each curated KO must be in ko02000; one that is not is rejected and
listed in the provenance.

KEGG gene IDs for `crg` are NCBI GeneIDs. Each one is checked against the
species' AREE crosswalk; a retired GeneID is remapped through the
*_retired_gene_ids.tsv sidecar and counted, and one that resolves to nothing is
dropped and listed in the provenance, never silently.

The output (metabolite_gene_links_mgigas.tsv + .provenance.json) is NOT
committed: KEGG REST is provided for academic use and KEGG does not clearly
grant redistribution of derived tables. Every raw response is cached under
data/mappings/_sources/kegg/ (also git-ignored) with its sha256 recorded, so a
build can be reproduced offline and audited. See docs/metabolite_gene_links.md.
"""
from __future__ import annotations

import csv
import hashlib
import json
import re
import urllib.request
from collections import defaultdict
from datetime import date
from pathlib import Path
from typing import Callable

from common import MAPPINGS_DIR, REPO_ROOT, load_assemblies, sha256_file

from .metabolite_links import (
    COMPOUND_MAP_PATH,
    CURRENCY_PATH,
    LINK_COLUMNS,
    LINK_TYPE_ENZYME,
    LINK_TYPE_TRANSPORTER,
    MAX_GENES_PER_COMPOUND,
    TRANSPORTER_SUBSTRATES_PATH,
    load_compound_map,
    load_currency_compounds,
    load_transporter_substrates,
    provenance_path,
)

BUILDER_VERSION = "1.1.0"
KEGG_REST = "https://rest.kegg.jp"
KEGG_ORGANISM = {29159: "crg"}
DEFAULT_SLUG = {29159: "mgigas"}
SOURCES_DIR = MAPPINGS_DIR / "_sources" / "kegg"
# KEGG REST accepts at most 10 entries per get/link call.
BATCH = 10

Fetcher = Callable[[str], str]


def output_path(taxid: int) -> Path:
    return MAPPINGS_DIR / f"metabolite_gene_links_{DEFAULT_SLUG[taxid]}.tsv"


class CachedKegg:
    """Fetch KEGG REST operations, caching each raw response and its sha256."""

    def __init__(self, cache_dir: Path = SOURCES_DIR, offline: bool = False, fetch: Fetcher | None = None):
        self.cache_dir = cache_dir
        self.offline = offline
        self._fetch = fetch or self._http_get
        self.log = []

    @staticmethod
    def _http_get(operation: str) -> str:
        with urllib.request.urlopen(f"{KEGG_REST}/{operation}", timeout=60) as response:  # noqa: S310
            return response.read().decode("utf-8")

    def get(self, operation: str) -> str:
        path = self.cache_dir / (re.sub(r"[^A-Za-z0-9_.+-]+", "_", operation) + ".txt")
        if path.exists():
            text = path.read_text()
        elif self.offline:
            raise FileNotFoundError(f"--offline, and no cached KEGG response for {operation!r} at {path}")
        else:
            text = self._fetch(operation)
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(text)
        self.log.append({
            "url": f"{KEGG_REST}/{operation}",
            "local_copy": str(path.relative_to(REPO_ROOT)) if path.is_relative_to(REPO_ROOT) else str(path),
            "sha256": hashlib.sha256(text.encode("utf-8")).hexdigest(),
        })
        return text


def _batches(items: list, size: int = BATCH):
    for i in range(0, len(items), size):
        yield items[i:i + size]


def _pairs(text: str) -> list:
    """Parse a KEGG `link` response into (source, target) pairs without db prefixes."""
    out = []
    for line in text.splitlines():
        if "\t" in line:
            a, b = line.split("\t", 1)
            out.append((a.split(":", 1)[-1].strip(), b.split(":", 1)[-1].strip()))
    return out


def parse_entries(text: str, fields: tuple) -> dict:
    """Parse a flat-file `get` response into {ENTRY id: {field: first-line value}}."""
    entries, current = {}, None
    for line in text.splitlines():
        if line.startswith("ENTRY"):
            current = line.split()[1]
            entries[current] = {}
        elif current and line[:12].strip() in fields and line[:12].strip() not in entries[current]:
            entries[current][line[:12].strip()] = line[12:].strip()
    return entries


def parse_brite(text: str) -> dict:
    """Parse a KEGG BRITE hierarchy (`get/br:ko02000`) into {KO: family path}.

    A KO listed under several families gets every path, joined by '; '. The
    path is the chain of headings above the KO line, e.g.
    'Solute carrier family (SLC) > SLC16: Monocarboxylate transporter'.
    """
    headings, families = {}, defaultdict(list)
    for line in text.splitlines():
        if not line or line[0] not in "ABCDEFGH":
            continue
        level, body = line[0], re.sub(r"<[^>]+>", "", line[1:]).strip()
        match = re.match(r"(K\d{5})\s", body)
        if match:
            path = " > ".join(headings[k] for k in sorted(headings) if k < level and headings[k])
            if path not in families[match.group(1)]:
                families[match.group(1)].append(path)
            continue
        headings[level] = body
        for deeper in [k for k in headings if k > level]:
            del headings[deeper]
    return {ko: "; ".join(paths) for ko, paths in families.items()}


def equation_side(equation: str, compound: str) -> str | None:
    """'left' / 'right' / 'both' for a compound in a KEGG equation ('A + B <=> C')."""
    if "<=>" not in equation:
        return None
    left, right = equation.split("<=>", 1)
    on_left = re.search(rf"\b{compound}\b", left) is not None
    on_right = re.search(rf"\b{compound}\b", right) is not None
    if on_left and on_right:
        return "both"
    return "left" if on_left else "right" if on_right else None


def load_gene_resolver(taxid: int) -> tuple:
    """(current GeneIDs, retired -> current map, crosswalk path) for a taxid."""
    crosswalk = next(
        (a["identifier_crosswalk"] for a in load_assemblies().values()
         if a.get("ncbi_taxid") == taxid and a.get("is_ncbi_reference") and a.get("identifier_crosswalk")),
        None,
    )
    if crosswalk is None:
        raise ValueError(f"no identifier crosswalk registered for the reference assembly of taxid {taxid}")
    path = REPO_ROOT / crosswalk
    with open(path) as fh:
        current = {r["ncbi_gene_id"] for r in csv.DictReader(
            (line for line in fh if not line.startswith("#")), delimiter="\t")}
    retired = {}
    retired_path = path.with_name(path.name.replace("_gene_id_crosswalk.tsv", "_retired_gene_ids.tsv"))
    if retired_path.exists():
        with open(retired_path) as fh:
            for r in csv.DictReader((line for line in fh if not line.startswith("#")), delimiter="\t"):
                if r.get("current_gene_id") and r["current_gene_id"] != "-":
                    retired[r["discontinued_gene_id"]] = r["current_gene_id"]
    return current, retired, path


def build(
    taxid: int = 29159,
    out_path: Path | None = None,
    offline: bool = False,
    fetch: Fetcher | None = None,
    cache_dir: Path = SOURCES_DIR,
    gene_resolver: tuple | None = None,
) -> Path:
    """Build the link table and its provenance sidecar; return the table path."""
    if taxid not in KEGG_ORGANISM:
        raise ValueError(f"no KEGG organism code recorded for taxid {taxid} (known: {sorted(KEGG_ORGANISM)})")
    org = KEGG_ORGANISM[taxid]
    out_path = out_path or output_path(taxid)
    kegg = CachedKegg(cache_dir=cache_dir, offline=offline, fetch=fetch)
    current, retired, crosswalk_path = gene_resolver or load_gene_resolver(taxid)

    release = kegg.get("info/kegg").strip().splitlines()
    genes_by_ko = defaultdict(set)
    raw_genes = set()
    for gene, ko in _pairs(kegg.get(f"link/ko/{org}")):
        raw_genes.add(gene)
        genes_by_ko[ko].add(gene)

    resolved, n_retired_remapped, unresolved = {}, 0, []
    for gene in sorted(raw_genes):
        if gene in current:
            resolved[gene] = gene
        elif gene in retired and retired[gene] in current:
            resolved[gene] = retired[gene]
            n_retired_remapped += 1
        else:
            unresolved.append(gene)

    currency = load_currency_compounds(CURRENCY_PATH)
    compound_map = load_compound_map(COMPOUND_MAP_PATH)
    curated = sorted({c["kegg_compound"] for cs in compound_map.values() for c in cs})
    queried = [c for c in curated if c not in currency]

    reactions_by_compound = {
        cpd: sorted({r for _, r in _pairs(kegg.get(f"link/reaction/cpd:{cpd}"))}) for cpd in queried
    }
    all_reactions = sorted({r for rs in reactions_by_compound.values() for r in rs})
    kos_by_reaction = defaultdict(set)
    equations = {}
    for batch in _batches(all_reactions):
        for reaction, ko in _pairs(kegg.get("link/ko/" + "+".join(f"rn:{r}" for r in batch))):
            kos_by_reaction[reaction].add(ko)
        for entry, fields in parse_entries(kegg.get("get/" + "+".join(f"rn:{r}" for r in batch)),
                                           ("EQUATION",)).items():
            equations[entry] = fields.get("EQUATION", "")
    names = {}
    for batch in _batches(queried):
        for entry, fields in parse_entries(kegg.get("get/" + "+".join(f"cpd:{c}" for c in batch)),
                                           ("NAME",)).items():
            names[entry] = fields.get("NAME", "").rstrip(";")

    rows = []
    for cpd in queried:
        for reaction in reactions_by_compound[cpd]:
            side = equation_side(equations.get(reaction, ""), cpd)
            for ko in sorted(kos_by_reaction.get(reaction, ())):
                for gene in sorted(genes_by_ko.get(ko, ())):
                    if gene not in resolved:
                        continue
                    rows.append({
                        "kegg_compound": cpd, "compound_name": names.get(cpd, ""),
                        "kegg_reaction": reaction, "equation_side": side or "unknown", "ko": ko,
                        "ncbi_gene_id": resolved[gene], "link_type": LINK_TYPE_ENZYME,
                        "link_basis": "kegg_reaction", "brite_family": "",
                    })

    # Transporter links: curated KO -> compound, KO verified in BRITE ko02000.
    transporter_families = parse_brite(kegg.get("get/br:ko02000"))
    transporter_rows = load_transporter_substrates(TRANSPORTER_SUBSTRATES_PATH)
    rejected, without_genes, used_kos = [], set(), set()
    for entry in transporter_rows:
        cpd, ko = entry["kegg_compound"], entry["ko"]
        if cpd not in queried:
            continue
        if ko not in transporter_families:
            rejected.append({"ko": ko, "kegg_compound": cpd, "reason": "not in KEGG BRITE ko02000"})
            continue
        used_kos.add(ko)
        genes = [g for g in sorted(genes_by_ko.get(ko, ())) if g in resolved]
        if not genes:
            without_genes.add(ko)
        for gene in genes:
            rows.append({
                "kegg_compound": cpd, "compound_name": names.get(cpd, ""),
                "kegg_reaction": "", "equation_side": "", "ko": ko,
                "ncbi_gene_id": resolved[gene], "link_type": LINK_TYPE_TRANSPORTER,
                "link_basis": entry["substrate_basis"], "brite_family": transporter_families[ko],
            })
    rows = [dict(t) for t in sorted({tuple(r.items()) for r in rows},
                                    key=lambda t: tuple(v for _, v in t))]

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as fh:
        fh.write(f"# AREE metabolite -> gene link table, taxid {taxid} (KEGG organism {org}).\n")
        fh.write(f"# Built by src/mappings/build_metabolite_links.py v{BUILDER_VERSION} on "
                 f"{date.today().isoformat()} from the KEGG REST API. NOT for redistribution; see\n")
        fh.write("# the .provenance.json sidecar and docs/metabolite_gene_links.md.\n")
        fh.write("\t".join(LINK_COLUMNS) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in LINK_COLUMNS) + "\n")

    genes_per_compound = defaultdict(set)
    genes_per_type = defaultdict(set)
    for r in rows:
        genes_per_compound[r["kegg_compound"]].add(r["ncbi_gene_id"])
        genes_per_type[(r["kegg_compound"], r["link_type"])].add(r["ncbi_gene_id"])
    release_line = next((line.strip() for line in release if line.strip().startswith("kegg")), release[0])
    ko_line = next((" ".join(line.split()) for line in release if line.strip().startswith("ko ")), "")
    provenance = {
        "artifact": out_path.name,
        "artifact_sha256": sha256_file(out_path),
        "simulated": False,
        "builder": "src/mappings/build_metabolite_links.py",
        "builder_version": BUILDER_VERSION,
        "date_generated": date.today().isoformat(),
        "ncbi_taxid": taxid,
        "kegg_organism": org,
        "kegg_release": {"summary": release_line, "ko": ko_line},
        "source_summary": f"KEGG REST {org}, accessed {date.today().isoformat()} ({ko_line or release_line})",
        "licence_note": ("KEGG REST is provided for academic use by academic users. This derived table "
                         "is not committed or redistributed by AREE; rebuild it locally."),
        "identifier_crosswalk": str(crosswalk_path.relative_to(REPO_ROOT))
        if crosswalk_path.is_relative_to(REPO_ROOT) else str(crosswalk_path),
        "compound_map": {"path": str(COMPOUND_MAP_PATH.relative_to(REPO_ROOT)),
                         "sha256": sha256_file(COMPOUND_MAP_PATH)},
        "transporter_substrates": {"path": str(TRANSPORTER_SUBSTRATES_PATH.relative_to(REPO_ROOT)),
                                   "sha256": sha256_file(TRANSPORTER_SUBSTRATES_PATH),
                                   "rejected_not_in_ko02000": rejected,
                                   "kos_without_crg_genes": sorted(without_genes),
                                   "kos_used": sorted(used_kos)},
        "currency_metabolites": {"path": str(CURRENCY_PATH.relative_to(REPO_ROOT)),
                                 "sha256": sha256_file(CURRENCY_PATH),
                                 "skipped": sorted(c for c in curated if c in currency)},
        "max_genes_per_compound": MAX_GENES_PER_COMPOUND,
        "sources": kegg.log,
        "coverage": {
            "kegg_genes_with_ko": len(raw_genes),
            "kegg_genes_current": len(raw_genes) - n_retired_remapped - len(unresolved),
            "kegg_genes_retired_remapped": n_retired_remapped,
            "kegg_genes_unresolved": len(unresolved),
            "compounds_curated": len(curated),
            "compounds_queried": len(queried),
            "compounds_with_links": len(genes_per_compound),
            "hub_compounds": sorted(f"{c}:{t}" for (c, t), g in genes_per_type.items()
                                    if len(g) > MAX_GENES_PER_COMPOUND),
            "genes_per_compound": {c: len(g) for c, g in sorted(genes_per_compound.items())},
            "genes_per_compound_by_link_type": {
                f"{c}:{t}": len(g) for (c, t), g in sorted(genes_per_type.items())},
            "rows": len(rows),
            "rows_by_link_type": {t: sum(1 for r in rows if r["link_type"] == t)
                                  for t in (LINK_TYPE_ENZYME, LINK_TYPE_TRANSPORTER)},
        },
        "unresolved_kegg_gene_ids": unresolved,
        "known_limitations": [
            "KO assignments are KEGG's orthology transfer to C. gigas, not oyster biochemistry; every "
            "link is 'inferred' or 'one_to_many_ortholog', never 'exact'.",
            "KEGG equations are written reversibly; equation_side is not substrate/product.",
            "Transporter substrates are curated (transporter_substrates.tsv); KEGG BRITE ko02000 "
            "confirms only that a KO is a transporter and its family, not what it carries. Every "
            "transporter link is 'inferred'.",
            "ATP-driven pumps (e.g. SERCA, K05853; Na+/K+-ATPase) are not linked: their KEGG-named "
            "substrates are ions or ATP, which are not curated metabolites or are currency compounds.",
            "Only compounds in metabolite_compound_map.tsv are queried; add names there first.",
        ],
    }
    provenance_path(out_path).write_text(json.dumps(provenance, indent=2) + "\n")
    return out_path
