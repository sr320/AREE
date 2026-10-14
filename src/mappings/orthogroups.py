"""Orthogroups across the registered species, from OrthoFinder on RefSeq proteomes.

WP6 of docs/mussel_clam_integration_plan.md. An orthogroup run is inferred once
for every species with an identifier crosswalk, from the one-protein-per-gene
proteomes that `aree fetch-reference --proteins` writes, so every member of
every orthogroup is an NCBI GeneID the crosswalks already resolve.

A run is committed as a directory under data/reference/orthogroups/<run_id>/:

* ``gene_orthogroups.tsv`` - ncbi_gene_id, ncbi_taxid, orthogroup_id;
* ``orthogroups.tsv`` - per orthogroup: copy class, number of species and
  genes, and the copy number in each species;
* ``provenance.json`` - OrthoFinder version and command, every proteome's
  assembly and SHA-256, and per-species assignment rates.

``data/reference/orthogroups/active.yaml`` names the run that harmonization
uses. Crosswalks are not rewritten: they stay NCBI-derived, and an orthogroup
run is a separate inference with its own provenance.

Copy classes matter for cross-species pooling. Bivalve genomes carry large
lineage-specific expansions (HSP70, C1q, IAPs), and pooling "the" ortholog of a
gene with 40 paralogs in another species is not a like-for-like comparison:

* ``single_copy_all``  - exactly one gene in every species of the run;
* ``single_copy``      - at most one gene per species, in two or more species;
* ``multi_copy``       - two or more genes in at least one species;
* ``species_specific`` - genes from one species only.
"""
from __future__ import annotations

import json
import os
import re
import shutil
import subprocess
from datetime import date
from functools import lru_cache
from pathlib import Path

from common import DATA_DIR, assembly_crosswalk_path, load_assemblies, load_yaml, sha256_file

ORTHOGROUPS_DIR = DATA_DIR / "reference" / "orthogroups"
ACTIVE_PATH = ORTHOGROUPS_DIR / "active.yaml"
BUILDER_VERSION = "1.0.0"


class OrthogroupError(RuntimeError):
    """An orthogroup run could not be built, read or trusted."""


# --------------------------------------------------------------------------- #
# Inputs
# --------------------------------------------------------------------------- #


def run_species() -> list[dict]:
    """One entry per registered species with a crosswalk: its NCBI reference assembly.

    The slug is the crosswalk's (mgigas, mgall, ...), so OrthoFinder's species
    names match the crosswalk files.
    """
    species = []
    for record in load_assemblies().values():
        crosswalk = assembly_crosswalk_path(record)
        if not record.get("is_ncbi_reference") or crosswalk is None:
            continue
        species.append({
            "slug": crosswalk.name.replace("_gene_id_crosswalk.tsv", ""),
            "ncbi_taxid": int(record["ncbi_taxid"]),
            "species": record["species"],
            "assembly_id": record["assembly_id"],
            "ncbi_assembly_accession": record["ncbi_assembly_accession"],
            "annotation_release": record.get("annotation_release"),
            "proteome": DATA_DIR / "reference" / record["ncbi_assembly_accession"] / "primary_proteins.faa",
        })
    return sorted(species, key=lambda s: s["slug"])


def stage_proteomes(species: list[dict], dest: Path) -> list[dict]:
    """Copy each species' primary proteome to <dest>/<slug>.fa and record its checksum."""
    dest.mkdir(parents=True, exist_ok=True)
    staged = []
    for sp in species:
        src = Path(sp["proteome"])
        if not src.exists():
            raise OrthogroupError(
                f"{src} is missing; fetch it with "
                f"`aree fetch-reference --assembly {sp['ncbi_assembly_accession']} --proteins`")
        target = dest / f"{sp['slug']}.fa"
        shutil.copyfile(src, target)
        with open(target) as fh:
            n = sum(1 for line in fh if line.startswith(">"))
        staged.append({**{k: v for k, v in sp.items() if k != "proteome"},
                       "proteome_file": target.name, "proteome_sha256": sha256_file(target),
                       "proteins": n})
    return staged


# --------------------------------------------------------------------------- #
# OrthoFinder
# --------------------------------------------------------------------------- #


def _env_for(binary: str) -> dict:
    """OrthoFinder calls DIAMOND, MCL and FastME by name: put its own bin first."""
    env = dict(os.environ)
    bin_dir = str(Path(binary).expanduser().resolve().parent)
    env["PATH"] = bin_dir + os.pathsep + env.get("PATH", "")
    return env


def orthofinder_version(binary: str) -> str:
    out = subprocess.run([binary, "-v"], capture_output=True, text=True, env=_env_for(binary))
    text = re.sub(r"\x1b\[[0-9;]*m|\x1b\]8;[^\x1b]*\x1b\\\\", "", out.stdout + out.stderr)
    match = re.search(r"OrthoFinder:?\s*v?(\d+\.\d+\.\d+)", text)
    if not match:
        raise OrthogroupError(f"could not read the OrthoFinder version from `{binary} -v`")
    return f"OrthoFinder {match.group(1)}"


def run_orthofinder(proteome_dir: Path, binary: str, threads: int) -> tuple[Path, list[str]]:
    """Run OrthoFinder; returns (results dir, command).

    OrthoFinder 3 has no orthogroups-only mode, so gene trees are built too;
    `-M dendroblast` builds them from the DIAMOND scores instead of a multiple
    alignment per orthogroup, which on an 8 GB workstation is the difference
    between hours and days. Only Orthogroups/Orthogroups.tsv is used here.
    """
    command = [binary, "-f", str(proteome_dir), "-t", str(threads), "-a", str(threads),
               "-M", "dendroblast"]
    proc = subprocess.run(command, capture_output=True, text=True, env=_env_for(binary))
    if proc.returncode != 0:
        raise OrthogroupError(f"OrthoFinder failed ({proc.returncode}): {proc.stderr[-2000:]}")
    results = sorted((proteome_dir / "OrthoFinder").glob("Results_*"))
    if not results:
        raise OrthogroupError(f"no OrthoFinder results under {proteome_dir / 'OrthoFinder'}")
    return results[-1], command


# --------------------------------------------------------------------------- #
# Parsing and summaries
# --------------------------------------------------------------------------- #


def _read_orthogroup_table(path: Path, slug_taxid: dict[str, int]) -> list[tuple[str, int, str]]:
    rows = []
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        unknown = set(header[1:]) - set(slug_taxid)
        if unknown:
            raise OrthogroupError(f"{path}: species columns {sorted(unknown)} are not in this run")
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            og = fields[0]
            for slug, cell in zip(header[1:], fields[1:]):
                for gene in filter(None, (g.strip() for g in cell.split(","))):
                    rows.append((gene, slug_taxid[slug], og))
    return rows


def parse_orthogroups(results_dir: Path, slug_taxid: dict[str, int]) -> list[tuple[str, int, str]]:
    """(GeneID, taxid, orthogroup) for every gene, including OrthoFinder's unassigned
    genes, which it gives their own single-gene orthogroup IDs."""
    og_dir = Path(results_dir) / "Orthogroups"
    rows = _read_orthogroup_table(og_dir / "Orthogroups.tsv", slug_taxid)
    unassigned = og_dir / "Orthogroups_UnassignedGenes.tsv"
    if unassigned.exists():
        rows += _read_orthogroup_table(unassigned, slug_taxid)
    seen: dict[tuple[str, int], str] = {}
    for gene, taxid, og in rows:
        if (gene, taxid) in seen:
            raise OrthogroupError(f"GeneID {gene} (taxid {taxid}) is in {seen[(gene, taxid)]} and {og}")
        seen[(gene, taxid)] = og
    return sorted(rows, key=lambda r: (r[2], r[1], int(r[0]) if r[0].isdigit() else r[0]))


def copy_class(copies: dict[int, int], n_run_species: int) -> str:
    present = {t: n for t, n in copies.items() if n > 0}
    if len(present) == 1:
        return "species_specific"
    if any(n > 1 for n in present.values()):
        return "multi_copy"
    return "single_copy_all" if len(present) == n_run_species else "single_copy"


def summarize(rows: list[tuple[str, int, str]], taxids: list[int]) -> list[dict]:
    copies: dict[str, dict[int, int]] = {}
    for _, taxid, og in rows:
        copies.setdefault(og, dict.fromkeys(taxids, 0))[taxid] += 1
    out = []
    for og in sorted(copies):
        c = copies[og]
        out.append({"orthogroup_id": og, "copy_class": copy_class(c, len(taxids)),
                    "n_species": sum(1 for n in c.values() if n), "n_genes": sum(c.values()),
                    **{f"copies_{t}": c[t] for t in taxids}})
    return out


# --------------------------------------------------------------------------- #
# Writing a run
# --------------------------------------------------------------------------- #


def write_run(run_dir: Path, rows, summary, staged: list[dict], *, version: str,
              command: list[str] | None, results_dir: Path, today: str | None = None) -> dict:
    run_dir.mkdir(parents=True, exist_ok=True)
    genes_path = run_dir / "gene_orthogroups.tsv"
    with open(genes_path, "w") as fh:
        fh.write("ncbi_gene_id\tncbi_taxid\torthogroup_id\n")
        for gene, taxid, og in rows:
            fh.write(f"{gene}\t{taxid}\t{og}\n")
    taxids = [s["ncbi_taxid"] for s in staged]
    og_path = run_dir / "orthogroups.tsv"
    cols = ["orthogroup_id", "copy_class", "n_species", "n_genes", *[f"copies_{t}" for t in taxids]]
    with open(og_path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for row in summary:
            fh.write("\t".join(str(row[c]) for c in cols) + "\n")

    assigned = {}
    multi_species_ogs = {r["orthogroup_id"] for r in summary if r["n_species"] > 1}
    for sp in staged:
        genes = [r for r in rows if r[1] == sp["ncbi_taxid"]]
        assigned[sp["slug"]] = {
            "genes": len(genes),
            "in_multi_species_orthogroups": sum(1 for r in genes if r[2] in multi_species_ogs),
        }
    classes: dict[str, int] = {}
    for row in summary:
        classes[row["copy_class"]] = classes.get(row["copy_class"], 0) + 1

    provenance = {
        "run_id": run_dir.name,
        "builder": "src/mappings/orthogroups.py",
        "builder_version": BUILDER_VERSION,
        "date_generated": today or date.today().isoformat(),
        "method": version,
        "command": command,
        "results_source": str(results_dir),
        "inference_scope": "OrthoFinder orthogroups (Orthogroups/Orthogroups.tsv): MCL clusters "
                           "of the all-vs-all DIAMOND graph. Gene trees (dendroblast) and "
                           "hierarchical orthogroups were also produced by the run but are not used.",
        "species": staged,
        "per_species": assigned,
        "orthogroups": {"total": len(summary), "by_copy_class": dict(sorted(classes.items()))},
        "outputs": {genes_path.name: sha256_file(genes_path), og_path.name: sha256_file(og_path)},
    }
    (run_dir / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    return provenance


def activate(run_id: str, path: Path | None = None) -> None:
    path = path or ACTIVE_PATH
    path.write_text(
        "# The orthogroup run harmonization attaches to evidence records.\n"
        "# See src/mappings/orthogroups.py and docs/identifier_mapping.md.\n"
        f"run: {run_id}\n")
    gene_orthogroup.cache_clear()


# --------------------------------------------------------------------------- #
# Lookup (used by harmonization)
# --------------------------------------------------------------------------- #


def active_run_dir(path: Path | None = None) -> Path | None:
    """The active run's directory, read at call time (not import time)."""
    path = path or ACTIVE_PATH
    if not path.exists():
        return None
    run = (load_yaml(path) or {}).get("run")
    return ORTHOGROUPS_DIR / run if run else None


@lru_cache(maxsize=2)
def _gene_table(run_dir_str: str) -> dict[str, str]:
    table = {}
    with open(Path(run_dir_str) / "gene_orthogroups.tsv") as fh:
        next(fh)
        for line in fh:
            gene, _, og = line.rstrip("\n").split("\t")
            table[gene] = og
    return table


def gene_orthogroup(gene_id: str) -> str | None:
    """The active run's orthogroup for an NCBI GeneID, or None.

    NCBI GeneIDs are unique across taxa (tests/test_multi_species.py), so the
    GeneID alone identifies the gene.
    """
    run_dir = active_run_dir()
    if run_dir is None or gene_id is None:
        return None
    return _gene_table(str(run_dir)).get(str(gene_id))


gene_orthogroup.cache_clear = _gene_table.cache_clear  # type: ignore[attr-defined]


def active_run_ref() -> dict | None:
    """Run ID and checksum of the active run, for harmonization manifests."""
    run_dir = active_run_dir()
    if run_dir is None:
        return None
    return {"run": run_dir.name, "gene_orthogroups_sha256": sha256_file(run_dir / "gene_orthogroups.tsv")}
