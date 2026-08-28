"""Build a real, provenance-tracked gene identifier crosswalk from public reference sources.

Sources
-------
1. NCBI ``gene_info`` (Gene database), filtered to a single NCBI taxonomy id.
   Provides: GeneID, Symbol, Synonyms, Ensembl xrefs, description, gene type.
2. UniProtKB (REST), filtered to the same organism.
   Provides: accession, reviewed status, cross-reference back to NCBI GeneID.

Output
------
A one-row-per-NCBI-GeneID TSV plus a JSON provenance sidecar recording source
URLs, checksums, retrieval date, row counts, and per-column coverage.

Design decisions that matter for interpretation
-----------------------------------------------
* **One row per NCBI GeneID.** GeneID is the anchor of the crosswalk because it
  is the most stable identifier across annotation releases. Emitting one row per
  (gene, protein) pair would make a plain GeneID lookup return multiple rows,
  which ``harmonize.identifiers`` would (correctly, but misleadingly) downgrade
  to ``many_to_one_ortholog``.
* **UniProt is one-to-many.** A single gene often has several TrEMBL accessions.
  ``uniprot_accession`` holds a single representative (reviewed entries preferred);
  ``uniprot_accessions_all`` holds the full semicolon-delimited set so that a
  lookup by *any* accession still resolves exactly.
* **``locus_id`` is populated as ``LOC<GeneID>`` for every gene**, not only for
  genes that currently lack an official symbol. The LOC form unambiguously denotes
  that GeneID by NCBI convention, and studies analysed against an older annotation
  frequently report the LOC form for a gene that has since been given a name. This
  is the single largest source of recoverable matches across annotation versions.
* **``orthogroup_id`` is deliberately left empty.** Real orthogroups require an
  ortholog inference run (e.g. OrthoFinder/OrthoDB). Inventing them here would
  manufacture cross-species confidence that does not exist. See
  ``docs/identifier_mapping.md``.
* **Retired GeneIDs are captured in a sidecar table.** Studies analysed against
  an older annotation report GeneIDs that NCBI has since discontinued and
  replaced. Without ``<slug>_retired_gene_ids.tsv`` those identifiers resolve as
  ``unresolved`` even though NCBI records an authoritative replacement.
"""
from __future__ import annotations

import argparse
import gzip
import json
import re
import shutil
import subprocess
import sys
from collections import defaultdict
from datetime import date
from pathlib import Path

from common import REPO_ROOT, sha256_file

CROSSWALK_DIR = REPO_ROOT / "data" / "reference" / "crosswalk"
SOURCES_DIR = CROSSWALK_DIR / "_sources"

NCBI_GENE_INFO_URL = (
    "https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Invertebrates/All_Invertebrates.gene_info.gz"
)
NCBI_GENE_HISTORY_URL = "https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_history.gz"
UNIPROT_STREAM_URL = (
    "https://rest.uniprot.org/uniprotkb/stream"
    "?query=organism_id:{taxid}"
    "&fields=accession,reviewed,protein_name,gene_names,gene_primary,xref_geneid,xref_refseq"
    "&format=tsv&compressed=true"
)

BUILDER_VERSION = "1.0.0"

COLUMNS = [
    "ncbi_gene_id",
    "ensembl_gene_id",
    "uniprot_accession",
    "uniprot_accessions_all",
    "locus_id",
    "gene_symbol",
    "gene_synonyms",
    "gene_description",
    "gene_type",
    "orthogroup_id",
]

LOC_SYMBOL_RE = re.compile(r"^LOC\d+$")


def _blank(value: str) -> str:
    """NCBI uses '-' as its null token; normalise it (and whitespace) to ''."""
    value = (value or "").strip()
    return "" if value in {"-", "NA", "null"} else value


# --------------------------------------------------------------------------- #
# Source acquisition
# --------------------------------------------------------------------------- #


def download_gene_info(taxid: int, dest: Path) -> Path:
    """Stream the invertebrate gene_info archive and keep only rows for `taxid`.

    The upstream archive is ~230 MB compressed and covers every invertebrate, so
    it is filtered during streaming rather than stored. Records are grouped by
    tax_id upstream, so the scan stops once the block has been passed.
    """
    dest.parent.mkdir(parents=True, exist_ok=True)
    awk = (
        f'NR==1{{print; next}} $1=={taxid}{{print; seen=1; next}} seen{{exit}}'
    )
    cmd = f"curl -sS --fail --retry 3 {NCBI_GENE_INFO_URL!r} | gunzip -c | awk -F'\\t' {awk!r}"
    with open(dest, "w") as fh:
        proc = subprocess.run(cmd, shell=True, stdout=fh, stderr=subprocess.PIPE, text=True)
    # SIGPIPE (exit 141) is expected: awk exits early and closes the pipe.
    if proc.returncode not in (0, 141) and dest.stat().st_size == 0:
        raise RuntimeError(f"gene_info download failed: {proc.stderr.strip()}")
    return dest


def download_gene_history(taxid: int, dest: Path) -> Path:
    """Stream NCBI gene_history and keep only rows for `taxid` (same layout as gene_info)."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    awk = f'NR==1{{print; next}} $1=={taxid}{{print; seen=1; next}} seen{{exit}}'
    cmd = f"curl -sS --fail --retry 3 {NCBI_GENE_HISTORY_URL!r} | gunzip -c | awk -F'\\t' {awk!r}"
    with open(dest, "w") as fh:
        proc = subprocess.run(cmd, shell=True, stdout=fh, stderr=subprocess.PIPE, text=True)
    if proc.returncode not in (0, 141) and dest.stat().st_size == 0:
        raise RuntimeError(f"gene_history download failed: {proc.stderr.strip()}")
    return dest


def download_uniprot(taxid: int, dest: Path) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    url = UNIPROT_STREAM_URL.format(taxid=taxid)
    cmd = f"curl -sS --fail --retry 3 {url!r}"
    proc = subprocess.run(cmd, shell=True, capture_output=True, check=False)
    if proc.returncode != 0:
        raise RuntimeError(f"UniProt download failed: {proc.stderr.decode().strip()}")
    with open(dest, "wb") as fh:
        fh.write(gzip.decompress(proc.stdout))
    return dest


# --------------------------------------------------------------------------- #
# Parsing
# --------------------------------------------------------------------------- #


def parse_gene_info(path: Path, taxid: int) -> dict:
    """Parse a taxid-filtered gene_info file into {gene_id: record}."""
    genes = {}
    with open(path) as fh:
        header = fh.readline()
        if not header.startswith("#tax_id"):
            raise ValueError(f"{path} does not look like an NCBI gene_info file")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 10 or f[0] != str(taxid):
                continue
            gene_id = f[1].strip()
            symbol = _blank(f[2])
            synonyms = _blank(f[4])
            dbxrefs = _blank(f[5])

            ensembl = ""
            if dbxrefs:
                hits = [
                    x.split(":", 1)[1]
                    for x in dbxrefs.split("|")
                    if x.startswith("Ensembl:")
                ]
                ensembl = ";".join(sorted(set(hits)))

            genes[gene_id] = {
                "ncbi_gene_id": gene_id,
                "ensembl_gene_id": ensembl,
                "locus_id": f"LOC{gene_id}",
                "gene_symbol": symbol,
                "gene_synonyms": synonyms.replace("|", ";") if synonyms else "",
                "gene_description": _blank(f[8]),
                "gene_type": _blank(f[9]),
                "orthogroup_id": "",
            }
    return genes


def parse_uniprot(path: Path) -> dict:
    """Parse a UniProt TSV export into {gene_id: [accessions]}, reviewed entries first."""
    by_gene = defaultdict(list)
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        try:
            i_acc = header.index("Entry")
            i_rev = header.index("Reviewed")
            i_gid = header.index("GeneID")
        except ValueError as exc:  # pragma: no cover - guards a malformed export
            raise ValueError(f"unexpected UniProt column layout in {path}: {exc}") from exc

        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) <= i_gid:
                continue
            gene_ids = [g.strip() for g in f[i_gid].split(";") if g.strip()]
            if not gene_ids:
                continue
            reviewed = f[i_rev].strip().lower() == "reviewed"
            for gid in gene_ids:
                by_gene[gid].append((not reviewed, f[i_acc].strip()))

    # Sort so reviewed accessions lead; accession order is otherwise deterministic.
    return {gid: [acc for _, acc in sorted(set(pairs))] for gid, pairs in by_gene.items()}


def parse_gene_history(path: Path, taxid: int) -> tuple:
    """Parse NCBI gene_history into (replacements, dead).

    ``replacements`` maps a discontinued GeneID to the current GeneID that
    replaced it. ``dead`` lists discontinued GeneIDs with no replacement — NCBI
    writes '-' in the GeneID column for those, and they must stay unresolvable
    rather than being silently attached to some other gene.
    """
    replacements, dead = {}, []
    if not path or not Path(path).exists():
        return replacements, dead
    with open(path) as fh:
        header = fh.readline()
        if not header.startswith("#tax_id"):
            raise ValueError(f"{path} does not look like an NCBI gene_history file")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 4 or f[0] != str(taxid):
                continue
            current = _blank(f[1])
            discontinued = _blank(f[2])
            if not discontinued:
                continue
            if current:
                replacements[discontinued] = {
                    "discontinued_gene_id": discontinued,
                    "current_gene_id": current,
                    "discontinued_symbol": _blank(f[3]),
                    "discontinue_date": _blank(f[4]) if len(f) > 4 else "",
                }
            else:
                dead.append(discontinued)
    return replacements, dead


def write_retired_table(replacements: dict, dead: list, out_path: Path, taxid: int) -> None:
    """Write the retired-GeneID sidecar consulted by harmonize.identifiers."""
    cols = ["discontinued_gene_id", "current_gene_id", "discontinued_symbol", "discontinue_date"]
    with open(out_path, "w") as fh:
        fh.write(
            f"# Retired NCBI GeneIDs for taxid {taxid}, from NCBI gene_history.\n"
            f"# Built by src/mappings/build_crosswalk.py v{BUILDER_VERSION} on {date.today().isoformat()}\n"
            "# Studies analysed against an older annotation report these identifiers; each row\n"
            "# records the authoritative NCBI replacement. Resolution through this table is\n"
            "# labelled 'inferred', not 'exact', because it crosses an annotation version.\n"
            f"# {len(dead)} further GeneIDs were discontinued with NO replacement and are\n"
            "# deliberately absent, so that they stay 'unresolved'.\n"
        )
        fh.write("\t".join(cols) + "\n")
        for gid in sorted(replacements, key=int):
            fh.write("\t".join(replacements[gid][c] for c in cols) + "\n")


# --------------------------------------------------------------------------- #
# Assembly
# --------------------------------------------------------------------------- #


def build_rows(genes: dict, uniprot_by_gene: dict) -> list:
    rows = []
    for gene_id in sorted(genes, key=int):
        rec = dict(genes[gene_id])
        accessions = uniprot_by_gene.get(gene_id, [])
        rec["uniprot_accession"] = accessions[0] if accessions else ""
        rec["uniprot_accessions_all"] = ";".join(accessions)
        rows.append(rec)
    return rows


def coverage_stats(rows: list) -> dict:
    total = len(rows)

    def pct(n):
        return round(100.0 * n / total, 2) if total else 0.0

    n_ensembl = sum(1 for r in rows if r["ensembl_gene_id"])
    n_uniprot = sum(1 for r in rows if r["uniprot_accession"])
    n_named = sum(1 for r in rows if r["gene_symbol"] and not LOC_SYMBOL_RE.match(r["gene_symbol"]))
    n_multi_uniprot = sum(1 for r in rows if ";" in r["uniprot_accessions_all"])

    by_type = defaultdict(int)
    for r in rows:
        by_type[r["gene_type"] or "unspecified"] += 1

    return {
        "genes_total": total,
        "with_ensembl_xref": {"n": n_ensembl, "pct": pct(n_ensembl)},
        "with_uniprot_accession": {"n": n_uniprot, "pct": pct(n_uniprot)},
        "with_multiple_uniprot_accessions": {"n": n_multi_uniprot, "pct": pct(n_multi_uniprot)},
        "with_named_symbol": {"n": n_named, "pct": pct(n_named)},
        "with_orthogroup": {"n": 0, "pct": 0.0},
        "genes_by_type": dict(sorted(by_type.items(), key=lambda kv: -kv[1])),
    }


def write_crosswalk(rows: list, out_path: Path, taxid: int, organism: str, assembly: str) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as fh:
        fh.write(
            f"# AREE gene identifier crosswalk — {organism} (NCBI taxid {taxid})\n"
            f"# Built by src/mappings/build_crosswalk.py v{BUILDER_VERSION} on {date.today().isoformat()}\n"
            "# Sources: NCBI Gene (gene_info) and UniProtKB. See the .provenance.json\n"
            "# sidecar for source URLs, checksums, and coverage statistics.\n"
            f"# Current NCBI reference annotation at build time: {assembly}\n"
            "# NOTE: orthogroup_id is intentionally empty — real orthogroups require an\n"
            "# ortholog inference run and are not fabricated here.\n"
        )
        fh.write("\t".join(COLUMNS) + "\n")
        for r in rows:
            fh.write("\t".join(r.get(c, "") for c in COLUMNS) + "\n")


def write_provenance(
    out_path: Path,
    crosswalk_path: Path,
    taxid: int,
    organism: str,
    assembly: str,
    sources: list,
    stats: dict,
) -> None:
    doc = {
        "artifact": crosswalk_path.name,
        "artifact_sha256": sha256_file(crosswalk_path),
        "builder": "src/mappings/build_crosswalk.py",
        "builder_version": BUILDER_VERSION,
        "date_generated": date.today().isoformat(),
        "organism": organism,
        "ncbi_taxid": taxid,
        "ncbi_reference_annotation_at_build": assembly,
        "sources": sources,
        "coverage": stats,
        "known_limitations": [
            "UniProt cross-references to NCBI GeneID are sparse for this organism; most "
            "TrEMBL entries carry no GeneID or RefSeq xref and therefore cannot be linked.",
            "orthogroup_id is empty; cross-species evidence pooling requires an ortholog "
            "inference run that this builder does not perform.",
            "Ensembl xrefs come from NCBI dbXrefs and reflect Ensembl Metazoa's assembly, "
            "which may differ from the current NCBI reference assembly.",
            "gene_synonyms is emitted for curation and future symbol resolution but is not "
            "currently consulted by harmonize.identifiers.",
            "Retired GeneIDs resolve through the *_retired_gene_ids.tsv sidecar and are "
            "labelled 'inferred' because the remap crosses an annotation version; GeneIDs "
            "discontinued without a replacement stay 'unresolved' by design.",
        ],
    }
    out_path.write_text(json.dumps(doc, indent=2) + "\n")


# --------------------------------------------------------------------------- #
# Entry point
# --------------------------------------------------------------------------- #


def build(
    taxid: int = 29159,
    organism: str = "Magallana gigas (Crassostrea gigas), Pacific oyster",
    assembly: str = "GCF_963853765.1 (xbMagGiga1.1), annotation release GCF_963853765.1-RS_2024_06",
    gene_info_path: Path | None = None,
    uniprot_path: Path | None = None,
    gene_history_path: Path | None = None,
    out_dir: Path = CROSSWALK_DIR,
    download: bool = False,
    slug: str = "mgigas",
) -> Path:
    SOURCES_DIR.mkdir(parents=True, exist_ok=True)
    gene_info_path = Path(gene_info_path) if gene_info_path else SOURCES_DIR / f"{slug}_gene_info.tsv"
    uniprot_path = Path(uniprot_path) if uniprot_path else SOURCES_DIR / f"{slug}_uniprot.tsv"
    gene_history_path = (
        Path(gene_history_path) if gene_history_path else SOURCES_DIR / f"{slug}_gene_history.tsv"
    )

    if download or not gene_info_path.exists():
        print(f"Downloading NCBI gene_info for taxid {taxid} (streams ~230 MB, filters in place)...")
        download_gene_info(taxid, gene_info_path)
    if download or not uniprot_path.exists():
        print(f"Downloading UniProtKB entries for taxid {taxid}...")
        download_uniprot(taxid, uniprot_path)
    if download or not gene_history_path.exists():
        print(f"Downloading NCBI gene_history for taxid {taxid} (streams ~154 MB)...")
        download_gene_history(taxid, gene_history_path)

    genes = parse_gene_info(gene_info_path, taxid)
    if not genes:
        raise RuntimeError(f"no gene_info records found for taxid {taxid} in {gene_info_path}")
    uniprot_by_gene = parse_uniprot(uniprot_path)
    replacements, dead = parse_gene_history(gene_history_path, taxid)
    # A replacement is only useful if the replacing gene is actually in the crosswalk.
    replacements = {k: v for k, v in replacements.items() if v["current_gene_id"] in genes}

    rows = build_rows(genes, uniprot_by_gene)
    stats = coverage_stats(rows)

    out_dir.mkdir(parents=True, exist_ok=True)
    crosswalk_path = out_dir / f"{slug}_gene_id_crosswalk.tsv"
    write_crosswalk(rows, crosswalk_path, taxid, organism, assembly)
    retired_path = out_dir / f"{slug}_retired_gene_ids.tsv"
    write_retired_table(replacements, dead, retired_path, taxid)
    stats["retired_gene_ids_remappable"] = len(replacements)
    stats["retired_gene_ids_without_replacement"] = len(dead)

    sources = [
        {
            "name": "NCBI Gene gene_info",
            "url": NCBI_GENE_INFO_URL,
            "local_filtered_copy": str(gene_info_path.relative_to(REPO_ROOT))
            if gene_info_path.is_relative_to(REPO_ROOT)
            else str(gene_info_path),
            "sha256_filtered": sha256_file(gene_info_path),
            "records": len(genes),
        },
        {
            "name": "UniProtKB",
            "url": UNIPROT_STREAM_URL.format(taxid=taxid),
            "local_copy": str(uniprot_path.relative_to(REPO_ROOT))
            if uniprot_path.is_relative_to(REPO_ROOT)
            else str(uniprot_path),
            "sha256": sha256_file(uniprot_path),
            "entries_with_geneid_xref": sum(len(v) for v in uniprot_by_gene.values()),
            "genes_linked": len(uniprot_by_gene),
        },
        {
            "name": "NCBI Gene gene_history",
            "url": NCBI_GENE_HISTORY_URL,
            "local_filtered_copy": str(gene_history_path.relative_to(REPO_ROOT))
            if gene_history_path.is_relative_to(REPO_ROOT)
            else str(gene_history_path),
            "sha256_filtered": sha256_file(gene_history_path)
            if gene_history_path.exists()
            else None,
            "remappable_retired_ids": len(replacements),
            "discontinued_without_replacement": len(dead),
        },
    ]
    write_provenance(
        out_dir / f"{slug}_gene_id_crosswalk.provenance.json",
        crosswalk_path,
        taxid,
        organism,
        assembly,
        sources,
        stats,
    )
    return crosswalk_path


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description="Build a real gene identifier crosswalk for AREE.")
    p.add_argument("--taxid", type=int, default=29159)
    p.add_argument("--slug", default="mgigas", help="Output filename prefix.")
    p.add_argument("--organism", default="Magallana gigas (Crassostrea gigas), Pacific oyster")
    p.add_argument(
        "--assembly",
        default="GCF_963853765.1 (xbMagGiga1.1), annotation release GCF_963853765.1-RS_2024_06",
    )
    p.add_argument("--gene-info", default=None, help="Pre-downloaded taxid-filtered gene_info TSV.")
    p.add_argument("--uniprot", default=None, help="Pre-downloaded UniProt TSV export.")
    p.add_argument("--gene-history", default=None, help="Pre-downloaded taxid-filtered gene_history TSV.")
    p.add_argument("--out-dir", default=str(CROSSWALK_DIR))
    p.add_argument("--download", action="store_true", help="Force re-download of sources.")
    args = p.parse_args(argv)

    if not shutil.which("curl"):
        print("curl is required to download reference sources.", file=sys.stderr)
        return 1

    path = build(
        taxid=args.taxid,
        organism=args.organism,
        assembly=args.assembly,
        gene_info_path=args.gene_info,
        uniprot_path=args.uniprot,
        gene_history_path=args.gene_history,
        out_dir=Path(args.out_dir),
        download=args.download,
        slug=args.slug,
    )
    print(f"Wrote {path}")
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
