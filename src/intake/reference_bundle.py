"""Fetch and verify a reference bundle for a registered NCBI RefSeq assembly.

A *reference bundle* is what the raw-reanalysis workflows read from
``data/reference/<ACCESSION>/``:

* ``rna.fna``      RefSeq transcript FASTA (Salmon index input)
* ``genomic.gtf``  RefSeq annotation (Salmon, methylation region annotation)
* ``tx2gene.tsv``  transcript -> NCBI GeneID, headerless (DESeq2 / tximport)
* ``genomic.fna``  genome FASTA, only with ``--genome`` (Bismark)
* ``primary_proteins.faa``  one protein per NCBI GeneID (the longest), only
  with ``--proteins``; the input to orthogroup inference
  (``aree build-orthogroups``). Built from NCBI's ``protein.faa`` and
  ``feature_table.txt``, which are kept beside it.

plus the compressed NCBI originals, NCBI's ``md5checksums.txt``, and
``reference_provenance.json``, which records where every file came from and its
checksums, and two consistency checks that matter downstream:

* the annotation release named in the GTF header against the release
  registered for the assembly in ``data/reference/genome_assemblies.yaml``. NCBI
  can re-annotate an assembly in place, and evidence must not claim the
  registered release when the transcriptome came from a different one;
* the share of tx2gene GeneIDs found in the species' identifier crosswalk. The
  workflow's gene IDs are what ``aree harmonize`` resolves, so a bundle whose
  genes the crosswalk does not know would produce ``unresolved`` evidence.

This generalizes the curl/awk recipe that built the first oyster bundle (see
docs/first_raw_reanalysis.md); ``build_tx2gene`` reproduces that awk exactly.
Salmon and Bismark indexes are not built here; the workflows build them.
"""
from __future__ import annotations

import gzip
import hashlib
import json
import re
import shutil
import urllib.request
from datetime import date
from pathlib import Path

from common import (
    DATA_DIR,
    assembly_crosswalk_path,
    resolve_assembly,
    sha256_file,
)

NCBI_GENOMES_BASE = "https://ftp.ncbi.nlm.nih.gov/genomes/all"
TIMEOUT = 120
PROVENANCE_NAME = "reference_provenance.json"
BUNDLER_VERSION = "1.0.0"

# NCBI file suffix -> name of the decompressed file the workflows read.
CORE_FILES = {"_rna.fna.gz": "rna.fna", "_genomic.gtf.gz": "genomic.gtf"}
GENOME_FILE = {"_genomic.fna.gz": "genomic.fna"}
PROTEIN_FILES = {"_protein.faa.gz": "protein.faa", "_feature_table.txt.gz": "feature_table.txt"}
PRIMARY_PROTEOME_RULE = (
    "longest protein (feature_table product_length) per GeneID among CDS features; "
    "ties broken by the lowest product_accession"
)


class ReferenceBundleError(RuntimeError):
    """A bundle could not be fetched, or does not match what was recorded."""


# --------------------------------------------------------------------------- #
# NCBI genomes FTP layout
# --------------------------------------------------------------------------- #


def accession_dir_url(accession: str) -> str:
    """Parent directory holding one assembly's folder, e.g. .../GCF/965/363/235/."""
    m = re.fullmatch(r"(GC[AF])_(\d{3})(\d{3})(\d{3})\.\d+", accession)
    if not m:
        raise ReferenceBundleError(f"not an NCBI assembly accession: {accession!r}")
    prefix, a, b, c = m.groups()
    return f"{NCBI_GENOMES_BASE}/{prefix}/{a}/{b}/{c}/"


def find_assembly_folder(listing_html: str, accession: str) -> str:
    """The assembly's folder name (``<accession>_<assembly name>``) in a directory listing.

    The folder name embeds NCBI's assembly name with its own character
    substitutions, so it is read from the listing rather than reconstructed.
    """
    names = set(re.findall(rf'href="({re.escape(accession)}_[^"/]+)/?"', listing_html))
    if len(names) != 1:
        raise ReferenceBundleError(
            f"expected one folder for {accession} in the NCBI listing, found {sorted(names)}"
        )
    return names.pop()


def parse_md5_listing(text: str) -> dict[str, str]:
    """``md5checksums.txt`` -> {file name: md5}."""
    out = {}
    for line in text.splitlines():
        parts = line.split()
        if len(parts) == 2:
            md5, name = parts
            out[name.removeprefix("./")] = md5
    return out


def _http_get_text(url: str) -> str:
    with urllib.request.urlopen(url, timeout=TIMEOUT) as fh:  # noqa: S310 - fixed NCBI host
        return fh.read().decode()


def _http_download(url: str, dest: Path) -> None:
    tmp = dest.with_name(dest.name + ".part")
    with urllib.request.urlopen(url, timeout=TIMEOUT) as fh, open(tmp, "wb") as out:  # noqa: S310
        shutil.copyfileobj(fh, out, length=1 << 20)
    tmp.replace(dest)


def md5_file(path: Path) -> str:
    h = hashlib.md5()  # noqa: S324 - NCBI publishes MD5; used for integrity, not security
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


# --------------------------------------------------------------------------- #
# Derived files
# --------------------------------------------------------------------------- #

_TX_RE = re.compile(r'transcript_id "([^"]+)"')
_GENEID_RE = re.compile(r"GeneID:(\d+)")


def build_tx2gene(gtf_path: Path) -> list[tuple[str, str]]:
    """(transcript_id, NCBI GeneID) for every GTF ``transcript`` feature, sorted, unique.

    Same output as the awk in docs/first_raw_reanalysis.md, except that a
    transcript line with no GeneID is skipped rather than inheriting the
    previous line's gene (the awk never reset its variables).
    """
    pairs = set()
    with open(gtf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "transcript":
                continue
            tx, gene = _TX_RE.search(fields[8]), _GENEID_RE.search(fields[8])
            if tx and gene and tx.group(1):
                pairs.add((tx.group(1), gene.group(1)))
    return sorted(pairs)


def write_tx2gene(pairs: list[tuple[str, str]], path: Path) -> None:
    with open(path, "w") as fh:
        for tx, gene in pairs:
            fh.write(f"{tx}\t{gene}\n")


def gtf_header(gtf_path: Path) -> dict[str, str]:
    """The ``#!key value`` lines NCBI writes at the top of a GTF."""
    header = {}
    with open(gtf_path) as fh:
        for line in fh:
            if not line.startswith("#"):
                break
            if line.startswith("#!"):
                key, _, value = line[2:].strip().partition(" ")
                header[key] = value.strip()
    return header


def annotation_release_matches(registered: str | None, annotation_source: str | None) -> bool:
    """Whether the GTF's ``annotation-source`` names the registered release.

    NCBI writes e.g. ``NCBI RefSeq GCF_963853765.1-RS_2024_06`` for recent
    annotations and ``NCBI RefSeq <Species> Annotation Release 100`` for older
    ones; the registry records the release itself, so containment is the
    comparison, ignoring the word RefSeq.
    """
    if not registered or not annotation_source:
        return False

    # Older NCBI annotations are named "NCBI <Species> Annotation Release 100" by
    # the Datasets API but "NCBI RefSeq <Species> Annotation Release 100" in the
    # GTF header; the word RefSeq does not change which release it is.
    def norm(text: str) -> str:
        return " ".join(w for w in text.split() if w.lower() != "refseq")

    return norm(registered) in norm(annotation_source)


def crosswalk_coverage(gene_ids: set[str], crosswalk_path: Path | None) -> dict:
    """How many of the bundle's GeneIDs the species' crosswalk knows."""
    if crosswalk_path is None or not crosswalk_path.exists():
        return {"crosswalk": str(crosswalk_path) if crosswalk_path else None, "checked": False}
    known = set()
    with open(crosswalk_path) as fh:
        header = None
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if header is None:
                header = fields.index("ncbi_gene_id")
                continue
            known.add(fields[header])
    found = len(gene_ids & known)
    return {
        "crosswalk": crosswalk_path.name,
        "checked": True,
        "bundle_genes": len(gene_ids),
        "found_in_crosswalk": found,
        "pct": round(100 * found / len(gene_ids), 2) if gene_ids else 0.0,
        "missing_examples": sorted(gene_ids - known)[:10],
    }


def build_primary_proteome(protein_faa: Path, feature_table: Path) -> list[tuple[str, str, str]]:
    """(GeneID, protein accession, sequence) for one protein per gene, sorted by GeneID.

    Orthogroup inference wants one sequence per gene: isoforms of one gene would
    otherwise look like paralogs. The longest protein is the usual choice (it is
    what OrthoFinder's primary_transcript.py keeps). Lengths and the
    protein-to-gene link come from NCBI's feature table, not from FASTA headers.
    """
    best: dict[str, tuple[int, str]] = {}
    header = None
    with open(feature_table) as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if line.startswith("#"):
                header = [c.lstrip("# ").strip() for c in fields]
                continue
            if header is None:
                raise ReferenceBundleError(f"{feature_table}: no header line")
            row = dict(zip(header, fields))
            gene, acc, length = row.get("GeneID", ""), row.get("product_accession", ""), row.get("product_length", "")
            if row.get("feature") != "CDS" or not gene or not acc or not length.isdigit():
                continue
            current = best.get(gene)
            candidate = (int(length), acc)
            if current is None or candidate[0] > current[0] or (
                    candidate[0] == current[0] and candidate[1] < current[1]):
                best[gene] = candidate
    wanted = {acc: gene for gene, (_, acc) in best.items()}

    sequences: dict[str, list[str]] = {}
    accession = None
    with open(protein_faa) as fh:
        for line in fh:
            if line.startswith(">"):
                accession = line[1:].split()[0]
                if accession in wanted:
                    sequences[accession] = []
                else:
                    accession = None
            elif accession is not None:
                sequences[accession].append(line.strip())
    missing = sorted(set(wanted) - set(sequences))
    if missing:
        raise ReferenceBundleError(
            f"{len(missing)} chosen protein(s) are absent from {protein_faa}: {missing[:5]}")
    out = [(gene, acc, "".join(sequences[acc])) for acc, gene in wanted.items()]
    return sorted(out, key=lambda t: int(t[0]))


def write_primary_proteome(proteins: list[tuple[str, str, str]], path: Path) -> None:
    """FASTA headed by GeneID (OrthoFinder takes the first word as the gene ID)."""
    with open(path, "w") as fh:
        for gene, acc, seq in proteins:
            fh.write(f">{gene} {acc}\n")
            for i in range(0, len(seq), 80):
                fh.write(seq[i:i + 80] + "\n")


def _decompress(src: Path, dest: Path) -> None:
    tmp = dest.with_name(dest.name + ".part")
    with gzip.open(src, "rb") as fin, open(tmp, "wb") as fout:
        shutil.copyfileobj(fin, fout, length=1 << 20)
    tmp.replace(dest)


# --------------------------------------------------------------------------- #
# Fetch, record, verify
# --------------------------------------------------------------------------- #


def default_bundle_dir(accession: str) -> Path:
    return DATA_DIR / "reference" / accession


def _resolve_refseq_assembly(assembly_value: str) -> dict:
    record = resolve_assembly(assembly_value)
    if record is None:
        raise ReferenceBundleError(
            f"{assembly_value!r} is not in data/reference/genome_assemblies.yaml; "
            "register the assembly first (docs/adding_a_species.md)."
        )
    accession = record.get("ncbi_assembly_accession", "")
    if not accession.startswith("GCF_") or not record.get("annotation_release"):
        raise ReferenceBundleError(
            f"{record['assembly_id']} has no NCBI RefSeq annotation, so NCBI publishes no "
            "transcriptome or GTF for it. Raw reanalysis needs a RefSeq-annotated assembly."
        )
    return record


def plan_bundle(assembly_value: str, with_genome: bool = False,
                get_text=_http_get_text, with_proteins: bool = False) -> dict:
    """Resolve what would be fetched, without downloading any sequence data.

    Reads only the NCBI directory listing and ``md5checksums.txt``.
    """
    record = _resolve_refseq_assembly(assembly_value)
    accession = record["ncbi_assembly_accession"]
    parent = accession_dir_url(accession)
    folder = find_assembly_folder(get_text(parent), accession)
    base = f"{parent}{folder}/"
    md5s = parse_md5_listing(get_text(f"{base}md5checksums.txt"))
    wanted = {**CORE_FILES, **(GENOME_FILE if with_genome else {}),
              **(PROTEIN_FILES if with_proteins else {})}
    files = []
    for suffix, local_name in wanted.items():
        remote = f"{folder}{suffix}"
        if remote not in md5s:
            raise ReferenceBundleError(f"{remote} is not listed in {base}md5checksums.txt")
        files.append({"remote": remote, "url": f"{base}{remote}", "ncbi_md5": md5s[remote],
                      "local": local_name})
    return {"assembly": record, "base_url": base, "folder": folder, "files": files,
            "md5_listing": md5s}


def fetch_bundle(assembly_value: str, dest: Path | None = None, with_genome: bool = False,
                 get_text=_http_get_text, download=_http_download,
                 today: str | None = None, with_proteins: bool = False) -> dict:
    """Download, verify and record a reference bundle. Returns the provenance dict.

    Compressed files already present with the right MD5 are reused, so a
    re-run (or adding ``--genome`` later) downloads only what is missing.
    """
    plan = plan_bundle(assembly_value, with_genome=with_genome, get_text=get_text,
                       with_proteins=with_proteins)
    record = plan["assembly"]
    dest = Path(dest) if dest else default_bundle_dir(record["ncbi_assembly_accession"])
    dest.mkdir(parents=True, exist_ok=True)

    md5_path = dest / "md5checksums.txt"
    md5_path.write_text(get_text(f"{plan['base_url']}md5checksums.txt"))

    files = []
    for f in plan["files"]:
        gz = dest / f["remote"]
        if not (gz.exists() and md5_file(gz) == f["ncbi_md5"]):
            download(f["url"], gz)
        got = md5_file(gz)
        if got != f["ncbi_md5"]:
            gz.unlink(missing_ok=True)
            raise ReferenceBundleError(
                f"{f['remote']}: MD5 {got} does not match NCBI's {f['ncbi_md5']}; deleted, re-run"
            )
        local = dest / f["local"]
        _decompress(gz, local)
        files.append({
            "file": f["local"], "sha256": sha256_file(local),
            "source_file": f["remote"], "source_url": f["url"],
            "source_md5_ncbi": f["ncbi_md5"], "source_sha256": sha256_file(gz),
        })

    gtf = dest / "genomic.gtf"
    pairs = build_tx2gene(gtf)
    if not pairs:
        raise ReferenceBundleError(f"no transcript -> GeneID pairs found in {gtf}")
    tx2gene = dest / "tx2gene.tsv"
    write_tx2gene(pairs, tx2gene)
    files.append({"file": "tx2gene.tsv", "sha256": sha256_file(tx2gene),
                  "derived_from": "genomic.gtf", "derivation": "intake.reference_bundle.build_tx2gene"})

    primary = None
    if with_proteins:
        proteins = build_primary_proteome(dest / "protein.faa", dest / "feature_table.txt")
        primary_path = dest / "primary_proteins.faa"
        write_primary_proteome(proteins, primary_path)
        files.append({"file": "primary_proteins.faa", "sha256": sha256_file(primary_path),
                      "derived_from": "protein.faa, feature_table.txt",
                      "derivation": "intake.reference_bundle.build_primary_proteome"})
        primary = {"genes": len(proteins), "rule": PRIMARY_PROTEOME_RULE,
                   "crosswalk_coverage": crosswalk_coverage({g for g, _, _ in proteins},
                                                            assembly_crosswalk_path(record))}

    header = gtf_header(gtf)
    source = header.get("annotation-source")
    gene_ids = {g for _, g in pairs}
    provenance = {
        "bundle": str(dest.name),
        "bundler": "src/intake/reference_bundle.py",
        "bundler_version": BUNDLER_VERSION,
        "date_generated": today or date.today().isoformat(),
        "assembly_id": record["assembly_id"],
        "species": record.get("species"),
        "ncbi_taxid": record.get("ncbi_taxid"),
        "ncbi_assembly_accession": record["ncbi_assembly_accession"],
        "registered_annotation_release": record.get("annotation_release"),
        "gtf_header": header,
        "annotation_release_matches_registry": annotation_release_matches(
            record.get("annotation_release"), source),
        "source_base_url": plan["base_url"],
        "md5checksums_sha256": sha256_file(md5_path),
        "files": files,
        "tx2gene": {"transcripts": len(pairs), "genes": len(gene_ids)},
        "crosswalk_coverage": crosswalk_coverage(gene_ids, assembly_crosswalk_path(record)),
        **({"primary_proteome": primary} if primary else {}),
        "not_included": "Salmon and Bismark indexes; the workflows build them from these files.",
    }
    (dest / PROVENANCE_NAME).write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n")
    return provenance


def check_bundle(dest: Path) -> list[str]:
    """Problems with an existing bundle: missing files or checksum drift. Empty means OK."""
    dest = Path(dest)
    prov_path = dest / PROVENANCE_NAME
    if not prov_path.exists():
        return [f"{prov_path} is missing; fetch the bundle with `aree fetch-reference`"]
    provenance = json.loads(prov_path.read_text())
    problems = []
    for entry in provenance["files"]:
        path = dest / entry["file"]
        if not path.exists():
            problems.append(f"{entry['file']}: missing")
        elif sha256_file(path) != entry["sha256"]:
            problems.append(f"{entry['file']}: SHA-256 differs from {PROVENANCE_NAME}")
    if not provenance.get("annotation_release_matches_registry"):
        problems.append(
            "GTF annotation-source "
            f"{provenance.get('gtf_header', {}).get('annotation-source')!r} does not name the "
            f"registered release {provenance.get('registered_annotation_release')!r}"
        )
    return problems
