"""Build an identifier crosswalk for a genome with no NCBI Gene records. WP7.

The NCBI route (build_crosswalk.py) anchors every gene on its NCBI GeneID,
which only exists for genomes NCBI has annotated (RefSeq, GCF_). For a genome
annotated by its authors (a GenBank GCA_ annotation, or a GFF3 on figshare or
a lab server), this route builds the crosswalk from two author files and one
recorded homology search:

1. **The authors' GFF3.** Each gene becomes one row, anchored on ``gene_id``:
   the gene's GFF ID with any ``gene-`` prefix removed (for an NCBI GenBank
   GFF, that is the locus tag, e.g. ``MCOR_1``). Its mRNA IDs and protein IDs
   go into ``transcript_ids`` and ``protein_ids`` so a study reporting either
   still resolves to the gene. ``ncbi_gene_id`` stays empty unless the GFF
   carries a GeneID Dbxref.
2. **The authors' proteins.** The longest protein per gene is the query.
3. **DIAMOND blastp against a UniProtKB release** (Swiss-Prot by default),
   with explicit e-value, identity and coverage thresholds recorded in the
   provenance. The best hit is ``uniprot_accession``; every hit passing the
   thresholds is in ``uniprot_accessions_all``; ``uniprot_evidence`` is
   ``diamond_blastp``.

How confident a resolution is depends on what it went through, not on the
route that built the file. A study reporting the authors' gene, transcript or
protein ID names that gene in that annotation, so it resolves ``exact``. A
study reporting a UniProt accession reaches the gene only through our
homology search, so it resolves ``inferred`` (see harmonize.identifiers). The
authors' own homology names (GenBank ``standard_name``) become
``gene_symbol``, which resolves ``inferred`` as every symbol does.

The authors' own UniProt assignment, where the GFF records one (GenBank
``inference=...UniProtKB/Swiss-Prot:<acc>``), is kept in
``author_uniprot_accession``, and the provenance reports how often it agrees
with ours: a check on the recorded search that costs nothing.
"""
from __future__ import annotations

import gzip
import json
import re
import shutil
import subprocess
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from datetime import date
from pathlib import Path
from urllib.parse import unquote

from common import sha256_file

BUILDER_VERSION = "1.0.0"
ROUTE = "gff_homology"

COLUMNS = [
    "gene_id",
    "ncbi_gene_id",
    "ensembl_gene_id",
    "uniprot_accession",
    "uniprot_accessions_all",
    "uniprot_evidence",
    "locus_id",
    "gene_symbol",
    "gene_synonyms",
    "gene_description",
    "gene_type",
    "orthogroup_id",
    "legacy_locus_tags",
    "transcript_ids",
    "protein_ids",
    "author_uniprot_accession",
    "homology_identity",
    "homology_evalue",
    "homology_bitscore",
]

GENE_TYPES = {"gene", "pseudogene", "ncRNA_gene"}
TRANSCRIPT_TYPES = {"mRNA", "transcript", "lnc_RNA", "ncRNA", "tRNA", "rRNA", "snRNA", "snoRNA", "pseudogenic_transcript"}
_PREFIX_RE = re.compile(r"^(gene|rna|cds|id)-", re.I)
_UNIPROT_INFERENCE_RE = re.compile(r"UniProtKB/(?:Swiss-Prot|TrEMBL):([A-Z0-9]+)(?:\.\d+)?")
_UNIPROT_HEADER_RE = re.compile(r"^(?:sp|tr)\|([A-Z0-9]+)\|")


@dataclass(frozen=True)
class Thresholds:
    """Homology thresholds; every hit kept must pass all of them."""
    max_evalue: float = 1e-10
    min_identity: float = 30.0
    min_query_cover: float = 50.0
    min_subject_cover: float = 50.0
    max_target_seqs: int = 5

    def as_dict(self) -> dict:
        return {"max_evalue": self.max_evalue, "min_identity_pct": self.min_identity,
                "min_query_cover_pct": self.min_query_cover, "min_subject_cover_pct": self.min_subject_cover,
                "max_target_seqs": self.max_target_seqs}


@dataclass
class Annotation:
    genes: dict = field(default_factory=dict)          # gene_id -> record
    transcript_gene: dict = field(default_factory=dict)  # transcript id -> gene_id
    protein_gene: dict = field(default_factory=dict)     # protein id -> gene_id


def strip_prefix(value: str) -> str:
    return _PREFIX_RE.sub("", value)


def _attributes(text: str) -> dict:
    out = {}
    for part in text.strip().strip(";").split(";"):
        if "=" in part:
            key, value = part.split("=", 1)
            out[key.strip()] = unquote(value.strip())
    return out


def _open(path: Path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def parse_gff(path: Path) -> Annotation:
    """Genes, transcript -> gene and protein -> gene links from a GFF3.

    Works for NCBI GenBank GFFs (gene/mRNA/CDS with locus_tag and protein_id)
    and for author GFFs such as EVidenceModeler or BRAKER output (gene/mRNA/CDS
    linked by Parent, proteins named after the mRNA). A CDS's protein is its
    ``protein_id`` attribute, else its ``Name``, else the mRNA it belongs to.
    """
    ann = Annotation()
    raw_gene: dict[str, str] = {}                 # raw GFF ID -> gene_id
    raw_tx: dict[str, str] = {}                   # raw GFF ID -> gene_id
    names: dict[str, Counter] = defaultdict(Counter)
    products: dict[str, Counter] = defaultdict(Counter)
    author_uniprot: dict[str, Counter] = defaultdict(Counter)
    with _open(path) as fh:
        for line in fh:
            if line.startswith("##FASTA"):
                break
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            ftype, attrs = cols[2], _attributes(cols[8])
            raw_id = attrs.get("ID", "")
            if ftype in GENE_TYPES and raw_id:
                gene_id = strip_prefix(raw_id)
                if gene_id in ann.genes:
                    raise ValueError(f"{path.name}: gene ID {gene_id!r} occurs twice")
                dbx = dict(x.split(":", 1) for x in attrs.get("Dbxref", "").split(",") if ":" in x)
                name = attrs.get("Name", "")
                ann.genes[gene_id] = {
                    "gene_id": gene_id,
                    "ncbi_gene_id": dbx.get("GeneID", ""),
                    "locus_id": attrs.get("locus_tag", "") or gene_id,
                    "gene_type": attrs.get("gene_biotype", "") or ("pseudogene" if ftype == "pseudogene" else
                                                                   "protein_coding" if ftype == "gene" else ftype),
                    "_author_name": name if name and name not in (gene_id, raw_id, attrs.get("locus_tag")) else "",
                    "_author_description": attrs.get("description", "") or attrs.get("Note", ""),
                }
                raw_gene[raw_id] = gene_id
            elif ftype in TRANSCRIPT_TYPES and raw_id:
                parent = attrs.get("Parent", "").split(",")[0]
                gene_id = raw_gene.get(parent)
                if gene_id is None:
                    continue
                raw_tx[raw_id] = gene_id
                ann.transcript_gene[strip_prefix(raw_id)] = gene_id
                if attrs.get("transcript_id"):
                    ann.transcript_gene[attrs["transcript_id"]] = gene_id
            elif ftype == "CDS":
                parent = attrs.get("Parent", "").split(",")[0]
                gene_id = raw_tx.get(parent) or raw_gene.get(parent)
                if gene_id is None:
                    continue
                protein = attrs.get("protein_id") or attrs.get("Name") or strip_prefix(parent)
                ann.protein_gene[protein] = gene_id
                if "." in protein and attrs.get("protein_id"):
                    ann.protein_gene.setdefault(protein.rsplit(".", 1)[0], gene_id)
                for key in ("standard_name", "gene"):
                    if attrs.get(key):
                        names[gene_id][attrs[key]] += 1
                product = attrs.get("product", "")
                if product and product != protein:
                    products[gene_id][product] += 1
                for acc in _UNIPROT_INFERENCE_RE.findall(attrs.get("inference", "")):
                    author_uniprot[gene_id][acc] += 1
    for gene_id, rec in ann.genes.items():
        ranked = [n for n, _ in names[gene_id].most_common()]
        if rec["_author_name"]:
            ranked = [rec["_author_name"]] + [n for n in ranked if n != rec["_author_name"]]
        rec["gene_symbol"] = ranked[0] if ranked else ""
        rec["gene_synonyms"] = ";".join(ranked[1:])
        rec["_author_product"] = products[gene_id].most_common(1)[0][0] if products[gene_id] else ""
        rec["author_uniprot_accession"] = (author_uniprot[gene_id].most_common(1)[0][0]
                                           if author_uniprot[gene_id] else "")
    return ann


def read_fasta(path: Path) -> dict[str, str]:
    seqs: dict[str, list[str]] = {}
    current = None
    with _open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                current = line[1:].split()[0]
                seqs[current] = []
            elif current is not None:
                seqs[current].append(line.strip())
    return {k: "".join(v).rstrip("*") for k, v in seqs.items()}


def longest_protein_per_gene(proteins: dict[str, str], ann: Annotation) -> tuple[dict, list[str]]:
    """gene_id -> (protein id, sequence) for the longest protein, and headers matched to no gene."""
    best: dict[str, tuple[str, str]] = {}
    unmatched = []
    for header, seq in proteins.items():
        gene_id = (ann.protein_gene.get(header) or ann.transcript_gene.get(strip_prefix(header))
                   or ann.protein_gene.get(header.rsplit(".", 1)[0])
                   or (strip_prefix(header) if strip_prefix(header) in ann.genes else None))
        if gene_id is None:
            unmatched.append(header)
            continue
        ann.protein_gene.setdefault(header, gene_id)
        if gene_id not in best or len(seq) > len(best[gene_id][1]):
            best[gene_id] = (header, seq)
    return best, unmatched


def write_query(best: dict, out: Path) -> Path:
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w") as fh:
        for gene_id in sorted(best):
            fh.write(f">{gene_id}\n{best[gene_id][1]}\n")
    return out


def diamond_version(diamond: str) -> str:
    return subprocess.run([diamond, "version"], capture_output=True, text=True, check=True).stdout.strip()


DIAMOND_FIELDS = ["qseqid", "sseqid", "pident", "evalue", "bitscore", "qcovhsp", "scovhsp", "stitle"]


def run_diamond(diamond: str, query: Path, database_fasta: Path, workdir: Path, thresholds: Thresholds,
                threads: int = 4, sensitivity: str = "--sensitive") -> tuple[Path, list[str]]:
    """makedb then blastp; returns (hits TSV, the blastp command as run)."""
    workdir.mkdir(parents=True, exist_ok=True)
    db = workdir / "uniprot.dmnd"
    if not db.exists():
        subprocess.run([diamond, "makedb", "--in", str(database_fasta), "--db", str(db), "--threads", str(threads)],
                       check=True, capture_output=True)
    out = workdir / "hits.tsv"
    command = [diamond, "blastp", sensitivity, "--query", str(query), "--db", str(db), "--out", str(out),
               "--outfmt", "6", *DIAMOND_FIELDS, "--evalue", str(thresholds.max_evalue),
               "--id", str(thresholds.min_identity), "--query-cover", str(thresholds.min_query_cover),
               "--subject-cover", str(thresholds.min_subject_cover),
               "--max-target-seqs", str(thresholds.max_target_seqs), "--max-hsps", "1",
               "--threads", str(threads)]
    subprocess.run(command, check=True, capture_output=True)
    return out, command


def parse_hits(path: Path, thresholds: Thresholds) -> dict[str, list[dict]]:
    """gene_id -> hits passing every threshold, best bitscore first.

    The thresholds are re-applied here, so a hits file produced elsewhere (or
    with looser settings) cannot slip weaker hits into the crosswalk.
    """
    hits: dict[str, list[dict]] = defaultdict(list)
    with open(path) as fh:
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            if len(cols) < len(DIAMOND_FIELDS):
                continue
            rec = dict(zip(DIAMOND_FIELDS, cols))
            pident, evalue = float(rec["pident"]), float(rec["evalue"])
            qcov, scov = float(rec["qcovhsp"]), float(rec["scovhsp"])
            if (evalue > thresholds.max_evalue or pident < thresholds.min_identity
                    or qcov < thresholds.min_query_cover or scov < thresholds.min_subject_cover):
                continue
            m = _UNIPROT_HEADER_RE.match(rec["sseqid"])
            accession = m.group(1) if m else rec["sseqid"]
            title = rec["stitle"].split(" ", 1)[1] if " " in rec["stitle"] else ""
            hits[rec["qseqid"]].append({
                "accession": accession, "identity": pident, "evalue": evalue,
                "bitscore": float(rec["bitscore"]),
                "protein_name": title.split(" OS=")[0].strip(),
                "organism": (title.split(" OS=")[1].split(" OX=")[0].strip() if " OS=" in title else ""),
            })
    for gene_hits in hits.values():
        gene_hits.sort(key=lambda h: (-h["bitscore"], h["evalue"], h["accession"]))
        seen = set()
        gene_hits[:] = [h for h in gene_hits if not (h["accession"] in seen or seen.add(h["accession"]))]
    return dict(hits)


def build_rows(ann: Annotation, hits: dict[str, list[dict]]) -> list[dict]:
    transcripts: dict[str, set] = defaultdict(set)
    for tx, gene_id in ann.transcript_gene.items():
        if tx != gene_id:
            transcripts[gene_id].add(tx)
    proteins: dict[str, set] = defaultdict(set)
    for prot, gene_id in ann.protein_gene.items():
        proteins[gene_id].add(prot)
    rows = []
    for gene_id in sorted(ann.genes, key=_natural_key):
        rec = ann.genes[gene_id]
        gene_hits = hits.get(gene_id, [])
        best = gene_hits[0] if gene_hits else None
        description = rec["_author_description"] or rec["_author_product"]
        if not description and best:
            description = f"similar to {best['protein_name']} ({best['organism']})".replace(" ()", "")
        rows.append({
            "gene_id": gene_id,
            "ncbi_gene_id": rec["ncbi_gene_id"],
            "ensembl_gene_id": "",
            "uniprot_accession": best["accession"] if best else "",
            "uniprot_accessions_all": ";".join(h["accession"] for h in gene_hits),
            "uniprot_evidence": "diamond_blastp" if best else "",
            "locus_id": rec["locus_id"],
            "gene_symbol": rec["gene_symbol"],
            "gene_synonyms": rec["gene_synonyms"],
            "gene_description": description,
            "gene_type": rec["gene_type"],
            "orthogroup_id": "",
            "legacy_locus_tags": "",
            "transcript_ids": ";".join(sorted(transcripts[gene_id], key=_natural_key)),
            "protein_ids": ";".join(sorted(proteins[gene_id], key=_natural_key)),
            "author_uniprot_accession": rec["author_uniprot_accession"],
            "homology_identity": f"{best['identity']:.1f}" if best else "",
            "homology_evalue": f"{best['evalue']:.2g}" if best else "",
            "homology_bitscore": f"{best['bitscore']:.1f}" if best else "",
        })
    return rows


def _natural_key(value: str):
    return [int(t) if t.isdigit() else t for t in re.split(r"(\d+)", value)]


def coverage_stats(rows: list[dict], n_proteins_unmatched: int) -> dict:
    total = len(rows)

    def pct(n):
        return round(100.0 * n / total, 2) if total else 0.0

    coding = [r for r in rows if r["protein_ids"]]
    with_hit = [r for r in rows if r["uniprot_accession"]]
    both = [r for r in rows if r["uniprot_accession"] and r["author_uniprot_accession"]]
    agree = [r for r in both if r["author_uniprot_accession"] in r["uniprot_accessions_all"].split(";")]
    return {
        "genes_total": total,
        "genes_with_protein": {"n": len(coding), "pct": pct(len(coding))},
        "with_uniprot_accession": {"n": len(with_hit), "pct": pct(len(with_hit)),
                                   "pct_of_genes_with_protein": round(100.0 * len(with_hit) / len(coding), 2)
                                   if coding else 0.0},
        "with_multiple_uniprot_accessions": sum(1 for r in rows if ";" in r["uniprot_accessions_all"]),
        "with_author_symbol": {"n": sum(1 for r in rows if r["gene_symbol"]),
                               "pct": pct(sum(1 for r in rows if r["gene_symbol"]))},
        "with_ncbi_gene_id": sum(1 for r in rows if r["ncbi_gene_id"]),
        "with_orthogroup": {"n": 0, "pct": 0.0},
        "author_uniprot_agreement": {
            "genes_with_both": len(both),
            "author_accession_among_our_hits": len(agree),
            "pct": round(100.0 * len(agree) / len(both), 2) if both else None,
        },
        "proteins_matched_to_no_gene": n_proteins_unmatched,
        "genes_by_type": dict(Counter(r["gene_type"] or "unspecified" for r in rows).most_common()),
    }


def write_crosswalk(rows: list[dict], out_path: Path, taxid: int, organism: str, annotation: str) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as fh:
        fh.write(
            f"# AREE gene identifier crosswalk — {organism} (NCBI taxid {taxid})\n"
            f"# Built by src/mappings/gff_crosswalk.py v{BUILDER_VERSION} on {date.today().isoformat()}\n"
            "# Route: gff_homology. Genes from the authors' annotation (anchor: gene_id);\n"
            "# UniProt accessions from a recorded DIAMOND blastp search, so a UniProt lookup\n"
            "# resolves 'inferred'. See the .provenance.json sidecar.\n"
            f"# Annotation: {annotation}\n"
            "# NOTE: orthogroup_id is intentionally empty; orthogroups come from an orthogroup run.\n"
        )
        fh.write("\t".join(COLUMNS) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(c, "")).replace("\t", " ") for c in COLUMNS) + "\n")


def build(*, gff: Path, proteins: Path, uniprot_fasta: Path, uniprot_release: str, uniprot_url: str | None,
          taxid: int, organism: str, annotation: str, out_dir: Path, slug: str, workdir: Path,
          diamond: str = "diamond", threads: int = 4, thresholds: Thresholds | None = None,
          hits_path: Path | None = None, today: str | None = None) -> Path:
    """Build `<slug>_gene_id_crosswalk.tsv` and its provenance sidecar.

    `hits_path` reuses a DIAMOND hits table instead of running DIAMOND (tests,
    or a rebuild); the thresholds are re-applied either way, and the provenance
    records which happened.
    """
    thresholds = thresholds or Thresholds()
    ann = parse_gff(gff)
    if not ann.genes:
        raise ValueError(f"{gff} contains no gene features")
    best, unmatched = longest_protein_per_gene(read_fasta(proteins), ann)
    if not best:
        raise ValueError(f"no protein in {proteins} matched a gene in {gff}; check the FASTA headers")
    query = write_query(best, workdir / "longest_proteins.faa")

    if hits_path is None:
        if shutil.which(diamond) is None and not Path(diamond).exists():
            raise FileNotFoundError(f"DIAMOND not found: {diamond!r}. Pass --diamond <path>.")
        hits_path, command = run_diamond(diamond, query, uniprot_fasta, workdir, thresholds, threads)
        search = {"tool": diamond_version(diamond), "command": command, "run_by_builder": True}
    else:
        search = {"tool": "precomputed hits table", "command": None, "run_by_builder": False}
    hits = parse_hits(hits_path, thresholds)

    rows = build_rows(ann, hits)
    out_path = out_dir / f"{slug}_gene_id_crosswalk.tsv"
    write_crosswalk(rows, out_path, taxid, organism, annotation)
    provenance = {
        "artifact": out_path.name,
        "artifact_sha256": sha256_file(out_path),
        "builder": "src/mappings/gff_crosswalk.py",
        "builder_version": BUILDER_VERSION,
        "route": ROUTE,
        "anchor_column": "gene_id",
        "date_generated": today or date.today().isoformat(),
        "organism": organism,
        "ncbi_taxid": taxid,
        "ncbi_reference_annotation_at_build": annotation,
        "sources": [
            {"name": "annotation_gff3", "file": gff.name, "sha256": sha256_file(gff)},
            {"name": "annotation_proteins", "file": proteins.name, "sha256": sha256_file(proteins),
             "query": "longest protein per gene", "n_queries": len(best)},
            {"name": "uniprot_database", "file": uniprot_fasta.name, "sha256": sha256_file(uniprot_fasta),
             "release": uniprot_release, "url": uniprot_url},
            {"name": "homology_hits", "file": hits_path.name, "sha256": sha256_file(hits_path)},
        ],
        "homology_search": {**search, "thresholds": thresholds.as_dict(),
                            "selection": "best bitscore per gene; all hits passing thresholds kept"},
        "coverage": coverage_stats(rows, len(unmatched)),
        "resolution": {
            "gene_id, locus_id, transcript_ids, protein_ids": "exact (identity in the authors' annotation)",
            "uniprot_accession": "inferred (homology), many_to_one_ortholog when several genes share it",
            "gene_symbol": "inferred, as for every route",
        },
        "known_limitations": [
            "UniProt accessions are homology assignments from one recorded search, not curated "
            "cross-references; a gene's best Swiss-Prot hit is usually another species' protein.",
            "gene_symbol is the authors' name for the gene where the GFF gives one (for GenBank "
            "annotations, the CDS standard_name, itself homology-derived); it is not unique.",
            "orthogroup_id is empty; cross-species comparison needs this species in an orthogroup run.",
            "Genes without a protein (non-coding, or proteins whose FASTA headers match no GFF ID) "
            "have no homology fields.",
        ],
    }
    out_path.with_name(f"{out_path.stem}.provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    return out_path
