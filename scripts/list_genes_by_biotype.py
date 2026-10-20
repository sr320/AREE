#!/usr/bin/env python3
"""Write the gene_ids of chosen biotypes from an NCBI GTF, for --exclude_genes.

Total-RNA (not rRNA-depleted, not poly-A selected) libraries can be dominated
by rRNA: 65-97% of gene counts in ZHONG2025_PONV1_MORTALITY. The RNA-seq
workflow's --exclude_genes drops listed genes before DESeq2, so they never
appear as results. This script builds such a list from the annotation the run uses, so
the ids match the tx2gene map: NCBI GeneID numbers, the same ids as
tx2gene.tsv's second column.

Usage:
  python scripts/list_genes_by_biotype.py \\
      --gtf data/reference/GCF_963853765.1/genomic.gtf --biotype rRNA \\
      [--tx2gene data/reference/GCF_963853765.1/tx2gene.tsv] \\
      --out data/reference/gene_lists/GCF_963853765.1-RS_2024_06_rRNA.txt

The output starts with '#' lines recording the GTF, its sha256, the biotypes
and the date, then one gene_id per line.
"""
from __future__ import annotations

import argparse
import datetime
import hashlib
import re
import sys
from pathlib import Path

GENE_ID = re.compile(r'db_xref "GeneID:(\d+)"')
BIOTYPE = re.compile(r'gene_biotype "([^"]+)"')
SYMBOL = re.compile(r'gene "([^"]+)"')


def genes_by_biotype(gtf: Path, biotypes: set[str]) -> dict[str, str]:
    """GeneID -> gene symbol for every 'gene' record whose biotype is listed."""
    out: dict[str, str] = {}
    with open(gtf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue
            bt = BIOTYPE.search(fields[8])
            gid = GENE_ID.search(fields[8])
            if bt and gid and bt.group(1) in biotypes:
                sym = SYMBOL.search(fields[8])
                out[gid.group(1)] = sym.group(1) if sym else ""
    return out


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--gtf", type=Path, required=True)
    p.add_argument("--biotype", action="append", required=True, help="gene_biotype to list (repeatable)")
    p.add_argument("--tx2gene", type=Path, help="report how many listed genes the tx2gene map contains")
    p.add_argument("--out", type=Path, required=True)
    args = p.parse_args()

    genes = genes_by_biotype(args.gtf, set(args.biotype))
    if not genes:
        print(f"ERROR: no genes with biotype {args.biotype} in {args.gtf}")
        return 1
    header = [
        f"# genes with gene_biotype {', '.join(sorted(set(args.biotype)))}",
        f"# source GTF: {args.gtf.name} sha256 {sha256(args.gtf)}",
        f"# generated {datetime.date.today().isoformat()} by scripts/list_genes_by_biotype.py",
        f"# {len(genes)} genes; ids are NCBI GeneIDs as in tx2gene.tsv",
    ]
    if args.tx2gene:
        with open(args.tx2gene) as fh:
            in_map = {line.split("\t")[1].strip() for line in fh if "\t" in line}
        n = sum(g in in_map for g in genes)
        header.append(f"# {n} of {len(genes)} are in {args.tx2gene.name}")
        print(f"{n} of {len(genes)} listed genes are in {args.tx2gene}")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text("\n".join(header) + "\n"
                        + "\n".join(f"{g}\t# {genes[g]}" for g in sorted(genes, key=int)) + "\n")
    print(f"wrote {len(genes)} gene_ids to {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
