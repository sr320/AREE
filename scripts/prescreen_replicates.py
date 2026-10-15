#!/usr/bin/env python3
"""Pre-screen a study's replicates on a subsample of reads, before downloading it.

Two real studies (IOCAS2022_OA_ENERGY, WANG2024_HEAT_GILL) passed every
metadata check and failed `check_replicate_dispersion.py` only after a full
download and reanalysis: their "biological replicates" varied like technical
replicates. The same statistic can be computed on the first ~1 M reads of
each library, a few percent of the download.

For every run in the chosen condition levels this script fetches only the
first --n-reads reads (pairs, for paired-end), quantifies them with Salmon
against an existing index, and computes the same per-group median dispersion
as `check_replicate_dispersion.py` (whose functions it imports). Nothing is
written outside --workdir and --out.

Sources:
  ena   stream the start of each FASTQ in fastq_manifest.tsv (curl | head).
        Same files the full download would use, but ENA can be very slow.
  ncbi  `fastq-dump -X N` (sra-tools) reads the first N spots from NCBI's
        copy, usually much faster. Screening only: no ENA MD5 provenance.

Usage:
  python scripts/prescreen_replicates.py \\
      --samplesheet data/studies/STUDY/samplesheet.tsv \\
      --manifest data/studies/STUDY/fastq_manifest.tsv \\
      --salmon-index data/reference/GCF_963853765.1/salmon_index \\
      [--level control_level --level treatment_level] [--n-reads 1000000] \\
      [--source ncbi] [--jobs 3] [--workdir /tmp/prescreen_STUDY] [--out prescreen.tsv]

Exits 1 if any group's dispersion is below --min-dispersion. A pass is not a
substitute for running `check_replicate_dispersion.py` on the full reanalysis.
Reads from the start of a run are not a random sample (early flow-cell tiles),
which matters little for a within-group dispersion comparison but is a caveat.
"""
from __future__ import annotations

import argparse
import importlib.util
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
_spec = importlib.util.spec_from_file_location("check_replicate_dispersion", HERE / "check_replicate_dispersion.py")
crd = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(crd)

DEFAULT_N_READS = 1_000_000
# check_replicate_dispersion keeps genes with mean normalized count > 500 at
# full depth (20-30 M reads). At 1 M reads the same genes have about a
# twentieth of the counts; 25 keeps roughly the same set of well-expressed genes.
DEFAULT_MIN_MEAN = 25


def https(url: str) -> str:
    return "https://" + url.split("://", 1)[1] if url.startswith("ftp://") else url


def _count_records(path: Path) -> int:
    cat = "gzip -dc" if path.suffix == ".gz" else "cat"
    out = subprocess.run(["bash", "-c", f"{cat} '{path}' | wc -l"], capture_output=True, text=True).stdout
    return int(out.strip() or 0) // 4


def stream_head(url: str, n_reads: int, dest: Path, attempts: int = 3) -> int:
    """Write the first n_reads FASTQ records of a gzipped URL to dest (gzipped).

    Returns the number of records written. Fewer than n_reads is accepted only
    if the stream ended cleanly (the whole file is shorter than n_reads).
    """
    for attempt in range(1, attempts + 1):
        cmd = (f"set -o pipefail; curl -fsSL --retry 3 --connect-timeout 60 '{https(url)}' "
               f"| gzip -dc 2>/dev/null | head -n {4 * n_reads} | gzip -1 > '{dest}'")
        proc = subprocess.run(["bash", "-c", cmd], capture_output=True, text=True)
        got = _count_records(dest)
        # curl and gzip exit non-zero (SIGPIPE) when head stops reading early: expected.
        if got == n_reads or (proc.returncode == 0 and got > 0):
            return got
        print(f"  {dest.name}: attempt {attempt} got {got} reads; retrying", file=sys.stderr)
    raise RuntimeError(f"could not stream {n_reads} reads from {url}")


def fetch_ena(run: str, files: pd.DataFrame, n_reads: int, reads_dir: Path) -> list[tuple[Path, int]]:
    out = []
    for _, f in files.sort_values("mate").iterrows():
        dest = reads_dir / f"{run}_{f['mate']}.fastq.gz"
        out.append((dest, stream_head(f["url"], n_reads, dest)))
    return out


def fetch_ncbi(run: str, n_mates: int, n_reads: int, reads_dir: Path, attempts: int = 3) -> list[tuple[Path, int]]:
    split = ["--split-files"] if n_mates == 2 else []
    for attempt in range(1, attempts + 1):
        proc = subprocess.run(["fastq-dump", "-X", str(n_reads), *split, "--skip-technical",
                               "-O", str(reads_dir), run], capture_output=True, text=True)
        names = [f"{run}_1.fastq", f"{run}_2.fastq"] if n_mates == 2 else [f"{run}.fastq"]
        paths = [reads_dir / n for n in names]
        if proc.returncode == 0 and all(p.exists() for p in paths):
            return [(p, _count_records(p)) for p in paths]
        print(f"  {run}: fastq-dump attempt {attempt} failed: {proc.stderr.strip()[-200:]}", file=sys.stderr)
    raise RuntimeError(f"fastq-dump could not fetch {run}")


def quantify(run: str, fastqs: list[Path], index: Path, out: Path, threads: int) -> float:
    reads = ["-1", str(fastqs[0]), "-2", str(fastqs[1])] if len(fastqs) == 2 else ["-r", str(fastqs[0])]
    subprocess.run(["salmon", "quant", "-i", str(index), "-l", "A", *reads, "-p", str(threads),
                    "--validateMappings", "-q", "-o", str(out / run)],
                   check=True, capture_output=True)
    meta = pd.read_json(out / run / "aux_info" / "meta_info.json", typ="series")
    return float(meta["percent_mapped"])


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--samplesheet", type=Path, required=True)
    p.add_argument("--manifest", type=Path, required=True, help="fastq_manifest.tsv from aree fetch-samplesheet")
    p.add_argument("--salmon-index", type=Path, required=True)
    p.add_argument("--tx2gene", type=Path, default=crd.DEFAULT_TX2GENE)
    p.add_argument("--level", action="append", default=[],
                   help="condition level to screen (repeatable; default: every level with >= 2 runs)")
    p.add_argument("--n-reads", type=int, default=DEFAULT_N_READS, help="reads (pairs) per library")
    p.add_argument("--source", choices=["ena", "ncbi"], default="ena")
    p.add_argument("--jobs", type=int, default=3, help="runs fetched in parallel")
    p.add_argument("--min-mean", type=float, default=DEFAULT_MIN_MEAN)
    p.add_argument("--min-dispersion", type=float, default=crd.DEFAULT_MIN_DISPERSION)
    p.add_argument("--threads", type=int, default=4, help="Salmon threads")
    p.add_argument("--workdir", type=Path, default=Path("prescreen_work"))
    p.add_argument("--keep-reads", action="store_true", help="keep the subsampled FASTQs")
    p.add_argument("--out", type=Path)
    args = p.parse_args()

    needed = ["salmon"] + (["fastq-dump"] if args.source == "ncbi" else ["curl"])
    missing = [t for t in needed if not shutil.which(t)]
    if missing:
        print(f"ERROR: not on PATH: {', '.join(missing)} (activate the aree environment).")
        return 2
    sheet = pd.read_csv(args.samplesheet, sep="\t", dtype=str).set_index("run_accession")
    manifest = pd.read_csv(args.manifest, sep="\t", dtype=str)
    levels = args.level or [lv for lv, n in sheet["condition"].value_counts().items() if n >= 2]
    runs = list(sheet.index[sheet["condition"].isin(levels)])
    if not runs:
        print(f"ERROR: no runs in levels {levels}")
        return 2

    reads_dir, salmon_dir = args.workdir / "reads", args.workdir / "salmon"
    reads_dir.mkdir(parents=True, exist_ok=True)
    salmon_dir.mkdir(parents=True, exist_ok=True)
    todo = [r for r in runs if not (salmon_dir / r / "quant.sf").exists()]
    for r in sorted(set(runs) - set(todo)):
        print(f"{r}: already quantified, reusing")

    def fetch(run: str) -> list[tuple[Path, int]]:
        files = manifest[manifest["run_accession"] == run]
        if args.source == "ncbi":
            return fetch_ncbi(run, len(files), args.n_reads, reads_dir)
        return fetch_ena(run, files, args.n_reads, reads_dir)

    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        for run, fetched in zip(todo, pool.map(fetch, todo)):
            n = {count for _, count in fetched}
            if len(n) != 1:
                print(f"ERROR: {run} mates gave different read counts {sorted(n)}")
                return 2
            fastqs = [path for path, _ in fetched]
            mapped = quantify(run, fastqs, args.salmon_index, salmon_dir, args.threads)
            print(f"{run} ({sheet.loc[run, 'condition']}): {n.pop():,} reads, {mapped:.1f}% mapped")
            if not args.keep_reads:
                for fq in fastqs:
                    fq.unlink()

    tx2gene = pd.read_csv(args.tx2gene, sep="\t", header=None, names=["tx", "gene"], dtype=str)
    sub = crd.library_ids(sheet.loc[runs])
    counts = crd.gene_counts([salmon_dir], sub, tx2gene)
    groups = sub.set_index("sample_id").loc[counts.columns, "condition"]
    table = crd.group_dispersion(counts, groups, min_mean=args.min_mean)
    table["passes"] = table["median_dispersion"] >= args.min_dispersion
    table.insert(1, "reads_per_library", args.n_reads)
    table.insert(2, "source", args.source)
    print(f"\nPRE-SCREEN (first {args.n_reads:,} reads per library, from {args.source})")
    print(table.to_string(index=False))
    if args.out:
        table.to_csv(args.out, sep="\t", index=False)
    if not table["passes"].all():
        print(f"\nFAIL: dispersion below {args.min_dispersion} in "
              f"{', '.join(table.loc[~table['passes'], 'condition'])}. "
              "Do not download this study as is: its replicates look technical, not biological.")
        return 1
    print("\nPASS on the subsample. Run check_replicate_dispersion.py again after the full reanalysis.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
