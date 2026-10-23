#!/usr/bin/env python3
"""Per-segment PONV1 read counts for ZHONG2025_PONV1_MORTALITY.

The workflow's VIRAL_SUMMARY sums both PONV1 segments. The prespecified rule
(study YAML header) also needs both segments detected, so this reads each
library's VIRAL_QUANT quant.sf from the Nextflow work dir and joins it to the
sample sheet, the summed viral-load table and the flowcell lane in each
library's first read header.

Usage: ponv1_per_segment.py <nextflow log name/workdir TSV> <samplesheet.tsv>
                            <viral_load.tsv> <lanes.tsv: run<TAB>flowcell:lane>
The first argument comes from:
  nextflow log <session> -f name,workdir,status -filter 'process == "VIRAL_QUANT"'
"""
import csv
import os
import re
import sys

FARM = {"StnI", "StnL", "StnM"}
SEG1, SEG2 = "PQ030832.1", "PV579160.1"
THRESHOLD_PER_MILLION = 100  # prespecified: --viral_detection_per_million default

dirs = {}
for line in open(sys.argv[1]):
    if line.strip():
        name, workdir = line.split("\t")[:2]
        dirs[re.search(r"\((\w+)\)", name).group(1)] = workdir
sheet = {r["run_accession"]: r for r in csv.DictReader(open(sys.argv[2]), delimiter="\t")}
summed = {r["library_id"]: r for r in csv.DictReader(open(sys.argv[3]), delimiter="\t")}
lanes = dict(line.rstrip("\n").split("\t") for line in open(sys.argv[4]) if line.strip())

w = csv.writer(sys.stdout, delimiter="\t", lineterminator="\n")
w.writerow(["run_accession", "animal", "station", "arm", "flowcell_lane", "pairs_processed",
            "seg1_PQ030832_pairs", "seg2_PV579160_pairs", "ponv1_pairs", "ponv1_per_million",
            "both_segments_detected", "prespecified_call", "workflow_summed_call"])
for run in sorted(dirs):
    quant = {r["Name"]: float(r["NumReads"]) for r in
             csv.DictReader(open(os.path.join(dirs[run], "viral_quant", "quant.sf")), delimiter="\t")}
    s = sheet[run]
    v = summed[run]
    n = int(v["pairs_processed"])
    a, b = quant[SEG1], quant[SEG2]
    pm = (a + b) / n * 1e6
    both = a > 0 and b > 0
    call = "positive" if pm >= THRESHOLD_PER_MILLION and both else "negative"
    w.writerow([run, s["oyster"], s["station"], "farm" if s["station"] in FARM else "wild",
                lanes.get(run, ""), n, f"{a:.0f}", f"{b:.0f}", f"{a + b:.0f}", f"{pm:.2f}",
                both, call, v["virus_detected"]])
