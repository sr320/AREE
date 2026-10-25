#!/usr/bin/env bash
# Run the three CALLA2026_OSHV comparisons that have not been run at full depth,
# then copy their outputs next to the study and register them.
# See docs/first_raw_reanalysis.md ("Running the remaining comparisons").
#
# Written for Chlo, the workstation that ran the first three comparisons and
# still holds the 42 FASTQs (verified against ENA's MD5s on 2026-09-30):
#
#   bash scripts/run_calla2026_remaining.sh                 # all three, in order
#   bash scripts/run_calla2026_remaining.sh miyagi_oshv1_france_vs_control
#
# Run it from the repository root of the study/PRJNA1329250/CALLA2026_OSHV
# checkout. For each comparison it:
#   1. checks, against fastq_manifest.tsv, every FASTQ the comparison uses
#      (size and MD5; SKIP_MD5=1 checks size only);
#   2. runs workflows/rnaseq with -profile local,workstation and an --outdir
#      under results/ (not /private/tmp, which lost a run's outputs before);
#   3. copies the standardized table, workflow manifest, viral-load table and
#      HTML report into data/studies/CALLA2026_OSHV/;
#   4. sets that comparison's results_file in the study YAML;
#   5. harmonizes the table with `aree harmonize --comparison` as a check.
# A comparison whose outputs are already committed is skipped.
#
# Environment overrides:
#   RAW_DIR       FASTQ directory (default data/raw/CALLA2026_OSHV)
#   AREE_WORKDIR  Nextflow work directory (default from -profile workstation:
#                 /Volumes/ST5/aree/work). Must be APFS or HFS+, not exFAT,
#                 and no path may contain a space.
#   NXF_EXTRA     extra arguments for nextflow run, e.g. "-c my_machine.config"
#   SKIP_MD5=1    check sizes only (MD5 of ~140 GB takes a while on an HDD)
set -euo pipefail

STUDY=CALLA2026_OSHV
D=data/studies/$STUDY
YAML=registry/studies/$STUDY.yaml
RAW_DIR=${RAW_DIR:-data/raw/$STUDY}
ALL="midori_oshv1_usa_vs_control miyagi_oshv1_australia_vs_control miyagi_oshv1_france_vs_control"

[ -f "$YAML" ] && [ -d workflows/rnaseq ] || { echo "run this from the AREE repository root"; exit 1; }
[ -d "$RAW_DIR" ] || { echo "no FASTQ directory at $RAW_DIR (set RAW_DIR, or: ln -sfn /Volumes/blue-block/aree/raw data/raw)"; exit 1; }
[ -f data/reference/GCF_963853765.1/rna.fna ] || { echo "no reference at data/reference/GCF_963853765.1 (aree fetch-reference --assembly GCF_963853765.1)"; exit 1; }
case "$PWD${AREE_WORKDIR:-}" in *" "*) echo "a path contains a space; Nextflow, Trim Galore and Quarto break on it"; exit 1 ;; esac

# The editable `aree` install may point at another checkout; use this one.
aree_here() { PYTHONPATH="$PWD/src${PYTHONPATH:+:$PYTHONPATH}" python -m aree.cli "$@"; }

# comparison_id -> sample-sheet levels and exclusions. The second population is
# spelled "Myagi" in the deposited BioSamples and so in the sample sheet.
levels() {
  case $1 in
    midori_oshv1_usa_vs_control)       echo "Midori_Control Midori_USA DO_CT_4" ;;
    miyagi_oshv1_australia_vs_control) echo "Myagi_Control Myagi_Australia -" ;;
    miyagi_oshv1_france_vs_control)    echo "Myagi_Control Myagi_France -" ;;
    *) echo "unknown comparison: $1" >&2; return 1 ;;
  esac
}

file_size() { wc -c < "$1" | tr -d ' '; }
file_md5() { if command -v md5sum >/dev/null; then md5sum "$1" | cut -d' ' -f1; else md5 -q "$1"; fi; }

verify_reads() {
  local ctrl=$1 trt=$2 bad=0 n=0
  local runs
  runs=$(awk -F'\t' -v c="$ctrl" -v t="$trt" 'NR>1 && ($5==c || $5==t) {print $2}' "$D/samplesheet.tsv")
  for r in $runs; do
    while IFS=$'\t' read -r run mate _url md5 bytes; do
      f=$RAW_DIR/${run}_${mate}.fastq.gz; n=$((n+1))
      if [ ! -s "$f" ]; then echo "  MISSING $f"; bad=1; continue; fi
      [ "$(file_size "$f")" = "$bytes" ] || { echo "  WRONG SIZE $f"; bad=1; continue; }
      if [ "${SKIP_MD5:-0}" != 1 ] && [ "$(file_md5 "$f")" != "$md5" ]; then echo "  WRONG MD5 $f"; bad=1; fi
    done < <(awk -F'\t' -v r="$r" 'NR>1 && $1==r' "$D/fastq_manifest.tsv")
  done
  echo "  checked $n FASTQs for $ctrl + $trt"
  return $bad
}

set_results_file() {
  # Replace the first `results_file: null` after this comparison's id line.
  python - "$YAML" "$1" "$2" <<'EOF'
import re, sys
path, cid, rel = sys.argv[1:]
text = open(path).read()
start = text.index(f'comparison_id: "{cid}"')
m = re.compile(r"results_file: null").search(text, start)
nxt = text.find("comparison_id:", start + 1)
if m is None or (nxt != -1 and m.start() > nxt):
    sys.exit(f"{cid}: results_file is not null in {path}; leaving it")
open(path, "w").write(text[:m.start()] + f'results_file: "{rel}"' + text[m.end():])
print(f"  set results_file for {cid}")
EOF
}

run_one() {
  local cid=$1 ctrl trt excl
  read -r ctrl trt excl < <(levels "$cid")
  local table=$D/${STUDY}_${cid}_dge_standardized.tsv
  echo "== $cid ($ctrl vs $trt, excluding: $excl)"
  if [ -f "$table" ]; then echo "  $table exists; skipping"; return 0; fi

  verify_reads "$ctrl" "$trt" || { echo "  FASTQs not ready; not running $cid"; return 1; }

  local out=results/rnaseq/$STUDY/$cid
  local excl_arg=()
  [ "$excl" != - ] && excl_arg=(--exclude_samples "$excl")
  mkdir -p "$out/pipeline_info"
  # RENDER_REPORT has failed on every earlier full-depth run of this study
  # (Quarto cleanup on an exFAT work directory), after every other task
  # finished. So a non-zero exit is checked against the outputs below rather
  # than stopping the script.
  set +e
  # shellcheck disable=SC2086
  nextflow run workflows/rnaseq \
    -c config/$STUDY.config -profile local,workstation -resume ${NXF_EXTRA:-} \
    --comparison_id "$cid" --control_level "$ctrl" --treatment_level "$trt" \
    "${excl_arg[@]+"${excl_arg[@]}"}" \
    --reads "$RAW_DIR/SRR*_{1,2}.fastq.gz" \
    --outdir "$out" \
    -with-report "$out/pipeline_info/execution_report.html" \
    -with-trace "$out/pipeline_info/trace.tsv"
  local status=$?
  set -e
  [ $status = 0 ] || echo "  WARN: nextflow exited $status; checking outputs"

  local src_table=$out/rnaseq/standardized/${STUDY}_${cid}_dge_standardized.tsv
  local src_manifest=$out/rnaseq/manifests/${STUDY}_${cid}_manifest.json
  local src_viral=$out/rnaseq/viral_load/${STUDY}_${cid}_viral_load.tsv
  local src_report=$out/rnaseq/reports/${STUDY}_${cid}_rnaseq_report.html
  for f in "$src_table" "$src_manifest" "$src_viral"; do
    [ -s "$f" ] || { echo "  FAILED: $f was not produced"; return 1; }
  done
  cp "$src_table" "$table"
  cp "$src_manifest" "$D/${STUDY}_${cid}_workflow_manifest.json"
  cp "$src_viral" "$D/${STUDY}_${cid}_viral_load.tsv"
  if [ -s "$src_report" ]; then
    cp "$src_report" "$D/${STUDY}_${cid}_rnaseq_report.html"
  else
    echo "  WARN: no HTML report; find it in the failed RENDER_REPORT task directory"
  fi
  set_results_file "$cid" "$table"
  aree_here harmonize --study "$STUDY" --comparison "$cid" --input "$table"
  echo "  DONE $cid"
}

for cid in ${*:-$ALL}; do
  levels "$cid" >/dev/null
  run_one "$cid"
done

echo
echo "Next: review the viral-load tables, update the CALLA2026 status in"
echo "docs/first_raw_reanalysis.md, docs/implementation_status.md and"
echo "tests/test_calla2026_study.py, run pytest, then commit data/studies/$STUDY/ and $YAML."
