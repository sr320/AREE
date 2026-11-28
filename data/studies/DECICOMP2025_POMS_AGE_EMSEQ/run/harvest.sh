#!/bin/bash
# Copy one comparison's kept outputs off scrubbed (dereferencing the
# symlinks published into --outdir) into the study's data directory in this
# checkout, under the names the other methylation studies commit, plus the
# run's trace into run/pipeline_info/. Run from the repository root.
#
#   harvest.sh <outdir> <comparison_id>
#
# Fails if any required file is missing, so a purge or a failed comparison
# is noticed at once rather than at commit time.
set -euo pipefail
OUT="${1:?outdir}"; CMP="${2:?comparison_id}"
STUDY=DECICOMP2025_POMS_AGE_EMSEQ
M="$OUT/methylation"; P="${STUDY}_${CMP}"
DEST="data/studies/$STUDY"
mkdir -p "$DEST/run/pipeline_info/$CMP"

cp -L "$M/standardized/${P}_dmr_standardized.tsv" "$DEST/${P}_dmr_standardized.tsv"
gzip -f -n "$DEST/${P}_dmr_standardized.tsv"
cp -L "$M/standardized/${P}_standardize_qc.tsv"   "$DEST/${P}_standardize_qc.tsv"
cp -L "$M/manifest/${P}_manifest.json"            "$DEST/${P}_workflow_manifest.json"
cp -L "$M/report/${P}_report.html"                "$DEST/${P}_methylation_report.html"
cp -L "$M/dmr/${P}_methylkit_qc.tsv"              "$DEST/${P}_dmr_qc.tsv"
cp -L "$M/annotated/${P}_annotation_qc.tsv"       "$DEST/${P}_annotation_qc.tsv"
cp -L "$M/manifest/versions.yml"                  "$DEST/run/pipeline_info/$CMP/versions.yml"
cp -L "$OUT/pipeline_info_${CMP}"/* "$DEST/run/pipeline_info/$CMP/" 2>/dev/null || true
echo "harvested $CMP -> $DEST"
