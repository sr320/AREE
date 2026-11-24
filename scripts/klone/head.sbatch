#!/bin/bash -l
# Nextflow head job for one AREE study on Klone, as submitted by
# `aree autopilot launch` (which sets --job-name=aree-<STUDY_ID>, the
# account, partition, time, --chdir and the exported AREE_* variables).
#
#   head.sbatch <workflow dir> <study config> <outdir> <autopilot override> [run config ...]
#
# Always runs with -resume: on the first submission there is nothing to
# resume and Nextflow starts fresh; on a resubmission it continues from the
# cache in the work dir. Run configs are passed with -c in order, and the
# autopilot's override (memory escalations) last, when it exists.
set -euo pipefail

WORKFLOW="${1:?workflow dir}"
CONFIG="${2:?study config}"
OUTDIR="${3:?outdir}"
OVERRIDE="${4:?autopilot override path}"
shift 4

: "${AREE_MINIFORGE:?set by aree autopilot launch}"
: "${AREE_NF_ENV:=nextflow}"
: "${AREE_NXF_ROOT:=/gscratch/scrubbed/$USER/aree}"

# shellcheck disable=SC1091
source "$AREE_MINIFORGE/etc/profile.d/conda.sh"
conda activate "$AREE_NF_ENV"

export NXF_HOME="$AREE_NXF_ROOT/nxf-home"
# NXF_TEMP, never TMPDIR: see docs/running_on_klone.md#chunked-alignment.
export NXF_TEMP="$AREE_NXF_ROOT/tmp"
unset TMPDIR
mkdir -p "$NXF_HOME" "$NXF_TEMP" "$OUTDIR/pipeline_info"

extra=()
for cfg in "$@"; do
  extra+=(-c "$cfg")
done
if [[ -f "$OVERRIDE" ]]; then
  extra+=(-c "$OVERRIDE")
fi

nextflow run "$WORKFLOW" \
  -config "$CONFIG" \
  "${extra[@]}" \
  -profile apptainer,klone \
  -resume \
  --outdir "$OUTDIR" \
  --tracedir "$OUTDIR/pipeline_info"
