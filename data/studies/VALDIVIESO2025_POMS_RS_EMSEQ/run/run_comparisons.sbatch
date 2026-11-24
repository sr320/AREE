#!/bin/bash -l
#SBATCH --job-name=aree-VALDIVIESO2025_POMS_RS_EMSEQ
#SBATCH --account=coenv
#SBATCH --partition=cpu-g2
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=4-00:00:00
#
# Nextflow head job: run the eight comparisons in comparisons.tsv in order,
# each with -resume so later comparisons reuse the alignments, and harvest
# each comparison's outputs into the checkout (group storage) as soon as it
# finishes. Tasks go to ckpt-all through -profile apptainer,klone.
#
#   sbatch --chdir=<repo> [--dependency=afterok:<fetch job>] \
#     data/studies/VALDIVIESO2025_POMS_RS_EMSEQ/run/run_comparisons.sbatch [comparison_id ...]
#
# With comparison ids, runs only those (in comparisons.tsv order).
set -euo pipefail

STUDY=VALDIVIESO2025_POMS_RS_EMSEQ
RUNDIR="data/studies/$STUDY/run"
ROOT="${AREE_NXF_ROOT:-/gscratch/scrubbed/$USER/aree}"
RAW="${AREE_RAW:-$ROOT/raw/$STUDY}"
OUT="$ROOT/results/$STUDY"
LOGS="$ROOT/runs/$STUDY"
CFGDIR="$LOGS/configs"

[[ -f "$RAW/FASTQ_OK" ]] || { echo "No $RAW/FASTQ_OK: FASTQs not verified"; exit 1; }

source "${AREE_MINIFORGE:-/gscratch/srlab/$USER/miniforge3}/etc/profile.d/conda.sh"
conda activate "${AREE_NF_ENV:-myflow}"
# The myflow env was moved from ~/miniforge3 to group storage; its launcher
# still points NXF_DIST at the old path, so name the bundled jar directly.
jar=$(ls "$CONDA_PREFIX"/share/nextflow/dist/*/nextflow-*-one.jar 2>/dev/null | tail -1)
[[ -n "$jar" ]] && export NXF_BIN="$jar"
export NXF_HOME="$ROOT/nxf-home" NXF_TEMP="$ROOT/tmp"
unset TMPDIR
export AREE_SLURM_ACCOUNT="${AREE_SLURM_ACCOUNT:-srlab}"
mkdir -p "$NXF_HOME" "$NXF_TEMP" "$OUT" "$CFGDIR"

# One directory of symlinks per tissue, so each comparison aligns only its
# own 20 libraries (as in the first run) with the plain ERR*_{1,2} glob.
for tissue in G M; do
  td="$ROOT/raw/${STUDY}_$tissue"
  mkdir -p "$td"
  awk -F'\t' -v t="$tissue" 'NR==1{for(i=1;i<=NF;i++)h[$i]=i; next}
       $h["tissue_code"]==t{print $h["run_accession"]}' "data/studies/$STUDY/samplesheet.tsv" \
    | while read -r run; do
        for m in 1 2; do ln -sfn "$RAW/${run}_$m.fastq.gz" "$td/${run}_$m.fastq.gz"; done
      done
done

mapfile -t rows < <(tail -n +2 "$RUNDIR/comparisons.tsv")
for row in "${rows[@]}"; do
  IFS=$'\t' read -r cmp tissue sheet <<< "$row"
  if [[ $# -gt 0 && ! " $* " =~ " $cmp " ]]; then continue; fi
  cfg="$CFGDIR/comparison_$cmp.config"
  cat > "$cfg" <<CFG
params.reads = '$ROOT/raw/${STUDY}_$tissue/ERR*_{1,2}.fastq.gz'
params.methylation.comparison_id = '$cmp'
params.methylation.sample_sheet  = '$sheet'
CFG
  echo "== $(date -Is) $cmp"
  nextflow -log "$LOGS/nextflow_${cmp}_$(date +%Y%m%d_%H%M%S).log" run workflows/methylation \
    -config "config/$STUDY.config" \
    -c "$RUNDIR/klone.config" -c "$cfg" \
    -profile apptainer,klone -resume \
    --outdir "$OUT" --tracedir "$OUT/pipeline_info_$cmp" < /dev/null
  bash "$RUNDIR/harvest.sh" "$OUT" "$cmp"
done
echo "== $(date -Is) all requested comparisons done and harvested"
