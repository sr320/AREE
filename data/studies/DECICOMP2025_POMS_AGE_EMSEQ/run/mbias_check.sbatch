#!/bin/bash -l
#SBATCH --job-name=aree-mbias-DECICOMP2025
#SBATCH --account=srlab
#SBATCH --partition=ckpt-all
#SBATCH --requeue
#SBATCH --open-mode=append
#SBATCH --cpus-per-task=12
#SBATCH --mem=48G
#SBATCH --time=4:30:00
#
# M-bias and mapping check on the first 1,000,000 read pairs of one library,
# to set methylation.extractor_args in klone.config before the full run.
# Same tools and versions as the workflow (Trim Galore 0.6.10, Bismark
# 0.24.2 containers) and the study's bismark_align_args.
#
# Everything is written to group storage, not scrubbed: the subset is small,
# and the Bismark index it builds is kept for reuse. A requeued job skips
# the steps whose outputs exist.
#
#   sbatch --chdir=<repo> data/studies/DECICOMP2025_POMS_AGE_EMSEQ/run/mbias_check.sbatch [RUN]
set -euo pipefail

RUN="${1:-ERR15996532}"
PAIRS=1000000
MIRROR="${AREE_MIRROR_URL:-https://gannet.fish.washington.edu/v1_web/aree/raw/PRJEB105019}"
BASE="/mmfs1/gscratch/srlab/$USER/aree"
REF="$BASE/reference/GCF_963853765.1"
IDX="$REF/bismark"
OUT="$BASE/mbias/DECICOMP2025_POMS_AGE_EMSEQ/$RUN"
IMG="$BASE/apptainer-images"
export APPTAINER_CACHEDIR="$BASE/apptainer-cache" APPTAINER_TMPDIR="${TMPDIR:-/tmp}"
mkdir -p "$OUT" "$IMG" "$IDX" "$APPTAINER_CACHEDIR"
cd "$OUT"

pull() { [[ -f "$IMG/$1.sif" ]] || apptainer pull "$IMG/$1.sif" "docker://quay.io/biocontainers/$2"; }
pull bismark     bismark:0.24.2--hdfd78af_0
pull trim-galore trim-galore:0.6.10--hdfd78af_0
run() { local img="$1"; shift; apptainer exec --bind /mmfs1 "$IMG/$img.sif" "$@"; }

echo "== $(date -Is) subset $RUN, first $PAIRS pairs"
for m in 1 2; do
  f="${RUN}_sub_$m.fastq.gz"
  if [[ ! -s "$f" ]]; then
    # head closes the pipe once it has its lines; curl then exits non-zero.
    (curl -sS -L --fail "$MIRROR/${RUN}_$m.fastq.gz" || true) | zcat 2>/dev/null \
      | head -n $((PAIRS * 4)) | gzip > "$f.tmp" || true
    [[ $(zcat "$f.tmp" | wc -l) -eq $((PAIRS * 4)) ]] || { echo "subset of mate $m is short"; exit 1; }
    mv "$f.tmp" "$f"
  fi
done

if [[ ! -f "$IDX/INDEX_OK" ]]; then
  echo "== $(date -Is) bismark_genome_preparation"
  ln -sfn "$REF/genomic.fna" "$IDX/genomic.fa"
  run bismark bismark_genome_preparation --bowtie2 --parallel 6 "$IDX"
  date -Is > "$IDX/INDEX_OK"
fi

echo "== $(date -Is) trim_galore"
[[ -s "${RUN}_sub_1_val_1.fq.gz" ]] || run trim-galore trim_galore --paired --cores 4 "${RUN}_sub_1.fastq.gz" "${RUN}_sub_2.fastq.gz"

echo "== $(date -Is) bismark (directional, --score_min L,0,-0.6)"
BAM="${RUN}_sub_1_val_1_bismark_bt2_pe.bam"
[[ -s "$BAM" && -s "${BAM%.bam}_PE_report.txt" ]] || run bismark bismark --genome "$IDX" --score_min L,0,-0.6 \
  --multicore 3 --temp_dir "$OUT/tmp" -1 "${RUN}_sub_1_val_1.fq.gz" -2 "${RUN}_sub_2_val_2.fq.gz"

echo "== $(date -Is) M-bias"
run bismark bismark_methylation_extractor --paired-end --mbias_only --multicore 4 "$BAM"

grep -E "Sequence pairs analysed|unique best hit|Mapping efficiency|C methylated in" "${BAM%.bam}_PE_report.txt"
echo "== $(date -Is) done: $OUT/${BAM%.bam}.M-bias.txt"
