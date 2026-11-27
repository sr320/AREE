#!/bin/bash -l
#SBATCH --job-name=aree-fetch-DECICOMP2025
#SBATCH --account=srlab
#SBATCH --partition=ckpt-all
#SBATCH --requeue
#SBATCH --open-mode=append
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=2-00:00:00
#
# Fetch the 120 FASTQ files (60 runs, 710.5 GB) listed in fastq_manifest.tsv
# to scrubbed, and accept each one only at the manifest's MD5 and size.
# Writes $RAW/FASTQ_OK when every file verifies, so the run job can be
# chained with --dependency=afterok.
#
# SOURCE. The lab mirror on gannet (aree mirror-fastq, docs/fastq_mirror.md)
# holds all 120 files as ENA served them: on 2026-11-26 every MD5 and size in
# its ledger (raw/PRJEB105019/mirror_manifest.tsv) equalled this manifest's.
# So each file is tried from the mirror first, and from ENA over HTTPS after
# 20 failed attempts. Either way the check is the manifest's ENA MD5.
#
# ERR15996548 has no ENA-generated FASTQ; the manifest lists the submitter's
# two files (F14R_C3_T0_S3_DNA_R{1,2}.fastq.gz). They are fetched under those
# names and linked as ERR15996548_{1,2}.fastq.gz for the config's reads glob.
#
# Runs on the checkpoint partition: a preempted job is requeued and starts
# again from the top, which is safe because every file is re-verified and a
# partial file is resumed. Retry handling follows
# VALDIVIESO2025_POMS_RS_EMSEQ/run/fetch_fastq.sbatch: a write error (curl
# 23, seen when scrubbed was 99% full) logs free space and waits 10 minutes,
# and a file longer than the manifest size is deleted rather than resumed.
#
#   sbatch --chdir=<repo> data/studies/DECICOMP2025_POMS_AGE_EMSEQ/run/fetch_fastq.sbatch
set -euo pipefail

STUDY=DECICOMP2025_POMS_AGE_EMSEQ
MANIFEST="data/studies/$STUDY/fastq_manifest.tsv"
RAW="${AREE_RAW:-/gscratch/scrubbed/$USER/aree/raw/$STUDY}"
MIRROR="${AREE_MIRROR_URL:-https://gannet.fish.washington.edu/v1_web/aree/raw/PRJEB105019}"
JOBS="${FETCH_JOBS:-8}"
mkdir -p "$RAW"
rm -f "$RAW/FASTQ_OK"

fetch_one() {
  local url="$1" md5="$2" bytes="$3" raw="$4" mirror="$5"
  local name; name="$(basename "$url")"
  local f="$raw/$name"
  local src rc
  for attempt in $(seq 1 40); do
    if [[ -f "$f" && $(stat -c %s "$f") -eq $bytes ]]; then
      if [[ $(md5sum "$f" | cut -d' ' -f1) == "$md5" ]]; then
        echo "OK $name"; return 0
      fi
      echo "BADMD5 $name (attempt $attempt); refetching"; rm -f "$f"
    elif [[ -f "$f" && $(stat -c %s "$f") -gt $bytes ]]; then
      echo "TOOLONG $name (attempt $attempt); refetching"; rm -f "$f"
    fi
    if [[ $attempt -le 20 ]]; then src="$mirror/$name"; else src="${url/ftp:\/\//https://}"; fi
    rc=0
    curl -sS -L --fail --retry 5 --retry-delay 30 -C - -o "$f" "$src" || rc=$?
    if [[ $rc -ne 0 ]]; then
      echo "$(date -Is) curl exit $rc $name (attempt $attempt, $src)"
      if [[ $rc -eq 23 ]]; then
        df -h "$raw" | tail -1
        sleep 600
      else
        sleep 30
      fi
    fi
  done
  echo "FAILED $name"; return 1
}
export -f fetch_one

echo "== $(date -Is) fetch start on $(hostname), job ${SLURM_JOB_ID:-none}, $JOBS at a time"
tail -n +2 "$MANIFEST" | cut -f3,4,5 \
  | xargs -P "$JOBS" -n 3 bash -c 'fetch_one "$0" "$1" "$2" "'"$RAW"'" "'"$MIRROR"'"' \
  || echo "some files failed; verifying what is on disk"

# Final check of every file against the manifest, then the run-named links
# for files the manifest lists under another name (source = submitted_ftp).
n=0; ok=0
while IFS=$'\t' read -r run mate url md5 bytes source; do
  n=$((n + 1))
  f="$RAW/$(basename "$url")"
  if [[ -f "$f" && $(stat -c %s "$f") -eq $bytes && $(md5sum "$f" | cut -d' ' -f1) == "$md5" ]]; then
    ok=$((ok + 1))
    want="${run}_${mate}.fastq.gz"
    [[ "$(basename "$f")" == "$want" ]] || ln -sfn "$(basename "$f")" "$RAW/$want"
  else
    echo "NOT VERIFIED $(basename "$f")"
  fi
done < <(tail -n +2 "$MANIFEST")
echo "== $(date -Is) verified $ok of $n files"
[[ $ok -eq $n ]] || exit 1
date -Is > "$RAW/FASTQ_OK"
