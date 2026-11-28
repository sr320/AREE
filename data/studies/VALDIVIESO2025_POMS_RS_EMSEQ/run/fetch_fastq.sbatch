#!/bin/bash -l
#SBATCH --job-name=aree-fetch-VALDIVIESO2025
#SBATCH --account=srlab
#SBATCH --partition=ckpt-all
#SBATCH --requeue
#SBATCH --open-mode=append
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=2-00:00:00
#
# Fetch the 80 FASTQ files (40 runs, 624.9 GB) listed in fastq_manifest.tsv
# from ENA to scrubbed, and accept each one only at the manifest's MD5 and
# size. Writes $RAW/FASTQ_OK when every file verifies, so the run job can be
# chained with --dependency=afterok.
#
# Runs on the checkpoint partition: a preempted job is requeued and starts
# again from the top, which is safe because every file is re-verified and a
# partial file is resumed. The first attempt (job 43496224, coenv/compute,
# 2026-11-23..25) ended with 61 of 80 files failed: ENA dropped transfers
# (curl 18), then hours of write errors on scrubbed (curl 23, cause not
# established; scrubbed was 97% full) used up the 20 attempts per file. So a
# write error now logs free space and waits 10 minutes, a file longer than
# the manifest size is deleted rather than resumed (curl 22 / HTTP 416
# otherwise repeats), and there are 40 attempts per file.
#
#   sbatch --chdir=<repo> data/studies/VALDIVIESO2025_POMS_RS_EMSEQ/run/fetch_fastq.sbatch
set -euo pipefail

STUDY=VALDIVIESO2025_POMS_RS_EMSEQ
MANIFEST="data/studies/$STUDY/fastq_manifest.tsv"
RAW="${AREE_RAW:-/gscratch/scrubbed/$USER/aree/raw/$STUDY}"
JOBS="${FETCH_JOBS:-8}"
mkdir -p "$RAW"
rm -f "$RAW/FASTQ_OK"

fetch_one() {
  local url="$1" md5="$2" bytes="$3" raw="$4"
  local f="$raw/$(basename "$url")"
  local https="${url/ftp:\/\//https://}"
  local rc
  for attempt in $(seq 1 40); do
    if [[ -f "$f" && $(stat -c %s "$f") -eq $bytes ]]; then
      if [[ $(md5sum "$f" | cut -d' ' -f1) == "$md5" ]]; then
        echo "OK $(basename "$f")"; return 0
      fi
      echo "BADMD5 $(basename "$f") (attempt $attempt); refetching"; rm -f "$f"
    elif [[ -f "$f" && $(stat -c %s "$f") -gt $bytes ]]; then
      echo "TOOLONG $(basename "$f") (attempt $attempt); refetching"; rm -f "$f"
    fi
    rc=0
    curl -sS -L --fail --retry 5 --retry-delay 30 -C - -o "$f" "$https" || rc=$?
    if [[ $rc -ne 0 ]]; then
      echo "$(date -Is) curl exit $rc $(basename "$f") (attempt $attempt)"
      if [[ $rc -eq 23 ]]; then
        df -h "$raw" | tail -1
        sleep 600
      else
        sleep 30
      fi
    fi
  done
  echo "FAILED $(basename "$f")"; return 1
}
export -f fetch_one

tail -n +2 "$MANIFEST" | cut -f3-5 \
  | xargs -P "$JOBS" -L 1 bash -c 'fetch_one "$0" "$1" "$2" "'"$RAW"'"' \
  || echo "some files failed; verifying what is on disk"

n=$(tail -n +2 "$MANIFEST" | wc -l)
ok=0
while IFS=$'\t' read -r _ _ url md5 bytes _; do
  f="$RAW/$(basename "$url")"
  [[ -f "$f" && $(stat -c %s "$f") -eq $bytes && $(md5sum "$f" | cut -d' ' -f1) == "$md5" ]] && ok=$((ok+1))
done < <(tail -n +2 "$MANIFEST")
echo "verified $ok of $n files"
[[ $ok -eq $n ]] || exit 1
date -Is > "$RAW/FASTQ_OK"
