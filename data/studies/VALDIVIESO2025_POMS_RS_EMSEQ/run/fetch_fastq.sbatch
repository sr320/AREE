#!/bin/bash -l
#SBATCH --job-name=aree-fetch-VALDIVIESO2025
#SBATCH --account=coenv
#SBATCH --partition=compute
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=2-00:00:00
#
# Fetch the 80 FASTQ files (40 runs, 624.9 GB) listed in fastq_manifest.tsv
# from ENA to scrubbed, and accept each one only at the manifest's MD5 and
# size. Writes $RAW/FASTQ_OK when every file verifies, so the run job can be
# chained with --dependency=afterok.
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
  for attempt in $(seq 1 20); do
    if [[ -f "$f" && $(stat -c %s "$f") -eq $bytes ]]; then
      if [[ $(md5sum "$f" | cut -d' ' -f1) == "$md5" ]]; then
        echo "OK $(basename "$f")"; return 0
      fi
      echo "BADMD5 $(basename "$f") (attempt $attempt); refetching"; rm -f "$f"
    fi
    curl -sS -L --fail --retry 5 --retry-delay 30 -C - -o "$f" "$https" \
      || { echo "curl failed $(basename "$f") (attempt $attempt)"; sleep 30; }
  done
  echo "FAILED $(basename "$f")"; return 1
}
export -f fetch_one

tail -n +2 "$MANIFEST" | cut -f3-5 \
  | xargs -P "$JOBS" -L 1 bash -c 'fetch_one "$0" "$1" "$2" "'"$RAW"'"'

n=$(tail -n +2 "$MANIFEST" | wc -l)
ok=0
while IFS=$'\t' read -r _ _ url md5 bytes _; do
  f="$RAW/$(basename "$url")"
  [[ -f "$f" && $(stat -c %s "$f") -eq $bytes && $(md5sum "$f" | cut -d' ' -f1) == "$md5" ]] && ok=$((ok+1))
done < <(tail -n +2 "$MANIFEST")
echo "verified $ok of $n files"
[[ $ok -eq $n ]] && date -Is > "$RAW/FASTQ_OK"
