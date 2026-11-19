#!/usr/bin/env bash
# Fallback watchdog for when scrontab is not available: run one tick, then
# queue the next one 30 minutes out. Start it once with
#
#   sbatch scripts/klone/autopilot.sbatch <miniforge prefix> [conda env]
#
# and stop it with `scancel -n aree-autopilot`.
#SBATCH --job-name=aree-autopilot
#SBATCH --account=srlab
#SBATCH --partition=cpu-g2
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:15:00
#SBATCH --output=/gscratch/srlab/%u/aree-autopilot/tick-%j.log
set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Under sbatch, BASH_SOURCE is Slurm's copy of this script; use the submit dir.
REPO="${SLURM_SUBMIT_DIR:-$SCRIPT_DIR/../..}"

# Queue the next tick first, so a failing tick does not end the chain.
sbatch --begin=now+30minutes --chdir="$REPO" "$REPO/scripts/klone/autopilot.sbatch" "$@"
"$REPO/scripts/klone/autopilot_tick.sh" "$@"
