#!/usr/bin/env bash
# One Klone autopilot tick (report-only). Called by scrontab or by
# autopilot.sbatch; never run it on a login node.
#
#   autopilot_tick.sh <miniforge prefix> [conda env, default aree]
set -euo pipefail

MINIFORGE="${1:?usage: autopilot_tick.sh <miniforge prefix> [conda env]}"
ENV_NAME="${2:-aree}"
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
STATE_DIR="${AREE_AUTOPILOT_STATE_DIR:-$HOME/.local/share/aree-autopilot}"
mkdir -p "$STATE_DIR"

# shellcheck disable=SC1091
source "$MINIFORGE/etc/profile.d/conda.sh"
conda activate "$ENV_NAME"

# Two ticks must never overlap (a slow sacct plus a short interval).
exec 9>"$STATE_DIR/tick.lock"
if ! flock -n 9; then
  echo "another tick is running; skipping"
  exit 0
fi

cd "$REPO"
aree autopilot tick
