"""Write a machine-readable provenance manifest for a workflow/analysis run."""
from __future__ import annotations

import json
from pathlib import Path


def write_manifest(output_path, **fields) -> Path:
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as fh:
        json.dump(fields, fh, indent=2, sort_keys=True, default=str)
    return output_path
