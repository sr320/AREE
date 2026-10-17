"""Record what PXD000905 holds, and why it is registered but not reanalyzed.

CORPOREAU2014_OSHV_LOAD is a Coomassie 2-DE study (one gel per oyster, 5
high-load vs 5 low-load OsHV-1 injected oysters). Spot volumes were quantified
in Progenesis SameSpots and never deposited; PRIDE PXD000905 holds only the
nanoLC-MS/MS identification runs of the 25 excised spots (Bruker .yep, named
<spot number><source gel>, e.g. 337CA61) with their Mascot results. The paper
(Corporeau et al. 2014, J Proteomics 109:176-187; accepted manuscript on
Ifremer Archimer) reports an average fold change for those 25 spots only.

    python -I data/studies/CORPOREAU2014_OSHV_LOAD/build_provenance.py --download STAGING_DIR

fetches the PRIDE records and the Archimer manuscript (checksummed, not
committed: Elsevier copyright) and writes
data/studies/CORPOREAU2014_OSHV_LOAD/pride_provenance.json.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
import urllib.request
from collections import Counter
from datetime import date
from pathlib import Path

ACCESSION = "PXD000905"
STUDY_ID = "CORPOREAU2014_OSHV_LOAD"
PRIDE_API = f"https://www.ebi.ac.uk/pride/ws/archive/v3/projects/{ACCESSION}"
ARCHIMER_PDF = "https://archimer.ifremer.fr/doc/00199/31035/29439.pdf"

REPO_ROOT = Path(__file__).resolve().parents[3]
OUT = REPO_ROOT / "data" / "studies" / STUDY_ID


def fetch(url: str) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent": "AREE-curation"})
    with urllib.request.urlopen(req, timeout=300) as resp:
        return resp.read()


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--download", type=Path, required=True, help="staging directory for fetched files")
    args = ap.parse_args()
    args.download.mkdir(parents=True, exist_ok=True)

    project = json.loads(fetch(PRIDE_API))
    files = json.loads(fetch(f"{PRIDE_API}/files/all"))
    pdf = fetch(ARCHIMER_PDF)
    if not pdf.startswith(b"%PDF"):
        sys.exit("Archimer did not return a PDF")
    (args.download / "archimer_29439.pdf").write_bytes(pdf)

    spots = []
    for f in sorted(x["fileName"] for x in files if x["fileCategory"]["value"] == "RAW"):
        m = re.fullmatch(r"Analysis_Corporeau_(\d+)([CI]A)(\d+)\.yep", f)
        if not m:
            sys.exit(f"unrecognised spot run {f}")
        spots.append({"raw_file": f, "spot": int(m.group(1)),
                      "source_gel_group": {"CA": "low_load_control", "IA": "high_load_challenged"}[m.group(2)],
                      "source_oyster": int(m.group(3))})
    if len({s["spot"] for s in spots}) != len(spots):
        sys.exit("a spot number appears in more than one run")

    provenance = {
        "study_id": STUDY_ID,
        "accession": ACCESSION,
        "repository": "PRIDE",
        "license": project["license"],
        "submission_type": project["submissionType"],
        "date_retrieved": date.today().isoformat(),
        "generated_by": f"data/studies/{STUDY_ID}/build_provenance.py",
        "source_urls": [PRIDE_API, f"{PRIDE_API}/files/all", ARCHIMER_PDF],
        "pride_reference": project["references"][0]["referenceLine"],
        "pride_instrument": [i["name"] for i in project["instruments"]],
        "pride_files_by_category": dict(sorted(Counter(f["fileCategory"]["value"] for f in files).items())),
        "n_spot_runs": len(spots),
        "spot_runs": spots,
        "source_oysters_of_excised_spots": sorted({(s["source_gel_group"], s["source_oyster"]) for s in spots}),
        "quantitative_data_deposited": False,
        "input_sha256": {"archimer_29439.pdf (accepted manuscript)": hashlib.sha256(pdf).hexdigest()},
    }
    (OUT / "pride_provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    print(f"wrote pride_provenance.json: {len(spots)} spot runs from oysters "
          f"{provenance['source_oysters_of_excised_spots']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
