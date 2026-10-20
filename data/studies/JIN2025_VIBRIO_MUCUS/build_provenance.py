"""Record what PXD064651 holds, and why it is registered but cannot be reanalyzed.

PRIDE PXD064651 (CC0, PARTIAL, no publication): mucus of 3 control (C1-C3)
and 3 Vibrio ZJ-51-infected (T1-T3) Pacific oysters, IBT 16-plex isobaric
labelling, 20 high-pH fractions, Q Exactive HF-X; searched with Mascot and
quantified with IQuant. The deposit is three tarred .raw archives, a tarred
MGF and a tarred Mascot .dat. It names the six samples but says nothing about
which reporter channel carries which sample, and no quantification table is
deposited.

    python -I data/studies/JIN2025_VIBRIO_MUCUS/build_provenance.py --download STAGING_DIR

fetches the PRIDE records and the first 3 MB of dat.tar.gz, reads the
parameter header of the first Mascot .dat inside it, and writes
data/studies/JIN2025_VIBRIO_MUCUS/pride_provenance.json.
"""
from __future__ import annotations

import argparse
import json
import re
import sys
import urllib.request
import zlib
from collections import Counter
from datetime import date
from pathlib import Path

ACCESSION = "PXD064651"
STUDY_ID = "JIN2025_VIBRIO_MUCUS"
PRIDE_API = f"https://www.ebi.ac.uk/pride/ws/archive/v3/projects/{ACCESSION}"
DAT_TAR = "https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/07/PXD064651/dat.tar.gz"
DAT_KEYS = ("COM", "DB", "MODS", "IT_MODS", "QUANTITATION", "INSTRUMENT", "FILE", "sequences", "queries", "fastafile")

REPO_ROOT = Path(__file__).resolve().parents[3]
OUT = REPO_ROOT / "data" / "studies" / STUDY_ID


def fetch(url: str, byte_range: str | None = None) -> bytes:
    headers = {"User-Agent": "AREE-curation"}
    if byte_range:
        headers["Range"] = f"bytes={byte_range}"
    with urllib.request.urlopen(urllib.request.Request(url, headers=headers), timeout=600) as resp:
        return resp.read()


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--download", type=Path, required=True, help="staging directory for fetched files")
    args = ap.parse_args()
    args.download.mkdir(parents=True, exist_ok=True)

    project = json.loads(fetch(PRIDE_API))
    files = json.loads(fetch(f"{PRIDE_API}/files/all"))

    # First 3 MB of the gzip stream: the tar header of the first member plus
    # the start of that Mascot .dat, which carries the search parameters.
    head = zlib.decompressobj(16 + zlib.MAX_WBITS).decompress(fetch(DAT_TAR, "0-3000000"))
    member = head[:100].split(b"\0", 1)[0].decode()
    text = head[512:].decode("latin-1")
    params = {}
    for key in DAT_KEYS:
        m = re.search(rf"^{key}=(.*)$", text, flags=re.M)
        params[key] = m.group(1).strip() if m else None
    if params["QUANTITATION"] is None:
        sys.exit("could not read the Mascot .dat parameter header")

    provenance = {
        "study_id": STUDY_ID,
        "accession": ACCESSION,
        "repository": "PRIDE",
        "license": project["license"],
        "submission_type": project["submissionType"],
        "submission_date": project["submissionDate"],
        "publication_date": project["publicationDate"],
        "references": project.get("references") or [],
        "submitter_affiliation": [s.get("affiliation") for s in project.get("submitters", [])],
        "pride_instrument": [i["name"] for i in project["instruments"]],
        "pride_quantification_methods": [q["name"] for q in project.get("quantificationMethods", [])],
        "pride_files": {f["fileName"]: {"category": f["fileCategory"]["value"], "bytes": f["fileSizeBytes"]}
                        for f in sorted(files, key=lambda f: f["fileName"])},
        "pride_files_by_category": dict(sorted(Counter(f["fileCategory"]["value"] for f in files).items())),
        "samples_named_in_pride": ["C1", "C2", "C3", "T1", "T2", "T3"],
        "channel_to_sample_key_deposited": False,
        "quantification_table_deposited": False,
        "mascot_dat_first_member": member,
        "mascot_dat_parameters": params,
        "date_retrieved": date.today().isoformat(),
        "generated_by": f"data/studies/{STUDY_ID}/build_provenance.py",
        "source_urls": [PRIDE_API, f"{PRIDE_API}/files/all", f"{DAT_TAR} (bytes 0-3000000)"],
    }
    (OUT / "pride_provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    print(json.dumps({"member": member, **params}, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
