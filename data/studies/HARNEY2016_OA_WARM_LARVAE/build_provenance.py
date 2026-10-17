"""Record what PXD002316 holds, and why it is registered but not harmonized.

HARNEY2016_OA_WARM_LARVAE is a Coomassie 2-DE study. Spot volumes were
quantified from gel images in Progenesis SameSpots; neither the images nor
the spot volumes were deposited. PRIDE PXD002316 holds only the LC-MS/MS
identification runs of excised spots. The only quantitative result is Table 2
of the paper (fold change and P for 12 significance-selected spots), and the
spots' UniProt accessions have mostly been deleted from UniProtKB since.

    python -I data/studies/HARNEY2016_OA_WARM_LARVAE/build_provenance.py --download STAGING_DIR

fetches the PRIDE project and file records, the authors' accepted manuscript
from Ifremer Archimer (checksummed, not committed) and the current UniProtKB
record for each spot accession, and writes

  data/studies/HARNEY2016_OA_WARM_LARVAE/pride_provenance.json

The spot -> accession list is transcribed from the paper (Table 1 and Results,
which give the C. gigas accessions behind spots 7 and 12; Table 1 itself lists
their nearest characterised homologs in other species). It is the only
transcribed content, and it is short enough to check by eye against the
manuscript; no fold changes or P values are transcribed.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
import urllib.error
import urllib.request
from collections import Counter
from datetime import date
from pathlib import Path

ACCESSION = "PXD002316"
STUDY_ID = "HARNEY2016_OA_WARM_LARVAE"
PRIDE_API = f"https://www.ebi.ac.uk/pride/ws/archive/v3/projects/{ACCESSION}"
ARCHIMER_PDF = "https://archimer.ifremer.fr/doc/00299/40985/40076.pdf"
UNIPROT_API = "https://rest.uniprot.org/uniprotkb/{}?format=json"

REPO_ROOT = Path(__file__).resolve().parents[3]
PROVENANCE = REPO_ROOT / "data" / "studies" / STUDY_ID / "pride_provenance.json"

# Paper spot number -> (annotation, C. gigas UniProt accession as published).
SPOTS = {
    1: ("Matrilin-2", "K1QDH1"),
    2: ("78 kDa glucose-regulated protein (GRP78)", "K1QIR8"),
    3: ("Protein disulfide-isomerase 2 (PDI 2)", "K1Q6X5"),
    4: ("ATP synthase subunit b", "K1RWW5"),
    5: ("Severin", "K1PE57"),
    6: ("Severin", "K1PE57"),
    7: ("Uncharacterised; nearest characterised: shell matrix protein G9MBW9 (Pinctada maxima)", "K1QWT8"),
    8: ("Actin-2", "K1QWP8"),
    9: ("Glyceraldehyde-3-phosphate dehydrogenase (GAPDH)", "K1Q350"),
    10: ("Superoxide dismutase (SOD)", "K1QGQ5"),
    11: ("Vitellogenin-6", "K1QNA2"),
    12: ("Uncharacterised; nearest characterised: USP-like protein isoform 2 H9LJ92 (Crassostrea ariakensis)", "K1QVF8"),
}


def fetch(url: str) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent": "AREE-curation"})
    with urllib.request.urlopen(req, timeout=300) as resp:
        return resp.read()


def uniprot_status(acc: str) -> dict:
    try:
        rec = json.loads(fetch(UNIPROT_API.format(acc)))
    except urllib.error.HTTPError as e:
        return {"status": f"http_{e.code}", "ncbi_gene_ids": []}
    if rec.get("entryType") == "Inactive":
        return {"status": rec["inactiveReason"]["inactiveReasonType"].lower(),
                "reason": rec["inactiveReason"].get("deletedReason"), "ncbi_gene_ids": []}
    return {"status": "active", "entry_type": rec["entryType"],
            "ncbi_gene_ids": [x["id"] for x in rec.get("uniProtKBCrossReferences", []) if x["database"] == "GeneID"]}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--download", type=Path, required=True, help="staging directory for fetched files")
    args = ap.parse_args()
    stage = args.download
    stage.mkdir(parents=True, exist_ok=True)

    project_raw = fetch(PRIDE_API)
    files_raw = fetch(f"{PRIDE_API}/files/all")
    pdf_raw = fetch(ARCHIMER_PDF)
    for name, data in (("pride_project.json", project_raw), ("pride_files.json", files_raw),
                       ("archimer_40076.pdf", pdf_raw)):
        (stage / name).write_bytes(data)
    if not pdf_raw.startswith(b"%PDF"):
        sys.exit("Archimer did not return a PDF")
    project, files = json.loads(project_raw), json.loads(files_raw)

    raw_runs = sorted(f["fileName"].removesuffix(".RAW") for f in files if f["fileCategory"]["value"] == "RAW")
    accessions = sorted({acc for _, acc in SPOTS.values()})
    status = {acc: uniprot_status(acc) for acc in accessions}

    provenance = {
        "study_id": STUDY_ID,
        "accession": ACCESSION,
        "repository": "PRIDE",
        "license": project["license"],
        "submission_type": project["submissionType"],
        "date_retrieved": date.today().isoformat(),
        "generated_by": f"data/studies/{STUDY_ID}/build_provenance.py",
        "source_urls": [PRIDE_API, f"{PRIDE_API}/files/all", ARCHIMER_PDF, UNIPROT_API.format("<accession>")],
        "pride_reference": project["references"][0]["referenceLine"],
        "pride_instrument": [i["name"] for i in project["instruments"]],
        "pride_description_spot_count": "22 differentially expressed spots, 18 proteins (pre-revision; the paper reports 12 spots at P < 0.02)",
        "pride_files_by_category": dict(sorted(Counter(f["fileCategory"]["value"] for f in files).items())),
        "pride_spot_runs": raw_runs,
        "n_pride_spot_runs": len(raw_runs),
        "quantitative_data_deposited": False,
        "input_sha256": {"archimer_40076.pdf (accepted manuscript)": hashlib.sha256(pdf_raw).hexdigest()},
        "paper_spots": {str(n): {"annotation": a, "uniprot_accession": acc, **status[acc]} for n, (a, acc) in SPOTS.items()},
        "uniprot_summary": dict(Counter(s["status"] for s in status.values())),
    }
    PROVENANCE.write_text(json.dumps(provenance, indent=2, ensure_ascii=False) + "\n")
    print(f"wrote {PROVENANCE.relative_to(REPO_ROOT)}: {len(raw_runs)} PRIDE spot runs, "
          f"UniProt {provenance['uniprot_summary']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
