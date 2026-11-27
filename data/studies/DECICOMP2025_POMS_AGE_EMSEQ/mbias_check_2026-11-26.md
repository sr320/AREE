# DECICOMP2025_POMS_AGE_EMSEQ: mapping and M-bias check (2026-11-26)

Run on Klone before the full run, to set `methylation.extractor_args` in
`run/klone.config`. Script: `run/mbias_check.sbatch` (jobs 43570502, 43571349,
43571417; the first two stopped for want of samtools, see below).

**Input.** First 1,000,000 read pairs of ERR15996532 (F14R_C3_T0_S4), streamed
from the gannet mirror of ENA's FASTQ. 2 x 150 bp, NovaSeq 6000, EM-seq.
Trim Galore 0.6.10 `--paired` (Illumina adapter, auto-detected) kept 999,814
pairs. Bismark 0.24.2, directional, `--score_min L,0,-0.6`, genome
GCF_963853765.1. One library only: the other 59 were not checked.

## Mapping

| Pairs analysed | Unique best hit | Mapping efficiency | CpG | CHG | CHH |
|---|---|---|---|---|---|
| 999,814 | 485,469 | 48.6% | 12.8% | 0.3% | 0.2% |

CHH and CHG methylation of 0.2-0.3% indicate near-complete enzymatic
conversion. Non-directional mode was not tested.

## M-bias (CpG context)

About 9,000 CpG calls per position in R1 and 7,500 in R2, so one position's
value has a standard error of about 0.35 points.

- Plateau (positions 30-120): R1 12.95%, R2 12.75%.
- R1 positions 1-2: 15.5%, 14.3% (high). From position 3, on the plateau.
- R1 position 150: 11.7% (about 1 point low; positions 141-149 are 12.3-13.1%).
- R2 positions 1-3: 0.4%, 4.5%, 10.4% (end-repair fill-in). Positions 4-12:
  11.2-12.2%, about 1 point low. From position 13, on the plateau.
- R2 position 150: 100% on 892 calls (artifact).
- CHH: flat at 0.2% in both reads, except R2 position 150 (100% on 74 calls).

**Decision.** `extractor_args = '--ignore 2 --ignore_r2 12 --ignore_3prime 1 --ignore_3prime_r2 1'`.
VALDIVIESO2025_POMS_RS_EMSEQ (same laboratory and kit) used
`--ignore 2 --ignore_r2 20 --ignore_3prime 1 --ignore_3prime_r2 1`; this
library's R2 reaches the plateau sooner.

## Tooling note

The `quay.io/biocontainers/bismark:0.24.2--hdfd78af_0` image has no samtools.
Bismark then writes SAM and cannot merge `--multicore` outputs, and
`bismark_methylation_extractor` refuses to start. The check aligns with one
Bismark instance in the container and runs the extractor from the pinned
conda environment that `modules/methylation/bismark_chunked.nf` declares
(bismark 0.24.2, samtools 1.18).

Outputs on Klone group storage:
`/gscratch/srlab/sr320/aree/mbias/DECICOMP2025_POMS_AGE_EMSEQ/ERR15996532/`.
