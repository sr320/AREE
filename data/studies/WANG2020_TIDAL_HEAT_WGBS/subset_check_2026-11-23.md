# WANG2020_TIDAL_HEAT_WGBS: alignment-mode and M-bias check (2026-11-23)

Run on Chlo before the full run, to set `bismark_align_args` and
`extractor_args` in `config/WANG2020_TIDAL_HEAT_WGBS.config`.

**Input.** First 1,000,000 read pairs of SRR10065364, cut from the start of the
full ENA download (a direct ENA stream dropped mid-file). 2 x 150 bp,
HiSeq X Ten. Trim Galore 0.6.10 `--paired` (adapter AGATCGGAAGAGC) wrote
99.2% (R1) and 98.2% (R2) of bases. Bismark 0.24.2, bowtie2 2.5.5, genome
GCF_963853765.1.

## Alignment mode

| Mode | score_min | Unique best hit | Mapping efficiency | Unique pairs on CTOT / CTOB |
|---|---|---|---|---|
| directional | L,0,-0.2 (default) | 325,005 | 32.5% | (not searched) |
| directional | L,0,-0.6 | 588,460 | **58.9%** | (not searched) |
| non-directional | L,0,-0.6 | 588,488 | 58.9% | 120 / 89 |

These are **directional** libraries: 209 of 588,488 unique pairs (0.04%) came from
the complementary strands. The authors reported ~30% mapping at default
scoring, which matches the 32.5% here.

**Decision.** `--score_min L,0,-0.6`, directional (unchanged from the
registered config).

## M-bias (CpG context, directional L,0,-0.6 BAM)

- The plateau in both reads is about 12-13% CpG methylation.
- R1 positions 1-2: 14.9%, 17.5% (high).
- R2 positions 1-2: 0.2%, 1.8% (end-repair fill-in). R2 positions 3-15:
  9.6-12.3%, about 1.5 points low. From position 16, R2 is on the plateau.
- R2 position 150: 100% on 504 calls (artifact). R1 has no 3' bias.

**Decision.** `extractor_args = '--ignore 2 --ignore_r2 15 --ignore_3prime_r2 1'`.

## Throughput

The 1 M-pair directional alignment took about 20 min with `--multicore 3`, but
that includes loading six bowtie2 indexes from blue-block (USB) and it ran
alongside the 12-connection download. CIBNOR2021's full runs give the better
estimate, ~11 M pairs/h, or about 3.5 days of alignment for 925 M pairs.

Scripts: `~/aree-runs/WANG2020_TIDAL_HEAT_WGBS/subset_check_2026-11-23.sh` (not
committed). Scratch outputs in `/Volumes/OWC-Express/aree/tmp/subset_WANG2020`.
