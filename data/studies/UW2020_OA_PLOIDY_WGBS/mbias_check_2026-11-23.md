# UW2020_OA_PLOIDY_WGBS: alignment-mode and M-bias check (2026-11-23)

Run on Chlo before the full run, to set `bismark_align_args` and `extractor_args`
in `config/UW2020_OA_PLOIDY_WGBS.config`.

**Input.** First 2,000,000 read pairs of SRR13207057 (2N_LOW_03), streamed from
ENA. Raw reads are 135 bp (NovaSeq 6000). Trim Galore 0.6.10 `--paired`
(adapter AGATCGGAAGAGC) removed only 0.4% of bases. Bismark 0.24.2, bowtie2 2.5.5,
genome GCF_963853765.1.

## Alignment mode

| Mode | score_min | Pairs | Unique best hit | Mapping efficiency |
|---|---|---|---|---|
| directional | L,0,-0.2 (default) | 1,999,732 | 385,894 | 19.3% |
| directional | L,0,-0.6 | 1,999,732 | 597,035 | 29.9% |
| **non-directional** | L,0,-0.6 | 500,000 | 293,084 | **58.6%** |

Unique alignments by strand in non-directional mode: OT 73,401, CTOT 73,172,
CTOB 73,022, OB 73,489, about 25% each. In directional mode, 62% of pairs had no
alignment. The base composition agrees: both R1 and R2 are depleted of C (~9%) and
of G (~11%) at every position, as expected for a mixture of all four
bisulfite strands.

**Decision.** `--non_directional --score_min L,0,-0.6`. In directional mode
half of every library could not be placed. The ENA record has no library
construction protocol, so the library is identified as non-directional from
these data alone.

## M-bias (CpG context)

Methylation is flat at 10-12% from position 1 to 135 in both R1 and R2, in
directional and non-directional BAMs. There is no 5' bias, unlike the R2
positions 1-20 in CIBNOR2021_HEAT_RRSS. In the directional BAM the final
position (135) read 0% in R1 and 100% in R2 on very few calls, so one 3' base
is dropped per read.

**Decision.** `extractor_args = '--ignore_3prime 1 --ignore_3prime_r2 1'`.

Overall methylation in the directional subset was CpG 11.3%, CHG 0.9%, CHH 1.6%.
CHH methylation is not fully converted and is a ceiling on non-conversion, as
*C. gigas* has negligible non-CpG methylation.

## Cost

Non-directional alignment ran at ~5.4 M pairs/h with `--multicore 3` (12 bowtie2
processes). For the 1,240 M pairs across 24 libraries, alignment is about 9-10 days on Chlo.

Scripts: `~/aree-runs/UW2020_OA_PLOIDY_WGBS/mbias_subset_2026-11-23.sh` (not
committed). Scratch outputs in `/Volumes/OWC-Express/aree/tmp/mbias_UW2020`.
