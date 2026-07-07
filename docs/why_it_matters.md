# Why this resource matters

Aquaculture breeding and management decisions increasingly rely on molecular
evidence — which genes, proteins, or epigenetic marks track with an animal's
ability to survive heat events, disease outbreaks, acidification, or low
salinity. That evidence already exists across dozens of published studies on
Pacific oyster (*Crassostrea gigas*) and related shellfish. The problem is
not a lack of data. It is that the data are scattered across different
sequencing platforms, genome versions, phenotype definitions, and identifier
systems, so no one can easily ask "which candidate genes show up
consistently across independent thermal-tolerance studies?" without
re-doing the comparison from scratch every time.

The proposal this resource implements describes the need directly:

> "Develop standardized open-access, user-friendly, reproducible
> bioinformatics pipelines for resilience biomarker discovery through
> systematic reanalysis, data integration, and meta-analysis."

The Aquaculture Resilience Evidence Engine (AREE) is a working system built
to make that operational, not a literature summary or a static list of
papers.

## What this means for a breeder or researcher in practice

- **You can ask a specific question and get a traceable answer.** Instead of
  reading six papers and mentally reconciling their different effect-size
  scales, sample sizes, and phenotype definitions, you can query AREE's
  harmonized evidence table for a phenotype (say, thermal tolerance) and see
  every study's contribution on a common scale, with the original study
  context preserved alongside it.
- **You can see convergence, and you can see disagreement.** When two
  independent studies point at the same gene in the same direction, that is
  worth attention. When they point in opposite directions, AREE surfaces
  that conflict explicitly (see
  [interpreting_meta_analysis.md](interpreting_meta_analysis.md)) instead of
  quietly averaging it into a number that looks reassuring but hides real
  biological or methodological disagreement.
- **You can tell resilience evidence apart from stress-response noise.**
  A gene that merely responds to heat is not the same as a gene that
  predicts which animals survive heat. AREE enforces this distinction
  structurally (see
  [resilience_vs_exposure.md](resilience_vs_exposure.md)) so a candidate list
  isn't quietly inflated with genes that just react to being stressed.
- **You get a transparent, auditable priority list, not a black box.** Every
  candidate's score is broken into named components you can inspect —
  number of supporting studies, direction agreement, effect size, phenotype
  relevance, and more — and no candidate is promoted to a high-confidence
  tier without meeting explicit replication and consistency requirements
  (see [interpreting_candidate_scores.md](interpreting_candidate_scores.md)).
  This matters for breeding programs: a candidate marker that will inform
  selection decisions needs a paper trail a skeptical geneticist can follow
  back to source data, not a single p-value from one study.
- **The system tells you what it doesn't know.** Missing metadata, uncertain
  identifier mappings, low replication, and processed-only data availability
  are recorded and surfaced, not hidden (see
  [governance_and_provenance.md](governance_and_provenance.md)). A candidate
  list you can trust is one that is honest about its own weak points.

## Why this is worth building now

Reanalyzing public aquaculture omics data one paper at a time does not scale,
and every ad hoc reanalysis effort has to re-solve the same problems: which
identifier system to standardize on, how to compare effect sizes across
assay types, how to avoid conflating "changed under stress" with "predicts
resilience." AREE solves those problems once, as reusable infrastructure, so
that adding the next public *C. gigas* dataset — or eventually a second
shellfish species — is a matter of filling in a registration template and
running a handful of commands, not re-inventing a harmonization pipeline.

## Related documentation

- [architecture.md](architecture.md) — how the system is built
- [roadmap.md](roadmap.md) — what is functional today vs. planned
- [design.md](design.md) — the full data model and design rationale
