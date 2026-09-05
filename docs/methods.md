# Methods

*This section is written in a manuscript-appropriate style for reuse in a
paper's Methods, describing the AREE framework as implemented in this
repository. Numeric examples reference the demonstration dataset described
below; production use requires substituting real registered studies.*

## Data model

The Aquaculture Resilience Evidence Engine (AREE) organizes biomarker
evidence into three related tiers. A **study** record captures
registration-time metadata describing what was done — species, strain or
population, reference genome assembly and annotation version, assay type,
analysis mode (raw reanalysis versus processed-results harmonization), and
data availability — without asserting what was found. Because a single study
frequently reports multiple phenotypes or stressor exposures (for example, an
acute and a chronic thermal challenge from the same experiment), phenotype,
stressor, tissue, life stage, and sample-size fields are attached to one or
more **comparisons** nested within the study record, rather than to the study
as a whole. Each comparison is validated against a JSON Schema
(`schemas/study.schema.json`) and against controlled vocabularies for
phenotype, stressor, tissue, life stage, and assay type, so that terms such
as "thermal tolerance" or "ocean acidification" carry a fixed, shared
definition across independently curated studies. Each phenotype term is
additionally tagged with a `resilience_relevance` classification —
resilience, stress_response, disease, or exposure_only — that is propagated
into every downstream comparison and evidence record, so that a measured
change under stressor exposure is never structurally conflated with a
validated measure of organismal resilience.

Assay-specific results (differential expression tables, differentially
methylated region tables, protein abundance tables, metabolite feature
tables) are converted into a single assay-agnostic **evidence record**
schema (`schemas/evidence.schema.json`): one row per molecular feature per
comparison, carrying the original and standardized feature identifiers,
effect size and effect-size type, standard error or p-value, sample size,
the controlled-vocabulary phenotype/stressor/tissue/life-stage context, a
mapping-confidence label, quality flags, and full provenance (source file
and checksum, workflow version, tool versions, reference genome/annotation,
generation date, and generating agent). This harmonization step is identical
in downstream schema regardless of whether the source study underwent raw
reanalysis or contributed only processed summary results, differing only in
which quality flags are attached (processed-results-only studies are flagged
as such because upstream QC cannot be independently verified).

## Identifier harmonization

Molecular feature identifiers are standardized through a documented
precedence hierarchy: NCBI Gene ID, then Ensembl gene ID, then UniProt
accession, then a reference-annotation locus identifier, then gene symbol,
then orthogroup. The original identifier as reported by the source study is
always retained (`feature_id_original`); the harmonized identifier
(`feature_id_standardized`) is resolved against a crosswalk table and
assigned one of six mapping-confidence levels — exact, one-to-one ortholog,
one-to-many ortholog, many-to-one ortholog, inferred, or unresolved — ranked
from most to least trustworthy. Ambiguous legacy gene symbols that do not
resolve cleanly against the primary crosswalk are checked against a curated
exception table before falling back to unresolved. Records with unresolved
mapping remain in the evidence table for transparency but are excluded from
identifier-level pooling, since no stable shared identity exists to group on.

## Meta-analysis

Evidence records are grouped for pooling by the combination of standardized
feature identifier, phenotype, and feature type; pooling is never performed
across feature types, since effect-size scales (for example, log2 fold
change for expression versus percent methylation difference) are not
directly comparable. Within each group, study-level effect sizes are
combined using a random-effects model with the DerSimonian–Laird estimator
of between-study variance (tau²), implemented directly in Python rather than
through a compiled R meta-analysis package, to keep the estimator fully
auditable and dependency-light. For a group of *k* independent effect
estimates *y_i* with standard errors *se_i*, fixed-effect weights
*w_i = 1/se_i²* are used to compute Cochran's Q and a method-of-moments
estimate of tau²; random-effects weights *w*_i = 1/(se_i² + tau²)* then
produce the pooled effect, its standard error, a 95% confidence interval, and
a two-sided z-test p-value. Higgins' I² (the proportion of total variation
attributable to between-study heterogeneity rather than sampling error) is
reported alongside every pooled estimate. When a contributing record lacks a
directly reported standard error, one is approximated from the reported
effect size and p-value under a normal approximation
(*se = |effect| / Φ⁻¹(1 − p/2)*), a documented fallback rather than a silent
assumption. A single-study group falls back to that study's own estimate
(tau² = 0, I² = 0) rather than a degenerate random-effects computation.
Because a genome-wide reanalysis contributes one pooled test per gene, pooled
p-values are adjusted for multiple testing with the Benjamini–Hochberg
procedure within each test family, defined as all features pooled for one
phenotype, feature type, data origin, and species; both the raw and adjusted
values and the family size are reported. Direction consistency is reported
separately as the fraction of contributing
records whose effect-size sign agrees with the majority sign across the
group; low direction consistency combined with high I² is treated as
evidence of genuine cross-study disagreement rather than as noise to be
averaged away, and is preserved rather than suppressed in downstream outputs.

## Candidate scoring and tiering

Each pooled feature/phenotype/feature-type group constitutes a candidate.
Candidates are assigned a transparent score computed as a weighted sum of
ten independently inspectable 0–1 components — number of independent
studies, total biological sample size, effect magnitude, adjusted
significance, direction consistency, phenotype relevance (derived from the
resilience-relevance classification of the phenotype), breadth across
tissues and life stages, assay-type diversity, identifier mapping
confidence, and data-quality flag burden — with weights fixed and documented
in source (`src/prioritize/scoring.py`), less a penalty proportional to I².
The scoring function is pure, so identical inputs deterministically produce
an identical score.

Critically, the numeric score does not by itself determine tier membership.
Three tiers are defined with explicit, score-independent gating conditions.
A candidate is assigned to the high-priority cross-study tier only if it is
supported by at least two independent studies, has a phenotype relevance
above the exposure-only floor, has a direction consistency of at least 0.7,
and clears a minimum data-quality threshold. A candidate lacking sufficient
independent-study support but showing convergent evidence across at least
two distinct molecular assay types (for example, transcriptomic and
methylation evidence for the same standardized feature) is assigned to a
multi-omics convergence tier, with the specific cross-assay linkage shown
explicitly. All remaining candidates — including those with a strong effect
in a single study — are assigned to an emerging tier and explicitly labeled
as requiring independent replication. This gating structure is the
mechanism by which the framework prevents a single study's statistical
significance from being read as validation: no combination of effect size,
significance, or sample size in one study can promote a candidate past the
replication and consistency gates.

## Demonstration dataset

The framework is exercised against six synthetic study records representing
*Crassostrea gigas* resilience contexts (acute and chronic heat challenge,
ocean acidification, pathogen challenge, salinity stress, nutritional
limitation affecting larval viability, and hypoxia affecting growth),
spanning RNA-seq, DNA methylation, proteomics, and metabolomics assays, and
including both raw-reanalysis and processed-results-only registration modes,
one study with deliberately imperfect legacy gene-symbol identifiers, and one
pair of studies reporting opposite-direction effects for the same gene under
the same phenotype. All demonstration data are explicitly labeled as
simulated and are not derived from real sequencing reads or a verified
reference genome accession; they establish that the pipeline runs end to end
and produces the intended range of mapping-confidence and heterogeneity
outcomes, not that any specific candidate identified in the demonstration
run is a real biological finding.

## Related documentation

- [design.md](design.md)
- [identifier_mapping.md](identifier_mapping.md)
- [interpreting_meta_analysis.md](interpreting_meta_analysis.md)
- [interpreting_candidate_scores.md](interpreting_candidate_scores.md)
