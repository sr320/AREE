You are building an open, reproducible research resource called the:



Aquaculture Resilience Evidence Engine

Short name: AREE



Its purpose is to make Objective 1 of the attached proposal operational:



“Develop standardized open-access, user-friendly, reproducible bioinformatics pipelines for resilience biomarker discovery through systematic reanalysis, data integration, and meta-analysis.”



The resource must initially focus on Pacific oyster (Crassostrea gigas), while being deliberately extensible to additional shellfish species and other aquaculture organisms.



Do not build a generic website or a static resource list. Build a functioning, versioned evidence-generation system that converts public omics datasets into harmonized, comparable biomarker evidence.



## Core scientific problem



Public aquaculture datasets are fragmented across studies, assays, genome versions, phenotype definitions, environmental treatments, and identifier systems. The resource must make it possible to:



1. Register and characterize a public dataset.

2. Standardize study metadata, phenotype labels, exposure descriptions, and sample annotations.

3. Reanalyze datasets through transparent, reusable workflows.

4. Convert assay-specific results into a shared evidence schema.

5. Compare effects across studies despite different scales and experimental designs.

6. Identify genes, pathways, proteins, genomic regions, or molecular features repeatedly associated with resilience-related phenotypes.

7. Produce outputs that can populate a future resilience biomarker database.



The system must distinguish clearly between:

- resilience-associated evidence,

- stress-response evidence,

- disease-associated evidence,

- environmental exposure evidence,

- and evidence that is only suggestive because phenotype definitions or study quality are weak.



Never imply that a candidate biomarker is validated merely because it is statistically significant in one study.



## Required deliverable



Create a GitHub-ready repository that includes:



1. A working data model.

2. Dataset-registration templates.

3. Reproducible workflow scaffolds.

4. Standardized result schemas.

5. Cross-study harmonization and meta-analysis code.

6. A candidate biomarker prioritization framework.

7. A lightweight user-facing interface or report generator.

8. Documentation sufficient for another lab to add a public study without direct assistance.



The deliverable should be useful even before all raw public data are downloaded or processed. Use representative synthetic/demo data where needed, but structure the repository so real public datasets can be added with minimal modification.



## Product definition



Build the resource as a modular platform with five connected layers:



### Layer 1: Study Registry and Dataset Intake



Create a machine-readable registry for public datasets.



Each study must have:

- study_id

- DOI, BioProject, GEO, SRA, ENA, ProteomeXchange, or other accession

- citation

- species

- strain, population, family, or breeding line where available

- genome assembly and annotation version

- assay type

- tissue

- life stage

- experimental treatment

- control condition

- stressor class

- exposure intensity, duration, and timing

- phenotype measured

- phenotype direction and units

- resilience classification

- sample size

- biological replication

- sequencing or instrument platform

- data availability status

- raw versus processed data availability

- quality-control status

- analysis status

- provenance links

- limitations and caveats



Use controlled vocabularies where possible.



Implement a phenotype ontology with terms such as:

- survival

- mortality

- thermal tolerance

- growth under stress

- pathogen load

- disease susceptibility

- disease resistance

- reproductive performance

- larval viability

- metabolic resilience

- recovery following stress

- immune responsiveness

- acidification tolerance

- salinity tolerance

- hypoxia tolerance



Implement an environmental-stressor ontology with:

- temperature

- ocean acidification / pH

- salinity

- hypoxia

- pathogen challenge

- harmful algal bloom exposure

- nutritional limitation

- freshwater exposure

- pollutant exposure

- handling or transport stress

- multi-stressor exposure



Create templates in:

- YAML for individual study registration

- CSV for batch study entry

- JSON schema for validation



### Layer 2: Standardized Reanalysis Workflows



Implement reproducible workflow scaffolds using Nextflow.



Use Docker or Apptainer containers, environment lock files, and explicit software versioning.



Build separate modular workflows or workflow templates for:



A. RNA-seq

- FASTQ quality control

- adapter trimming

- alignment or pseudoalignment

- gene-level quantification

- sample QC

- differential expression

- effect-size calculation

- enrichment-ready output

- standardized QC and analysis report



B. DNA methylation / WGBS / EM-seq

- read QC

- alignment

- methylation calling

- coverage filtering

- DML / DMR analysis

- annotation to gene, promoter, exon, intron, and intergenic contexts

- standardized genomic-region result table



C. Proteomics

- input harmonization for protein abundance tables

- normalization

- missingness report

- differential abundance

- protein-to-gene identifier translation

- standardized protein-level evidence output



D. Metabolomics

- metabolite feature table intake

- annotation confidence tracking

- normalization and QC

- differential abundance

- pathway or metabolite-class mapping

- standardized feature-level evidence output



Do not assume all studies have raw data available. Support two modes:

- raw-data reanalysis mode

- processed-results harmonization mode



Every workflow must output:

- a machine-readable manifest

- parameters used

- software versions

- input checksums

- QC metrics

- warnings or failure states

- standardized result tables

- a human-readable HTML or Quarto report



### Layer 3: Cross-study Harmonization



Create a shared evidence table that makes results comparable across assay types and studies.



At minimum, each evidence record must include:

- evidence_id

- study_id

- sample comparison

- feature_id_original

- feature_id_standardized

- feature_type

- species

- genome/annotation version

- ortholog/reference gene identifier where applicable

- molecular direction

- effect size

- effect-size type

- standard error or confidence interval where available

- p-value

- adjusted p-value

- sample size

- tissue

- life stage

- stressor

- phenotype

- phenotype direction

- analysis method

- quality flags

- source file

- workflow version

- date generated



Create identifier-harmonization functions that can support:

- NCBI Gene IDs

- Ensembl IDs

- UniProt accessions

- gene symbols

- transcript IDs

- protein IDs

- locus IDs from oyster reference annotations

- orthogroups



Use a clearly documented identifier hierarchy. Preserve original identifiers rather than overwriting them.



Create a mapping-confidence field:

- exact

- one-to-one ortholog

- one-to-many ortholog

- many-to-one ortholog

- inferred

- unresolved



### Layer 4: Meta-analysis and Candidate Prioritization



Develop a reproducible evidence-synthesis framework.



The framework should be able to:

- pool effect sizes across comparable studies;

- use random-effects meta-analysis when appropriate;

- report heterogeneity statistics;

- retain study-level contributions;

- separate evidence by tissue, life stage, stressor, phenotype, and assay type;

- identify candidates with repeated directional support;

- flag inconsistent or context-dependent candidates;

- identify candidates supported across multiple omics layers;

- prioritize candidates for future validation.



Create at least three candidate-ranking categories:



1. High-priority cross-study candidates

   - evidence from at least two independent studies

   - interpretable phenotype relationship

   - effect direction reasonably consistent

   - acceptable metadata and quality scores



2. Multi-omics convergence candidates

   - support from at least two molecular layers, such as transcriptomics plus methylation or proteomics

   - include a transparent explanation of how molecular layers were linked



3. Emerging candidates

   - biologically plausible or statistically strong in one study

   - clearly labeled as requiring replication



Build a transparent candidate score, not a black-box score.



The score should incorporate:

- number of independent studies

- total biological sample size

- effect magnitude

- adjusted significance

- direction consistency

- phenotype relevance

- evidence across tissues or life stages

- assay diversity

- identifier-mapping confidence

- data quality

- heterogeneity

- known limitations



Every candidate must have an “evidence card” that shows:

- candidate identifier(s)

- species and ortholog context

- evidence summary

- studies supporting it

- forest plot or effect summary where feasible

- assay types represented

- phenotype/stressor contexts

- direction of association

- limitations

- recommended next validation step



### Layer 5: User-facing Outputs



Build a lightweight but real interface.



Choose the most practical combination:

- Quarto-generated static site plus searchable data tables, or

- Streamlit interface plus static documentation site.



The interface must allow a user to:

- browse registered studies;

- filter by stressor, phenotype, assay, tissue, life stage, and species;

- inspect study-quality and metadata completeness;

- view harmonized molecular evidence;

- search for a gene, protein, orthogroup, or pathway;

- view candidate biomarker evidence cards;

- download filtered evidence as CSV or TSV;

- download workflow manifests and reproducibility metadata;

- view pipeline status for each study.



Do not make user login or cloud deployment a dependency for the first release.



## Demonstration dataset



Create a realistic demo dataset that demonstrates the entire system.



Use at least six simulated or openly available C. gigas study records covering different resilience contexts, for example:

- thermal tolerance / heat survival

- ocean acidification response

- pathogen or disease challenge

- salinity stress

- larval performance

- growth or survival under environmental stress



The demo must include at least:

- two RNA-seq studies,

- one methylation study,

- one proteomics or metabolomics study,

- one study with only processed results,

- one example where identifier mapping is imperfect,

- one example with conflicting direction across studies.



Use synthetic numerical data where required, but make the study metadata and file structures realistic. Clearly label all simulated data as simulated.



## Repository structure



Use a clear repository structure similar to:



AREE/

├── README.md

├── LICENSE

├── CITATION.cff

├── CONTRIBUTING.md

├── CODE_OF_CONDUCT.md

├── docs/

├── schemas/

├── registry/

│   ├── studies/

│   ├── controlled_vocabularies/

│   └── study_registry.csv

├── workflows/

│   ├── rnaseq/

│   ├── methylation/

│   ├── proteomics/

│   └── metabolomics/

├── modules/

├── containers/

├── config/

├── data/

│   ├── demo/

│   ├── reference/

│   └── mappings/

├── src/

│   ├── intake/

│   ├── harmonize/

│   ├── meta_analysis/

│   ├── prioritize/

│   ├── reporting/

│   └── validation/

├── app/

├── reports/

├── tests/

├── notebooks/

└── .github/

    └── workflows/



Use Python for core data processing and validation. R is appropriate for DESeq2/edgeR/limma-style workflows, meta-analysis, and publication-quality plots. Prefer Quarto for integrated reports. Use Make, Just, or a simple task runner for common commands.



## Required command-line workflows



Ensure the README documents commands equivalent to:



1. Validate a study registration file:

   aree validate-study registry/studies/STUDY_ID.yaml



2. Add a study to the registry:

   aree register-study registry/studies/STUDY_ID.yaml



3. Harmonize a processed study result table:

   aree harmonize --study STUDY_ID --input path/to/results.tsv



4. Run a demo meta-analysis:

   aree meta-analyze --phenotype thermal_tolerance --feature-type gene



5. Generate biomarker evidence cards:

   aree build-evidence-cards --phenotype survival



6. Build a static report or launch the interface:

   quarto render docs/

   streamlit run app/main.py



These commands may initially invoke scaffolded workflows and demo data, but they must run successfully on a clean machine after installation.



## Scientific and reproducibility requirements



Every output must preserve provenance.



For each result, retain:

- source accession

- input file name and checksum

- parameter set

- workflow version

- tool versions

- reference genome and annotation version

- date generated

- analyst or automated workflow identity

- any manual curation decisions



Do not hide missing metadata. Missingness must be explicit and reportable.



Do not collapse all “stress” treatments into one category. Preserve the original treatment and separately map it to controlled ontology terms.



Do not overstate causality. The resource identifies associations and evidence convergence, not confirmed mechanistic causation.



Do not remove contradictory findings. Surface them in candidate evidence cards and meta-analysis summaries.



## Validation and quality requirements



Include automated tests for:

- schema validation;

- malformed study metadata;

- duplicate study IDs;

- required provenance fields;

- identifier mapping confidence assignment;

- effect-size calculations;

- meta-analysis on demo inputs;

- candidate-score reproducibility;

- generation of an evidence card;

- successful build of the demo report.



Include GitHub Actions for:

- Python tests;

- R checks where feasible;

- schema validation;

- rendering the demo report;

- linting;

- checking that demo commands complete successfully.



## Documentation requirements



Write documentation for:

- installing the resource;

- adding a study;

- defining a phenotype;

- defining resilience versus exposure-only studies;

- selecting raw reanalysis versus processed-results harmonization;

- mapping identifiers;

- interpreting meta-analysis outputs;

- interpreting candidate biomarker scores;

- adding a new species;

- handling genome-version changes;

- contributing code or curated datasets.



Also write:

1. A concise “Why this resource matters” page for aquaculture researchers and breeders.

2. A technical architecture page.

3. A data-governance and provenance page.

4. A roadmap distinguishing MVP features from later production features.

5. A short manuscript-style Methods section describing the framework.



## Development sequence



Work in phases.



Phase 1: Design

- Inspect the proposal.

- Produce a concise design document.

- Define the data model, controlled vocabularies, provenance model, ranking framework, and repository architecture.

- Identify assumptions and explicitly list them.



Phase 2: Implement MVP

- Create the repository structure.

- Implement schemas and validation.

- Implement registry ingestion.

- Create demo studies.

- Implement harmonization for processed RNA-seq, methylation, proteomics, and metabolomics result tables.

- Implement a working meta-analysis and candidate prioritization pipeline.

- Generate evidence cards and a static Quarto report or minimal Streamlit interface.



Phase 3: Reproducibility infrastructure

- Add workflow scaffolds, containers, tests, and CI.

- Add example raw-data workflow configurations, even if full large-data execution is not performed during this build.



Phase 4: Review

- Run all demo commands.

- Identify incomplete or mocked components.

- Produce a candid implementation-status table with:

  - complete and runnable,

  - scaffolded but not yet production-ready,

  - planned future work.



## Final response format



When finished, provide:



1. A concise summary of the resource built.

2. The repository tree.

3. The most important files created.

4. Commands to run the demo.

5. A table mapping each component to Objective 1.

6. A candid list of remaining work before production deployment.

7. A suggested first set of 8–12 public C. gigas datasets to curate, described by dataset type and resilience context rather than inventing accession numbers.



Prioritize a credible, functioning, reproducible MVP over visual polish.