# AREE container strategy

STATUS: documentation only. No images referenced below have been built or
pushed as part of this build. This describes the *intended* strategy so a lab
with real data and a container registry can make the scaffolds runnable.

## Approach

Each Nextflow process in `modules/<assay>/*.nf` declares a `container`
directive pointing at a versioned public image (mostly existing
`biocontainers`/`quay.io/biocontainers` or `bioconductor` images so nothing
bespoke needs to be built for standard steps). Where AREE-specific glue code
runs (e.g. schema-conformant TSV emission, manifest writing), the process uses
a small shared AREE utility image instead of installing ad hoc packages inline.

This repository does **not** vendor Dockerfiles for third-party tools — it
references existing, versioned upstream images by tag. This keeps the
scaffold honest: we are declaring *which* image a step is designed to run in,
not asserting that we built/tested it here.

## Declared images by workflow

### RNA-seq (`workflows/rnaseq/`, `modules/rnaseq/`)

| Step | Image (illustrative tag) |
|---|---|
| FastQC | `biocontainers/fastqc:0.12.1--hdfd78af_0` |
| fastp | `biocontainers/fastp:0.23.4--h5f740d0_0` |
| Salmon (pseudoalignment + quant) | `combinelab/salmon:1.10.3` |
| STAR (alternative alignment) | `biocontainers/star:2.7.11b--h43eeafb_0` |
| Sample QC / MultiQC | `biocontainers/multiqc:1.21--pyhdfd78af_0` |
| DESeq2 differential expression | `bioconductor/bioconductor_docker:RELEASE_3_18` |
| Report render (Quarto) | `ghcr.io/quarto-dev/quarto:1.5.57` |

### Methylation / WGBS / EM-seq (`workflows/methylation/`, `modules/methylation/`)

| Step | Image (illustrative tag) |
|---|---|
| FastQC | `biocontainers/fastqc:0.12.1--hdfd78af_0` |
| Trim Galore | `biocontainers/trim-galore:0.6.10--hdfd78af_0` |
| Bismark (genome prep + align + dedup) | `biocontainers/bismark:0.24.2--hdfd78af_1` |
| Bismark methylation extractor | `biocontainers/bismark:0.24.2--hdfd78af_1` |
| methylKit DML/DMR + annotation (R/Bioconductor) | `bioconductor/bioconductor_docker:RELEASE_3_18` |
| Report render (Quarto) | `ghcr.io/quarto-dev/quarto:1.5.57` |

### Proteomics (`workflows/proteomics/`, `modules/proteomics/`)

| Step | Image (illustrative tag) |
|---|---|
| Table harmonization / normalization (Python) | `python:3.11-slim` (AREE utility layer, see below) |
| Differential abundance (limma, R/Bioconductor) | `bioconductor/bioconductor_docker:RELEASE_3_18` |
| Report render (Quarto) | `ghcr.io/quarto-dev/quarto:1.5.57` |

### Metabolomics (`workflows/metabolomics/`, `modules/metabolomics/`)

| Step | Image (illustrative tag) |
|---|---|
| Feature table intake / QC / normalization (Python) | `python:3.11-slim` (AREE utility layer, see below) |
| Differential abundance (R, limma/Wilcoxon) | `bioconductor/bioconductor_docker:RELEASE_3_18` |
| Report render (Quarto) | `ghcr.io/quarto-dev/quarto:1.5.57` |

## AREE utility layer

Several lightweight Python steps (harmonization-shaping, manifest emission,
QC-metric summarization) are declared against a plain `python:3.11-slim` base
plus `pandas`/`pyyaml` rather than a bespoke image, since no Dockerfile is
built in this MVP. A production deployment should replace this with a pinned,
built-and-pushed `ghcr.io/<org>/aree-utils:<version>` image containing an
exact `requirements.txt` lock — tracked as future work, not implemented here.

## Apptainer / Singularity

Every `nextflow.config` in `workflows/*/` includes an `apptainer` profile that
reuses the same image references (Apptainer can pull Docker Hub / `quay.io`
images directly via the `docker://` URI scheme), for HPC environments where
Docker is unavailable. This has not been tested against a real Apptainer
installation in this build.

## What is real vs. aspirational here

- Real: the process-level `container` directives are syntactically valid
  Nextflow and reference images that exist upstream (tags may need bumping to
  whatever is current when actually run).
- Aspirational: nobody has pulled, run, or verified compatibility of these
  images against the exact process command lines in this repository. Version
  pins should be re-verified before first real execution.
