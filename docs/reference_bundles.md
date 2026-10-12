# Reference bundles

A raw reanalysis needs the reference its reads are mapped to: a transcriptome,
an annotation, and a transcript-to-gene map for RNA-seq, plus the genome for
methylation. AREE calls that set a **reference bundle** and keeps one per NCBI
RefSeq assembly under `data/reference/<ACCESSION>/`.

Bundles are fetched, never hand-assembled:

```bash
aree fetch-reference --assembly GCF_965363235.1            # RNA-seq bundle
aree fetch-reference --assembly GCF_963853765.1 --genome   # + genome FASTA for Bismark
```

`--assembly` takes any registered assembly id or accession from
`data/reference/genome_assemblies.yaml`. The command refuses an assembly that
is not registered, and one without an NCBI RefSeq annotation (NCBI publishes no
transcriptome or GTF for those; for example the Roslin oyster assembly).

## What a bundle contains

| File | Used by | Origin |
|---|---|---|
| `rna.fna` | RNA-seq `--transcript_fasta` (Salmon index) | NCBI `*_rna.fna.gz`, decompressed |
| `genomic.gtf` | RNA-seq `--gtf`, methylation `annotation_gtf` | NCBI `*_genomic.gtf.gz`, decompressed |
| `tx2gene.tsv` | RNA-seq `--tx2gene` (DESeq2/tximport) | Built from `genomic.gtf`: one `transcript_id`, NCBI GeneID pair per transcript feature, headerless, sorted |
| `genomic.fna` | methylation `genome_fasta` (Bismark) | NCBI `*_genomic.fna.gz`, only with `--genome` |
| `md5checksums.txt` | verification | NCBI's own checksum list for the folder |
| `reference_provenance.json` | provenance | Written by `aree fetch-reference` |

The compressed NCBI originals are kept next to the decompressed files, so a
re-run (or adding `--genome` later) reuses every file whose MD5 still matches
NCBI's and downloads only what is missing.

Salmon and Bismark indexes are **not** part of the bundle. The workflows build
them on first use. A Salmon index built once can be reused by setting
`salmon_index` in the study config, as the oyster studies do.

## What gets checked

Every downloaded file is compared with NCBI's published MD5 before it is used;
a mismatch deletes the file and stops. `reference_provenance.json` then records:

- source URL, NCBI MD5 and SHA-256 of every compressed file, and SHA-256 of
  every file the workflows read;
- the GTF header (`genome-build`, `annotation-source`, `annotation-date`);
- **whether the GTF's annotation release is the one registered** for the
  assembly. NCBI can re-annotate an assembly in place. If it has, the command
  exits non-zero: update `genome_assemblies.yaml` and rebuild the species'
  crosswalk before using the bundle, or evidence would claim an annotation it
  did not come from;
- **how many of the bundle's GeneIDs the species' identifier crosswalk knows.**
  The workflow's gene IDs are what `aree harmonize` resolves, so anything below
  100% becomes `unresolved` evidence. It should be 100% when the crosswalk and
  the bundle come from the same annotation release.

Re-verify a bundle at any time, for example on a new machine:

```bash
aree fetch-reference --assembly GCF_965363235.1 --check
```

`--plan` resolves the NCBI folder and checksums without downloading sequence
data, which is a quick way to confirm an assembly is fetchable.

## Where bundles live

Bundles are large (the oyster genomic GTF alone is 459 MB decompressed) and are
not committed: `data/reference/GCF_*` is gitignored. On a workstation with
external storage, make `data/reference/<ACCESSION>` a symlink before fetching,
and the bundle is written through it:

```bash
mkdir -p ~/aree-data/reference/GCF_965363235.1
ln -sfn ~/aree-data/reference/GCF_965363235.1 data/reference/GCF_965363235.1
aree fetch-reference --assembly GCF_965363235.1
```

`aree start-study` links every `data/reference/GCF_*` symlink into a study's
worktree, so one bundle serves every study branch on that machine. See
*Storage layout* in [first_raw_reanalysis.md](first_raw_reanalysis.md).

## What the workflow manifest records

The RNA-seq workflow stages the transcriptome, GTF and tx2gene it reads, plus
the bundle's `reference_provenance.json` when it sits beside them, into the
manifest step. The manifest's `reference` block then carries each file's
SHA-256 and whether it still matches the bundle provenance, and names the
bundle, its accession and its annotation source. Warnings are raised when a
file differs from what was fetched, when the bundle is for a different assembly
than the run declares, when its annotation release is not the registered one,
or when a raw run used reference files with no bundle provenance at all.

## Using a bundle in a study config

```groovy
params {
    rnaseq {
        genome_assembly    = 'GCF_965363235.1'
        annotation_version = 'GCF_965363235.1-RS_2025_08'
    }
    transcript_fasta = 'data/reference/GCF_965363235.1/rna.fna'
    gtf              = 'data/reference/GCF_965363235.1/genomic.gtf'
    tx2gene          = 'data/reference/GCF_965363235.1/tx2gene.tsv'
    salmon_index     = null   // built on first run; point here afterwards to reuse it
}
```

`genome_assembly` and `annotation_version` go into the workflow manifest, and the
study YAML's `genome_assembly` selects the identifier crosswalk at harmonization
(see [identifier_mapping.md](identifier_mapping.md#which-crosswalk-is-in-force)),
so the two must name the same assembly.

## Related documentation

- [first_raw_reanalysis.md](first_raw_reanalysis.md)
- [handling_genome_versions.md](handling_genome_versions.md)
- [adding_a_species.md](adding_a_species.md)
- [mussel_clam_integration_plan.md](mussel_clam_integration_plan.md)
