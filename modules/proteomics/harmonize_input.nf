// Harmonize a raw/wide protein abundance matrix (proteins x samples) into a
// long-format table plus a validated sample sheet join, for raw_reanalysis
// mode only.
//
// STATUS: structurally complete DSL2 process with a real pandas script. The
// script implements actual wide-to-long reshape, dtype coercion, and
// structural validation (duplicate accessions, sample-sheet coverage). It has
// NOT been executed against a real raw abundance matrix in this build — no
// such file exists in this repository (see workflows/proteomics/README.md).
// "Raw" here means a raw/unfiltered wide protein-abundance matrix (e.g. a
// MaxQuant/Skyline/DIA-NN export), not raw mass-spec spectra — spectral
// reprocessing is explicitly out of scope for this scaffold.

process HARMONIZE_INPUT {
    tag "${study_id}:${comparison_id}"
    label 'process_low'
    container 'python:3.11-slim'
    publishDir "${params.outdir}/proteomics/harmonize_input", mode: params.publish_mode

    input:
    val study_id
    val comparison_id
    path raw_abundance_matrix
    path sample_sheet

    output:
    tuple val(study_id), val(comparison_id), path("${study_id}_${comparison_id}_long_abundance.tsv"), emit: long_table
    path "${study_id}_${comparison_id}_harmonize_input.log", emit: log
    path "versions.yml", emit: versions

    script:
    """
    #!/usr/bin/env python3
    import sys
    import pandas as pd

    study_id = "${study_id}"
    comparison_id = "${comparison_id}"
    matrix_path = "${raw_abundance_matrix}"
    sample_sheet_path = "${sample_sheet}"

    log_lines = []

    def log(msg):
        print(msg)
        log_lines.append(str(msg))

    log(f"AREE proteomics harmonize_input: study={study_id} comparison={comparison_id}")
    log(f"Reading raw wide abundance matrix from {matrix_path}")

    # Wide matrix expected shape: first column = protein accession, remaining
    # columns = one per sample (raw or log-space intensities from a
    # MaxQuant/Skyline/DIA-NN-style export). Sep is sniffed between tab/comma.
    sep = "\\t"
    with open(matrix_path) as fh:
        header = fh.readline()
    if header.count(",") > header.count("\\t"):
        sep = ","

    wide = pd.read_csv(matrix_path, sep=sep)
    if wide.shape[1] < 2:
        log("ERROR: raw abundance matrix must have a protein accession column plus >=1 sample column")
        sys.exit(1)

    id_col = wide.columns[0]
    wide = wide.rename(columns={id_col: "protein_accession"})

    n_before = len(wide)
    dup_mask = wide["protein_accession"].duplicated(keep=False)
    n_dup = int(dup_mask.sum())
    if n_dup > 0:
        log(f"WARNING: {n_dup} duplicate protein_accession rows found; collapsing via mean")
        wide = wide.groupby("protein_accession", as_index=False).mean(numeric_only=True)
    log(f"Loaded {n_before} input rows -> {len(wide)} unique protein_accession rows")

    sample_cols = [c for c in wide.columns if c != "protein_accession"]
    log(f"Detected {len(sample_cols)} sample columns: {sample_cols}")

    # Validate against sample sheet (expects columns: sample_id, group; group
    # values typically 'treatment'/'control' or comparison-specific labels).
    samples = pd.read_csv(sample_sheet_path, sep=None, engine="python")
    required_cols = {"sample_id", "group"}
    missing_cols = required_cols - set(samples.columns)
    if missing_cols:
        log(f"ERROR: sample_sheet missing required columns: {sorted(missing_cols)}")
        sys.exit(1)

    unmatched = sorted(set(samples["sample_id"]) - set(sample_cols))
    if unmatched:
        log(f"WARNING: sample_sheet sample_id(s) not found in abundance matrix columns: {unmatched}")

    # Reshape wide -> long: one row per (protein_accession, sample_id, abundance)
    long_df = wide.melt(
        id_vars=["protein_accession"],
        value_vars=[c for c in sample_cols if c in set(samples["sample_id"]) or True],
        var_name="sample_id",
        value_name="abundance",
    )
    long_df = long_df.merge(samples[["sample_id", "group"]], on="sample_id", how="left")

    unmapped_samples = long_df.loc[long_df["group"].isna(), "sample_id"].unique().tolist()
    if unmapped_samples:
        log(f"WARNING: {len(unmapped_samples)} sample column(s) had no group in sample_sheet: {unmapped_samples}")

    out_path = f"{study_id}_{comparison_id}_long_abundance.tsv"
    long_df.to_csv(out_path, sep="\\t", index=False)
    log(f"Wrote long-format abundance table: {out_path} ({len(long_df)} rows)")

    with open(f"{study_id}_{comparison_id}_harmonize_input.log", "w") as fh:
        fh.write("\\n".join(log_lines) + "\\n")

    with open("versions.yml", "w") as fh:
        fh.write("HARMONIZE_INPUT:\\n")
        fh.write(f"    python: \\"{sys.version.split()[0]}\\"\\n")
        fh.write(f"    pandas: \\"{pd.__version__}\\"\\n")
    """

    stub:
    """
    touch ${study_id}_${comparison_id}_long_abundance.tsv
    touch ${study_id}_${comparison_id}_harmonize_input.log
    touch versions.yml
    """
}
