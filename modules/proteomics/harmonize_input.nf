// Harmonize a raw/wide protein abundance matrix (proteins x samples) into a
// long-format table joined to the sample sheet, for raw_reanalysis mode only.
//
// "Raw" here means a wide protein-abundance matrix (a MaxQuant/Skyline/DIA-NN
// intensity export, or a spectral-count/NSAF table), not raw mass-spec
// spectra; spectral reprocessing is out of scope.
//
// Only the samples named in the sample sheet enter the comparison: a matrix
// usually holds every sample of a study, and each comparison is run with its
// own sheet. The sheet's `group` column (or `condition`) must say `control` or
// `treatment` for every sample, so the reference level is never guessed.
//
// STATUS: executed by the raw-mode CI smoke test (tests/fixtures/proteomics_raw).
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

    def fail(msg):
        log(f"ERROR: {msg}")
        with open(f"{study_id}_{comparison_id}_harmonize_input.log", "w") as fh:
            fh.write("\\n".join(log_lines) + "\\n")
        sys.exit(1)

    log(f"AREE proteomics harmonize_input: study={study_id} comparison={comparison_id}")

    # First column = protein accession, remaining columns = one per sample.
    sep = "\\t"
    with open(matrix_path) as fh:
        header = fh.readline()
    if header.count(",") > header.count("\\t"):
        sep = ","
    wide = pd.read_csv(matrix_path, sep=sep)
    if wide.shape[1] < 2:
        fail("raw abundance matrix must have a protein accession column plus >=1 sample column")
    wide = wide.rename(columns={wide.columns[0]: "protein_accession"})
    if wide["protein_accession"].duplicated().any():
        dups = sorted(wide.loc[wide["protein_accession"].duplicated(), "protein_accession"].unique())
        fail(f"{len(dups)} duplicate protein_accession rows (e.g. {dups[:5]}); resolve them upstream")
    log(f"Loaded {len(wide)} proteins x {wide.shape[1] - 1} sample columns")

    samples = pd.read_csv(sample_sheet_path, sep=None, engine="python", dtype=str)
    if "group" not in samples.columns and "condition" in samples.columns:
        samples = samples.rename(columns={"condition": "group"})
    missing_cols = {"sample_id", "group"} - set(samples.columns)
    if missing_cols:
        fail(f"sample_sheet missing required columns: {sorted(missing_cols)}")
    bad = sorted(set(samples["group"]) - {"control", "treatment"})
    if bad:
        fail(f"sample_sheet group must be 'control' or 'treatment'; found {bad}")
    if samples["sample_id"].duplicated().any():
        fail("sample_sheet lists a sample_id more than once")

    absent = sorted(set(samples["sample_id"]) - set(wide.columns))
    if absent:
        fail(f"sample_sheet sample_id(s) not in the abundance matrix: {absent}")
    unused = [c for c in wide.columns[1:] if c not in set(samples["sample_id"])]
    if unused:
        log(f"{len(unused)} matrix column(s) are not in this comparison's sample sheet and are left out: {unused}")
    log(f"Comparison samples: {samples['group'].value_counts().to_dict()}")

    long_df = wide.melt(
        id_vars=["protein_accession"],
        value_vars=list(samples["sample_id"]),
        var_name="sample_id",
        value_name="abundance",
    ).merge(samples[["sample_id", "group"]], on="sample_id", how="left")

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
