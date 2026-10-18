"""Validate a study registration file against schemas/study.schema.json and
the controlled vocabularies referenced by its fields.

Errors are hard failures (schema violations, unknown vocabulary terms).
Warnings are soft signals (e.g. resilience_classification does not match the
ontology's default resilience_relevance for the chosen phenotype) that a
curator may have deliberately overridden with a documented rationale.
"""
from __future__ import annotations

import csv
from dataclasses import dataclass, field
from pathlib import Path

import jsonschema

from common import (
    DATA_DIR,
    DEMO_STUDIES_DIR,
    SCHEMAS_DIR,
    STUDIES_DIR,
    assembly_crosswalk_path,
    load_json,
    load_vocab,
    load_yaml,
    resolve_assembly,
    resolve_species,
)


@dataclass
class ValidationResult:
    path: str
    valid: bool
    errors: list = field(default_factory=list)
    warnings: list = field(default_factory=list)

    def __bool__(self):
        return self.valid


def _load_study_schema() -> dict:
    return load_json(SCHEMAS_DIR / "study.schema.json")


def validate_study_file(path) -> ValidationResult:
    path = Path(path)
    errors: list[str] = []
    warnings: list[str] = []

    try:
        study = load_yaml(path)
    except Exception as exc:  # malformed YAML
        return ValidationResult(str(path), False, errors=[f"Could not parse YAML: {exc}"])

    schema = _load_study_schema()
    validator = jsonschema.Draft7Validator(schema)
    for err in sorted(validator.iter_errors(study), key=lambda e: e.path):
        loc = "/".join(str(p) for p in err.path) or "<root>"
        errors.append(f"{loc}: {err.message}")

    if errors:
        return ValidationResult(str(path), False, errors=errors, warnings=warnings)

    stem = path.stem
    if stem not in ("_TEMPLATE",) and study.get("study_id") != stem:
        errors.append(
            f"study_id {study.get('study_id')!r} does not match filename stem {stem!r}"
        )

    # Real and simulated studies are kept in separate directories so the study
    # registry never lists demo fixtures as studies. A file in the wrong one is
    # an error rather than something to infer around.
    parent = path.resolve().parent
    simulated = bool(study.get("simulated", False))
    if parent == STUDIES_DIR.resolve() and simulated:
        errors.append(
            "simulated: true study is in registry/studies/; simulated demo studies "
            "belong in registry/demo_studies/"
        )
    elif parent == DEMO_STUDIES_DIR.resolve() and not simulated:
        errors.append(
            "study in registry/demo_studies/ is not simulated: true; real studies "
            "belong in registry/studies/"
        )

    phenotype_vocab = load_vocab("phenotype_ontology")
    stressor_ids = set(load_vocab("stressor_ontology").keys())
    tissue_ids = set(load_vocab("tissue_types").keys())
    life_stage_ids = set(load_vocab("life_stages").keys())
    assay_ids = set(load_vocab("assay_types").keys())
    quality_flag_ids = set(load_vocab("quality_flags").keys())

    # Species must resolve to a canonical taxon. A synonym is accepted (the
    # literature uses both Crassostrea and Magallana for this animal) but is
    # surfaced, because the registry row will read differently from the
    # canonical name that evidence is grouped by.
    species_name = study.get("species")
    species = resolve_species(species_name)
    if species is None:
        errors.append(
            f"species {species_name!r} is not in the species vocabulary "
            "(registry/controlled_vocabularies/species.yaml). Add the species there "
            "first — see docs/adding_a_species.md."
        )
    elif species.is_synonym:
        warnings.append(
            f"species {species_name!r} is an accepted synonym of "
            f"{species.scientific_name!r} (NCBI taxid {species.ncbi_taxid}); evidence "
            "will be grouped under the canonical name."
        )

    # The genome assembly must be one AREE has a record for, so that evidence can
    # be traced to a real, versioned reference rather than a free-text label.
    assembly_value = study.get("genome_assembly")
    assembly = resolve_assembly(assembly_value)
    if assembly is None:
        errors.append(
            f"genome_assembly {assembly_value!r} is not a known assembly "
            "(data/reference/genome_assemblies.yaml). Add it there with a verified "
            "NCBI accession — see docs/handling_genome_versions.md."
        )
    declared_flags = set(study.get("quality_flags") or [])

    if assembly is not None and species is not None and assembly.get("ncbi_taxid") != species.ncbi_taxid:
        # Mapping reads to a related species' genome is sometimes the only option
        # (e.g. Mytilus chilensis has no assembly). It is allowed only when the
        # study says so, because the gene identifiers then belong to the
        # reference species, not to the animal the evidence is about.
        mismatch = (
            f"genome_assembly {assembly['assembly_id']!r} belongs to taxid "
            f"{assembly.get('ncbi_taxid')} ({assembly.get('species')}), but species "
            f"{species_name!r} is taxid {species.ncbi_taxid}"
        )
        if "reference_species_mismatch" in declared_flags:
            warnings.append(
                f"{mismatch}. Accepted because the study declares "
                "reference_species_mismatch; identifiers will be standardized against "
                "the assembly's crosswalk, evidence partitioned under the study's species."
            )
        else:
            errors.append(
                f"{mismatch}. If the reads really were mapped to another species' genome, "
                "add the quality flag reference_species_mismatch and explain it in limitations."
            )
    elif "reference_species_mismatch" in declared_flags and assembly is not None and species is not None:
        errors.append(
            "quality flag reference_species_mismatch is declared, but the assembly and "
            "the species are the same taxon"
        )

    # A real study on an assembly with no crosswalk can be registered (so the
    # registry is honest about what exists) but cannot be harmonized yet.
    if assembly is not None and not study.get("simulated", False):
        crosswalk = assembly_crosswalk_path(assembly)
        if crosswalk is None:
            warnings.append(
                f"genome_assembly {assembly['assembly_id']!r} has no identifier_crosswalk "
                "in data/reference/genome_assemblies.yaml, so this study can be registered "
                "but not harmonized. Build one with `aree build-crosswalk --taxid "
                f"{assembly.get('ncbi_taxid')}` and register it there."
            )
        elif not crosswalk.exists():
            warnings.append(
                f"the identifier crosswalk registered for {assembly['assembly_id']!r} "
                f"({crosswalk}) does not exist on this machine; build it before harmonizing."
            )

    # Species in a hybridizing complex are routinely mislabelled in public
    # deposits. The study must say how its animals were assigned, and anything
    # short of genotyping is flagged rather than trusted.
    if species is not None and species.complex:
        identity = study.get("species_identity")
        if not identity:
            errors.append(
                f"species {species.scientific_name!r} is in the {species.complex} and is "
                "often mislabelled in public deposits; declare species_identity "
                "(basis: genotyped | morphology | collection_site | authors_statement | "
                "unstated) — see docs/adding_a_species.md."
            )
        elif identity.get("basis") != "genotyped" and "species_assignment_unverified" not in declared_flags:
            errors.append(
                f"species_identity.basis is {identity.get('basis')!r}, not genotyped, for a "
                f"species in the {species.complex}; add the quality flag "
                "species_assignment_unverified."
            )

    for assay in study.get("assay_type", []):
        if assay not in assay_ids:
            errors.append(f"assay_type {assay!r} is not in assay_types vocabulary")

    for flag in study.get("quality_flags", []) or []:
        if flag not in quality_flag_ids:
            errors.append(f"quality_flags entry {flag!r} is not in quality_flags vocabulary")

    seen_comparison_ids = set()
    for comp in study.get("comparisons", []):
        cid = comp.get("comparison_id")
        if cid in seen_comparison_ids:
            errors.append(f"duplicate comparison_id {cid!r} within study")
        seen_comparison_ids.add(cid)

        if comp.get("tissue") not in tissue_ids:
            errors.append(f"comparison {cid!r}: tissue {comp.get('tissue')!r} not in tissue_types vocabulary")
        if comp.get("life_stage") not in life_stage_ids:
            errors.append(f"comparison {cid!r}: life_stage {comp.get('life_stage')!r} not in life_stages vocabulary")
        if comp.get("stressor_standardized") not in stressor_ids:
            errors.append(
                f"comparison {cid!r}: stressor_standardized {comp.get('stressor_standardized')!r} "
                "not in stressor_ontology vocabulary"
            )
        phenotype_id = comp.get("phenotype")
        if phenotype_id not in phenotype_vocab:
            errors.append(f"comparison {cid!r}: phenotype {phenotype_id!r} not in phenotype_ontology vocabulary")
        else:
            expected_relevance = phenotype_vocab[phenotype_id]["resilience_relevance"]
            declared = comp.get("resilience_classification")
            if declared != expected_relevance:
                warnings.append(
                    f"comparison {cid!r}: resilience_classification={declared!r} differs from "
                    f"phenotype_ontology default ({expected_relevance!r}) for phenotype "
                    f"{phenotype_id!r}. Confirm this override is intentional and documented."
                )
        # A baseline contrast compares groups that differ in resilience under
        # common unstressed conditions; that is resilience evidence by
        # construction. Any other classification usually means the groups
        # differ by a treatment, which makes it an exposure_response.
        if (comp.get("contrast_type") == "baseline_phenotype_difference"
                and comp.get("resilience_classification") != "resilience"):
            warnings.append(
                f"comparison {cid!r}: contrast_type 'baseline_phenotype_difference' with "
                f"resilience_classification={comp.get('resilience_classification')!r}. A baseline "
                "contrast compares groups that differ in resilience; if the groups differ by a "
                "treatment, use 'exposure_response' (docs/contrast_types.md)."
            )

    # If a sample sheet was generated for this study, the BioProject it was
    # generated from must be the one the study claims. This catches a study
    # YAML drifting away from the data it actually describes, and it is the one
    # accession check that can be made offline.
    prov_path = (
        DATA_DIR / "studies" / str(study.get("study_id")) / "ena_provenance.json"
    )
    if prov_path.exists():
        declared = (study.get("accessions") or {}).get("bioproject")
        try:
            generated = load_json(prov_path).get("bioproject")
        except (ValueError, OSError) as exc:
            errors.append(f"could not read {prov_path.name}: {exc}")
        else:
            if declared != generated:
                errors.append(
                    f"accessions.bioproject is {declared!r} but the sample sheet in "
                    f"data/studies/{study.get('study_id')}/ was generated from "
                    f"{generated!r}. One of them is wrong."
                )

    # Excluded samples must name real samples, so a typo cannot silently
    # exclude nothing while the registry claims an exclusion.
    sheet_path = DATA_DIR / "studies" / str(study.get("study_id")) / "samplesheet.tsv"
    excluded = [
        (comp.get("comparison_id"), entry.get("sample_id"))
        for comp in study.get("comparisons", [])
        for entry in comp.get("excluded_samples") or []
    ]
    if excluded:
        if not sheet_path.exists():
            errors.append(
                "excluded_samples are declared but there is no sample sheet at "
                f"data/studies/{study.get('study_id')}/samplesheet.tsv to check them against"
            )
        else:
            with open(sheet_path, newline="") as fh:
                known = {row["sample_id"] for row in csv.DictReader(fh, delimiter="\t")}
            for cid, sample_id in excluded:
                if sample_id not in known:
                    errors.append(
                        f"comparison {cid!r}: excluded sample {sample_id!r} is not in the sample sheet"
                    )

    valid = len(errors) == 0
    return ValidationResult(str(path), valid, errors=errors, warnings=warnings)
