import json



def test_harmonize_emits_provenance_manifest(isolated_reports):
    from harmonize.core import harmonize_study
    harmonize_study("GIGAS_SAL04", date_generated="2026-01-01")

    manifest_path = isolated_reports["manifests_dir"] / "GIGAS_SAL04_low_salinity_vs_control_manifest.json"
    assert manifest_path.exists()

    manifest = json.loads(manifest_path.read_text())
    # Required provenance fields per docs/design.md section 4.
    for field in ["input_checksum", "workflow_version", "date_generated", "generated_by", "parameters"]:
        assert field in manifest
    assert manifest["input_checksum"].startswith("sha256:")


def test_processed_only_study_manifest_carries_warning(isolated_reports):
    from harmonize.core import harmonize_study
    harmonize_study("GIGAS_SAL04", date_generated="2026-01-01")
    manifest_path = isolated_reports["manifests_dir"] / "GIGAS_SAL04_low_salinity_vs_control_manifest.json"
    manifest = json.loads(manifest_path.read_text())
    assert any("raw_data_not_available" in w for w in manifest["warnings"])
