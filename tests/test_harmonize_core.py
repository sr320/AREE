"""Unit tests for harmonization dispatch helpers."""
from __future__ import annotations

import pytest

from harmonize.core import _comparison_assay_type


def test_comparison_assay_type_single_and_multi_assay():
    single = {"study_id": "S1", "assay_type": ["rnaseq"]}
    assert _comparison_assay_type(single, {"comparison_id": "c1"}) == "rnaseq"

    multi = {"study_id": "S2", "assay_type": ["rnaseq", "proteomics"]}
    assert _comparison_assay_type(multi, {"comparison_id": "c1", "assay_type": "proteomics"}) == "proteomics"

    with pytest.raises(ValueError, match="multiple assay_types"):
        _comparison_assay_type(multi, {"comparison_id": "c1"})

    with pytest.raises(ValueError, match="not in study assay_type"):
        _comparison_assay_type(multi, {"comparison_id": "c1", "assay_type": "metabolomics"})
