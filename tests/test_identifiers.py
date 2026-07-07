from harmonize.identifiers import resolve_identifier


def test_exact_ncbi_gene_id_match():
    resolved = resolve_identifier("LOC105333935", "ncbi_gene_id")
    assert resolved.mapping_confidence == "exact"
    assert resolved.feature_id_standardized == "LOC105333935"
    assert resolved.orthogroup_id == "OG0000001"


def test_exact_uniprot_match():
    resolved = resolve_identifier("K1QN12", "uniprot_accession")
    assert resolved.mapping_confidence == "exact"
    assert resolved.feature_id_standardized == "LOC105333935"


def test_unresolved_stable_id():
    resolved = resolve_identifier("NOT_A_REAL_ID", "ncbi_gene_id")
    assert resolved.mapping_confidence == "unresolved"
    assert resolved.feature_id_standardized is None


def test_gene_symbol_exact_match_is_inferred_not_exact():
    resolved = resolve_identifier("hsp70", "gene_symbol")
    assert resolved.mapping_confidence == "inferred"
    assert resolved.feature_id_standardized == "LOC105333935"


def test_ambiguous_symbol_one_to_many():
    resolved = resolve_identifier("hsp70-like2", "gene_symbol")
    assert resolved.mapping_confidence == "one_to_many_ortholog"
    assert resolved.feature_id_standardized == "LOC105333935"


def test_ambiguous_symbol_many_to_one():
    resolved = resolve_identifier("cgi-actin", "gene_symbol")
    assert resolved.mapping_confidence == "many_to_one_ortholog"


def test_unknown_symbol_is_unresolved():
    resolved = resolve_identifier("unknown_EST_223", "gene_symbol")
    assert resolved.mapping_confidence == "unresolved"
    assert resolved.feature_id_standardized is None


def test_null_id_is_unresolved():
    resolved = resolve_identifier(None, "gene_symbol")
    assert resolved.mapping_confidence == "unresolved"
