"""AREE as a Model Context Protocol (MCP) server.

`aree mcp` starts a stdio MCP server that lets an MCP client (Claude Desktop,
Claude Code, or any other) answer questions about AREE's studies, harmonized
evidence and candidate biomarkers by calling read-only tools. The tools are thin
wrappers over `query`, so the model's numbers come from the pipeline's outputs,
never from its own recall. See docs/querying_aree.md.

Requires the optional `mcp` extra (Python >= 3.10): pip install -e ".[mcp]"
"""
from __future__ import annotations

import query

INSTRUCTIONS = f"""\
You are answering questions from the Aquaculture Resilience Evidence Engine (AREE):
harmonized, provenance-tracked omics evidence on resilience-related phenotypes in
Pacific oyster (Crassostrea / Magallana gigas) and other aquaculture species.

How to answer:
- Get numbers only from these tools. Never invent genes, studies, effect sizes or accessions.
- Cite the study_id (and evidence_id where you quote a single record) behind every claim.
- Always state when evidence is simulated (simulated=true): it is synthetic demo data, not a finding.
- Map free-text phenotypes/stressors to controlled ids with list_vocabulary before filtering.
- Surface conflicting directions and heterogeneity; do not average them away.
- Keep the evidence classes distinct: resilience-associated, stress-response, disease-associated,
  exposure-only. A stress_response candidate is evidence of a response, not of a resilient animal.
- Never call a candidate validated, and never imply causation. Emerging-tier candidates require
  replication. {query.INTERPRETATION_NOTE}
- If a tool reports a missing artifact, say which `aree` command produces it.
"""


def _server_class():
    try:  # mcp >= 2
        from mcp.server.mcpserver import MCPServer
        return MCPServer
    except ImportError:
        pass
    try:  # mcp 1.x
        from mcp.server.fastmcp import FastMCP
        return FastMCP
    except ImportError as exc:
        raise SystemExit(
            "The AREE MCP server needs the 'mcp' package (Python >= 3.10). "
            'Install it with: pip install -e ".[mcp]"'
        ) from exc


def build_server(sources: query.Sources = query.DEFAULT_SOURCES):
    """Construct the MCP server with every AREE tool registered."""
    server = _server_class()(name="aree", instructions=INSTRUCTIONS)

    @server.tool()
    def aree_overview() -> dict:
        """Start here. Counts of studies, evidence records and candidates (real vs simulated),
        the phenotypes, stressors, feature types and species that have evidence, and which
        pipeline outputs are present."""
        return query.overview(sources)

    @server.tool()
    def list_vocabulary(name: str) -> dict:
        """Controlled-vocabulary term ids and labels, to turn free text into filter values.
        name: phenotype, stressor, assay_type, feature_type, tissue, life_stage,
        mapping_confidence, quality_flag, or species."""
        return query.list_vocabulary(name)

    @server.tool()
    def list_studies(species: str | None = None, stressor: str | None = None,
                     phenotype: str | None = None, assay_type: str | None = None,
                     include_simulated: bool = True) -> dict:
        """Registered studies, optionally filtered. stressor/phenotype/assay_type are controlled
        ids (e.g. temperature, survival, rnaseq); species matches any part of the name."""
        return query.list_studies(species, stressor, phenotype, assay_type, include_simulated, sources)

    @server.tool()
    def get_study(study_id: str) -> dict:
        """Full registration record for one study: accessions, citation, every comparison's
        tissue, life stage, original and standardized stressor, phenotype, sample size,
        limitations, QC and provenance."""
        return query.get_study(study_id, sources)

    @server.tool()
    def search_features(query_text: str, limit: int = 25, include_simulated: bool = True) -> dict:
        """Find genes/proteins/methylation regions/metabolites whose identifier (standardized,
        original, or orthogroup) contains query_text, e.g. 'LOC105333935', 'hsp70', 'OG0001234'.
        One row per standardized feature with its studies, phenotypes and directions."""
        return query.search_features(query_text, limit, include_simulated, sources)

    @server.tool()
    def get_evidence(feature: str | None = None, study_id: str | None = None,
                     phenotype: str | None = None, stressor: str | None = None,
                     species: str | None = None, feature_type: str | None = None,
                     tissue: str | None = None, max_adjusted_p: float | None = None,
                     include_simulated: bool = True, exact_feature: bool = True,
                     limit: int = 50) -> dict:
        """Harmonized evidence records (effect size, SE/CI, p, adjusted p, direction, context,
        mapping confidence, quality flags, source file, workflow version), most significant
        first. Reports direction counts and whether up/down directions conflict."""
        return query.get_evidence(feature, study_id, phenotype, stressor, species, feature_type,
                                  tissue, max_adjusted_p, include_simulated, exact_feature, limit, sources)

    @server.tool()
    def list_candidates(phenotype: str | None = None, tier: str | None = None,
                        evidence_class: str | None = None, feature: str | None = None,
                        species_taxid: int | None = None, feature_type: str | None = None,
                        include_simulated: bool = True, limit: int = 25) -> dict:
        """Ranked candidate biomarkers with transparent score components (k_studies, direction
        consistency, I-squared, pooled effect, adjusted p). tier: high_priority_cross_study,
        multi_omics_convergence, metabolite_linked_convergence, emerging. evidence_class:
        resilience_associated, stress_response, disease_associated, exposure_only."""
        return query.list_candidates(phenotype, tier, evidence_class, feature, species_taxid,
                                     feature_type, include_simulated, limit, sources)

    @server.tool()
    def get_evidence_card(feature: str, phenotype: str | None = None,
                          include_simulated: bool = True) -> dict:
        """The rendered evidence card for a candidate: supporting studies, effect summary,
        contexts, limitations, what the evidence does and does not show, and the recommended
        next validation step."""
        return query.get_evidence_card(feature, phenotype, include_simulated, sources)

    return server


def run(transport: str = "stdio") -> None:
    build_server().run(transport=transport)


if __name__ == "__main__":
    run()
