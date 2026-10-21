# Querying AREE: the Explore page and the MCP server

AREE's evidence is generated as tables: the harmonized evidence table, ranked
candidates, evidence cards and the study registry. There are two ways to query
them without writing pipeline code. A third, a chat box on the website, is being
built on the same query layer ([Ask AREE](ask_aree.md)). Both read the same outputs through the same
read-only query layer (`src/query/`), so they always agree with each other and
with the evidence cards.

| | Explore page | MCP server |
|---|---|---|
| Who it's for | Anyone with a browser | Anyone using Claude (Desktop, Code) or another MCP client |
| How you ask | Filters, gene search, SQL | Plain-language questions |
| Runs where | In your browser (DuckDB-WASM) | On your machine, beside a checkout of this repository |
| Needs | Nothing | `pip install -e ".[mcp]"` (Python ≥ 3.10) |
| Writes anything? | No | No. Every tool is read-only |

## The Explore page

<https://sr320.github.io/AREE/explore.html> ([source](explore.qmd)) has five tabs:

- **Gene / feature lookup.** Search a standardized or original identifier, or an
  orthogroup. The tab shows every evidence record across studies, a summary that
  flags significant records disagreeing in direction, and an effect-size plot with
  confidence intervals once the search narrows to one feature.
- **Candidates.** Ranked candidates, filterable by phenotype, tier, evidence class and contrast
  type, with a "Predictive baseline only" toggle and a search over gene IDs, symbols and
  descriptions ([contrast types](contrast_types.md)).
- **Filter evidence.** Filter by phenotype, stressor, species, feature type, tissue and
  an adjusted-p cutoff.
- **Studies.** Every registered study.
- **SQL.** Free DuckDB SQL over the `evidence`, `candidates` and `studies` tables, with
  worked examples (genes significant in three or more studies, conflicting directions,
  evidence volume by context).

Every result downloads as CSV. The published page holds real studies only:
`aree export-site-data` leaves the simulated demo studies (`registry/demo_studies/`)
and their evidence out unless run with `--include-simulated`, and the Pages build
never harmonizes them. For a locally built page that includes them, the page's
toggle still separates them, and rows keep their `simulated` column either way.

### How it works

`aree export-site-data` writes `docs/explore/data/`:

| File | Contents |
|---|---|
| `evidence.parquet` | Every harmonized evidence record, provenance columns included, sorted by standardized then original feature id |
| `features.parquet` | One row per distinct identifier combination (standardized, original, orthogroup), with its record count |
| `candidates.parquet` | Every ranked candidate from `aree build-evidence-cards`, with `gene_symbol` and `gene_description` for real genes and `is_predictive_baseline` as a true/false column |
| `studies.parquet` | One row per study YAML |
| `cards.parquet` | The markdown of every exported candidate's evidence card, keyed by file name. The page doesn't load it; the MCP server reads it when it answers from published data |
| `meta.json` | Export time, row counts (real / simulated), filter values, vocabulary labels |

The page loads [DuckDB-WASM](https://duckdb.org/docs/api/wasm/overview) and queries
the Parquet files in place over HTTP range requests, fetching only the columns and row
groups a query needs. No query leaves the browser. With ~840k evidence records (≈750 MB
as TSV, ≈32 MB as Parquet), these are the measured costs:

| Action | Transferred |
|---|---|
| Open the page | ≈0.25 MB (file footers, candidates, studies) |
| First gene lookup | ≈1.3 MB feature index (once), plus ≈1.5 MB of evidence |
| Later lookups | ≈1.5–4 MB each |
| A broad evidence filter or SQL scan | Up to the whole evidence file |

A lookup first resolves what you typed against the small feature index (exact match,
else substring). It then fetches evidence by exact id, and DuckDB skips every row group
of the sorted file whose id range excludes it. Broad filters and SQL run only when you
press their button, never on page load.

The Pages workflow runs `aree export-site-data` after the pipeline. To build the page
locally:

```bash
pip install -e ".[site]"
make demo                 # and/or make real-pool
aree export-site-data
quarto preview docs       # then open /explore.html
```

`quarto preview` (or any HTTP server) is required. Opening the HTML from disk does not
work, because browsers block `fetch` from `file://` URLs.

## The MCP server

The [Model Context Protocol](https://modelcontextprotocol.io) lets an AI assistant call
tools you provide. `aree mcp` exposes AREE's query layer as tools, so you can ask Claude
questions such as:

- "Which genes are significant for thermal tolerance in more than one real study? Do any
  of them disagree in direction?"
- "Show me the evidence card for LOC105333935 under thermal tolerance and tell me what
  validation it needs."
- "Which registered studies measured salinity tolerance in larvae, and are their raw
  data available?"
- "Are there multi-omics convergence candidates for pathogen challenge?"

The assistant answers by calling the tools, so its numbers come from AREE's outputs,
not from its memory. The server's instructions tell it to cite `study_id`s (and
`evidence_id`s for single records), to say when evidence is simulated, to report
conflicting directions rather than averaging them away, to keep the evidence classes
apart, and never to call a candidate validated or causal.

### Tools

| Tool | Returns |
|---|---|
| `aree_overview` | Counts (real vs. simulated), the phenotypes, stressors, feature types and species with evidence, and which pipeline outputs are present |
| `list_vocabulary` | Controlled-vocabulary ids and labels (phenotype, stressor, tissue, …), used to map free text onto filter values |
| `list_studies` | Registered studies, filterable by species, stressor, phenotype and assay |
| `get_study` | One study's full registration record, plus how many evidence records it has contributed |
| `search_features` | Features whose identifier or orthogroup matches, one summary row per feature |
| `get_evidence` | Evidence records with effect, SE/CI, p, adjusted p, context, mapping confidence, quality flags and provenance; flags conflicting directions |
| `find_replicated_features` | Features significant in at least *k* studies (default 2, adjusted p ≤ 0.05), with how many studies went up versus down and whether directions conflict |
| `list_candidates` | Ranked candidates with their transparent score components |
| `get_evidence_card` | The rendered evidence card markdown for a candidate |

Results are capped (50 records by default, 500 at most) and report when they are
truncated. Use the Explore page or the TSVs in `reports/` for bulk export.

### Setup

```bash
pip install -e ".[mcp]"      # Python >= 3.10
```

The `mcp` extra holds the MCP SDK below 2.0, because an `mcp` 2.x server rejects the
initialize handshake that Claude Desktop and Claude Code send.

The server answers from one of two places:

- **The published data, with no pipeline run.** `aree fetch-site-data` downloads the
  tables behind the Explore page, evidence cards included (about 70 MB), into
  `docs/explore/data/`. Then serve them with `aree mcp --data-dir docs/explore/data`.
  Like the site, they hold real studies only. Fetch again to pick up a newer build.
- **Your checkout's pipeline outputs.** Plain `aree mcp` reads `reports/`. Run the
  pipeline first (`make demo`, `make real-pool`, then `aree build-evidence-cards`). If
  an output is missing, the tools report which command produces it.

Either way the queries run in DuckDB. Published Parquet is read in place, and a
`reports/` evidence table is loaded once per server start (about 20 s for a million
records).

**Claude Code**, from the repository root:

```bash
claude mcp add aree -- aree mcp
```

For the published data, give the directory as an absolute path:
`claude mcp add aree -- aree mcp --data-dir /absolute/path/to/AREE/docs/explore/data`.

**Claude Desktop**: quit Claude Desktop first (⌘Q). A running Desktop rewrites
`claude_desktop_config.json` from its own copy and discards edits made meanwhile.
Then add this to the file (`~/Library/Application Support/Claude/` on macOS), using the
absolute path printed by `which aree`:

```json
{
  "mcpServers": {
    "aree": {
      "command": "/absolute/path/to/env/bin/aree",
      "args": ["mcp"]
    }
  }
}
```

For the published data, use `"args": ["mcp", "--data-dir", "/absolute/path/to/AREE/docs/explore/data"]`.
Reopen Claude Desktop. The AREE tools then appear in the tools menu. If they don't,
`~/Library/Logs/Claude/mcp-server-aree.log` shows why.

**Other MCP clients** work the same way: run `aree mcp` over stdio. For a shared
deployment, `aree mcp --transport streamable-http` serves over HTTP. Put
authentication in front of it before exposing it beyond your machine.

## Interpreting answers

Both interfaces return the same caveat with every evidence result, and it applies
to whatever you conclude from them. Evidence records are statistical associations
from individual studies. A significant record or a high candidate score is not a
validated biomarker and does not establish causation. Emerging-tier candidates
require replication, and rows with `simulated = true` are synthetic demo data. See
[Interpreting candidate scores](interpreting_candidate_scores.md) and
[Interpreting meta-analysis outputs](interpreting_meta_analysis.md).

A language model can still misread a table, so check any claim you intend to rely on
against the cited `study_id`/`evidence_id` on the Explore page or in
`reports/evidence/evidence_table.tsv`.
