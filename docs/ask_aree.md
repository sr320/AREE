# Ask AREE: the public chat service

Ask AREE answers plain-language questions about AREE's evidence on the website, for
visitors who don't use the [Explore page](explore.qmd) or an MCP client. Claude
writes each answer by calling AREE's read-only query layer as tools, the same nine
queries the [MCP server](querying_aree.md#the-mcp-server) offers. Its numbers
therefore come from AREE's published tables, and it cites the `study_id`s behind them.

Status: the service (`src/chat/`) and the chat page ([`ask.qmd`](ask.qmd)) run locally
and are tested. The Cloud Run deployment is next (see [Roadmap](#roadmap)). Until it is
live, the page isn't linked from the site's navigation, and on the public site it says
the chat isn't available yet.

## How it works

```
GitHub Pages (static)              Ask AREE service (one container)          Anthropic API
chat box  ── POST /api/chat ──▶    chat.server  (limits, key)        ──▶     Claude
          ◀── event stream ──      chat.conversation  (tool loop)
                                   chat.tools → query  (DuckDB over the published Parquet)
```

- **`chat.server`** is a FastAPI app. `POST /api/chat` takes a question and streams
  the answer back as server-sent events; `GET /api/health` reports status, model and
  data build. The service holds the Anthropic API key, so the static site never does.
- **`chat.conversation`** runs one turn. Claude may call tools for up to six rounds,
  and the round after that must answer. Answer text streams as it is written. Each
  tool call and its result are streamed too, so the page can show how an answer was
  reached. Conversations live in memory for an hour, up to ten questions each.
- **`chat.tools`** adapts the query layer for a public audience:
  - only real studies are served, never simulated demo data;
  - evidence records carry a compact column set (an `evidence_id` leads back to full
    provenance) and results are capped at about 10k tokens;
  - every tool input is validated against its JSON schema before it runs.
- **`chat.prompt`** combines the MCP server's rules (cite study ids, report conflicting
  directions, keep evidence classes apart, never call a candidate validated or causal)
  with website rules: answer only about AREE, lead with the answer, say plainly when
  AREE holds no evidence.

The model is Claude Opus 5.5 at medium effort by default. If its safety classifiers
decline a question, a server-side fallback answers on another model rather than
returning a refusal. The system prompt and tool definitions are cached, so the cost
of a question is mostly its tool results and the answer.

## Run it locally

```bash
pip install -e ".[chat]"          # Python >= 3.10
aree fetch-site-data              # the published tables, about 70 MB
export ANTHROPIC_API_KEY=...      # your key; never commit it
aree chat-serve --data-dir docs/explore/data
```

The service listens on <http://127.0.0.1:8000>. `chat-serve` also accepts calls from
pages on `localhost`, so a locally previewed site can use it. To try it from a
terminal:

```bash
curl -N -X POST http://127.0.0.1:8000/api/chat -H 'Content-Type: application/json' -d '{"message": "Which phenotypes does AREE have evidence for?"}'
```

### The chat page

With `aree chat-serve` running, preview the site and open the Ask page:

```bash
quarto preview docs
```

Then go to `/ask.html` on the preview's address. A page served from `localhost` sends
questions to `http://127.0.0.1:8000`; add `?api=http://127.0.0.1:PORT` to the page's
address for a service on another port. Only local addresses are accepted there, so a
link can't point the page at someone else's server.

The page (`docs/ask/ask.js`) renders each answer as markdown, sanitized with
DOMPurify because model output is untrusted. Each study id that a query returned is
linked to its registration record. **How this answer was found** under each answer
lists every query and the exact result Claude received.

Deployment sets two values at the top of `docs/ask/ask.js`: `PRODUCTION_API`, the
service's address, and `TURNSTILE_SITE_KEY` when the service requires a bot check.

### The event stream

Each event is a line `data: {...}` with a `type`:

| `type` | Meaning |
|---|---|
| `start` | The turn began; carries the `conversation_id` to send with follow-up questions |
| `text` | A fragment of the answer |
| `tool_call` | A query Claude is running: tool `name` and `input` |
| `tool_result` | That query's result, exactly as Claude received it |
| `round_reset` | Discard answer text since the last tool call (a retried round) |
| `notice` | Something the reader should know, such as an answer cut off at the length limit |
| `error` | The question could not be answered (`code`, `message`) |
| `done` | The turn finished: tools used, token usage, estimated cost |

## Settings

Everything is an environment variable, so a deployment needs no code change.

| Variable | Default | Meaning |
|---|---|---|
| `ANTHROPIC_API_KEY` | none | The API key (in production, from a secret manager) |
| `AREE_DATA_DIR` | `reports/` outputs | Published tables to answer from (`aree fetch-site-data`) |
| `AREE_CHAT_MODEL` | `claude-opus-5-5` | Claude model |
| `AREE_CHAT_EFFORT` | `medium` | Reasoning effort: `low`, `medium`, `high` |
| `AREE_CHAT_DAILY_BUDGET_USD` | `10` | Estimated spend per UTC day before the chat pauses |
| `AREE_CHAT_PER_IP_HOURLY` / `_DAILY` | `10` / `30` | Questions per visitor |
| `AREE_CHAT_MAX_QUESTION_CHARS` | `1000` | Longest question accepted |
| `AREE_CHAT_MAX_TOOL_ROUNDS` | `6` | Tool rounds per question |
| `AREE_CHAT_MAX_TURNS` | `10` | Questions per conversation |
| `AREE_CHAT_ALLOWED_ORIGINS` | `https://sr320.github.io` | Pages allowed to call the API (comma-separated) |
| `AREE_CHAT_ALLOW_LOCALHOST` | `false` | Also allow `localhost` pages (`chat-serve` turns this on) |
| `AREE_CHAT_TURNSTILE_SECRET` | none | When set, every question needs a valid Cloudflare Turnstile token |
| `AREE_CHAT_TRUST_PROXY` | `false` | Take the visitor's address from `X-Forwarded-For` (set on Cloud Run) |
| `AREE_CHAT_LOG_QUESTIONS` | `false` | Include question text in the logs, for quality review |

## Cost and abuse controls

- **Daily budget.** Each response's cost is estimated from the token usage the API
  reports. Once the day's total passes `AREE_CHAT_DAILY_BUDGET_USD`, the chat answers
  with a pointer to the Explore page until UTC midnight. A turn abandoned midway still
  counts.
- **Per-visitor limits** on questions per hour and per day, and a per-question cap
  on tool rounds and answer length.
- **Bot check.** With `AREE_CHAT_TURNSTILE_SECRET` set, the service verifies a
  Cloudflare Turnstile token on every question.
- **Origins.** Only the AREE site (and, locally, `localhost`) may call the API from a browser.
- **The hard backstop** is a monthly spend limit on the API key's workspace in the
  Claude Console.

These limits are held in memory by one process. Run a single instance
(`--max-instances=1` on Cloud Run) to keep them exact. One instance serves many
visitors at once, because a turn mostly waits on the API.

Logs are one JSON line per turn: whether it completed, latency, tools used, token
usage and estimated cost. Question text is logged only with `AREE_CHAT_LOG_QUESTIONS`,
and visitor addresses never are.

## Interpreting answers

The chat inherits every caveat of the evidence it reads. An evidence record is a
statistical association from one study, and a high candidate score is not a
validated biomarker. A language model can also misread a table, so check any claim
you intend to rely on against the cited `study_id` or `evidence_id` on the Explore page.

## Roadmap

1. **Answer-quality evals**: about 25 questions with checkable expected facts, run
   before launch and whenever the prompt or model changes.
2. **Cloud Run deployment**: a container with the published tables baked in, deployed
   by GitHub Actions after each site build, with the API key in Secret Manager. Then set
   `PRODUCTION_API` in `docs/ask/ask.js` and add the page to the site's navigation.
