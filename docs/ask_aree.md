# Ask AREE: the public chat service

Ask AREE answers plain-language questions about AREE's evidence on the website, for
visitors who don't use the [Explore page](explore.qmd) or an MCP client. Claude
writes each answer by calling AREE's read-only query layer as tools, the same nine
queries the [MCP server](querying_aree.md#the-mcp-server) offers. Its numbers
therefore come from AREE's published tables, and it cites the `study_id`s behind them.

Status: live at [Ask AREE](ask.qmd). The service runs on Cloud Run
([Deploy to Cloud Run](#deploy-to-cloud-run)) and is redeployed with the fresh tables after
each site build.

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

## Answer-quality eval

Before launch, and whenever the system prompt, the tools, the model or its settings
change, run the eval: 26 questions in
[`evals/ask_aree/cases.yaml`](https://github.com/sr320/AREE/blob/main/evals/ask_aree/cases.yaml),
each sent through the chat's real turn loop and graded. The questions cover every tool,
questions AREE holds no evidence for, attempts to draw out an overclaim ("is this gene
a validated biomarker?"), off-topic requests and a follow-up conversation.

An answer passes only when every check passes:

| Check | How | Passes when |
|---|---|---|
| facts | program | The answer states the numbers and names the ids AREE's data say it should. Expected values are recomputed from the data each run, so they follow the published build. Large numbers may be rounded ("1.24 million", within 5%). |
| grounded ids | program | Every study id cited is a real registered study: none invented, none simulated. |
| answers | judge | The answer addresses the question, or says plainly that AREE holds nothing on it, or declines an off-topic request. |
| faithful | judge | Every claim about AREE's data is supported by the tool results the chat received. |
| no overclaim | judge | Nothing is called validated, proven or causal. |
| case criterion | judge | The question's own requirement in `cases.yaml` is met. |

The judge is Claude Sonnet 5.5, a different model from the one under test. It reads the
question, every tool result as the chat received it, and the answer, and gives a reason
for each verdict.

```bash
python scripts/eval_ask_aree.py --list           # the questions and expected facts; free
python scripts/eval_ask_aree.py --reps 2         # the full eval (needs ANTHROPIC_API_KEY)
python scripts/eval_ask_aree.py --regrade        # re-run only the judge on stored answers
```

The eval refuses to run when the eval script, the questions, or the chat and query code
it exercises changed since a person last reviewed them and ran
`python scripts/eval_ask_aree.py --approve-harness`. Results go to
`.claude/hillclimb/ask-aree/baseline/` (or `v1/`, `v2/`, ... with `--variant`):
`results.jsonl` with every grade and its reason, full transcripts in `traces/`, failed
attempts in `errors.jsonl` (never scored as wrong answers), and `summary.json`. A full
run (26 questions × 2 repeats) costs about $5 and takes about 5 minutes; at that size,
differences in pass rate smaller than about 14 points are within noise.

**Baseline** (2026-10-21; Claude Opus 5.5 at medium effort; data build of 2026-10-21):
67% of answers pass every check (±16 points). Facts, grounded ids, answers, no overclaim
and case criterion are all at 100%: no answer invented a study, stated a wrong count,
leaked simulated data or called a candidate validated. Every failure is on faithfulness.
Most failures are answers that generalize from the top rows of a shortened result to
the whole set ("every replicated gene comes from these two studies" after seeing 15 of
30); about five are outright errors, such as describing a study the chat never looked
up. Cost: $3.92 for 52 answers (median 16 s each), plus $1.23 for grading. The judge
also varies: re-grading the same answers flipped two verdicts, so treat single-case
differences with caution.

**Variants tried** (same questions, model, judge and rubric; paired against the baseline):

| Variant | Change | Pass | Paired change (95% CI) | Median length | Cost per answer |
|---|---|---|---|---|---|
| baseline | | 67% | | 292 words | $0.075 |
| v1 | Results say when they're partial; ordering documented; replication counts by feature type | 60% | −8 (−23 to +7) | 292 words | $0.075 |
| v2 (shipped) | v1, plus answers under 200 words with a few directly sourced facts | 69% | +2 (−14 to +18) | 220 words | $0.058 |

Neither change measurably moved the pass rate. In both the baseline and v1,
faithfulness fell as answers grew longer: 100% under 200 words, about 55% over 300.
v2 was kept because it matches the baseline's quality with shorter answers, 23% lower
cost and faster replies (median 13 s). One cost of brevity: answers about replication
sometimes name no conflicting example, because results list consistent features first.
At 26 questions × 2 runs this eval can't see differences much under 15 points, so
further tuning needs more questions and runs, and a check of the judge against human
grading.

## Deploy to Cloud Run

The service runs as one container on Google Cloud Run
([`containers/chat/Dockerfile`](https://github.com/sr320/AREE/blob/main/containers/chat/Dockerfile)).
The image holds the code, the study registry and the tables the site publishes,
fetched when the image is built. Dependency versions are locked in
`containers/chat/requirements.lock`. CI builds the image and checks that it starts on
every pull request.

After each successful site build, the **Deploy Ask AREE chat service** workflow
(`.github/workflows/chat-deploy.yml`) builds a fresh image, so the chat answers from the
same build as the Explore page. It then deploys the image and checks that it answers.
GitHub signs in to Google through Workload Identity Federation, so no Google key file
exists anywhere. The Anthropic API key lives only in Google Secret Manager.

The service runs on one instance, which keeps the rate limits and daily budget exact
and serves many visitors at once. It scales to zero when idle, so the first question
after a quiet spell waits a few seconds for it to start.

### One-time setup

1. **A Google Cloud project with billing.** If your Google account belongs to an
   organization (for example UW), its policies may forbid public (`allUsers`) access to
   Cloud Run services. In that case, ask its administrators, or use a project outside
   the organization.
2. **An Anthropic API key in a dedicated workspace**, with a monthly spend limit set on
   that workspace in the Claude Console. That limit is the hard cap behind the chat's
   own daily budget. The key must belong to a workspace: an organization-level key is
   refused.
3. **Run the setup script in Google Cloud Shell** (the terminal button in the Cloud
   console). Upload [`scripts/setup_cloud_run.sh`](https://github.com/sr320/AREE/blob/main/scripts/setup_cloud_run.sh)
   or paste it, then run:

   ```bash
   bash setup_cloud_run.sh PROJECT_ID
   ```

   It enables the APIs and creates the image repository and the secret (you paste the
   key; it isn't shown). It also creates a runtime account that can read only that
   secret, and a deploy account that only GitHub Actions in `sr320/AREE` can use. It is
   safe to run again.
4. **Set the GitHub repository variables** the script prints (`GCP_PROJECT_ID`,
   `GCP_REGION`, `GCP_WIF_PROVIDER`, `GCP_DEPLOY_SA`, `GCP_RUNTIME_SA`; none is a
   secret). Optionally set `AREE_CHAT_DAILY_BUDGET_USD`; the default is 10.
5. **Run the workflow** from the Actions tab. Its summary shows the service's address.

### Going live

1. Set `PRODUCTION_API` at the top of `docs/ask/ask.js` to the service's address.
2. Add the page to the site's navigation in `docs/_quarto.yml`.
3. Optionally, turn on the Cloudflare Turnstile bot check: put its secret key in Secret
   Manager as `aree-turnstile-secret`, pass it to the service as
   `AREE_CHAT_TURNSTILE_SECRET` (the `secrets` list in the workflow), and set
   `TURNSTILE_SITE_KEY` in `docs/ask/ask.js`.

To change the API key later, add a new secret version
(`gcloud secrets versions add aree-anthropic-api-key --data-file=-`) and re-run the
workflow. To stop the chat at once, set `AREE_CHAT_DAILY_BUDGET_USD` to 0 and re-run the
workflow, or delete the Cloud Run service. The page then reports that the chat isn't
reachable, and the Explore page is unaffected.
