// Ask AREE: the chat box on docs/ask.qmd. It streams answers from the Ask AREE
// service (src/chat, docs/ask_aree.md) as server-sent events over fetch.
(() => {
  "use strict";

  // The deployed service's address (Cloud Run, docs/ask_aree.md#deploy-to-cloud-run).
  // Empty means the public page says the chat isn't live.
  const PRODUCTION_API = "https://aree-chat-4l635c2jtq-uw.a.run.app";
  // Cloudflare Turnstile site key. Normally empty: the service reports its own key in
  // /api/health when it requires a bot check. Set it only to override that.
  const TURNSTILE_SITE_KEY = "";

  const LOCAL_ADDRESS = /^https?:\/\/(localhost|127\.0\.0\.1)(:\d+)?$/;
  const studyUrl = id => `https://github.com/sr320/AREE/blob/main/registry/studies/${encodeURIComponent(id)}.yaml`;
  const STUDY_TOKEN = /\b[A-Z][A-Z0-9]*(?:_[A-Z0-9]+)+\b/g;
  const MAX_TRACE_CHARS = 20000;

  const TOOL_LABELS = {
    aree_overview: "an overview of what AREE holds",
    list_vocabulary: "the controlled vocabularies",
    list_studies: "the registered studies",
    get_study: "a study record",
    search_features: "matching genes and features",
    get_evidence: "evidence records",
    find_replicated_features: "features significant in several studies",
    list_candidates: "ranked candidate biomarkers",
    get_evidence_card: "an evidence card",
  };

  const $ = id => document.getElementById(id);
  const root = $("ask-aree");
  if (!root) return;
  const statusLine = $("aa-status");
  const messages = $("aa-messages");
  const form = $("aa-form");
  const input = $("aa-input");
  const send = $("aa-send");

  // A page served from localhost talks to a local `aree chat-serve`; ?api= may
  // point it elsewhere, but only at another local address.
  const override = new URLSearchParams(location.search).get("api");
  const onLocalPage = ["localhost", "127.0.0.1"].includes(location.hostname);
  const API = override && LOCAL_ADDRESS.test(override) ? override
    : onLocalPage ? "http://127.0.0.1:8000" : PRODUCTION_API;

  let conversationId = null;
  let busy = false;
  let available = false;
  let turnstileRequired = false;
  let turnstileWidget = null;
  let turnstileToken = null;
  const knownStudies = new Set();

  function el(tag, className, text) {
    const node = document.createElement(tag);
    if (className) node.className = className;
    if (text !== undefined) node.textContent = text;
    return node;
  }

  function setStatus(text, kind = "muted") {
    statusLine.className = `aa-status small text-${kind}`;
    statusLine.textContent = text;
  }

  function setBusy(value) {
    busy = value;
    send.disabled = value || !available;
    input.disabled = !available;
    root.querySelectorAll(".aa-example").forEach(b => { b.disabled = value || !available; });
  }

  // ------------------------------------------------------------ rendering

  // Every study id a tool returned in this conversation; only these are linked.
  // Ids appear as `study_id` values, in `studies` lists, and in pipe-joined `studies`
  // strings; list items keep the key of the list they are in.
  function collectStudies(value, key) {
    if (typeof value === "string") {
      if (key === "study_id" || key === "known_study_ids") knownStudies.add(value);
      else if (key === "studies") value.split("|").forEach(s => knownStudies.add(s.trim()));
    } else if (Array.isArray(value)) {
      value.forEach(item => collectStudies(item, key));
    } else if (value && typeof value === "object") {
      for (const [k, v] of Object.entries(value)) collectStudies(v, k);
    }
  }

  function linkStudies(container) {
    const walker = document.createTreeWalker(container, NodeFilter.SHOW_TEXT, {
      acceptNode: node => node.parentElement.closest("a, pre") ? NodeFilter.FILTER_REJECT : NodeFilter.FILTER_ACCEPT,
    });
    const nodes = [];
    while (walker.nextNode()) nodes.push(walker.currentNode);
    for (const node of nodes) {
      const text = node.nodeValue;
      const matches = [...text.matchAll(STUDY_TOKEN)].filter(m => knownStudies.has(m[0]));
      if (!matches.length) continue;
      const fragment = document.createDocumentFragment();
      let last = 0;
      for (const m of matches) {
        fragment.append(text.slice(last, m.index));
        const link = el("a", "aa-study", m[0]);
        link.href = studyUrl(m[0]);
        link.target = "_blank";
        link.rel = "noopener";
        link.title = "Open this study's registration record";
        fragment.append(link);
        last = m.index + m[0].length;
      }
      fragment.append(text.slice(last));
      node.replaceWith(fragment);
    }
  }

  function renderMarkdown(target, markdown) {
    // Claude's text is untrusted: parse, then sanitize before it reaches the page.
    target.innerHTML = window.DOMPurify.sanitize(window.marked.parse(markdown));
    target.querySelectorAll("table").forEach(t => t.classList.add("table", "table-sm"));
    target.querySelectorAll("a[href]").forEach(a => {
      if (/^https?:/.test(a.getAttribute("href"))) { a.target = "_blank"; a.rel = "noopener"; }
    });
    linkStudies(target);
  }

  function addMessage(role) {
    const message = el("div", `aa-msg aa-${role}`);
    messages.append(message);
    return message;
  }

  function newAnswer() {
    const message = addMessage("assistant");
    const answer = el("div", "aa-answer");
    const progress = el("div", "aa-progress small text-muted", "Thinking…");
    const trace = el("details", "aa-trace");
    const summary = el("summary", "small");
    trace.append(summary);
    trace.hidden = true;
    message.append(answer, progress, trace);
    return { message, answer, progress, trace, summary, segment: null, calls: {}, nCalls: 0,
             started: performance.now(), finished: false };
  }

  function appendText(state, text) {
    if (!state.segment) {
      state.segment = { text: "", node: el("div", "aa-segment"), pending: false };
      state.answer.append(state.segment.node);
    }
    const segment = state.segment;
    segment.text += text;
    if (!segment.pending) {
      segment.pending = true;
      requestAnimationFrame(() => { segment.pending = false; renderMarkdown(segment.node, segment.text); });
    }
  }

  function addTraceEntry(state, event) {
    state.nCalls += 1;
    state.segment = null;  // text after a query starts a new paragraph
    state.trace.hidden = false;
    state.summary.textContent = `How this answer was found (${state.nCalls} ${state.nCalls === 1 ? "query" : "queries"})`;
    state.progress.textContent = `Looking up ${TOOL_LABELS[event.name] || event.name}…`;
    const entry = el("div", "aa-call");
    entry.append(el("div", "aa-call-name", event.name));
    entry.append(el("pre", "aa-call-input", JSON.stringify(event.input || {}, null, 2)));
    const result = el("pre", "aa-call-result", "running…");
    entry.append(result);
    state.trace.append(entry);
    state.calls[event.id] = result;
  }

  function fillTraceEntry(state, event) {
    const target = state.calls[event.id];
    let text = event.content;
    try {
      const parsed = JSON.parse(event.content);
      collectStudies(parsed);
      text = JSON.stringify(parsed, null, 2);
    } catch (_) { /* shown as sent */ }
    if (text.length > MAX_TRACE_CHARS) text = text.slice(0, MAX_TRACE_CHARS) + "\n… (shortened here; the model received the full result)";
    if (target) {
      target.textContent = text;
      target.classList.toggle("aa-call-error", Boolean(event.is_error));
    }
    state.progress.textContent = "Reading the results…";
  }

  function showNotice(state, text, kind) {
    const notice = el("div", `alert alert-${kind} py-2 small mb-2`, text);
    state.message.insertBefore(notice, state.progress);
  }

  function finish(state, interrupted) {
    if (state.finished) return;
    state.finished = true;
    if (state.segment) renderMarkdown(state.segment.node, state.segment.text);
    state.answer.querySelectorAll(".aa-segment").forEach(n => linkStudies(n));
    state.progress.remove();
    if (interrupted) showNotice(state, "The answer was interrupted. Try asking again.", "warning");
    if (state.nCalls) {
      const seconds = Math.round((performance.now() - state.started) / 1000);
      state.summary.textContent += ` · ${seconds} s`;
    }
  }

  // ------------------------------------------------------------ the service

  async function* readEvents(response) {
    const reader = response.body.getReader();
    const decoder = new TextDecoder();
    let buffer = "";
    for (;;) {
      const { value, done } = await reader.read();
      if (done) break;
      buffer += decoder.decode(value, { stream: true });
      let end;
      while ((end = buffer.indexOf("\n\n")) >= 0) {
        const block = buffer.slice(0, end);
        buffer = buffer.slice(end + 2);
        const data = block.split("\n").filter(l => l.startsWith("data: ")).map(l => l.slice(6)).join("\n");
        if (data) yield JSON.parse(data);
      }
    }
  }

  function handle(state, event) {
    switch (event.type) {
      case "start": conversationId = event.conversation_id; break;
      case "text": appendText(state, event.text); break;
      case "tool_call": addTraceEntry(state, event); break;
      case "tool_result": fillTraceEntry(state, event); break;
      case "round_reset":
        if (state.segment) { state.segment.node.remove(); state.segment = null; }
        break;
      case "notice": showNotice(state, event.text, "warning"); break;
      case "error":
        showNotice(state, event.message, "danger");
        if (event.code === "conversation_limit") conversationId = null;
        finish(state, false);
        break;
      case "done": finish(state, false); break;
      default: break;
    }
  }

  async function ask(question) {
    if (busy || !available) return;
    if (turnstileRequired && !turnstileToken) {
      setStatus("The verification check is still loading. Try again in a moment.", "warning");
      return;
    }
    setBusy(true);
    addMessage("user").textContent = question;
    const state = newAnswer();
    state.message.scrollIntoView({ behavior: "smooth", block: "start" });
    try {
      const response = await fetch(`${API}/api/chat`, {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ message: question, conversation_id: conversationId, turnstile_token: turnstileToken }),
      });
      resetTurnstile();
      if (!response.ok) {
        const body = await response.json().catch(() => ({}));
        showNotice(state, body.message || `Ask AREE answered with an error (${response.status}).`, "danger");
        finish(state, false);
        return;
      }
      for await (const event of readEvents(response)) handle(state, event);
      finish(state, !state.finished);
    } catch (_) {
      showNotice(state, "Couldn't reach Ask AREE. Check your connection and try again.", "danger");
      finish(state, false);
    } finally {
      setBusy(false);
      input.focus();
    }
  }

  // ------------------------------------------------------------ bot check

  function resetTurnstile() {
    turnstileToken = null;
    if (turnstileWidget !== null && window.turnstile) window.turnstile.reset(turnstileWidget);
  }

  function loadTurnstile(siteKey) {
    const script = document.createElement("script");
    script.src = "https://challenges.cloudflare.com/turnstile/v0/api.js?render=explicit";
    script.async = true;
    script.onload = () => {
      turnstileWidget = window.turnstile.render("#aa-turnstile", {
        sitekey: siteKey,
        callback: token => { turnstileToken = token; },
        "expired-callback": () => { turnstileToken = null; },
      });
    };
    document.head.append(script);
  }

  // ------------------------------------------------------------ start

  async function connect() {
    if (!API) {
      setStatus("Ask AREE isn't live yet. Until it is, the Explore page answers the same questions with filters and SQL.", "secondary");
      setBusy(false);
      return;
    }
    try {
      const response = await fetch(`${API}/api/health`);
      const health = await response.json();
      if (!health.accepting_questions) {
        setStatus("Ask AREE has reached today's limit. The Explore page still works.", "warning");
      } else if (health.turnstile && !(TURNSTILE_SITE_KEY || health.turnstile_site_key)) {
        setStatus("Ask AREE requires a verification check this page isn't configured for.", "danger");
      } else {
        available = true;
        turnstileRequired = Boolean(health.turnstile);
        if (turnstileRequired) loadTurnstile(TURNSTILE_SITE_KEY || health.turnstile_site_key);
        const build = health.data_build ? ` published ${health.data_build.slice(0, 10)}` : "";
        setStatus(`Answers come from AREE's real studies (data${build}).`);
      }
    } catch (_) {
      setStatus("Ask AREE isn't reachable right now. The Explore page still works.", "warning");
    }
    setBusy(false);
  }

  form.addEventListener("submit", e => {
    e.preventDefault();
    const question = input.value.trim();
    if (!question) return;
    input.value = "";
    ask(question);
  });
  input.addEventListener("keydown", e => {
    if (e.key === "Enter" && !e.shiftKey && !e.isComposing) { e.preventDefault(); form.requestSubmit(); }
  });
  root.querySelectorAll(".aa-example").forEach(button => {
    button.addEventListener("click", () => ask(button.textContent.trim()));
  });
  $("aa-new").addEventListener("click", () => {
    if (busy) return;
    conversationId = null;
    knownStudies.clear();
    messages.replaceChildren();
    input.focus();
  });

  connect();
})();
