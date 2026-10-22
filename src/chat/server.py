"""The "Ask AREE" chat service: a small web API in front of Claude and the query layer.

    POST /api/chat     {"message": ..., "conversation_id": ..., "turnstile_token": ...}
                       -> a text/event-stream of chat.conversation events
    GET  /api/health   -> status, model and data build

It holds the Anthropic API key, so the static site never does. Every setting is an
environment variable (`Settings.from_env`); run it with `aree chat-serve`, or in a
container with `uvicorn --factory chat.server:app_from_env`.
"""
from __future__ import annotations

import json
import logging
import os
import time
import urllib.parse
import urllib.request
from dataclasses import dataclass, field
from pathlib import Path

import anthropic
from fastapi import FastAPI, Request
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse, StreamingResponse
from pydantic import BaseModel

import query
from chat.conversation import ChatConfig, Conversations, run_turn
from chat.limits import DailyBudget, RateLimiter
from chat.prompt import system_prompt

log = logging.getLogger("aree.chat")

TURNSTILE_VERIFY_URL = "https://challenges.cloudflare.com/turnstile/v0/siteverify"
LOCALHOST_ORIGINS = r"https?://(localhost|127\.0\.0\.1)(:\d+)?"


def _env_bool(name: str, default: bool) -> bool:
    value = os.environ.get(name)
    return default if value is None else value.strip().lower() in {"1", "true", "yes"}


@dataclass
class Settings:
    data_dir: Path | None = None
    model: str = "claude-opus-5-5"
    effort: str = "medium"
    daily_budget_usd: float = 10.0
    per_ip_hourly: int = 10
    per_ip_daily: int = 30
    max_question_chars: int = 1000
    max_tool_rounds: int = 6
    max_turns: int = 10
    allowed_origins: list[str] = field(default_factory=lambda: ["https://sr320.github.io"])
    allow_localhost: bool = False
    turnstile_secret: str | None = None
    # Behind Cloud Run's proxy the client address is the first X-Forwarded-For entry.
    trust_proxy: bool = False
    log_questions: bool = False

    @classmethod
    def from_env(cls) -> Settings:
        env = os.environ.get
        data_dir = env("AREE_DATA_DIR")
        return cls(
            data_dir=Path(data_dir) if data_dir else None,
            model=env("AREE_CHAT_MODEL", cls.model),
            effort=env("AREE_CHAT_EFFORT", cls.effort),
            daily_budget_usd=float(env("AREE_CHAT_DAILY_BUDGET_USD", cls.daily_budget_usd)),
            per_ip_hourly=int(env("AREE_CHAT_PER_IP_HOURLY", cls.per_ip_hourly)),
            per_ip_daily=int(env("AREE_CHAT_PER_IP_DAILY", cls.per_ip_daily)),
            max_question_chars=int(env("AREE_CHAT_MAX_QUESTION_CHARS", cls.max_question_chars)),
            max_tool_rounds=int(env("AREE_CHAT_MAX_TOOL_ROUNDS", cls.max_tool_rounds)),
            max_turns=int(env("AREE_CHAT_MAX_TURNS", cls.max_turns)),
            allowed_origins=[o.strip() for o in env("AREE_CHAT_ALLOWED_ORIGINS", "https://sr320.github.io")
                             .split(",") if o.strip()],
            allow_localhost=_env_bool("AREE_CHAT_ALLOW_LOCALHOST", False),
            turnstile_secret=env("AREE_CHAT_TURNSTILE_SECRET") or None,
            trust_proxy=_env_bool("AREE_CHAT_TRUST_PROXY", False),
            log_questions=_env_bool("AREE_CHAT_LOG_QUESTIONS", False),
        )

    def chat_config(self) -> ChatConfig:
        return ChatConfig(model=self.model, effort=self.effort, max_tool_rounds=self.max_tool_rounds,
                          max_turns=self.max_turns, fallbacks=not self.model.startswith("claude-haiku"))


class ChatRequest(BaseModel):
    message: str
    conversation_id: str | None = None
    turnstile_token: str | None = None


def _client_ip(request: Request, trust_proxy: bool) -> str:
    if trust_proxy:
        forwarded = request.headers.get("x-forwarded-for", "")
        if forwarded:
            return forwarded.split(",")[0].strip()
    return request.client.host if request.client else "unknown"


def verify_turnstile(secret: str, token: str | None, ip: str) -> bool:
    """Ask Cloudflare whether the page's Turnstile token is genuine."""
    if not token:
        return False
    body = urllib.parse.urlencode({"secret": secret, "response": token, "remoteip": ip}).encode()
    try:
        with urllib.request.urlopen(TURNSTILE_VERIFY_URL, data=body, timeout=10) as response:
            return bool(json.loads(response.read()).get("success"))
    except OSError:
        log.warning("Turnstile verification unreachable")
        return False


def _refuse(status: int, code: str, message: str) -> JSONResponse:
    return JSONResponse({"error": code, "message": message}, status_code=status)


def _sse(event: dict) -> str:
    return f"event: {event['type']}\ndata: {json.dumps(event, default=str)}\n\n"


def create_app(settings: Settings | None = None, client=None, verify=verify_turnstile) -> FastAPI:
    """The web app. `client` and `verify` are replaceable for tests."""
    settings = settings or Settings.from_env()
    sources = query.Sources.from_site_data(settings.data_dir) if settings.data_dir else query.DEFAULT_SOURCES
    client = client or anthropic.Anthropic()
    config = settings.chat_config()
    system = system_prompt(settings.data_dir)
    conversations = Conversations()
    limiter = RateLimiter(settings.per_ip_hourly, settings.per_ip_daily)
    budget = DailyBudget(settings.daily_budget_usd)
    data_meta = settings.data_dir / "meta.json" if settings.data_dir else None
    data_build = json.loads(data_meta.read_text()).get("generated") if data_meta and data_meta.exists() else None

    app = FastAPI(title="Ask AREE", docs_url=None, redoc_url=None, openapi_url=None)
    app.add_middleware(CORSMiddleware, allow_origins=settings.allowed_origins,
                       allow_origin_regex=LOCALHOST_ORIGINS if settings.allow_localhost else None,
                       allow_methods=["GET", "POST"], allow_headers=["Content-Type"])

    @app.get("/api/health")
    def health():
        return {"ok": True, "model": settings.model, "data_build": data_build,
                "accepting_questions": not budget.exhausted(),
                "turnstile": settings.turnstile_secret is not None}

    @app.post("/api/chat")
    def chat(body: ChatRequest, request: Request):
        ip = _client_ip(request, settings.trust_proxy)
        message = body.message.strip()
        if not message:
            return _refuse(400, "empty", "Type a question first.")
        if len(message) > settings.max_question_chars:
            return _refuse(400, "too_long", f"Questions are limited to {settings.max_question_chars} "
                                            "characters.")
        if settings.turnstile_secret and not verify(settings.turnstile_secret, body.turnstile_token, ip):
            return _refuse(403, "verification", "The page could not verify this request. Reload and try again.")
        if budget.exhausted():
            return _refuse(503, "daily_budget", "Ask AREE has reached today's limit. The Explore page "
                                                "still answers queries: https://sr320.github.io/AREE/explore.html")
        refused = limiter.check(ip)
        if refused:
            return _refuse(429, "rate_limit", refused + " The Explore page has no limit.")
        conversation = conversations.get_or_create(body.conversation_id)

        def events():
            # Taken here, not in the handler, so a response that is never streamed holds nothing.
            if not conversation.lock.acquire(blocking=False):
                yield _sse({"type": "error", "code": "busy",
                            "message": "The previous answer in this conversation is still being written."})
                return
            started = time.monotonic()
            done = None
            try:
                yield _sse({"type": "start", "conversation_id": conversation.id})
                for event in run_turn(client, conversation, message, sources, config, system, budget.add):
                    if event["type"] == "done":
                        done = event
                    yield _sse(event)
            except anthropic.RateLimitError:
                yield _sse({"type": "error", "code": "busy",
                            "message": "Ask AREE is busy right now. Try again in a minute."})
            except anthropic.APIError:
                log.exception("Claude API error")
                yield _sse({"type": "error", "code": "api",
                            "message": "Ask AREE could not reach its language model. Try again shortly."})
            except Exception:
                # Anything else (a missing API key surfaces as a TypeError) still ends the
                # stream with an event the page can show, rather than a silent cut.
                log.exception("chat turn failed")
                yield _sse({"type": "error", "code": "internal",
                            "message": "Ask AREE hit an internal error. Try again shortly."})
            finally:
                conversation.lock.release()
                record = {"event": "chat_turn", "conversation": conversation.id,
                          "completed": done is not None, "latency_s": round(time.monotonic() - started, 2)}
                if done:
                    record.update(model=done["model"], tools=done["tools"], usage=done["usage"],
                                  cost_usd=done["cost_usd"])
                if settings.log_questions:
                    record["question"] = message
                log.info(json.dumps(record))

        return StreamingResponse(events(), media_type="text/event-stream",
                                 headers={"Cache-Control": "no-cache", "X-Accel-Buffering": "no"})

    return app


def app_from_env() -> FastAPI:
    """The container's entry point (`uvicorn --factory chat.server:app_from_env`): settings
    from the environment, and one JSON line per chat turn on stdout for the platform's logs."""
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    return create_app()
